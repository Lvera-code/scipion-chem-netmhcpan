# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Enzo Sierra (enzogael57@gmail.com)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
This protocol is used to evaluate MHC-I (HLA-A/B/C) cytotoxic T-cell
promiscuity of a set of peptide candidates with a local NetMHCpan-4.2
installation.
"""

import os

import pandas as pd
from pwchem.objects import Sequence, SequenceROI, SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.object import Float, Integer, String
from pyworkflow.protocol import params

from .. import Plugin as netmhcpanPlugin
from ..constants import MAX_PEPTIDE_MODE_LENGTH, MIN_PEPTIDE_LENGTH_MHCI, NETMHCPAN_REFERENCE_PANEL, PEPTIDE_LENGTHS
from ..utils.exceptions import NetMHCpanExecutionError
from ..utils.netmhcpan import build_traceback_report, parse_xls


class ProtNetMHCpanPromiscuity(EMProtocol):
    """
    AI Generated:

    Evaluates MHC-I (HLA-A/B/C) cytotoxic T-cell promiscuity of a set of
    peptide candidates using a local NetMHCpan-4.2 installation, and traces
    every resulting candidate back to the region of the parent sequence it
    came from.

    Deliberately kept as an independent protocol, not merged with
    ``ProtNetMHCIIpanPromiscuity`` (MHC-II/T-helper): these are biologically
    distinct antigen presentation pathways (professional antigen-presenting
    cell vs. any nucleated cell; CD8+ cytotoxic vs. CD4+ helper; different
    binding-core length rules), so they are reported as independent
    verdicts, never merged into a single promiscuity figure.

    Overview
    --------
    1. Take an arbitrary SetOfSequenceROIs as input (e.g. the output of any
       upstream epitope-prediction or filtering protocol).
    2. Discard ROIs shorter than the MHC-I binding footprint
       (MIN_PEPTIDE_LENGTH_MHCI, 8 aa).
    3. Route each remaining peptide to NetMHCpan-4.2's exact-peptide mode
       ('-p') if it is short enough, or to protein mode (internal sliding
       window over PEPTIDE_LENGTHS, via '-l') otherwise, since the binary
       crashes on long inputs in exact mode (MAX_PEPTIDE_MODE_LENGTH).
    4. Run NetMHCpan-4.2 against a user-configurable panel of HLA-A/B/C
       alleles (NETMHCPAN_REFERENCE_PANEL by default, 23 alleles: 12 HLA-A/B
       supertypes of Sidney et al. 2008 + 11 HLA-C alleles).
    5. Keep peptides binding at least a minimum number of alleles at or below
       a weak-binder %Rank_EL threshold ("promiscuous" candidates).
    6. Trace every valid candidate back to its parent ROI's position on the
       original sequence, and collapse redundant protein-mode sliding
       windows (same binding core and same allele count keep only the
       lowest %Rank_EL). Steps 1-6 run in netmhcpanStep, which persists the
       traceback report to extra/traceback_report.csv.
    7. createOutputStep reads that report back and builds the output ROIs,
       kept separate so a failure here can be retried without rerunning
       NetMHCpan-4.2.

    Output
    ------
    outputROIs: SetOfSequenceROIs with one ROI per deduplicated, promiscuous
    cytotoxic T-cell candidate, annotated with its binding core, allele
    counts and minimum %Rank_EL.
    """

    _label = 'netmhcpan mhc-i promiscuity'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       label='Sequence ROIs: ',
                       help='Peptide candidates to evaluate for MHC-I cytotoxic T-cell promiscuity.')
        form.addParam('allelePanel', params.StringParam, default=NETMHCPAN_REFERENCE_PANEL,
                       label='HLA-A/B/C allele panel: ',
                       help='Comma-separated alleles, NO spaces (NetMHCpan format). Defaults to '
                            'a 23-allele reference panel (12 HLA-A/B supertypes of Sidney et al. 2008 '
                            '+ 11 HLA-C alleles). To append a specific allele, edit this field '
                            'directly (e.g. add ",HLA-A*11:01" at the end).')

        gGroup = form.addGroup('Promiscuity')
        gGroup.addParam('rankWeak', params.FloatParam, label='Weak binder %Rank threshold: ', default=2.0,
                         help='NetMHCpan %Rank_EL threshold at or below which an allele counts as a '
                              'binder (SB or WB) for promiscuity. Note this is on a different scale '
                              'than NetMHCIIpan (MHC-I and MHC-II %Rank are not comparable 1:1).')
        gGroup.addParam('minPromiscuousAlleles', params.IntParam, label='Min. promiscuous alleles: ', default=3,
                         help='Minimum number of panel alleles that must classify as binders to report '
                              'a peptide as a valid candidate.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.netmhcpanStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------- Steps -----------------------------------

    def _getTracebackReportPath(self):
        return self._getExtraPath('traceback_report.csv')

    def _getPeptideRecords(self):
        # Iterating a Scipion SetOfXXX reuses the same Python object per row
        # (the underlying sqlite cursor): each item must be cloned when
        # materialized into a list, or all N references end up pointing to
        # the cursor's last state.
        rois = [roi.clone() for roi in self.inputROIs.get()]

        parentRecords = []
        peptides = []
        for roi in rois:
            seq = roi.getROISequence()
            if len(seq) < MIN_PEPTIDE_LENGTH_MHCI:
                continue
            peptides.append(seq)
            parentRecords.append({'sequence': seq, 'start': roi.getROIIdx(), 'roi_id': roi.getROIId()})

        return rois, peptides, parentRecords

    def netmhcpanStep(self):
        _, peptides, parentRecords = self._getPeptideRecords()
        if not peptides:
            return

        shortPeptides = [p for p in peptides if len(p) <= MAX_PEPTIDE_MODE_LENGTH]
        longPeptides = [p for p in peptides if len(p) > MAX_PEPTIDE_MODE_LENGTH]

        allelePanel = self.allelePanel.get()
        alleleNames = [a for a in allelePanel.split(',') if a]
        nAlleles = len(alleleNames)
        binary = netmhcpanPlugin.getNetMHCpanBinaryPath()

        reportFrames = []
        if shortPeptides:
            reportFrames.append(self._runMode(
                binary, allelePanel, nAlleles, alleleNames,
                modeArgs=['-p', '-f', self._writePeptideFile(shortPeptides)],
                xlsPath=self._getExtraPath('peptide_mode_output.xls'),
                modeDesc='exact peptide mode',
            ))
        if longPeptides:
            # Unlike NetMHCIIpan's protein mode, NetMHCpan needs the
            # candidate binding-core lengths passed explicitly via '-l':
            # MHC-I cores are 8-11 aa, not a single fixed window.
            reportFrames.append(self._runMode(
                binary, allelePanel, nAlleles, alleleNames,
                modeArgs=['-f', self._writeFragmentsFasta(longPeptides), '-l', PEPTIDE_LENGTHS],
                xlsPath=self._getExtraPath('protein_mode_output.xls'),
                modeDesc='protein mode (sliding window)',
            ))

        reportDf = pd.concat(reportFrames, ignore_index=True) if reportFrames else pd.DataFrame()
        tracebackDf = build_traceback_report(reportDf, parentRecords)
        tracebackDf.to_csv(self._getTracebackReportPath(), index=False)

    def createOutputStep(self):
        tracebackPath = self._getTracebackReportPath()
        if not os.path.isfile(tracebackPath):
            # No peptide met MIN_PEPTIDE_LENGTH_MHCI: netmhcpanStep never ran the program.
            return

        tracebackDf = pd.read_csv(tracebackPath, dtype={'parent_roi_id': str})
        if tracebackDf.empty:
            return

        rois, _, _ = self._getPeptideRecords()
        parentSeqById = {roi.getROIId(): roi._sequence for roi in rois}

        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for row in tracebackDf.itertuples(index=False):
            parentSeq = parentSeqById.get(row.parent_roi_id)
            roiId = f'ROI_{row.start}-{row.end}'
            roiSeq = Sequence(sequence=row.sequence_f5, name=roiId, id=roiId,
                               description='NetMHCpan cytotoxic T-cell candidate')
            newRoi = SequenceROI(sequence=parentSeq, seqROI=roiSeq, roiIdx=row.start, roiIdx2=row.end)
            newRoi._core9aa = String(row.core_9aa)
            newRoi._nPromiscuousAlleles = Integer(int(row.n_promiscuous_alleles))
            newRoi._promiscuousAlleles = String(row.promiscuous_alleles if pd.notna(row.promiscuous_alleles) else '')
            newRoi._nAllelesEvaluated = Integer(int(row.n_alleles_evaluated))
            newRoi._minRankEl = Float(row.min_rank_el)
            outROIs.append(newRoi)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            self._defineSourceRelation(self.inputROIs, outROIs)

    def _writePeptideFile(self, peptides):
        pepPath = self._getExtraPath('peptides.pep')
        with open(pepPath, 'w') as fh:
            fh.write('\n'.join(peptides) + '\n')
        return pepPath

    def _writeFragmentsFasta(self, peptides):
        fastaPath = self._getExtraPath('fragments.fasta')
        with open(fastaPath, 'w') as fh:
            for i, seq in enumerate(peptides):
                fh.write(f'>candidate_{i}\n{seq}\n')
        return fastaPath

    def _runMode(self, binary, allelePanel, nAlleles, alleleNames, modeArgs, xlsPath, modeDesc):
        args = ' '.join(modeArgs) + f' -a {allelePanel} -xls -xlsfile {xlsPath}'
        # Invoked as 'tcsh <script>' (via netmhcpanPlugin's own tcsh conda
        # env), not by executing the wrapper script directly: its shebang
        # ('#!/bin/tcsh -f') hardcodes an absolute system path that a
        # conda-installed tcsh does NOT satisfy even while the env is
        # active -- see TCSH_DIC's docstring in constants.py. Real bug
        # Blanca hit 2026-07-31 on a machine without system-wide tcsh.
        self.runJob(f'{netmhcpanPlugin.getTcshRunPrefix()} {binary}', args)

        if not os.path.isfile(xlsPath):
            raise NetMHCpanExecutionError(
                f"NetMHCpan-4.2 ({modeDesc}) finished without generating the expected output file "
                f"at '{xlsPath}'. Known causes: an input peptide exceeds the limit of the mode used, or "
                "the NMHOME line inside the wrapper script points to an outdated path."
            )
        return parse_xls(xlsPath, nAlleles, self.rankWeak.get(), self.minPromiscuousAlleles.get(), alleleNames)

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        return netmhcpanPlugin.validateNetMHCpanInstallation()

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            n = len(outROIs) if outROIs is not None else 0
            nAlleles = len([a for a in self.allelePanel.get().split(',') if a])
            summary.append(f'{n} promiscuous cytotoxic T-cell candidate(s) (>= '
                            f'{self.minPromiscuousAlleles.get()} of {nAlleles} alleles) after traceback '
                            'and deduplication.')
        return summary
