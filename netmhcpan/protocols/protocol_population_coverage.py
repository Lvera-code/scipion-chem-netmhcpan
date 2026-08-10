# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Enzo Sierra (enzogael57@gmail.com)
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
This protocol annotates each input ROI with a REAL, allele-frequency-
weighted population coverage estimate, using the same methodology as the
IEDB Population Coverage tool (Bui et al. 2006, BMC Bioinformatics).

Ported from the standalone B-Cell-Epitope-Prediction repo's
src/engines/population_coverage.py (Fase 4-bis, Carmen Elena Gomez
feedback 2026-07-30): 'n_promiscuous_alleles' (ProtNetMHCpanPromiscuity/
ProtNetMHCIIpanPromiscuity) counts how many reference-panel alleles a
candidate hits SB/WB, treating every allele equally -- but real HLA allele
frequencies vary enormously: 3 very common alleles is not the same as 3
rare ones, even though both give the same count. This protocol computes
the real fraction of the population likely to present a candidate via AT
LEAST one of its promiscuous alleles, assuming Hardy-Weinberg equilibrium
per locus and independence between loci (same simplification IEDB's own
basic tool makes).

Requires '_promiscuousAlleles' on each input ROI (set by
ProtNetMHCpanPromiscuity/ProtNetMHCIIpanPromiscuity's own
createOutputStep) -- the comma-joined list of alleles that ROI hit SB/WB.

Reference allele-frequency table (world_pooled_afnd.csv, bundled in this
plugin under netmhcpan/data/): derived from the Allele Frequency Net
Database via the MIT-licensed github.com/slowkow/allelefrequencies
mirror -- unlike LANL/CATNAP/IEDB's reference data, this one IS safe to
redistribute, so no scipion.conf variable is needed. 2 known data gaps
(documented in the CSV itself, see the standalone module's docstring for
detail): DRB3/4/5 secondary genes have no frequency (excluded from the
calculation, never assumed 0); DQ/DP combo alleles approximate their
frequency as the product of the two chain frequencies.

Purely informative: never filters (same treatment as
ProtLANLCATNAPCrossref/ProtIEDBCrossref/ProtBLASTPanelConservation).
"""

import os
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd
from pwchem.objects import SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.object import Float
from pyworkflow.protocol import params
from pyworkflow.utils import Message

_DEFAULT_FREQ_TABLE_PATH = Path(__file__).resolve().parent.parent / 'data' / 'world_pooled_afnd.csv'


def phenotypic_frequency(allele_frequencies: List[float]) -> float:
    """P(an individual carries AT LEAST ONE of these alleles in one locus), Hardy-Weinberg.

    Sum is capped at 1.0 before applying the formula (a sum >1 can only
    come from rounding/overlap in the source data, never a real
    probability >1).
    """
    total = min(sum(allele_frequencies), 1.0)
    return 1 - (1 - total) ** 2


def combined_population_coverage(phenotypic_frequencies_by_locus: List[float]) -> float:
    """Population coverage combining several loci, assuming independence between them."""
    coverage = 1.0
    for pf in phenotypic_frequencies_by_locus:
        coverage *= (1 - pf)
    return 1 - coverage


def load_allele_frequencies(path: Optional[str] = None) -> pd.DataFrame:
    """Load the (panel_allele -> locus, frequency) reference table.

    Rows with an empty frequency (known data gaps, see module docstring)
    are kept as NaN -- population_coverage_for_alleles excludes them
    explicitly rather than treating them as 0.
    """
    resolved = path or str(_DEFAULT_FREQ_TABLE_PATH)
    if not resolved or not Path(resolved).is_file():
        return pd.DataFrame(columns=['locus', 'frequency']).rename_axis('panel_allele')
    table = pd.read_csv(resolved)
    return table.set_index('panel_allele')[['locus', 'frequency']]


def population_coverage_for_alleles(alleles: List[str], freq_table: pd.DataFrame) -> Tuple[Optional[float], List[str]]:
    """Population coverage for one candidate, given the alleles it hit SB/WB.

    Returns:
        (coverage_pct, excluded_alleles): coverage_pct in [0,100], None if
        NONE of 'alleles' has a known frequency. excluded_alleles are
        alleles missing from freq_table (known data gaps) -- excluded from
        the calculation, never assumed to have frequency 0.
    """
    if not alleles or freq_table.empty:
        return None, list(alleles)

    known = freq_table.reindex(alleles).dropna(subset=['frequency'])
    excluded = [a for a in alleles if a not in known.index]
    if known.empty:
        return None, excluded

    phenotypic_by_locus = [
        phenotypic_frequency(group['frequency'].tolist())
        for _, group in known.groupby('locus')
    ]
    coverage = combined_population_coverage(phenotypic_by_locus)
    return round(coverage * 100, 2), excluded


class ProtPopulationCoverage(EMProtocol):
    """
    AI Generated:

    Annotates (does NOT filter) every input ROI with
    'population_coverage_pct' (Bui et al. 2006 methodology), computed from
    its '_promiscuousAlleles' attribute (set by ProtNetMHCpanPromiscuity/
    ProtNetMHCIIpanPromiscuity) and a bundled world-pooled allele-frequency
    reference table.

    Output
    ------
    outputROIs: the same SetOfSequenceROIs as the input, each ROI
    annotated with '_populationCoveragePct' (Float, NaN if none of its
    promiscuous alleles has a known frequency -- e.g. all DRB3/4/5).
    """

    _label = 'population coverage'

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       label='Sequence ROIs: ',
                       help="Candidates already annotated with '_promiscuousAlleles' "
                            '(output of ProtNetMHCpanPromiscuity/ProtNetMHCIIpanPromiscuity).')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.annotateStep)

    # ---------------------------------- Steps -----------------------------------

    def _getRois(self):
        # Iterating a Scipion SetOfXXX reuses the same Python object per row
        # (the underlying sqlite cursor): each item must be cloned when
        # materialized into a list, or all N references end up pointing to
        # the cursor's last state.
        return [roi.clone() for roi in self.inputROIs.get()]

    def annotateStep(self):
        rois = self._getRois()
        freqTable = load_allele_frequencies()

        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for roi in rois:
            allelesStr = getattr(roi, '_promiscuousAlleles', None)
            alleles = [a for a in (allelesStr.get() if allelesStr is not None else '').split(',') if a]
            coverage, _excluded = population_coverage_for_alleles(alleles, freqTable)
            roi._populationCoveragePct = Float(coverage) if coverage is not None else Float(float('nan'))
            outROIs.append(roi)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            self._defineSourceRelation(self.inputROIs, outROIs)

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        errors = []
        if not os.path.isfile(_DEFAULT_FREQ_TABLE_PATH):
            errors.append(f"Bundled allele frequency table not found at '{_DEFAULT_FREQ_TABLE_PATH}'.")
        return errors

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            if outROIs is not None:
                covered = [roi._populationCoveragePct.get() for roi in outROIs
                           if roi._populationCoveragePct.get() == roi._populationCoveragePct.get()]  # skip NaN
                if covered:
                    summary.append(f'{len(covered)}/{len(outROIs)} candidate(s) with a computable population '
                                    f'coverage estimate (mean {sum(covered) / len(covered):.2f}%).')
        return summary
