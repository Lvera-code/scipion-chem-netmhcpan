from pyworkflow.tests import setupTestProject, BaseTest

from pwem.protocols import ProtImportSequence
from pwchem.protocols import ProtDefineSeqROI

from ..protocols import ProtNetMHCpanPromiscuity


class TestNetMHCpanPromiscuity(BaseTest):
    NAME = 'GP120_P03377'
    DESCRIPTION = 'ENV_HV1BR (GP120), UniProt P03377'
    AMINOACIDSSEQ = (
        'MRVKEKYQHLWRWGWKWGTMLLGILMICSATEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVW'
        'ATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTPLCVSLKCTDLGNAT'
        'NTNSSNTNSSSGEMMMEKGEIKNCSFNISTSIRGKVQKEYAFFYKLDIIPIDNDTTSYTLTSCNTSVIT'
        'QACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVV'
        'IRSANFTDNAKTIIVQLNQSVEINCTRPNNNTRKSIRIQRGPGRAFVTIGKIGNMRQAHCNISRAKWNA'
        'TLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFYCNSTQLFNSTWFNSTWSTEGSNNTE'
        'GSDTITLPCRIKQFINMWQEVGKAMYAPPISGQIRCSSNITGLLLTRDGGNNNNGSEIFRPGGGDMRDN'
        'WRSELYKYKVVKIEPLGVAPTKAKRRVVQREKRAVGIGALFLGFLGAAGSTMGARSMTLTVQARQLLSG'
        'IVQQQNNLLRAIEAQQHLLQLTVWGIKQLQARILAVERYLKDQQLLGIWGCSGKLICTTAVPWNASWSN'
        'KSLEQIWNNMTWMEWDREINNYTSLIHSLIEESQNQQEKNEQELLELDKWASLWNWFNITNWLWYIKI'
        'FIMIVGGLVGLRIVFAVLSIVNRVRQGYSPLSFQTHLPTPRGPDRPEGIEEEGGERDRDRSIRLVNGSL'
        'ALIWDDLRSLCLFSYHRLRDLLLIVTRIVELLGRRGWEALKYWWNLLQYWSQELKNSAVSLLNATAIAV'
        'AEGTDRVIEVVQGACRAIRHIPRRIRQGLERILL'
    )

    # Windows fed to the protocol: 4 real strong HLA-A/B/C binders found by a
    # full 8-11mer sliding-window scan of this sequence against the 23-allele
    # reference panel (>= 3 promiscuous alleles required, %Rank_EL <= 2.0),
    # PLUS 2 real T-helper (MHC-II-optimized) windows reused from
    # test_netmhciipan.py's own fixture that score ZERO promiscuous MHC-I
    # alleles -- a genuine negative control confirming these two antigen
    # presentation pathways are evaluated independently, not that a filter
    # was tuned to pass. All 6 values (core_9aa, n_promiscuous_alleles) come
    # from a real local run of the netMHCpan-4.2 binary (not estimated),
    # since NetMHCpan is a static predictor: feeding the exact same submitted
    # peptide windows directly reproduces the exact same output.
    EXPECTED = sorted([
        (2, 10, 'RVKEKYQHL', 14, 23),
        (103, 111, 'QMHEDIISL', 15, 23),
        (285, 293, 'NAKTIIVQL', 16, 23),
        (562, 570, 'RAIEAQQHL', 17, 23),
    ])
    # Rejected (0 promiscuous MHC-I alleles): must NOT appear in the output.
    REJECTED_WINDOWS = [(639, 653), (656, 670)]
    WINDOWS = [(start, end) for start, end, *_ in EXPECTED] + REJECTED_WINDOWS

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        setupTestProject(cls)

        cls._runImportSeq()
        cls._waitOutput(cls.protImportSeq, 'outputSequence', sleepTime=5)

        cls.protSeedROIs = cls._runDefSeqROIs(cls.protImportSeq)
        cls._waitOutput(cls.protSeedROIs, 'outputROIs', sleepTime=5)

    @classmethod
    def _runImportSeq(cls):
        kwargs = {
            'inputSequenceName': cls.NAME,
            'inputSequenceDescription': cls.DESCRIPTION,
            'inputRawSequence': cls.AMINOACIDSSEQ,
        }
        cls.protImportSeq = cls.newProtocol(ProtImportSequence, **kwargs)
        cls.proj.launchProtocol(cls.protImportSeq, wait=False)

    @classmethod
    def _runDefSeqROIs(cls, inProt):
        inROIs = '\n'.join(
            '{}) Residues: {{"index": "{}-{}", "residues": "{}", "desc": "None"}}'.format(
                i, start, end, cls.AMINOACIDSSEQ[start - 1:end]
            )
            for i, (start, end) in enumerate(cls.WINDOWS, 1)
        )
        protDefSeqROIs = cls.newProtocol(ProtDefineSeqROI, chooseInput=0, inROIs=inROIs)
        protDefSeqROIs.inputSequence.set(inProt)
        protDefSeqROIs.inputSequence.setExtended('outputSequence')

        cls.proj.launchProtocol(protDefSeqROIs, wait=False)
        return protDefSeqROIs

    def test(self):
        protNetMHCpan = self.newProtocol(ProtNetMHCpanPromiscuity)
        protNetMHCpan.inputROIs.set(self.protSeedROIs)
        protNetMHCpan.inputROIs.setExtended('outputROIs')
        self.launchProtocol(protNetMHCpan, wait=True)

        outROIs = getattr(protNetMHCpan, 'outputROIs', None)
        self.assertIsNotNone(outROIs)
        got = sorted(
            (roi.getROIIdx(), roi.getROIIdx2(), roi._core9aa.get(),
             roi._nPromiscuousAlleles.get(), roi._nAllelesEvaluated.get())
            for roi in outROIs
        )
        self.assertEqual(got, self.EXPECTED)

        # The 2 rejected T-helper windows must not have produced any output ROI.
        gotWindows = {(idx, idx2) for idx, idx2, *_ in got}
        for rejectedWindow in self.REJECTED_WINDOWS:
            self.assertNotIn(rejectedWindow, gotWindows)
