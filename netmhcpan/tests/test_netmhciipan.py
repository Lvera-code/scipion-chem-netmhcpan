from pyworkflow.tests import setupTestProject, BaseTest

from pwem.protocols import ProtImportSequence
from pwchem.protocols import ProtDefineSeqROI

from ..protocols import ProtNetMHCIIpanPromiscuity


class TestNetMHCIIpanPromiscuity(BaseTest):
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

    # Real reference values: same (start, end) windows and NetMHCIIpan-4.3
    # results as the original end-to-end pipeline test (see
    # B-Cell-Epitope-Prediction / scipion-chem-bcellepitope history). Since
    # NetMHCIIpan is a static predictor, feeding the exact same submitted
    # peptide windows directly (bypassing the retired BepiPred/EpiDope/
    # consensus/BLAST phases) reproduces the exact same output.
    EXPECTED = sorted([
        (88, 123, 'WKNDMVEQM', 3, 27),
        (309, 323, 'IRIQRGPGR', 7, 27),
        (316, 330, 'RAFVTIGKI', 4, 27),
        (317, 331, 'FVTIGKIGN', 4, 27),
        (318, 332, 'FVTIGKIGN', 6, 27),
        (319, 333, 'FVTIGKIGN', 5, 27),
        (338, 352, 'KWNATLKQI', 3, 27),
        (339, 353, 'KWNATLKQI', 5, 27),
        (340, 354, 'WNATLKQIA', 4, 27),
        (349, 363, 'KLREQFGNN', 3, 27),
        (350, 364, 'LREQFGNNK', 3, 27),
        (639, 653, 'NNYTSLIHS', 3, 27),
        (640, 654, 'YTSLIHSLI', 4, 27),
        (642, 656, 'LIHSLIEES', 5, 27),
        (643, 657, 'LIHSLIEES', 8, 27),
        (644, 658, 'IHSLIEESQ', 7, 27),
        (656, 670, 'KNEQELLEL', 3, 27),
    ])
    WINDOWS = [(start, end) for start, end, *_ in EXPECTED]

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
        protNetMHCIIpan = self.newProtocol(ProtNetMHCIIpanPromiscuity)
        protNetMHCIIpan.inputROIs.set(self.protSeedROIs)
        protNetMHCIIpan.inputROIs.setExtended('outputROIs')
        self.launchProtocol(protNetMHCIIpan, wait=True)

        outROIs = getattr(protNetMHCIIpan, 'outputROIs', None)
        self.assertIsNotNone(outROIs)
        got = sorted(
            (roi.getROIIdx(), roi.getROIIdx2(), roi._core9aa.get(),
             roi._nPromiscuousAlleles.get(), roi._nAllelesEvaluated.get())
            for roi in outROIs
        )
        self.assertEqual(got, self.EXPECTED)
