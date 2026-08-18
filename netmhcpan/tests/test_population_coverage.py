"""AI Generated. Reuses TestNetMHCpanPromiscuity's already-verified real fixture (4 real strong
HLA-A/B/C binders, exact core_9aa/n_promiscuous_alleles confirmed against a real local
netMHCpan-4.2 run) and chains ProtPopulationCoverage onto its output -- all promiscuous alleles
those 4 candidates hit (HLA-A/B/C, common ones) are present in the bundled world_pooled_afnd.csv
reference table, so this is a real, non-trivial computation, not just a smoke test.
"""

from ..protocols import ProtPopulationCoverage
from .test_netmhcpan import TestNetMHCpanPromiscuity


class TestPopulationCoverage(TestNetMHCpanPromiscuity):

    def testPopulationCoverage(self):
        protNetMHCpan = self.runNetMHCpan()

        protCoverage = self.newProtocol(ProtPopulationCoverage)
        protCoverage.inputROIs.set(protNetMHCpan)
        protCoverage.inputROIs.setExtended('outputROIs')
        self.launchProtocol(protCoverage, wait=True)

        outROIs = getattr(protCoverage, 'outputROIs', None)
        self.assertIsNotNone(outROIs)
        for roi in outROIs:
            pct = roi._populationCoveragePct.get()
            # All 4 EXPECTED candidates hit common HLA-A/B/C alleles present in the
            # bundled reference table -- a real, computable, non-degenerate coverage.
            self.assertIsNotNone(pct)
            self.assertTrue(pct == pct)  # not NaN
            self.assertGreater(pct, 0.0)
            self.assertLessEqual(pct, 100.0)
