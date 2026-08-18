"""Deterministic tests for netmhcpan.utils.netmhciipan against synthetic .xls
data (no live NetMHCIIpan-4.3 binary involved).

Same rationale as test_netmhcpan_parsing.py: test_netmhciipan.py's end-to-end
test depends on the real binary's numeric output, which is not bit-exact
across independently built/downloaded installations -- borderline fixture
windows can drop below the promiscuity floor on one install while every
other window's core_9aa still reproduces exactly. These tests instead
exercise parse_xls/build_traceback_report/the Inverted-allele exclusion
logic against fixed, hand-built input, so they reproduce bit-for-bit
identically everywhere.
"""

import tempfile
import unittest

from ..utils.netmhciipan import (
    NetMHCIIpanParseError,
    REJECTED,
    VALID_CANDIDATE,
    build_traceback_report,
    parse_xls,
)

_HEADER_2ALLELES = (
    "#dummy command line\n"
    "\t\t\tDRB1_0101\t\t\t\tDRB1_0301\t\t\t\n"
    "Pos\tPeptide\tID\tTarget\tCore\tInverted\tScore_EL\tRank_EL\t"
    "Core\tInverted\tScore_EL\tRank_EL\tAve\tNB\n"
)


def _write_xls(rows):
    """rows: list of (peptide, core1, inverted1, rank1, core2, inverted2, rank2)."""
    fh = tempfile.NamedTemporaryFile(mode="w", suffix=".xls", delete=False)
    fh.write(_HEADER_2ALLELES)
    for i, (peptide, core1, inv1, rank1, core2, inv2, rank2) in enumerate(rows):
        fh.write(f"{i}\t{peptide}\t{peptide}\t0\t{core1}\t{inv1}\t0.1\t{rank1}\t"
                  f"{core2}\t{inv2}\t0.1\t{rank2}\t0\t0\n")
    fh.close()
    return fh.name


class TestParseXls(unittest.TestCase):
    ALLELE_NAMES = ["DRB1_0101", "DRB1_0301"]

    def test_promiscuity_and_verdict(self):
        xls_path = _write_xls([
            ("AAAAAAAAAAAAAAA", "AAAAAAAAA", "0", "1.0", "AAAAAAAAA", "0", "10.0"),
            ("BBBBBBBBBBBBBBB", "BBBBBBBBB", "0", "10.0", "BBBBBBBBB", "0", "10.0"),
        ])
        df = parse_xls(xls_path, n_alleles=2, rank_weak=5.0, min_promiscuous_alleles=1,
                        allele_names=self.ALLELE_NAMES)

        aRow = df[df["sequence"] == "AAAAAAAAAAAAAAA"].iloc[0]
        self.assertEqual(aRow["n_promiscuous_alleles"], 1)
        self.assertEqual(aRow["promiscuous_alleles"], "DRB1_0101")
        self.assertEqual(aRow["verdict"], VALID_CANDIDATE)

        bRow = df[df["sequence"] == "BBBBBBBBBBBBBBB"].iloc[0]
        self.assertEqual(bRow["n_promiscuous_alleles"], 0)
        self.assertEqual(bRow["verdict"], REJECTED)

    def test_inverted_allele_excluded_from_count_and_core(self):
        # Allele 2 has the best (lowest) rank but is Inverted: it must NOT
        # count towards promiscuity, and must NOT be able to win core_9aa,
        # even though its raw rank beats allele 1's.
        xls_path = _write_xls([
            ("CCCCCCCCCCCCCCC", "NORMALAAA", "0", "3.0", "INVERTEDB", "1", "0.1"),
        ])
        df = parse_xls(xls_path, n_alleles=2, rank_weak=5.0, min_promiscuous_alleles=1,
                        allele_names=self.ALLELE_NAMES)
        row = df.iloc[0]
        self.assertEqual(row["n_promiscuous_alleles"], 1)
        self.assertEqual(row["promiscuous_alleles"], "DRB1_0101")
        self.assertEqual(row["core_9aa"], "NORMALAAA")
        self.assertAlmostEqual(row["min_rank_el"], 3.0)

    def test_all_inverted_yields_zero_promiscuity(self):
        xls_path = _write_xls([
            ("DDDDDDDDDDDDDDD", "D1", "1", "0.1", "D2", "1", "0.2"),
        ])
        df = parse_xls(xls_path, n_alleles=2, rank_weak=5.0, min_promiscuous_alleles=1,
                        allele_names=self.ALLELE_NAMES)
        row = df.iloc[0]
        self.assertEqual(row["n_promiscuous_alleles"], 0)
        self.assertEqual(row["verdict"], REJECTED)
        self.assertTrue(row["min_rank_el"] == float("inf") or row["min_rank_el"] > 5.0)

    def test_locale_comma_decimal_separator(self):
        xls_path = _write_xls([
            ("EEEEEEEEEEEEEEE", "EEEEEEEEE", "0", "0,301", "EEEEEEEEE", "0", "1,5"),
        ])
        df = parse_xls(xls_path, n_alleles=2, rank_weak=2.0, min_promiscuous_alleles=1,
                        allele_names=self.ALLELE_NAMES)
        self.assertAlmostEqual(df.iloc[0]["min_rank_el"], 0.301)

    def test_malformed_numeric_value_raises(self):
        xls_path = _write_xls([
            ("FFFFFFFFFFFFFFF", "FFFFFFFFF", "0", "not_a_number", "FFFFFFFFF", "0", "1.0"),
        ])
        with self.assertRaises(NetMHCIIpanParseError):
            parse_xls(xls_path, n_alleles=2, rank_weak=2.0, min_promiscuous_alleles=1,
                       allele_names=self.ALLELE_NAMES)


class TestBuildTracebackReport(unittest.TestCase):
    ALLELE_NAMES = ["DRB1_0101", "DRB1_0301"]

    def test_offset_recovered_from_parent_sequence(self):
        xls_path = _write_xls([
            ("KWNATLKQIASKLRE", "KWNATLKQI", "0", "0.5", "KWNATLKQI", "0", "0.5"),
        ])
        reportDf = parse_xls(xls_path, n_alleles=2, rank_weak=5.0, min_promiscuous_alleles=1,
                              allele_names=self.ALLELE_NAMES)
        parentRecords = [{"sequence": "SRAKWNATLKQIASKLREQFGNN", "start": 339, "roi_id": "roi1"}]

        tracebackDf = build_traceback_report(reportDf, parentRecords)
        self.assertEqual(len(tracebackDf), 1)
        row = tracebackDf.iloc[0]
        # 'KWNATLKQIASKLRE' starts at offset 3 (0-indexed) inside the parent.
        self.assertEqual(row["start"], 342)
        self.assertEqual(row["end"], 356)

    def test_dedup_keeps_lowest_rank_among_equal_core_and_count(self):
        xls_path = _write_xls([
            ("HHHHHHHHHHHHHHH", "HHHHHHHHH", "0", "1.0", "HHHHHHHHH", "0", "10.0"),
            ("XHHHHHHHHHHHHHHX", "HHHHHHHHH", "0", "0.4", "HHHHHHHHH", "0", "10.0"),
        ])
        reportDf = parse_xls(xls_path, n_alleles=2, rank_weak=5.0, min_promiscuous_alleles=1,
                              allele_names=self.ALLELE_NAMES)
        parentRecords = [{"sequence": "AAHHHHHHHHHHHHHHHAA", "start": 1, "roi_id": "roi1"},
                          {"sequence": "AXHHHHHHHHHHHHHHXAA", "start": 1, "roi_id": "roi2"}]

        tracebackDf = build_traceback_report(reportDf, parentRecords)
        self.assertEqual(len(tracebackDf), 1)
        self.assertAlmostEqual(tracebackDf.iloc[0]["min_rank_el"], 0.4)
        self.assertEqual(tracebackDf.iloc[0]["parent_roi_id"], "roi2")
