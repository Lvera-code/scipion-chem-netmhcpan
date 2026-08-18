"""Deterministic tests for netmhcpan.utils.netmhcpan against synthetic .xls
data (no live NetMHCpan-4.2 binary involved).

test_netmhcpan.py's end-to-end test necessarily depends on the real
NetMHCpan-4.2 binary's numeric output, which is not bit-exact across
independently built/downloaded installations: core_9aa and window positions
reproduce exactly, but n_promiscuous_alleles can differ for alleles whose
%Rank_EL sits close to the classification cutoff. That means the one thing
an end-to-end test cannot fully guarantee is that this exact code is free of
bugs on any machine. These tests close that gap: they exercise parse_xls,
build_traceback_report and the protein-mode dedup logic against fixed,
hand-built input, so they reproduce bit-for-bit identically everywhere,
independent of any NetMHCpan build or version.
"""

import tempfile
import unittest

from ..utils.netmhcpan import (
    NetMHCpanParseError,
    REJECTED,
    VALID_CANDIDATE,
    build_traceback_report,
    parse_xls,
)

_HEADER_2ALLELES = (
    "#dummy command line\n"
    "\tHLA-A01:01\t\t\t\tHLA-A02:01\t\t\t\n"
    "Pos\tPeptide\tID\tcore\ticore\tEL_score\tEL_rank\tcore\ticore\tEL_score\tEL_rank\tAve\tNB\n"
)


def _write_xls(rows):
    """rows: list of (peptide, core1, rank1, core2, rank2) tuples."""
    fh = tempfile.NamedTemporaryFile(mode="w", suffix=".xls", delete=False)
    fh.write(_HEADER_2ALLELES)
    for i, (peptide, core1, rank1, core2, rank2) in enumerate(rows):
        fh.write(f"{i}\t{peptide}\t{peptide}\t{core1}\t{peptide}\t0.1\t{rank1}\t"
                  f"{core2}\t{peptide}\t0.1\t{rank2}\t0\t0\n")
    fh.close()
    return fh.name


class TestParseXls(unittest.TestCase):
    ALLELE_NAMES = ["HLA-A01:01", "HLA-A02:01"]

    def test_promiscuity_and_verdict(self):
        xls_path = _write_xls([
            ("AAAAAAAAA", "AAAAAAAAA", "0.5", "AAAAAAAAA", "10.0"),
            ("BBBBBBBBB", "BBBBBBBBB", "10.0", "BBBBBBBBB", "10.0"),
        ])
        df = parse_xls(xls_path, n_alleles=2, rank_weak=2.0, min_promiscuous_alleles=1,
                        allele_names=self.ALLELE_NAMES)

        aRow = df[df["sequence"] == "AAAAAAAAA"].iloc[0]
        self.assertEqual(aRow["n_promiscuous_alleles"], 1)
        self.assertEqual(aRow["promiscuous_alleles"], "HLA-A01:01")
        self.assertEqual(aRow["core_9aa"], "AAAAAAAAA")
        self.assertAlmostEqual(aRow["min_rank_el"], 0.5)
        self.assertEqual(aRow["verdict"], VALID_CANDIDATE)

        bRow = df[df["sequence"] == "BBBBBBBBB"].iloc[0]
        self.assertEqual(bRow["n_promiscuous_alleles"], 0)
        self.assertEqual(bRow["promiscuous_alleles"], "")
        self.assertEqual(bRow["verdict"], REJECTED)

    def test_best_core_is_lowest_rank_allele(self):
        # Allele 2 has the lower (better) rank, so its core must win, not allele 1's.
        xls_path = _write_xls([
            ("CCCCCCCCC", "CORE1AAAA", "5.0", "CORE2BBBB", "0.2"),
        ])
        df = parse_xls(xls_path, n_alleles=2, rank_weak=6.0, min_promiscuous_alleles=1,
                        allele_names=self.ALLELE_NAMES)
        row = df.iloc[0]
        self.assertEqual(row["core_9aa"], "CORE2BBBB")
        self.assertAlmostEqual(row["min_rank_el"], 0.2)

    def test_locale_comma_decimal_separator(self):
        # Real bug: on an es_ES-locale machine, NetMHCpan's own awk-based .xls
        # formatting emits '0,301' instead of '0.301'.
        xls_path = _write_xls([
            ("DDDDDDDDD", "DDDDDDDDD", "0,301", "DDDDDDDDD", "1,5"),
        ])
        df = parse_xls(xls_path, n_alleles=2, rank_weak=2.0, min_promiscuous_alleles=1,
                        allele_names=self.ALLELE_NAMES)
        self.assertAlmostEqual(df.iloc[0]["min_rank_el"], 0.301)

    def test_malformed_numeric_value_raises(self):
        xls_path = _write_xls([
            ("EEEEEEEEE", "EEEEEEEEE", "not_a_number", "EEEEEEEEE", "1.0"),
        ])
        with self.assertRaises(NetMHCpanParseError):
            parse_xls(xls_path, n_alleles=2, rank_weak=2.0, min_promiscuous_alleles=1,
                       allele_names=self.ALLELE_NAMES)

    def test_wrong_allele_count_raises(self):
        xls_path = _write_xls([("FFFFFFFFF", "FFFFFFFFF", "1.0", "FFFFFFFFF", "1.0")])
        with self.assertRaises(NetMHCpanParseError):
            parse_xls(xls_path, n_alleles=3, rank_weak=2.0, min_promiscuous_alleles=1,
                       allele_names=["A", "B", "C"])


class TestBuildTracebackReport(unittest.TestCase):
    ALLELE_NAMES = ["HLA-A01:01", "HLA-A02:01"]

    def test_offset_recovered_from_parent_sequence(self):
        xls_path = _write_xls([("VKEKYQHLW", "VKEKYQHLW", "0.5", "VKEKYQHLW", "0.5")])
        reportDf = parse_xls(xls_path, n_alleles=2, rank_weak=2.0, min_promiscuous_alleles=1,
                              allele_names=self.ALLELE_NAMES)
        parentRecords = [{"sequence": "MRVKEKYQHLWRWGW", "start": 1, "roi_id": "roi1"}]

        tracebackDf = build_traceback_report(reportDf, parentRecords)
        self.assertEqual(len(tracebackDf), 1)
        row = tracebackDf.iloc[0]
        # 'VKEKYQHLW' starts at offset 2 (0-indexed) inside 'MRVKEKYQHLWRWGW'.
        self.assertEqual(row["start"], 3)
        self.assertEqual(row["end"], 11)
        self.assertEqual(row["parent_roi_id"], "roi1")

    def test_rejected_rows_excluded(self):
        xls_path = _write_xls([("GGGGGGGGG", "GGGGGGGGG", "50.0", "GGGGGGGGG", "50.0")])
        reportDf = parse_xls(xls_path, n_alleles=2, rank_weak=2.0, min_promiscuous_alleles=1,
                              allele_names=self.ALLELE_NAMES)
        parentRecords = [{"sequence": "AAAGGGGGGGGGAAA", "start": 1, "roi_id": "roi1"}]

        tracebackDf = build_traceback_report(reportDf, parentRecords)
        self.assertTrue(tracebackDf.empty)

    def test_dedup_keeps_lowest_rank_among_equal_core_and_count(self):
        # Two protein-mode windows that both resolve to the same (core, count)
        # -- only the one with the lowest min_rank_el should survive.
        xls_path = _write_xls([
            ("HHHHHHHHH", "HHHHHHHHH", "1.0", "HHHHHHHHH", "10.0"),
            ("IHHHHHHHHJ", "HHHHHHHHH", "0.3", "HHHHHHHHH", "10.0"),
        ])
        reportDf = parse_xls(xls_path, n_alleles=2, rank_weak=2.0, min_promiscuous_alleles=1,
                              allele_names=self.ALLELE_NAMES)
        parentRecords = [{"sequence": "XXHHHHHHHHHXX", "start": 1, "roi_id": "roi1"},
                          {"sequence": "XIHHHHHHHHJXX", "start": 1, "roi_id": "roi2"}]

        tracebackDf = build_traceback_report(reportDf, parentRecords)
        self.assertEqual(len(tracebackDf), 1)
        self.assertAlmostEqual(tracebackDf.iloc[0]["min_rank_el"], 0.3)
        self.assertEqual(tracebackDf.iloc[0]["parent_roi_id"], "roi2")

    def test_empty_report_returns_empty_frame(self):
        import pandas as pd
        self.assertTrue(build_traceback_report(pd.DataFrame(), []).empty)
