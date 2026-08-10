"""Parsing of NetMHCIIpan-4.3 output and traceback to the parent ROI.

See ``ProtNetMHCIIpanPromiscuity`` for how this is used: peptides are read
from an arbitrary ``SetOfSequenceROIs``, run through NetMHCIIpan-4.3, and
promiscuous T-helper candidates are traced back to their parent ROI's
position on the original sequence.
"""

from typing import List

import numpy as np
import pandas as pd

_OUTPUT_COLUMNS = [
    "sequence", "core_9aa", "n_alleles_evaluated", "n_promiscuous_alleles", "promiscuous_alleles",
    "min_rank_el", "verdict",
]

_TRACEBACK_COLUMNS = [
    "sequence_f5", "core_9aa", "start", "end", "parent_roi_id",
    "n_promiscuous_alleles", "promiscuous_alleles", "n_alleles_evaluated", "min_rank_el",
]

VALID_CANDIDATE = "Valid candidate"
REJECTED = "Rejected"


class NetMHCIIpanParseError(Exception):
    """The NetMHCIIpan output .xls does not match the expected format."""


def _to_float_matrix(df: pd.DataFrame, xls_path: str) -> np.ndarray:
    """Coerce a block of NetMHCIIpan .xls numeric columns to float64.

    NetMHCIIpan's own scripts format these numbers with awk, which is
    locale-sensitive: on a machine where LC_NUMERIC uses ',' as the decimal
    separator (e.g. es_ES), the .xls output ends up with values like
    '14,841316' instead of '14.841316', and pandas silently leaves the whole
    column as strings instead of floats.
    """
    cleaned = df.apply(lambda col: col.astype(str).str.replace(",", ".", regex=False) if col.dtype == object else col)
    numeric = cleaned.apply(pd.to_numeric, errors="coerce")
    bad = numeric.isna() & df.notna()
    if bad.to_numpy().any():
        bad_values = sorted(set(df.to_numpy()[bad.to_numpy()].tolist()))[:5]
        raise NetMHCIIpanParseError(
            f"Could not parse numeric values in '{xls_path}': {bad_values}. Expected '.' or ',' as the "
            "decimal separator."
        )
    return numeric.to_numpy()


def parse_xls(
    xls_path: str, n_alleles: int, rank_weak: float, min_promiscuous_alleles: int,
    allele_names: List[str],
) -> pd.DataFrame:
    """Parse a NetMHCIIpan .xls output file and score the promiscuity of each peptide.

    Inverted alleles (the ``Inverted`` column) are excluded entirely before
    any computation: they neither count towards promiscuity nor can be the
    "winning" allele that determines ``core_9aa``/``min_rank_el``. This
    guarantees ``core_9aa`` is always a literal substring of the input
    peptide.

    Args:
        allele_names: Panel alleles in the SAME ORDER passed to '-a'
            (``allele_panel.split(",")``) -- ``rank_cols[i]``/``core_cols[i]``/
            ``inverted_cols[i]`` correspond to ``allele_names[i]``. Used to
            name ``promiscuous_alleles`` (needed for real population
            coverage downstream -- a raw count doesn't distinguish a
            candidate hitting 3 very rare alleles from one hitting 3
            globally common ones).

    Returns:
        DataFrame with columns ``sequence``, ``core_9aa``,
        ``n_alleles_evaluated``, ``n_promiscuous_alleles``,
        ``promiscuous_alleles`` (comma-joined allele names), ``min_rank_el``
        and ``verdict`` (``VALID_CANDIDATE`` / ``REJECTED``).
    """
    try:
        raw = pd.read_csv(xls_path, sep="\t", skiprows=2)
    except Exception as exc:
        raise NetMHCIIpanParseError(f"Could not parse NetMHCIIpan output at '{xls_path}': {exc}") from exc

    rank_cols = [c for c in raw.columns if c == "Rank_EL" or c.startswith("Rank_EL.")]
    core_cols = [c for c in raw.columns if c == "Core" or c.startswith("Core.")]
    inverted_cols = [c for c in raw.columns if c == "Inverted" or c.startswith("Inverted.")]
    if (
        len(rank_cols) != n_alleles
        or len(core_cols) != n_alleles
        or len(inverted_cols) != n_alleles
        or "Peptide" not in raw.columns
    ):
        raise NetMHCIIpanParseError(
            f"NetMHCIIpan .xls output format does not match what was expected: "
            f"found {len(rank_cols)} 'Rank_EL' column(s), {len(core_cols)} "
            f"'Core' column(s) and {len(inverted_cols)} 'Inverted' column(s) for "
            f"{n_alleles} evaluated allele(s). Columns found: {list(raw.columns)}."
        )

    rank_matrix = _to_float_matrix(raw[rank_cols], xls_path)
    core_matrix = raw[core_cols].to_numpy()
    is_inverted = _to_float_matrix(raw[inverted_cols], xls_path).astype(bool)
    row_idx = np.arange(len(raw))

    rank_matrix_normal = np.where(is_inverted, np.inf, rank_matrix)
    best_allele_idx = rank_matrix_normal.argmin(axis=1)
    best_core = core_matrix[row_idx, best_allele_idx]

    is_binder_normal = (rank_matrix_normal <= rank_weak)
    n_promiscuous_alleles = is_binder_normal.sum(axis=1)
    promiscuous_alleles = [
        ",".join(allele_names[j] for j in range(n_alleles) if is_binder_normal[i, j])
        for i in range(len(raw))
    ]

    result = pd.DataFrame(
        {
            "sequence": raw["Peptide"],
            "core_9aa": best_core,
            "n_alleles_evaluated": n_alleles,
            "n_promiscuous_alleles": n_promiscuous_alleles,
            "promiscuous_alleles": promiscuous_alleles,
            "min_rank_el": rank_matrix_normal.min(axis=1),
        }
    )
    result["verdict"] = result["n_promiscuous_alleles"].apply(
        lambda n: VALID_CANDIDATE if n >= min_promiscuous_alleles else REJECTED
    )
    return result[_OUTPUT_COLUMNS]


def build_traceback_report(report_df: pd.DataFrame, parent_records: List[dict]) -> pd.DataFrame:
    """Cross ``VALID_CANDIDATE`` peptides with their parent ROI and deduplicate windows.

    Args:
        report_df: Output of :func:`parse_xls` (includes rejected rows; only
            ``verdict == VALID_CANDIDATE`` rows are kept here).
        parent_records: One dict per evaluated ROI, with ``sequence``,
            ``start`` and ``roi_id`` keys.

    Returns:
        DataFrame with columns ``_TRACEBACK_COLUMNS``. Redundant protein-mode
        sliding windows (same ``core_9aa`` and same
        ``n_promiscuous_alleles``) are collapsed, keeping only the one with
        the lowest ``min_rank_el``.
    """
    if report_df.empty or not parent_records:
        return pd.DataFrame(columns=_TRACEBACK_COLUMNS)

    valid_df = report_df[report_df["verdict"] == VALID_CANDIDATE]
    if valid_df.empty:
        return pd.DataFrame(columns=_TRACEBACK_COLUMNS)

    records = []
    for candidate in valid_df.itertuples(index=False):
        for parent in parent_records:
            if candidate.sequence not in parent["sequence"]:
                continue
            offset = parent["sequence"].find(candidate.sequence)
            start_real = parent["start"] + offset
            end_real = start_real + len(candidate.sequence) - 1
            records.append(
                {
                    "sequence_f5": candidate.sequence,
                    "core_9aa": candidate.core_9aa,
                    "start": start_real,
                    "end": end_real,
                    "parent_roi_id": parent["roi_id"],
                    "n_promiscuous_alleles": candidate.n_promiscuous_alleles,
                    "promiscuous_alleles": candidate.promiscuous_alleles,
                    "n_alleles_evaluated": candidate.n_alleles_evaluated,
                    "min_rank_el": candidate.min_rank_el,
                }
            )

    traceback_df = pd.DataFrame.from_records(records, columns=_TRACEBACK_COLUMNS)
    return _deduplicate_protein_mode_windows(traceback_df)


def _deduplicate_protein_mode_windows(traceback_df: pd.DataFrame) -> pd.DataFrame:
    """Collapse redundant protein-mode windows: for equal (core_9aa,
    n_promiscuous_alleles), keep only the one with the lowest min_rank_el."""
    if traceback_df.empty:
        return traceback_df

    best_idx = traceback_df.groupby(["core_9aa", "n_promiscuous_alleles"], sort=False)["min_rank_el"].idxmin()
    return traceback_df.loc[best_idx].sort_index().reset_index(drop=True)
