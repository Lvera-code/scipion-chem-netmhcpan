"""Parsing of NetMHCpan-4.2 output and traceback to the parent ROI.

See ``ProtNetMHCpanPromiscuity`` for how this is used: peptides are read
from an arbitrary ``SetOfSequenceROIs``, run through NetMHCpan-4.2, and
promiscuous cytotoxic T-cell candidates are traced back to their parent
ROI's position on the original sequence.

Unlike NetMHCIIpan, NetMHCpan-4.2's .xls output has NO ``Inverted`` column:
the MHC-I binding core does not suffer the "inverted register" alignment
artifact of NetMHCIIpan's Gibbs-sampling training (a closed groove at both
ends leaves no room to invert the register), so no equivalent filter is
needed here -- confirmed by reading the real .xls columns, not assumed by
analogy with NetMHCIIpan. Column names also differ: ``EL_rank``/``core``
(lowercase), not ``Rank_EL``/``Core``.
"""

from typing import List

import numpy as np
import pandas as pd

_OUTPUT_COLUMNS = [
    "sequence", "core_9aa", "n_alleles_evaluated", "n_promiscuous_alleles", "min_rank_el", "verdict",
]

_TRACEBACK_COLUMNS = [
    "sequence_f5", "core_9aa", "start", "end", "parent_roi_id",
    "n_promiscuous_alleles", "n_alleles_evaluated", "min_rank_el",
]

VALID_CANDIDATE = "Valid candidate"
REJECTED = "Rejected"


class NetMHCpanParseError(Exception):
    """The NetMHCpan output .xls does not match the expected format."""


def _to_float_matrix(df: pd.DataFrame, xls_path: str) -> np.ndarray:
    """Coerce a block of NetMHCpan .xls numeric columns to float64.

    Same real bug as ``netmhciipan._to_float_matrix`` (NetMHCpan's own
    scripts format these numbers with awk, which is locale-sensitive): on a
    machine where LC_NUMERIC uses ',' as the decimal separator (e.g. es_ES),
    the .xls output ends up with values like '0,301' instead of '0.301',
    and pandas silently leaves the whole column as strings instead of
    floats -- confirmed real 2026-07-31 (Blanca hit this running the real
    test on her own machine: ``TypeError: '<=' not supported between
    instances of 'str' and 'float'`` in the ``is_binder`` comparison below).
    """
    cleaned = df.apply(lambda col: col.astype(str).str.replace(",", ".", regex=False) if col.dtype == object else col)
    numeric = cleaned.apply(pd.to_numeric, errors="coerce")
    bad = numeric.isna() & df.notna()
    if bad.to_numpy().any():
        bad_values = sorted(set(df.to_numpy()[bad.to_numpy()].tolist()))[:5]
        raise NetMHCpanParseError(
            f"Could not parse numeric values in '{xls_path}': {bad_values}. Expected '.' or ',' as the "
            "decimal separator."
        )
    return numeric.to_numpy()


def parse_xls(xls_path: str, n_alleles: int, rank_weak: float, min_promiscuous_alleles: int) -> pd.DataFrame:
    """Parse a NetMHCpan-4.2 .xls output file and score the promiscuity of each peptide.

    No ``Inverted`` column exists in this format (see module docstring), so
    every allele's ``EL_rank``/``core`` is used directly, unlike
    ``netmhciipan.parse_xls`` which must exclude inverted alleles first.

    Returns:
        DataFrame with columns ``sequence``, ``core_9aa``,
        ``n_alleles_evaluated``, ``n_promiscuous_alleles``, ``min_rank_el``
        and ``verdict`` (``VALID_CANDIDATE`` / ``REJECTED``).
    """
    try:
        raw = pd.read_csv(xls_path, sep="\t", skiprows=2)
    except Exception as exc:
        raise NetMHCpanParseError(f"Could not parse NetMHCpan output at '{xls_path}': {exc}") from exc

    rank_cols = [c for c in raw.columns if c == "EL_rank" or c.startswith("EL_rank.")]
    core_cols = [c for c in raw.columns if c == "core" or c.startswith("core.")]
    if len(rank_cols) != n_alleles or len(core_cols) != n_alleles or "Peptide" not in raw.columns:
        raise NetMHCpanParseError(
            f"NetMHCpan .xls output format does not match what was expected: "
            f"found {len(rank_cols)} 'EL_rank' column(s) and {len(core_cols)} "
            f"'core' column(s) for {n_alleles} evaluated allele(s). "
            f"Columns found: {list(raw.columns)}."
        )

    rank_matrix = _to_float_matrix(raw[rank_cols], xls_path)
    core_matrix = raw[core_cols].to_numpy()
    row_idx = np.arange(len(raw))

    best_allele_idx = rank_matrix.argmin(axis=1)
    best_core = core_matrix[row_idx, best_allele_idx]

    is_binder = rank_matrix <= rank_weak
    n_promiscuous_alleles = is_binder.sum(axis=1)

    result = pd.DataFrame(
        {
            "sequence": raw["Peptide"],
            "core_9aa": best_core,
            "n_alleles_evaluated": n_alleles,
            "n_promiscuous_alleles": n_promiscuous_alleles,
            "min_rank_el": rank_matrix.min(axis=1),
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
