#!/usr/bin/env python3
"""
Subject-level QC report for Homer3 wide GLM betas.

This script treats beta values of 0 and NaN as excluded/pruned channels
per repository data-integrity policy for `homer3_glm_betas_wide.csv`.

Expected input schema:
- Required ID column: Subject
- Beta columns: S##_D##_Cond##_HbO and S##_D##_Cond##_HbR

Outputs:
- Subject-level QC CSV (one row per subject)
- Cohort summary CSV (single row with aggregate diagnostics)
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple


BETA_COL_RE = re.compile(r"^S(\d+)_D(\d+)_Cond(\d+)_(HbO|HbR)$")
EXPECTED_CONDITIONS = ("01", "02", "03", "04")
EXPECTED_CHROMS = ("HbO", "HbR")
MISSING_TOKENS = {"", "na", "nan", "null"}


@dataclass(frozen=True)
class Schema:
    """Validated schema metadata for the Homer beta-wide table."""

    subject_column: str
    channels: Tuple[str, ...]
    conditions: Tuple[str, ...]
    # Maps (channel, condition, chromophore) -> original CSV column name.
    col_map: Mapping[Tuple[str, str, str], str]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate subject-level channel exclusion QC from "
            "data/tabular/homer3_glm_betas_wide.csv."
        )
    )
    parser.add_argument(
        "--input-csv",
        default="data/tabular/homer3_glm_betas_wide.csv",
        help="Path to Homer3 wide beta CSV.",
    )
    parser.add_argument(
        "--output-csv",
        default="data/results/homer3_betas_qc_subject_level.csv",
        help="Path to output subject-level QC CSV.",
    )
    parser.add_argument(
        "--summary-csv",
        default="data/results/homer3_betas_qc_cohort_summary.csv",
        help="Path to output cohort summary CSV.",
    )
    parser.add_argument(
        "--bad-channel-min-excluded-conds",
        type=int,
        default=2,
        help=(
            "Primary bad-channel rule: mark a channel bad if excluded in at least "
            "this many conditions (out of 4). Default=2."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow overwriting output files if they already exist.",
    )
    return parser.parse_args(argv)


def ensure_writable_output(path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Refusing to overwrite existing output: {path}. Pass --overwrite to allow."
        )
    path.parent.mkdir(parents=True, exist_ok=True)


def normalize_numeric_token(
    raw: str | None,
    *,
    column: str,
    subject: str,
) -> float:
    """
    Parse a numeric token into float with explicit failure on invalid strings.

    Allowed missing tokens are mapped to NaN. Non-numeric strings outside
    this set are rejected to avoid silent corruption of QC results.
    """
    token = "" if raw is None else str(raw).strip()
    if token.lower() in MISSING_TOKENS:
        return math.nan
    try:
        return float(token)
    except ValueError as exc:
        raise ValueError(
            "Non-numeric beta token encountered. "
            f"subject={subject}, column={column}, token={token!r}"
        ) from exc


def is_excluded_value(value: float) -> bool:
    """A beta is excluded if it is NaN/non-finite or exactly 0."""
    return (not math.isfinite(value)) or value == 0.0


def parse_schema(fieldnames: Sequence[str]) -> Schema:
    if not fieldnames:
        raise ValueError("Input CSV appears empty (missing header row).")
    if "Subject" not in fieldnames:
        raise ValueError("Expected required ID column 'Subject' in input CSV.")
    if fieldnames[0] != "Subject":
        # Fail-fast: the external export contract in this repo places Subject first.
        raise ValueError(
            "Expected first column to be 'Subject' for deterministic parsing, "
            f"but got {fieldnames[0]!r}."
        )

    col_map: Dict[Tuple[str, str, str], str] = {}
    channel_set = set()
    channel_sort_key: Dict[str, Tuple[int, int]] = {}
    condition_set = set()
    malformed: List[str] = []

    for col in fieldnames[1:]:
        match = BETA_COL_RE.match(col)
        if match is None:
            malformed.append(col)
            continue
        s_idx, d_idx, cond, chrom = match.groups()
        channel = f"S{int(s_idx):02d}_D{int(d_idx):02d}"
        key = (channel, cond, chrom)
        if key in col_map:
            raise ValueError(f"Duplicate beta column detected: {col!r}")
        col_map[key] = col
        channel_set.add(channel)
        channel_sort_key[channel] = (int(s_idx), int(d_idx))
        condition_set.add(cond)

    if malformed:
        sample = ", ".join(malformed[:10])
        raise ValueError(
            "Found non-matching beta column names. "
            "Expected format like S01_D01_Cond01_HbO. "
            f"Examples: {sample}"
        )

    if not col_map:
        raise ValueError("No beta columns found in input CSV.")

    observed_conditions = tuple(sorted(condition_set))
    if observed_conditions != EXPECTED_CONDITIONS:
        raise ValueError(
            "Unexpected condition set in beta columns. "
            f"Expected {EXPECTED_CONDITIONS}, observed {observed_conditions}."
        )

    channels = tuple(sorted(channel_set, key=lambda x: channel_sort_key[x]))
    for channel in channels:
        for cond in EXPECTED_CONDITIONS:
            for chrom in EXPECTED_CHROMS:
                key = (channel, cond, chrom)
                if key not in col_map:
                    raise ValueError(
                        "Missing required channel-condition-chromophore column for "
                        f"channel={channel}, condition=Cond{cond}, chrom={chrom}."
                    )

    return Schema(
        subject_column="Subject",
        channels=channels,
        conditions=EXPECTED_CONDITIONS,
        col_map=col_map,
    )


def channel_available_for_condition(
    row: Mapping[str, str],
    schema: Schema,
    *,
    subject: str,
    channel: str,
    condition: str,
) -> bool:
    """
    Availability rule: channel is available only when BOTH HbO and HbR are present.
    """
    hbo_col = schema.col_map[(channel, condition, "HbO")]
    hbr_col = schema.col_map[(channel, condition, "HbR")]
    hbo = normalize_numeric_token(row.get(hbo_col), column=hbo_col, subject=subject)
    hbr = normalize_numeric_token(row.get(hbr_col), column=hbr_col, subject=subject)
    return (not is_excluded_value(hbo)) and (not is_excluded_value(hbr))


def percent(numerator: int, denominator: int) -> float:
    if denominator <= 0:
        return math.nan
    return (100.0 * numerator) / float(denominator)


def compute_subject_metrics(
    row: Mapping[str, str],
    schema: Schema,
    *,
    bad_channel_min_excluded_conds: int,
) -> Dict[str, object]:
    subject = str(row[schema.subject_column]).strip()
    if subject == "":
        raise ValueError("Encountered empty Subject ID.")

    total_channels = len(schema.channels)
    if total_channels <= 0:
        raise ValueError("Schema contains zero channels; cannot compute QC.")

    per_condition_available_counts: Dict[str, int] = {}
    excluded_cond_count_by_channel: Dict[str, int] = {ch: 0 for ch in schema.channels}

    for cond in schema.conditions:
        available_count = 0
        for channel in schema.channels:
            available = channel_available_for_condition(
                row,
                schema,
                subject=subject,
                channel=channel,
                condition=cond,
            )
            if available:
                available_count += 1
            else:
                excluded_cond_count_by_channel[channel] += 1
        per_condition_available_counts[cond] = available_count

    bad_any = sorted(
        [ch for ch, n_exc in excluded_cond_count_by_channel.items() if n_exc >= 1]
    )
    bad_ge2 = sorted(
        [
            ch
            for ch, n_exc in excluded_cond_count_by_channel.items()
            if n_exc >= bad_channel_min_excluded_conds
        ]
    )
    bad_all = sorted(
        [ch for ch, n_exc in excluded_cond_count_by_channel.items() if n_exc == 4]
    )

    cond_passes = 0
    for cond in schema.conditions:
        p = per_condition_available_counts[cond] / float(total_channels)
        # Per requested rule: strict > 50% availability.
        # This thresholding style is aligned with fNIRS exclusion precedents using
        # "more than half of channels/trials usable" criteria (Novi 2023; Pinti 2024;
        # Dina 2025; see CITATIONS.md).
        if p > 0.5:
            cond_passes += 1

    out: Dict[str, object] = {
        "Subject": subject,
        "n_total_channels": total_channels,
        "n_bad_channels_any": len(bad_any),
        "pct_bad_channels_any": percent(len(bad_any), total_channels),
        "bad_channels_any": ";".join(bad_any),
        "n_bad_channels_ge2": len(bad_ge2),
        "pct_bad_channels_ge2": percent(len(bad_ge2), total_channels),
        "bad_channels_ge2": ";".join(bad_ge2),
        "n_bad_channels_all": len(bad_all),
        "pct_bad_channels_all": percent(len(bad_all), total_channels),
        "bad_channels_all": ";".join(bad_all),
    }
    for cond in schema.conditions:
        cond_label = f"cond{cond}"
        out[f"{cond_label}_available_channels"] = per_condition_available_counts[cond]
        out[f"{cond_label}_available_pct"] = percent(
            per_condition_available_counts[cond], total_channels
        )

    out["n_tasks_over_50pct_available"] = cond_passes
    out["pct_tasks_over_50pct_available"] = percent(cond_passes, len(schema.conditions))
    return out


def median(values: Sequence[float]) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    n = len(ordered)
    mid = n // 2
    if n % 2 == 1:
        return ordered[mid]
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def mean(values: Sequence[float]) -> float:
    if not values:
        return math.nan
    return sum(values) / float(len(values))


def build_cohort_summary(subject_rows: Sequence[Mapping[str, object]]) -> Dict[str, object]:
    if not subject_rows:
        raise ValueError("No subject rows were computed; cannot produce cohort summary.")

    ge2_counts = [int(r["n_bad_channels_ge2"]) for r in subject_rows]
    ge2_pct = [float(r["pct_bad_channels_ge2"]) for r in subject_rows]
    task_pass_pct = [float(r["pct_tasks_over_50pct_available"]) for r in subject_rows]
    n_tasks = [int(r["n_tasks_over_50pct_available"]) for r in subject_rows]

    return {
        "n_subjects": len(subject_rows),
        "mean_n_bad_channels_ge2": mean(ge2_counts),
        "median_n_bad_channels_ge2": median(ge2_counts),
        "min_n_bad_channels_ge2": min(ge2_counts),
        "max_n_bad_channels_ge2": max(ge2_counts),
        "mean_pct_bad_channels_ge2": mean(ge2_pct),
        "median_pct_bad_channels_ge2": median(ge2_pct),
        "min_pct_bad_channels_ge2": min(ge2_pct),
        "max_pct_bad_channels_ge2": max(ge2_pct),
        "mean_pct_tasks_over_50pct_available": mean(task_pass_pct),
        "median_pct_tasks_over_50pct_available": median(task_pass_pct),
        "min_pct_tasks_over_50pct_available": min(task_pass_pct),
        "max_pct_tasks_over_50pct_available": max(task_pass_pct),
        "n_subjects_with_any_bad_channels_ge2": sum(1 for x in ge2_counts if x > 0),
        "n_subjects_with_0_tasks_over_50pct": sum(1 for x in n_tasks if x == 0),
        "n_subjects_with_4_tasks_over_50pct": sum(1 for x in n_tasks if x == 4),
    }


def read_csv_rows(input_csv: Path) -> Tuple[Schema, List[Dict[str, str]]]:
    with input_csv.open("r", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"Missing CSV header row: {input_csv}")
        schema = parse_schema(reader.fieldnames)
        rows: List[Dict[str, str]] = []
        seen_subjects = set()
        for idx, row in enumerate(reader, start=2):
            subject = str(row.get(schema.subject_column, "")).strip()
            if subject == "":
                raise ValueError(f"Row {idx} has empty Subject value.")
            if subject in seen_subjects:
                raise ValueError(
                    "Duplicate Subject IDs are not allowed for subject-level QC. "
                    f"Found duplicate Subject={subject!r}."
                )
            seen_subjects.add(subject)
            rows.append(row)
    if not rows:
        raise ValueError(f"Input CSV has no data rows: {input_csv}")
    return schema, rows


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        raise ValueError(f"No rows to write for output: {path}")
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def print_summary(summary: Mapping[str, object], subject_output: Path, cohort_output: Path) -> None:
    print(f"[homer_betas_qc] Wrote subject-level QC: {subject_output}")
    print(f"[homer_betas_qc] Wrote cohort summary:  {cohort_output}")
    print("[homer_betas_qc] Cohort summary:")
    for key in (
        "n_subjects",
        "mean_n_bad_channels_ge2",
        "median_n_bad_channels_ge2",
        "min_n_bad_channels_ge2",
        "max_n_bad_channels_ge2",
        "mean_pct_bad_channels_ge2",
        "mean_pct_tasks_over_50pct_available",
        "n_subjects_with_any_bad_channels_ge2",
        "n_subjects_with_0_tasks_over_50pct",
        "n_subjects_with_4_tasks_over_50pct",
    ):
        print(f"  - {key}: {summary[key]}")


def run(
    input_csv: Path,
    output_csv: Path,
    summary_csv: Path,
    *,
    bad_channel_min_excluded_conds: int,
    overwrite: bool,
) -> int:
    if bad_channel_min_excluded_conds < 1 or bad_channel_min_excluded_conds > 4:
        raise ValueError(
            "--bad-channel-min-excluded-conds must be in [1, 4]. "
            f"Got {bad_channel_min_excluded_conds}."
        )

    ensure_writable_output(output_csv, overwrite)
    ensure_writable_output(summary_csv, overwrite)

    schema, rows = read_csv_rows(input_csv)

    per_subject = [
        compute_subject_metrics(
            row,
            schema,
            bad_channel_min_excluded_conds=bad_channel_min_excluded_conds,
        )
        for row in rows
    ]
    cohort_summary = build_cohort_summary(per_subject)

    write_csv(output_csv, per_subject)
    write_csv(summary_csv, [cohort_summary])
    print_summary(cohort_summary, output_csv, summary_csv)
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    return run(
        Path(args.input_csv),
        Path(args.output_csv),
        Path(args.summary_csv),
        bad_channel_min_excluded_conds=args.bad_channel_min_excluded_conds,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    raise SystemExit(main())
