"""Mask between-subject AUC outliers within each Homer beta column.

This stage operates on the derived wide AUC table after FIR-to-AUC conversion
and before tabular merge. For each exact channel x condition x chromophore
column, it computes the mean and sample standard deviation across subjects and
replaces values outside mean +/- 3 SD with NaN.

Scientific integrity note:
- This repo implements mean +/- 3 SD because it was explicitly requested for
  this project.
- The script also documents two limitations from the statistical literature:
  (1) mean/SD outlier rules are non-robust because the candidate outlier
      influences both the center and spread estimates (Leys et al., 2013),
  (2) with a sample z-score rule, thresholds this extreme are impossible to
      exceed for small n; for 3 SD, fewer than 11 observed subjects cannot
      yield a flagged value (Shiffler, 1988).
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

from homer_fir import compute_file_sha256


INPUT_CSV = "data/tabular/generated_data/homer3_glm_betas_wide_auc.csv"
OUTPUT_CSV = "data/tabular/generated_data/homer3_glm_betas_wide_auc_outliers_masked.csv"
OUT_AUDIT_CSV = "data/results/homer_auc_outlier_audit.csv"
OUT_SUMMARY_JSON = "data/results/homer_auc_outlier_summary.json"
BETA_COL_PATTERN = re.compile(r"^S\d+_D\d+_Cond\d{2}_(HbO|HbR)$")
ALLOWED_MISSING = {"", "NA", "NaN", "nan", "NAN", "NULL", "null"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Mask between-subject AUC outliers independently within each Homer "
            "channel x condition x chromophore beta column."
        )
    )
    parser.add_argument(
        "--input-csv",
        default=INPUT_CSV,
        help="Path to the raw derived AUC CSV.",
    )
    parser.add_argument(
        "--output-csv",
        default=OUTPUT_CSV,
        help="Path to the outlier-masked AUC CSV used downstream.",
    )
    parser.add_argument(
        "--out-audit-csv",
        default=OUT_AUDIT_CSV,
        help="Path to the row-level CSV audit of censored cells.",
    )
    parser.add_argument(
        "--out-summary-json",
        default=OUT_SUMMARY_JSON,
        help="Path to the JSON summary of between-subject outlier masking.",
    )
    parser.add_argument(
        "--sd-threshold",
        type=float,
        default=3.0,
        help="Two-sided standard-deviation threshold for masking outliers.",
    )
    parser.add_argument(
        "--min-n",
        type=int,
        default=11,
        help="Minimum observed subject count required to screen a beta column.",
    )
    return parser.parse_args()


def _minimum_detectable_n(sd_threshold: float) -> int:
    """Return the smallest n where a sample z-score can exceed `sd_threshold`.

    Shiffler (1988) showed that the largest attainable z-score in a sample of
    size n is bounded by (n - 1) / sqrt(n). This means mean +/- 3 SD screening
    cannot flag any point when n < 11, regardless of the observed values.
    """
    if not math.isfinite(sd_threshold) or sd_threshold <= 0:
        raise ValueError("--sd-threshold must be a finite value greater than 0.")

    n = 2
    while (n - 1) / math.sqrt(n) <= sd_threshold:
        n += 1
    return n


def _load_auc_table(path: Path) -> tuple[pd.DataFrame, list[str]]:
    if not path.exists():
        raise FileNotFoundError(f"Input CSV not found: {path}")

    df = pd.read_csv(path, dtype="string", keep_default_na=False)
    if df.empty:
        raise ValueError("Input AUC CSV is empty.")
    if "Subject" not in df.columns:
        raise ValueError("Expected column 'Subject' in AUC CSV.")
    if df["Subject"].str.strip().eq("").any():
        raise ValueError("Encountered blank Subject value in AUC CSV.")
    if df["Subject"].duplicated().any():
        dupes = df.loc[df["Subject"].duplicated(), "Subject"].unique().tolist()[:10]
        raise ValueError(
            "Duplicate subject rows detected in AUC CSV; expected one row per subject. "
            f"Examples: {dupes}"
        )

    beta_cols = [col for col in df.columns if BETA_COL_PATTERN.match(col)]
    if not beta_cols:
        raise ValueError("No Homer beta columns matched pattern like 'S01_D01_Cond01_HbO'.")

    for col in beta_cols:
        raw = df[col].astype("string").str.strip()
        numeric = pd.to_numeric(raw.replace(list(ALLOWED_MISSING), pd.NA), errors="coerce")
        bad_mask = raw.ne("") & ~raw.isin(ALLOWED_MISSING) & numeric.isna()
        if bad_mask.any():
            examples = raw[bad_mask].dropna().unique().tolist()[:10]
            raise ValueError(
                f"Non-numeric values found in Homer beta column '{col}'. Examples: {examples}"
            )
        df[col] = numeric.astype("Float64")

    return df, beta_cols


def mask_between_subject_outliers(
    input_csv: str = INPUT_CSV,
    output_csv: str = OUTPUT_CSV,
    out_audit_csv: str = OUT_AUDIT_CSV,
    out_summary_json: str = OUT_SUMMARY_JSON,
    *,
    sd_threshold: float = 3.0,
    min_n: int = 11,
) -> tuple[Path, Path, Path]:
    """Mask between-subject outliers within each Homer AUC beta column."""
    required_n = _minimum_detectable_n(sd_threshold)
    if min_n < required_n:
        raise ValueError(
            f"--min-n must be >= {required_n} for a {sd_threshold:g}-SD rule; "
            "smaller samples cannot yield a detectable outlier under sample mean/SD screening."
        )

    input_path = Path(input_csv)
    output_path = Path(output_csv)
    audit_path = Path(out_audit_csv)
    summary_path = Path(out_summary_json)
    df, beta_cols = _load_auc_table(input_path)

    output_df = df.copy(deep=True)
    audit_rows: list[dict[str, object]] = []
    per_column_summary: list[dict[str, object]] = []
    total_observed = 0

    for col in beta_cols:
        values = output_df[col].astype("Float64")
        observed_mask = values.notna()
        observed = values[observed_mask].astype(float)
        observed_n = int(observed.shape[0])
        total_observed += observed_n

        summary_row: dict[str, object] = {
            "column": col,
            "observed_subjects": observed_n,
            "masked_subjects": 0,
            "skip_reason": None,
        }
        if observed_n < min_n:
            summary_row["skip_reason"] = (
                f"observed_n<{min_n}; 3-SD screening is not attempted for undersized columns"
            )
            per_column_summary.append(summary_row)
            continue

        mean_value = float(observed.mean())
        sd_value = float(observed.std(ddof=1))
        summary_row["mean_across_subjects"] = mean_value
        summary_row["sd_across_subjects"] = sd_value

        if not math.isfinite(sd_value) or np.isclose(sd_value, 0.0):
            summary_row["skip_reason"] = "sd_non_finite_or_zero"
            per_column_summary.append(summary_row)
            continue

        lower_bound = mean_value - (sd_threshold * sd_value)
        upper_bound = mean_value + (sd_threshold * sd_value)
        summary_row["lower_bound"] = lower_bound
        summary_row["upper_bound"] = upper_bound

        flagged_mask = observed_mask & ((values.astype(float) < lower_bound) | (values.astype(float) > upper_bound))
        flagged_idx = output_df.index[flagged_mask]
        summary_row["masked_subjects"] = int(len(flagged_idx))

        for idx in flagged_idx:
            audit_rows.append(
                {
                    "Subject": str(output_df.at[idx, "Subject"]),
                    "channel_column": col,
                    "original_value": float(values.at[idx]),
                    "mean_across_subjects": mean_value,
                    "sd_across_subjects": sd_value,
                    "lower_bound": lower_bound,
                    "upper_bound": upper_bound,
                    "sd_threshold": sd_threshold,
                }
            )
        if len(flagged_idx) > 0:
            output_df.loc[flagged_idx, col] = pd.NA

        per_column_summary.append(summary_row)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_df.to_csv(output_path, index=False, na_rep="NaN")

    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_df = pd.DataFrame(
        audit_rows,
        columns=[
            "Subject",
            "channel_column",
            "original_value",
            "mean_across_subjects",
            "sd_across_subjects",
            "lower_bound",
            "upper_bound",
            "sd_threshold",
        ],
    )
    audit_df.to_csv(audit_path, index=False)

    summary_payload = {
        "analysis_step": "mask_homer_auc_between_subject_outliers",
        "rule": {
            "type": "between_subject_mean_plus_minus_sd",
            "sd_threshold": sd_threshold,
            "min_n": min_n,
            "minimum_detectable_n_for_threshold": required_n,
        },
        "input_csv": str(input_path),
        "input_csv_sha256": compute_file_sha256(input_path),
        "output_csv": str(output_path),
        "output_csv_sha256": compute_file_sha256(output_path),
        "audit_csv": str(audit_path),
        "audit_csv_sha256": compute_file_sha256(audit_path),
        "subject_rows": int(output_df.shape[0]),
        "beta_columns": len(beta_cols),
        "observed_beta_cells": total_observed,
        "masked_beta_cells": int(len(audit_rows)),
        "per_column_summary": per_column_summary,
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(
        json.dumps(summary_payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    return output_path, audit_path, summary_path


def main() -> None:
    args = parse_args()
    output_path, audit_path, summary_path = mask_between_subject_outliers(
        input_csv=args.input_csv,
        output_csv=args.output_csv,
        out_audit_csv=args.out_audit_csv,
        out_summary_json=args.out_summary_json,
        sd_threshold=args.sd_threshold,
        min_n=args.min_n,
    )
    print(f"[INFO] Wrote outlier-masked AUC table: {output_path}")
    print(f"[INFO] Wrote outlier audit CSV: {audit_path}")
    print(f"[INFO] Wrote outlier summary JSON: {summary_path}")


if __name__ == "__main__":
    main()
