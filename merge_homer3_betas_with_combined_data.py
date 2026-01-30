"""
Merge Homer3 GLM betas (wide) with tabular combined dataset (demographics/behavior).

This script performs an INNER JOIN between:
  - `data/tabular/homer3_glm_betas_wide.csv` (ID column: `Subject`, e.g. `sub_0001`)
  - `data/tabular/combined_sfv_data.csv` (ID column: `subject_id`, sometimes zero-padded)

Key behavior
  - Extracts the numeric subject id from both ID columns (handles padding like `0017` vs `17`,
    and Homer-style IDs like `sub_0017`).
  - Produces a single merged row per subject containing all columns from both inputs.
  - Does NOT impute missingness. This script preserves beta columns as-is; downstream modeling
    must explicitly treat 0/NaN as pruned/missing channels per repo policy.

Scientific integrity note
  - `homer3_glm_betas_wide.csv` can contain both 0 and NaN as stand-ins for pruned channels.
    Do not interpret these as true zero activation; handle as missing downstream.

Usage:
  python merge_homer3_betas_with_combined_data.py \\
    --homer-csv data/tabular/homer3_glm_betas_wide.csv \\
    --combined-csv data/tabular/combined_sfv_data.csv \\
    --out-csv data/results/homer3_betas_plus_combined_sfv_data_inner_join.csv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


_DIGITS_RE = re.compile(r"(\d+)")
_BETA_COL_RE = re.compile(r"^S\d+_D\d+_Cond\d{2}_(HbO|HbR)$")


def _normalize_subject_id(series: pd.Series, *, column_name: str) -> pd.Series:
    """
    Convert a subject id series into a nullable integer ID by extracting digits.

    Examples accepted:
      - 17, "17", "0017"
      - "sub_0017"
    """
    as_str = series.astype("string")
    extracted = as_str.str.extract(_DIGITS_RE, expand=False)
    if extracted.isna().any():
        bad = as_str[extracted.isna()].dropna().unique().tolist()[:10]
        raise ValueError(
            f"Failed to parse numeric IDs from column '{column_name}'. "
            f"Examples of unparseable values: {bad}"
        )
    return extracted.astype("Int64")


def _coerce_beta_columns_to_numeric(df: pd.DataFrame) -> pd.DataFrame:
    """
    Coerce Homer3 beta columns to numeric.

    Homer3 exports sometimes contain numeric values stored as strings. We convert all beta columns
    matching `S##_D##_Cond##_HbO/HbR` to numeric and fail hard if unexpected non-numeric tokens
    are present (to avoid silently corrupting downstream inference).
    """
    beta_cols = [c for c in df.columns if _BETA_COL_RE.match(c)]
    if not beta_cols:
        raise ValueError("No beta columns matched expected pattern like 'S01_D01_Cond01_HbO'.")

    allowed_missing = {"", "NA", "NaN", "nan", "NAN", "NULL", "null"}
    for c in beta_cols:
        raw = df[c].astype("string")
        raw_trim = raw.str.strip()
        coerced = pd.to_numeric(raw_trim, errors="coerce")
        bad_mask = raw_trim.notna() & coerced.isna() & (~raw_trim.isin(allowed_missing))
        if bad_mask.any():
            examples = raw_trim[bad_mask].dropna().unique().tolist()[:10]
            raise ValueError(
                f"Non-numeric values found in Homer3 beta column '{c}'. "
                f"Examples: {examples}. Fix the upstream CSV export before merging."
            )
        df[c] = coerced
    return df


def merge_homer_with_combined(
    homer_csv: Path, combined_csv: Path, out_csv: Path
) -> pd.DataFrame:
    homer = pd.read_csv(homer_csv)
    combined = pd.read_csv(combined_csv)

    if "Subject" not in homer.columns:
        raise ValueError(f"Expected column 'Subject' in {homer_csv}")
    if "subject_id" not in combined.columns:
        raise ValueError(f"Expected column 'subject_id' in {combined_csv}")

    homer = homer.copy()
    combined = combined.copy()

    homer = _coerce_beta_columns_to_numeric(homer)

    homer["subject_id"] = _normalize_subject_id(homer["Subject"], column_name="Subject")
    combined["subject_id"] = _normalize_subject_id(combined["subject_id"], column_name="subject_id")

    if homer["subject_id"].isna().any() or combined["subject_id"].isna().any():
        raise ValueError("Unexpected NA subject_id after normalization.")

    # Fail hard if either input contains duplicate subject IDs. The merge spec assumes one row per subject.
    for df_name, df in [("combined", combined), ("homer", homer)]:
        dup_mask = df["subject_id"].duplicated(keep=False)
        if dup_mask.any():
            counts = df.loc[dup_mask, "subject_id"].astype(int).value_counts().sort_index()
            examples = counts.head(10).to_dict()
            raise ValueError(
                f"Duplicate subject_id values detected in {df_name} dataset after normalization. "
                f"Expected exactly one row per subject for an inner one-to-one join. "
                f"Example duplicate counts (subject_id: n_rows): {examples}. "
                f"Fix the upstream CSV before running the merge."
            )

    # Preserve Homer Subject as a separate column to avoid confusion after join
    homer = homer.rename(columns={"Subject": "homer_subject"})

    merged = combined.merge(homer, how="inner", on="subject_id", validate="one_to_one")

    combined_ids = set(combined["subject_id"].dropna().astype(int).tolist())
    homer_ids = set(homer["subject_id"].dropna().astype(int).tolist())
    merged_ids = set(merged["subject_id"].dropna().astype(int).tolist())

    print(f"[merge] combined subjects: {len(combined_ids)}")
    print(f"[merge] homer subjects:    {len(homer_ids)}")
    print(f"[merge] inner-joined:      {len(merged_ids)}")
    print(f"[merge] combined-only (dropped): {len(combined_ids - merged_ids)}")
    print(f"[merge] homer-only (dropped):    {len(homer_ids - merged_ids)}")

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_csv, index=False)
    print(f"[merge] wrote: {out_csv}")
    return merged


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--homer-csv",
        type=Path,
        default=Path("data/tabular/homer3_glm_betas_wide.csv"),
    )
    parser.add_argument(
        "--combined-csv",
        type=Path,
        default=Path("data/tabular/combined_sfv_data.csv"),
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path("data/results/homer3_betas_plus_combined_sfv_data_inner_join.csv"),
    )
    args = parser.parse_args()

    merge_homer_with_combined(args.homer_csv, args.combined_csv, args.out_csv)


if __name__ == "__main__":
    main()
