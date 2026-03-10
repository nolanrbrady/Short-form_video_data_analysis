"""Validation harness for between-subject Homer AUC outlier masking.

Run:
  python tests/validate_auc_outlier_masking_py.py
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mask_homer_auc_between_subject_outliers import mask_between_subject_outliers


def write_auc_csv(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, na_rep="NaN")


def make_base_auc_df() -> pd.DataFrame:
    subjects = [f"sub_{i:04d}" for i in range(1, 12)]
    return pd.DataFrame(
        {
            "Subject": subjects,
            "S01_D01_Cond01_HbO": [0.0] * 10 + [1.0],
            "S01_D01_Cond01_HbR": [0.25] * 11,
            "S01_D01_Cond02_HbO": [0.0, np.nan] + [0.0] * 9,
            "S01_D01_Cond02_HbR": np.linspace(-0.5, 0.5, num=11),
        }
    )


def test_masks_only_between_subject_outliers(tmp_path: Path) -> None:
    input_csv = tmp_path / "auc.csv"
    output_csv = tmp_path / "auc_masked.csv"
    audit_csv = tmp_path / "audit.csv"
    summary_json = tmp_path / "summary.json"

    write_auc_csv(input_csv, make_base_auc_df())
    mask_between_subject_outliers(
        input_csv=str(input_csv),
        output_csv=str(output_csv),
        out_audit_csv=str(audit_csv),
        out_summary_json=str(summary_json),
    )

    masked = pd.read_csv(output_csv)
    assert np.isnan(masked.loc[masked["Subject"] == "sub_0011", "S01_D01_Cond01_HbO"].iloc[0])
    assert np.isclose(
        masked.loc[masked["Subject"] == "sub_0001", "S01_D01_Cond02_HbR"].iloc[0],
        -0.5,
    )

    audit = pd.read_csv(audit_csv)
    assert audit.shape[0] == 1
    assert audit.loc[0, "Subject"] == "sub_0011"
    assert audit.loc[0, "channel_column"] == "S01_D01_Cond01_HbO"

    summary = json.loads(summary_json.read_text(encoding="utf-8"))
    assert summary["masked_beta_cells"] == 1


def test_existing_nan_is_preserved_and_ignored_in_threshold_estimation(tmp_path: Path) -> None:
    input_csv = tmp_path / "auc.csv"
    output_csv = tmp_path / "auc_masked.csv"
    audit_csv = tmp_path / "audit.csv"
    summary_json = tmp_path / "summary.json"

    df = make_base_auc_df()
    write_auc_csv(input_csv, df)
    mask_between_subject_outliers(
        input_csv=str(input_csv),
        output_csv=str(output_csv),
        out_audit_csv=str(audit_csv),
        out_summary_json=str(summary_json),
    )

    masked = pd.read_csv(output_csv)
    assert np.isnan(masked.loc[masked["Subject"] == "sub_0002", "S01_D01_Cond02_HbO"].iloc[0])
    assert pd.read_csv(audit_csv).loc[:, "channel_column"].tolist() == ["S01_D01_Cond01_HbO"]


def test_zero_sd_columns_are_skipped(tmp_path: Path) -> None:
    input_csv = tmp_path / "auc.csv"
    output_csv = tmp_path / "auc_masked.csv"
    audit_csv = tmp_path / "audit.csv"
    summary_json = tmp_path / "summary.json"

    write_auc_csv(input_csv, make_base_auc_df())
    mask_between_subject_outliers(
        input_csv=str(input_csv),
        output_csv=str(output_csv),
        out_audit_csv=str(audit_csv),
        out_summary_json=str(summary_json),
    )

    summary = json.loads(summary_json.read_text(encoding="utf-8"))
    zero_sd_row = next(
        row for row in summary["per_column_summary"] if row["column"] == "S01_D01_Cond01_HbR"
    )
    assert zero_sd_row["skip_reason"] == "sd_non_finite_or_zero"
    masked = pd.read_csv(output_csv)
    assert np.isclose(
        masked.loc[masked["Subject"] == "sub_0001", "S01_D01_Cond01_HbR"].iloc[0],
        0.25,
    )


def test_duplicate_subjects_fail_fast(tmp_path: Path) -> None:
    input_csv = tmp_path / "auc.csv"
    output_csv = tmp_path / "auc_masked.csv"
    audit_csv = tmp_path / "audit.csv"
    summary_json = tmp_path / "summary.json"

    df = make_base_auc_df()
    dup_row = df.iloc[[0]].copy()
    dup_row.loc[:, "Subject"] = "sub_0001"
    write_auc_csv(input_csv, pd.concat([df, dup_row], ignore_index=True))

    try:
        mask_between_subject_outliers(
            input_csv=str(input_csv),
            output_csv=str(output_csv),
            out_audit_csv=str(audit_csv),
            out_summary_json=str(summary_json),
        )
    except ValueError as exc:
        assert "Duplicate subject rows detected" in str(exc)
        return
    raise AssertionError("Expected duplicate Subject rows to fail fast.")


def test_non_numeric_beta_token_fails_fast(tmp_path: Path) -> None:
    input_csv = tmp_path / "auc.csv"
    output_csv = tmp_path / "auc_masked.csv"
    audit_csv = tmp_path / "audit.csv"
    summary_json = tmp_path / "summary.json"

    df = make_base_auc_df()
    df["S01_D01_Cond02_HbR"] = df["S01_D01_Cond02_HbR"].astype("object")
    df.loc[0, "S01_D01_Cond02_HbR"] = "not-a-number"
    write_auc_csv(input_csv, df)

    try:
        mask_between_subject_outliers(
            input_csv=str(input_csv),
            output_csv=str(output_csv),
            out_audit_csv=str(audit_csv),
            out_summary_json=str(summary_json),
        )
    except ValueError as exc:
        assert "Non-numeric values found in Homer beta column" in str(exc)
        return
    raise AssertionError("Expected a non-numeric beta token to fail fast.")


def test_invalid_min_n_fails_for_three_sd_rule(tmp_path: Path) -> None:
    input_csv = tmp_path / "auc.csv"
    output_csv = tmp_path / "auc_masked.csv"
    audit_csv = tmp_path / "audit.csv"
    summary_json = tmp_path / "summary.json"

    write_auc_csv(input_csv, make_base_auc_df())

    try:
        mask_between_subject_outliers(
            input_csv=str(input_csv),
            output_csv=str(output_csv),
            out_audit_csv=str(audit_csv),
            out_summary_json=str(summary_json),
            min_n=10,
        )
    except ValueError as exc:
        assert "--min-n must be >= 11" in str(exc)
        return
    raise AssertionError("Expected min_n < 11 to fail for the 3-SD rule.")


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        test_masks_only_between_subject_outliers(tmp_path)
        test_existing_nan_is_preserved_and_ignored_in_threshold_estimation(tmp_path)
        test_zero_sd_columns_are_skipped(tmp_path)
        test_duplicate_subjects_fail_fast(tmp_path)
        test_non_numeric_beta_token_fails_fast(tmp_path)
        test_invalid_min_n_fails_for_three_sd_rule(tmp_path)

    print("[PASS] validate_auc_outlier_masking_py")


if __name__ == "__main__":
    main()
