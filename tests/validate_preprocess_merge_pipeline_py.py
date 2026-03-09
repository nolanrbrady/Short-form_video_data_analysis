"""Synthetic end-to-end validation for collapse -> merge -> certification.

Run:
  python tests/validate_preprocess_merge_pipeline_py.py
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from collapse_homer_fir_to_auc import collapse_homer_fir_to_auc
from homer_fir import LatentHRFReconstructor, compute_file_sha256, load_preprocessing_settings
from validate_homer_fir_auc_conversion import validate_homer_fir_auc_conversion


def write_settings(path: Path) -> None:
    payload = {
        "latent_hrf_reconstruction": {
            "trange_start": -0.5,
            "trange_stop": 1.0,
            "basis_spacing": 0.5,
            "basis_sigma": 0.5,
        },
        "fir_auc_summary": {
            "baseline_start": -0.5,
            "baseline_stop": 0.0,
            "auc_start": 0.0,
            "auc_stop": 0.5,
        },
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def compute_basis_weights_for_target_auc(settings_path: Path, target_auc: float) -> tuple[float, float]:
    settings = load_preprocessing_settings(settings_path)
    reconstructor = LatentHRFReconstructor(settings.reconstruction, settings.auc_summary)
    _, hrf_one = reconstructor.reconstruct_from_beta_vector((1, 2), np.array([1.0, 0.0]))
    _, hrf_two = reconstructor.reconstruct_from_beta_vector((1, 2), np.array([0.0, 1.0]))
    coeff_one = reconstructor.summarize_auc(hrf_one)
    coeff_two = reconstructor.summarize_auc(hrf_two)
    if np.isclose(coeff_one, 0.0) and np.isclose(coeff_two, 0.0):
        raise AssertionError("Basis AUC coefficients unexpectedly both equal zero.")
    if np.isclose(target_auc, 0.0):
        if not np.isclose(coeff_two, 0.0):
            return (float(1.0), float(-coeff_one / coeff_two))
        return (0.0, 1.0)
    if not np.isclose(coeff_two, 0.0):
        return (0.0, float(target_auc / coeff_two))
    return (float(target_auc / coeff_one), 0.0)


def write_raw_fir_csv(path: Path, settings_path: Path) -> None:
    header = [
        "Subject",
        "S01_D01_Cond01_HbO_Basis001",
        "S01_D01_Cond01_HbO_Basis002",
        "S01_D01_Cond01_HbR_Basis001",
        "S01_D01_Cond01_HbR_Basis002",
        "S01_D01_Cond02_HbO_Basis001",
        "S01_D01_Cond02_HbO_Basis002",
        "S01_D01_Cond02_HbR_Basis001",
        "S01_D01_Cond02_HbR_Basis002",
        "S01_D01_Cond03_HbO_Basis001",
        "S01_D01_Cond03_HbO_Basis002",
        "S01_D01_Cond03_HbR_Basis001",
        "S01_D01_Cond03_HbR_Basis002",
        "S01_D01_Cond04_HbO_Basis001",
        "S01_D01_Cond04_HbO_Basis002",
        "S01_D01_Cond04_HbR_Basis001",
        "S01_D01_Cond04_HbR_Basis002",
    ]
    desired_auc = {
        "sub_0001": {
            "S01_D01_Cond01_HbO": 1.0,
            "S01_D01_Cond01_HbR": -0.5,
            "S01_D01_Cond02_HbO": 0.0,
            "S01_D01_Cond02_HbR": 0.25,
            "S01_D01_Cond03_HbO": 2.0,
            "S01_D01_Cond03_HbR": -1.0,
            "S01_D01_Cond04_HbO": 3.0,
            "S01_D01_Cond04_HbR": -1.5,
        },
        "sub_0002": {
            "S01_D01_Cond01_HbO": -2.0,
            "S01_D01_Cond01_HbR": 1.0,
            "S01_D01_Cond02_HbO": 4.0,
            "S01_D01_Cond02_HbR": -2.0,
            "S01_D01_Cond03_HbO": np.nan,
            "S01_D01_Cond03_HbR": 2.5,
            "S01_D01_Cond04_HbO": -3.0,
            "S01_D01_Cond04_HbR": 1.5,
        },
    }

    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for subject, values in desired_auc.items():
            row = [subject]
            for key in (
                "S01_D01_Cond01_HbO",
                "S01_D01_Cond01_HbR",
                "S01_D01_Cond02_HbO",
                "S01_D01_Cond02_HbR",
                "S01_D01_Cond03_HbO",
                "S01_D01_Cond03_HbR",
                "S01_D01_Cond04_HbO",
                "S01_D01_Cond04_HbR",
            ):
                value = values[key]
                if np.isnan(value):
                    row.extend(["NaN", "NaN"])
                else:
                    beta_one, beta_two = compute_basis_weights_for_target_auc(settings_path, float(value))
                    row.extend([f"{beta_one:.17g}", f"{beta_two:.17g}"])
            writer.writerow(row)


def write_combined_csv(path: Path) -> None:
    df = pd.DataFrame(
        {
            "subject_id": ["0001", "2"],
            "covariate_alpha": [10, 20],
        }
    )
    df.to_csv(path, index=False)


def run_command(args: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(args, cwd=cwd, text=True, capture_output=True, check=False)


def test_end_to_end_merge_and_certification(tmp_path: Path) -> None:
    settings_path = tmp_path / "settings.json"
    raw_fir_csv = tmp_path / "raw_fir.csv"
    auc_csv = tmp_path / "auc.csv"
    combined_csv = tmp_path / "combined.csv"
    merged_csv = tmp_path / "merged.csv"
    cert_json = tmp_path / "cert.json"
    audit_csv = tmp_path / "audit.csv"
    dropped_csv = tmp_path / "dropped.csv"

    write_settings(settings_path)
    write_raw_fir_csv(raw_fir_csv, settings_path)
    write_combined_csv(combined_csv)

    collapse_homer_fir_to_auc(
        input_csv=str(raw_fir_csv),
        output_csv=str(auc_csv),
        settings_json=str(settings_path),
    )
    validate_homer_fir_auc_conversion(
        str(raw_fir_csv),
        str(auc_csv),
        settings_json=str(settings_path),
    )
    auc_df = pd.read_csv(auc_csv)
    row1 = auc_df.loc[auc_df["Subject"] == "sub_0001"].iloc[0]
    row2 = auc_df.loc[auc_df["Subject"] == "sub_0002"].iloc[0]
    assert np.isclose(row1["S01_D01_Cond02_HbO"], 0.0)
    assert np.isfinite(row1["S01_D01_Cond02_HbO"])
    assert np.isnan(row2["S01_D01_Cond03_HbO"])

    merge = run_command(
        [
            "Rscript",
            "merge_homer3_betas_with_combined_data.R",
            "--homer_csv",
            str(auc_csv),
            "--combined_csv",
            str(combined_csv),
            "--out_csv",
            str(merged_csv),
        ],
        cwd=ROOT,
    )
    assert merge.returncode == 0, merge.stderr or merge.stdout

    merged_df = pd.read_csv(merged_csv)
    assert merged_df.shape[0] == 2
    assert "S01_D01_Cond02_HbO" in merged_df.columns
    assert np.isclose(
        merged_df.loc[merged_df["subject_id"] == 1, "S01_D01_Cond02_HbO"].iloc[0], 0.0
    )

    certify = run_command(
        [
            "python",
            "certify_preprocess_merge_integrity.py",
            "--combined_csv",
            str(combined_csv),
            "--homer_csv",
            str(auc_csv),
            "--merged_csv",
            str(merged_csv),
            "--out_json",
            str(cert_json),
            "--out_audit_csv",
            str(audit_csv),
            "--out_dropped_csv",
            str(dropped_csv),
        ],
        cwd=ROOT,
    )
    assert certify.returncode == 0, certify.stderr or certify.stdout
    payload = json.loads(cert_json.read_text(encoding="utf-8"))
    assert payload["passed"] is True
    assert payload["counts"]["beta_columns_in_merged"] == 8


def test_auc_lint_fails_if_excluded_channel_is_not_nan(tmp_path: Path) -> None:
    settings_path = tmp_path / "settings.json"
    raw_fir_csv = tmp_path / "raw_fir.csv"
    auc_csv = tmp_path / "auc.csv"

    write_settings(settings_path)
    write_raw_fir_csv(raw_fir_csv, settings_path)
    collapse_homer_fir_to_auc(
        input_csv=str(raw_fir_csv),
        output_csv=str(auc_csv),
        settings_json=str(settings_path),
    )

    auc_df = pd.read_csv(auc_csv)
    auc_df.loc[auc_df["Subject"] == "sub_0002", "S01_D01_Cond03_HbO"] = 123.0
    auc_df.to_csv(auc_csv, index=False)
    provenance_path = auc_csv.with_suffix(".provenance.json")
    payload = json.loads(provenance_path.read_text(encoding="utf-8"))
    payload["output_csv_sha256"] = compute_file_sha256(auc_csv)
    provenance_path.write_text(json.dumps(payload), encoding="utf-8")

    try:
        validate_homer_fir_auc_conversion(
            str(raw_fir_csv),
            str(auc_csv),
            settings_json=str(settings_path),
        )
    except ValueError as exc:
        assert "should map to NaN" in str(exc)
        return
    raise AssertionError("Expected lint failure when excluded FIR vector maps to non-NaN AUC.")


def test_auc_lint_fails_if_provenance_is_stale(tmp_path: Path) -> None:
    settings_path = tmp_path / "settings.json"
    raw_fir_csv = tmp_path / "raw_fir.csv"
    auc_csv = tmp_path / "auc.csv"

    write_settings(settings_path)
    write_raw_fir_csv(raw_fir_csv, settings_path)
    collapse_homer_fir_to_auc(
        input_csv=str(raw_fir_csv),
        output_csv=str(auc_csv),
        settings_json=str(settings_path),
    )

    payload = json.loads(auc_csv.with_suffix(".provenance.json").read_text(encoding="utf-8"))
    payload["settings"]["fir_auc_summary"]["auc_stop"] = 999.0
    auc_csv.with_suffix(".provenance.json").write_text(
        json.dumps(payload),
        encoding="utf-8",
    )

    try:
        validate_homer_fir_auc_conversion(
            str(raw_fir_csv),
            str(auc_csv),
            settings_json=str(settings_path),
        )
    except ValueError as exc:
        assert "Provenance mismatch" in str(exc)
        return
    raise AssertionError("Expected lint failure when provenance sidecar is stale.")


def test_duplicate_combined_ids_fail_merge(tmp_path: Path) -> None:
    settings_path = tmp_path / "settings.json"
    raw_fir_csv = tmp_path / "raw_fir.csv"
    auc_csv = tmp_path / "auc.csv"
    combined_csv = tmp_path / "combined_dup.csv"
    merged_csv = tmp_path / "merged.csv"

    write_settings(settings_path)
    write_raw_fir_csv(raw_fir_csv, settings_path)
    collapse_homer_fir_to_auc(
        input_csv=str(raw_fir_csv),
        output_csv=str(auc_csv),
        settings_json=str(settings_path),
    )

    df = pd.DataFrame(
        {
            "subject_id": ["0001", "1", "2"],
            "covariate_alpha": [10, 11, 20],
        }
    )
    df.to_csv(combined_csv, index=False)

    merge = run_command(
        [
            "Rscript",
            "merge_homer3_betas_with_combined_data.R",
            "--homer_csv",
            str(auc_csv),
            "--combined_csv",
            str(combined_csv),
            "--out_csv",
            str(merged_csv),
        ],
        cwd=ROOT,
    )
    assert merge.returncode != 0
    assert "Duplicate subject_id values detected in combined dataset" in (
        merge.stderr + merge.stdout
    )


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        test_end_to_end_merge_and_certification(tmp_path)
        test_auc_lint_fails_if_excluded_channel_is_not_nan(tmp_path)
        test_auc_lint_fails_if_provenance_is_stale(tmp_path)
        test_duplicate_combined_ids_fail_merge(tmp_path)

    print("[PASS] validate_preprocess_merge_pipeline_py")


if __name__ == "__main__":
    main()
