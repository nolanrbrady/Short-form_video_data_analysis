"""Validation harness for shared FIR reconstruction and AUC collapse.

Run:
  python tests/validate_fir_auc_adapter_py.py
"""

from __future__ import annotations

import csv
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import collapse_homer_fir_to_auc as fir_auc
from homer_fir import (
    LatentHRFReconstructor,
    build_condition_feature_map,
    compute_file_sha256,
    default_auc_provenance_path,
    extract_beta_vector,
    load_preprocessing_settings,
    parse_fir_header,
    settings_to_dict,
)


def write_settings(path: Path, *, auc_start: float = 0.0, auc_stop: float = 0.5) -> None:
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
            "auc_start": auc_start,
            "auc_stop": auc_stop,
        },
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def write_production_settings(path: Path) -> None:
    payload = {
        "latent_hrf_reconstruction": {
            "trange_start": -10.0,
            "trange_stop": 130.0,
            "basis_spacing": 0.5,
            "basis_sigma": 0.5,
        },
        "fir_auc_summary": {
            "baseline_start": -10.0,
            "baseline_stop": 0.0,
            "auc_start": 0.0,
            "auc_stop": 120.0,
        },
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def write_oscillatory_settings(path: Path) -> None:
    payload = {
        "latent_hrf_reconstruction": {
            "trange_start": -0.5,
            "trange_stop": 2.0,
            "basis_spacing": 0.5,
            "basis_sigma": 0.5,
        },
        "fir_auc_summary": {
            "baseline_start": -0.5,
            "baseline_stop": 0.0,
            "auc_start": 0.0,
            "auc_stop": 1.5,
        },
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def write_toy_csv(path: Path) -> None:
    header = [
        "Subject",
        "S01_D01_Cond01_HbO_Basis001",
        "S01_D01_Cond01_HbO_Basis002",
        "S01_D01_Cond01_HbR_Basis001",
        "S01_D01_Cond01_HbR_Basis002",
        "S01_D01_Cond02_HbO_Basis001",
        "S01_D01_Cond02_HbO_Basis002",
    ]
    rows = [
        ["sub_0001", "1.0", "2.0", "0", "0", "2.0", "4.0"],
        ["sub_0002", "3.0", "5.0", "NaN", "NaN", "6.0", "8.0"],
    ]
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def write_partial_nan_csv(path: Path) -> None:
    header = [
        "Subject",
        "S01_D01_Cond01_HbO_Basis001",
        "S01_D01_Cond01_HbO_Basis002",
    ]
    rows = [["sub_0001", "1.0", "NaN"]]
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def write_duplicate_subject_csv(path: Path) -> None:
    header = [
        "Subject",
        "S01_D01_Cond01_HbO_Basis001",
        "S01_D01_Cond01_HbO_Basis002",
    ]
    rows = [["sub_0001", "1.0", "2.0"], ["sub_0001", "3.0", "4.0"]]
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def write_noncontiguous_basis_csv(path: Path) -> None:
    header = [
        "Subject",
        "S01_D01_Cond01_HbO_Basis001",
        "S01_D01_Cond01_HbO_Basis003",
    ]
    rows = [["sub_0001", "1.0", "2.0"]]
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def write_basis_mismatch_header_csv(path: Path) -> None:
    header = [
        "Subject",
        "S01_D01_Cond01_HbO_Basis001",
        "S01_D01_Cond01_HbO_Basis002",
        "S01_D01_Cond01_HbR_Basis001",
        "S01_D01_Cond01_HbR_Basis002",
        "S01_D01_Cond01_HbR_Basis003",
    ]
    rows = [["sub_0001", "1.0", "2.0", "3.0", "4.0", "5.0"]]
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def write_production_shape_csv(path: Path, settings_path: Path) -> None:
    settings = load_preprocessing_settings(settings_path)
    reconstructor = LatentHRFReconstructor(settings.reconstruction, settings.auc_summary)
    basis_count = reconstructor.basis_matrix.shape[1]
    header = ["Subject", *[f"S01_D01_Cond01_HbO_Basis{i:03d}" for i in range(1, basis_count + 1)]]

    beta = np.zeros(basis_count, dtype=float)
    beta[0] = 0.75
    beta[139] = -0.25
    beta[-1] = 0.5

    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerow(["sub_0001", *[f"{value:.17g}" for value in beta]])


def manual_hrf(beta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    expected_t = np.array([-0.5, 0.0, 0.5], dtype=float)
    g1 = np.exp(-((expected_t - 0.0) ** 2) / (2 * 0.5**2))
    g1 = g1 / g1.max()
    g2 = np.exp(-((expected_t - 0.5) ** 2) / (2 * 0.5**2))
    g2 = g2 / g2.max()
    return expected_t, g1 * beta[0] + g2 * beta[1]


def test_auc_matches_manual_trapezoid(settings_path: Path) -> None:
    settings = load_preprocessing_settings(settings_path)
    reconstructor = LatentHRFReconstructor(settings.reconstruction, settings.auc_summary)
    _, reconstructed = reconstructor.reconstruct_from_beta_vector((1, 2), np.array([1.0, 2.0]))
    t, expected_hrf = manual_hrf(np.array([1.0, 2.0]))
    assert np.allclose(reconstructor.time_axis, t)
    assert np.allclose(reconstructed, expected_hrf)

    baseline_value = np.mean(expected_hrf[(t >= -0.5) & (t <= 0.0)])
    expected_auc = np.trapezoid(
        (expected_hrf - baseline_value)[(t >= 0.0) & (t <= 0.5)],
        t[(t >= 0.0) & (t <= 0.5)],
    )
    observed_auc = reconstructor.summarize_auc(reconstructed)
    assert np.isclose(observed_auc, expected_auc)


def test_baseline_correction_removes_constant_offset(settings_path: Path) -> None:
    settings = load_preprocessing_settings(settings_path)
    reconstructor = LatentHRFReconstructor(settings.reconstruction, settings.auc_summary)
    _, hrf = manual_hrf(np.array([1.0, 2.0]))
    auc_one = reconstructor.summarize_auc(hrf)
    auc_two = reconstructor.summarize_auc(hrf + 100.0)
    assert np.isclose(auc_one, auc_two)


def test_legitimate_zero_auc_is_finite(settings_path: Path) -> None:
    settings = load_preprocessing_settings(settings_path)
    reconstructor = LatentHRFReconstructor(settings.reconstruction, settings.auc_summary)
    hrf = np.array([1.0, 1.0, 1.0], dtype=float)
    auc = reconstructor.summarize_auc(hrf)
    assert np.isfinite(auc)
    assert np.isclose(auc, 0.0)


def test_oscillatory_mixed_sign_basis_can_yield_negative_auc(settings_path: Path) -> None:
    settings = load_preprocessing_settings(settings_path)
    reconstructor = LatentHRFReconstructor(settings.reconstruction, settings.auc_summary)
    beta = np.array([1.5, -4.0, 2.0, -5.0], dtype=float)
    time_axis, hrf = reconstructor.reconstruct_from_beta_vector((1, 2, 3, 4), beta)

    assert np.any(hrf > 0.0)
    assert np.any(hrf < 0.0)

    baseline_mask = (time_axis >= -0.5) & (time_axis <= 0.0)
    auc_mask = (time_axis >= 0.0) & (time_axis <= 1.5)
    baseline_value = float(np.mean(hrf[baseline_mask]))
    expected_auc = np.trapezoid(hrf[auc_mask] - baseline_value, time_axis[auc_mask])
    observed_auc = reconstructor.summarize_auc(hrf)

    assert np.isclose(observed_auc, expected_auc)
    assert observed_auc < 0.0


def test_non_grid_aligned_window_fails(settings_path: Path) -> None:
    payload = json.loads(settings_path.read_text(encoding="utf-8"))
    payload["fir_auc_summary"]["auc_stop"] = 0.4
    settings_path.write_text(json.dumps(payload), encoding="utf-8")
    settings = load_preprocessing_settings(settings_path)
    try:
        LatentHRFReconstructor(settings.reconstruction, settings.auc_summary)
    except ValueError as exc:
        assert "not represented on the HRF time axis" in str(exc)
        return
    raise AssertionError("Expected ValueError for non-grid-aligned AUC window.")


def test_missing_settings_section_fails(settings_path: Path) -> None:
    payload = {"latent_hrf_reconstruction": {"trange_start": -0.5}}
    settings_path.write_text(json.dumps(payload), encoding="utf-8")
    try:
        load_preprocessing_settings(settings_path)
    except ValueError as exc:
        assert "Missing required settings section" in str(exc)
        return
    raise AssertionError("Expected ValueError for missing settings section.")


def test_invalid_sigma_fails(settings_path: Path) -> None:
    payload = json.loads(settings_path.read_text(encoding="utf-8"))
    payload["latent_hrf_reconstruction"]["basis_sigma"] = 0.0
    settings_path.write_text(json.dumps(payload), encoding="utf-8")
    try:
        load_preprocessing_settings(settings_path)
    except ValueError as exc:
        assert "values must be > 0" in str(exc)
        return
    raise AssertionError("Expected ValueError for nonpositive basis sigma.")


def test_windows_outside_trange_fail(settings_path: Path) -> None:
    payload = json.loads(settings_path.read_text(encoding="utf-8"))
    payload["fir_auc_summary"]["baseline_start"] = -11.0
    settings_path.write_text(json.dumps(payload), encoding="utf-8")
    try:
        load_preprocessing_settings(settings_path)
    except ValueError as exc:
        assert "falls outside reconstruction trange" in str(exc)
        return
    raise AssertionError("Expected ValueError when summary windows fall outside trange.")


def test_noncontiguous_basis_numbering_fails(input_csv: Path) -> None:
    header = next(csv.reader(input_csv.open("r", newline="")))
    try:
        parse_fir_header(header)
    except ValueError as exc:
        assert "Expected contiguous basis numbering" in str(exc)
        return
    raise AssertionError("Expected ValueError for noncontiguous basis numbering.")


def test_hbo_hbr_basis_mismatch_fails(input_csv: Path) -> None:
    header = next(csv.reader(input_csv.open("r", newline="")))
    specs = parse_fir_header(header)
    try:
        build_condition_feature_map(specs, "01")
    except ValueError as exc:
        assert "basis mismatch between HbO/HbR" in str(exc)
        return
    raise AssertionError("Expected ValueError for HbO/HbR basis mismatch.")


def test_collapse_script_outputs_expected_values(
    input_csv: Path, output_csv: Path, settings_path: Path
) -> None:
    fir_auc.collapse_homer_fir_to_auc(
        input_csv=str(input_csv),
        output_csv=str(output_csv),
        settings_json=str(settings_path),
    )
    df = pd.read_csv(output_csv)
    assert list(df.columns) == [
        "Subject",
        "S01_D01_Cond01_HbO",
        "S01_D01_Cond01_HbR",
        "S01_D01_Cond02_HbO",
    ]

    settings = load_preprocessing_settings(settings_path)
    reconstructor = LatentHRFReconstructor(settings.reconstruction, settings.auc_summary)
    expected_hbo = reconstructor.summarize_auc(manual_hrf(np.array([1.0, 2.0]))[1])
    expected_cond02 = reconstructor.summarize_auc(manual_hrf(np.array([2.0, 4.0]))[1])

    row1 = df.loc[df["Subject"] == "sub_0001"].iloc[0]
    row2 = df.loc[df["Subject"] == "sub_0002"].iloc[0]

    assert np.isclose(row1["S01_D01_Cond01_HbO"], expected_hbo)
    assert np.isnan(row1["S01_D01_Cond01_HbR"])
    assert np.isclose(row1["S01_D01_Cond02_HbO"], expected_cond02)
    assert np.isnan(row2["S01_D01_Cond01_HbR"])


def test_provenance_sidecar_matches_generated_outputs(
    input_csv: Path, output_csv: Path, settings_path: Path
) -> None:
    fir_auc.collapse_homer_fir_to_auc(
        input_csv=str(input_csv),
        output_csv=str(output_csv),
        settings_json=str(settings_path),
    )
    provenance_path = default_auc_provenance_path(output_csv)
    payload = json.loads(provenance_path.read_text(encoding="utf-8"))

    settings = load_preprocessing_settings(settings_path)
    assert payload["analysis_step"] == "collapse_homer_fir_to_auc"
    assert payload["basis_family"] == fir_auc.BASIS_FAMILY
    assert payload["input_csv"] == str(input_csv)
    assert payload["input_csv_sha256"] == compute_file_sha256(input_csv)
    assert payload["settings_json"] == str(settings_path)
    assert payload["settings_json_sha256"] == compute_file_sha256(settings_path)
    assert payload["output_csv"] == str(output_csv)
    assert payload["output_csv_sha256"] == compute_file_sha256(output_csv)
    assert payload["row_count"] == 2
    assert payload["feature_count"] == 3
    assert payload["basis_count_per_feature"] == 2
    assert payload["settings"] == settings_to_dict(settings)


def test_production_geometry_matches_real_homer_configuration(
    input_csv: Path, output_csv: Path, settings_path: Path
) -> None:
    settings = load_preprocessing_settings(settings_path)
    reconstructor = LatentHRFReconstructor(settings.reconstruction, settings.auc_summary)
    assert reconstructor.time_axis.shape == (280,)
    assert reconstructor.basis_matrix.shape == (280, 279)
    assert np.isclose(reconstructor.time_axis[0], -10.0)
    assert np.isclose(reconstructor.time_axis[-1], 129.5)
    assert np.any(np.isclose(reconstructor.time_axis, -10.0))
    assert np.any(np.isclose(reconstructor.time_axis, 0.0))
    assert np.any(np.isclose(reconstructor.time_axis, 120.0))

    fir_auc.collapse_homer_fir_to_auc(
        input_csv=str(input_csv),
        output_csv=str(output_csv),
        settings_json=str(settings_path),
    )
    df = pd.read_csv(output_csv)
    assert list(df.columns) == ["Subject", "S01_D01_Cond01_HbO"]

    with input_csv.open("r", newline="") as f:
        header = next(csv.reader(f))
    spec = next(iter(parse_fir_header(header).values()))
    with input_csv.open("r", newline="") as f:
        reader = csv.reader(f)
        next(reader)
        raw_row = next(reader)
    beta = extract_beta_vector(raw_row, spec)
    _, expected_hrf = reconstructor.reconstruct_from_beta_vector(spec.basis_idx, beta)
    expected_auc = reconstructor.summarize_auc(expected_hrf)

    observed_auc = df.loc[df["Subject"] == "sub_0001", "S01_D01_Cond01_HbO"].iloc[0]
    assert np.isclose(observed_auc, expected_auc)


def test_partial_missing_basis_vector_fails(
    input_csv: Path, output_csv: Path, settings_path: Path
) -> None:
    try:
        fir_auc.collapse_homer_fir_to_auc(
            input_csv=str(input_csv),
            output_csv=str(output_csv),
            settings_json=str(settings_path),
        )
    except ValueError as exc:
        assert "partial NaNs" in str(exc)
        return
    raise AssertionError("Expected ValueError for partial missing basis vector.")


def test_duplicate_subject_fails(
    input_csv: Path, output_csv: Path, settings_path: Path
) -> None:
    try:
        fir_auc.collapse_homer_fir_to_auc(
            input_csv=str(input_csv),
            output_csv=str(output_csv),
            settings_json=str(settings_path),
        )
    except ValueError as exc:
        assert "Duplicate subject row detected" in str(exc)
        return
    raise AssertionError("Expected ValueError for duplicate Subject rows.")


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)

        settings_path = tmp_path / "settings.json"
        write_settings(settings_path)
        test_auc_matches_manual_trapezoid(settings_path)
        write_settings(settings_path)
        test_baseline_correction_removes_constant_offset(settings_path)
        write_settings(settings_path)
        test_legitimate_zero_auc_is_finite(settings_path)
        write_oscillatory_settings(settings_path)
        test_oscillatory_mixed_sign_basis_can_yield_negative_auc(settings_path)
        write_settings(settings_path)
        test_non_grid_aligned_window_fails(settings_path)
        write_settings(settings_path)
        test_missing_settings_section_fails(settings_path)
        write_settings(settings_path)
        test_invalid_sigma_fails(settings_path)
        write_settings(settings_path)
        test_windows_outside_trange_fail(settings_path)

        input_csv = tmp_path / "toy_homer.csv"
        output_csv = tmp_path / "toy_auc.csv"
        write_settings(settings_path)
        write_toy_csv(input_csv)
        test_collapse_script_outputs_expected_values(input_csv, output_csv, settings_path)
        write_settings(settings_path)
        test_provenance_sidecar_matches_generated_outputs(input_csv, output_csv, settings_path)

        noncontig_csv = tmp_path / "noncontiguous.csv"
        write_noncontiguous_basis_csv(noncontig_csv)
        test_noncontiguous_basis_numbering_fails(noncontig_csv)

        mismatch_csv = tmp_path / "basis_mismatch.csv"
        write_basis_mismatch_header_csv(mismatch_csv)
        test_hbo_hbr_basis_mismatch_fails(mismatch_csv)

        partial_csv = tmp_path / "partial_nan.csv"
        write_settings(settings_path)
        write_partial_nan_csv(partial_csv)
        test_partial_missing_basis_vector_fails(partial_csv, output_csv, settings_path)

        dup_csv = tmp_path / "duplicate_subjects.csv"
        write_settings(settings_path)
        write_duplicate_subject_csv(dup_csv)
        test_duplicate_subject_fails(dup_csv, output_csv, settings_path)

        production_settings_path = tmp_path / "production_settings.json"
        production_input_csv = tmp_path / "production_shape_homer.csv"
        production_output_csv = tmp_path / "production_shape_auc.csv"
        write_production_settings(production_settings_path)
        write_production_shape_csv(production_input_csv, production_settings_path)
        test_production_geometry_matches_real_homer_configuration(
            production_input_csv,
            production_output_csv,
            production_settings_path,
        )

    print("[PASS] validate_fir_auc_adapter_py")


if __name__ == "__main__":
    main()
