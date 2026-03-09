"""Validation harness for plot_fir_betas_subjects.py.

Run:
  python tests/validate_fir_plot_script_py.py
"""

from __future__ import annotations

import csv
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import plot_fir_betas_subjects as fir_plot


def write_toy_csv(path: Path) -> None:
    header = [
        "Subject",
        "S01_D01_Cond01_HbO_Basis001",
        "S01_D01_Cond01_HbR_Basis001",
        "S01_D01_Cond01_HbO_Basis002",
        "S01_D01_Cond01_HbR_Basis002",
        "S01_D02_Cond01_HbO_Basis001",
        "S01_D02_Cond01_HbR_Basis001",
        "S01_D02_Cond01_HbO_Basis002",
        "S01_D02_Cond01_HbR_Basis002",
        "S01_D01_Cond02_HbO_Basis001",
        "S01_D01_Cond02_HbR_Basis001",
    ]
    rows = [
        [
            "sub_0001",
            "0",  # excluded channel -> all-zero basis vector
            "1.0",
            "0",
            "2.0",
            "3.0",
            "4.0",
            "5.0",
            "6.0",
            "777.0",  # Cond02 (must be ignored for Cond01 test)
            "888.0",
        ],
        [
            "sub_0002",
            "10.0",
            "11.0",
            "12.0",
            "13.0",
            "14.0",
            "15.0",
            "16.0",
            "17.0",
            "999.0",
            "1000.0",
        ],
        [
            "sub_0003",
            "1.0",
            "2.0",
            "nan",  # partial NaN is invalid for exact reconstruction
            "3.0",
            "4.0",
            "5.0",
            "6.0",
            "7.0",
            "111.0",
            "222.0",
        ],
    ]

    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def write_toy_settings(path: Path) -> None:
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


def test_parse_condition_columns() -> None:
    header = [
        "Subject",
        "S01_D01_Cond01_HbO_Basis001",
        "S01_D01_Cond01_HbR_Basis001",
        "S01_D01_Cond02_HbO_Basis001",
        "S01_D01_Cond02_HbR_Basis001",
    ]
    parsed = fir_plot.parse_condition_columns(header, "01")
    assert sorted(parsed.keys()) == ["S01_D01"]
    assert list(parsed["S01_D01"]["HbO"].basis_idx) == [1]
    assert list(parsed["S01_D01"]["HbR"].basis_idx) == [1]

    parsed_target = fir_plot.parse_condition_columns(
        header, "01", target_channel="S01_D01"
    )
    assert sorted(parsed_target.keys()) == ["S01_D01"]


def test_extract_target_rows_and_missingness(csv_path: Path) -> None:
    header, rows = fir_plot.extract_target_rows(
        csv_path, ["sub_0001", "sub_0002", "sub_0003"]
    )
    cond_map = fir_plot.parse_condition_columns(header, "01")
    basis_idx, betas = fir_plot.extract_subject_betas(rows["sub_0001"], cond_map)

    assert basis_idx == [1, 2]
    assert np.isnan(betas["S01_D01"]["HbO"][0])
    assert np.isnan(betas["S01_D01"]["HbO"][1])
    assert np.isclose(betas["S01_D01"]["HbR"][0], 1.0)
    assert np.isclose(betas["S01_D01"]["HbR"][1], 2.0)


def test_reconstruct_hrf_from_betas() -> None:
    beta = np.array([1.0, 2.0], dtype=float)
    t_hrf, hrf = fir_plot.reconstruct_hrf_from_betas(
        basis_idx=[1, 2],
        beta=beta,
        trange=(-0.5, 1.0),
        basis_spacing=0.5,
        basis_sigma=0.5,
    )

    expected_t = np.array([-0.5, 0.0, 0.5], dtype=float)
    assert np.allclose(t_hrf, expected_t)

    g1 = np.exp(-((expected_t - 0.0) ** 2) / (2 * 0.5**2))
    g1 = g1 / g1.max()
    g2 = np.exp(-((expected_t - 0.5) ** 2) / (2 * 0.5**2))
    g2 = g2 / g2.max()
    expected_hrf = g1 * 1.0 + g2 * 2.0
    assert np.allclose(hrf, expected_hrf)


def test_partial_missing_hrf_reconstruction_fails() -> None:
    try:
        fir_plot.reconstruct_hrf_from_betas(
            basis_idx=[1, 2],
            beta=np.array([1.0, np.nan], dtype=float),
            trange=(-0.5, 1.0),
            basis_spacing=0.5,
            basis_sigma=0.5,
        )
    except ValueError as exc:
        assert "missing/pruned" in str(exc)
        return
    raise AssertionError("Expected ValueError for partial missing basis weights.")


def test_missing_subject_fails(csv_path: Path) -> None:
    try:
        fir_plot.extract_target_rows(csv_path, ["sub_0001", "sub_9999"])
    except ValueError as exc:
        assert "Missing requested subject(s)" in str(exc)
        return
    raise AssertionError("Expected ValueError for missing target subject.")


def test_schema_mismatch_fails() -> None:
    bad_header = [
        "Subject",
        "S01_D01_Cond01_HbO_Basis001",
    ]
    try:
        fir_plot.parse_condition_columns(bad_header, "01")
    except ValueError as exc:
        assert "missing HbR columns" in str(exc)
        return
    raise AssertionError("Expected ValueError for HbO/HbR schema mismatch.")


def test_target_channel_missing_fails(csv_path: Path) -> None:
    header, _ = fir_plot.extract_target_rows(
        csv_path, ["sub_0001", "sub_0002", "sub_0003"]
    )
    try:
        fir_plot.parse_condition_columns(header, "01", target_channel="S99_D99")
    except ValueError as exc:
        assert "Requested target channel S99_D99 not found in Cond01" in str(exc)
        return
    raise AssertionError("Expected ValueError for missing target channel.")


def test_extract_subject_betas_fails_on_partial_nan(csv_path: Path) -> None:
    header, rows = fir_plot.extract_target_rows(
        csv_path, ["sub_0001", "sub_0002", "sub_0003"]
    )
    cond_map = fir_plot.parse_condition_columns(header, "01")
    try:
        fir_plot.extract_subject_betas(rows["sub_0003"], cond_map)
    except ValueError as exc:
        assert "partial NaNs" in str(exc)
        return
    raise AssertionError("Expected ValueError for partial NaNs in one basis vector.")


def test_output_files_created(csv_path: Path, out_dir: Path, settings_json: Path) -> None:
    out_paths = fir_plot.run_plotting(
        input_csv=str(csv_path),
        output_dir=str(out_dir),
        target_subjects=["sub_0002"],
        target_condition="01",
        target_channel="S01_D01",
        settings_json=str(settings_json),
        dpi=120,
    )
    assert len(out_paths) == 1
    for p in out_paths:
        assert p.exists(), f"Expected output plot missing: {p}"
        assert "S01_D01" in p.name
        assert "fir_hrf" in p.name


def test_run_plotting_requires_target_channel(csv_path: Path, out_dir: Path, settings_json: Path) -> None:
    try:
        fir_plot.run_plotting(
            input_csv=str(csv_path),
            output_dir=str(out_dir),
            target_subjects=["sub_0001"],
            target_condition="01",
            target_channel=None,
            settings_json=str(settings_json),
            dpi=120,
        )
    except ValueError as exc:
        assert "target_channel must be explicitly set" in str(exc)
        return
    raise AssertionError("Expected ValueError when target_channel is None.")


def test_run_plotting_fails_on_partial_missing_basis_weights(
    csv_path: Path, out_dir: Path, settings_json: Path
) -> None:
    try:
        fir_plot.run_plotting(
            input_csv=str(csv_path),
            output_dir=str(out_dir),
            target_subjects=["sub_0003"],
            target_condition="01",
            target_channel="S01_D01",
            settings_json=str(settings_json),
            dpi=120,
        )
    except ValueError as exc:
        assert "partial NaNs" in str(exc)
        return
    raise AssertionError("Expected ValueError for partial missing basis weights.")


def main() -> None:
    test_parse_condition_columns()

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        csv_path = tmp_path / "toy_homer.csv"
        out_dir = tmp_path / "plots"
        settings_json = tmp_path / "preprocessing_settings.json"
        write_toy_csv(csv_path)
        write_toy_settings(settings_json)

        test_extract_target_rows_and_missingness(csv_path)
        test_reconstruct_hrf_from_betas()
        test_partial_missing_hrf_reconstruction_fails()
        test_missing_subject_fails(csv_path)
        test_schema_mismatch_fails()
        test_target_channel_missing_fails(csv_path)
        test_extract_subject_betas_fails_on_partial_nan(csv_path)
        test_output_files_created(csv_path, out_dir, settings_json)
        test_run_plotting_requires_target_channel(csv_path, out_dir, settings_json)
        test_run_plotting_fails_on_partial_missing_basis_weights(
            csv_path, out_dir, settings_json
        )

    print("[PASS] validate_fir_plot_script_py")


if __name__ == "__main__":
    main()
