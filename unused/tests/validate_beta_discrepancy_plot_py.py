"""Validation harness for plot_beta_discrepancy_dynamics.py.

Run:
  python tests/validate_beta_discrepancy_plot_py.py
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import plot_beta_discrepancy_dynamics as beta_plot


def write_toy_merged_csv(path: Path) -> None:
    rows = [
        {
            "subject_id": "0001",
            "S04_D02_Cond01_HbR": -5.0,
            "S04_D02_Cond02_HbR": -7.0,
            "S04_D02_Cond03_HbR": -1.0,
            "S04_D02_Cond04_HbR": -3.0,
            "S04_D02_Cond01_HbO": 1.0,
            "S04_D02_Cond02_HbO": 2.0,
            "S04_D02_Cond03_HbO": 3.0,
            "S04_D02_Cond04_HbO": 4.0,
            "S03_D02_Cond01_HbO": 1.0,
            "S03_D02_Cond02_HbO": 2.0,
            "S03_D02_Cond03_HbO": 3.0,
            "S03_D02_Cond04_HbO": 4.0,
            "S03_D04_Cond01_HbO": 2.0,
            "S03_D04_Cond02_HbO": 3.0,
            "S03_D04_Cond03_HbO": 4.0,
            "S03_D04_Cond04_HbO": 5.0,
        },
        {
            "subject_id": "0002",
            "S04_D02_Cond01_HbR": -4.0,
            "S04_D02_Cond02_HbR": -6.0,
            "S04_D02_Cond03_HbR": -2.0,
            "S04_D02_Cond04_HbR": -3.5,
            "S04_D02_Cond01_HbO": 2.0,
            "S04_D02_Cond02_HbO": 0.0,  # zero should be treated as pruned/missing
            "S04_D02_Cond03_HbO": 4.0,
            "S04_D02_Cond04_HbO": 5.0,
            "S03_D02_Cond01_HbO": 2.0,
            "S03_D02_Cond02_HbO": 2.0,
            "S03_D02_Cond03_HbO": 4.0,
            "S03_D02_Cond04_HbO": 5.0,
            "S03_D04_Cond01_HbO": 4.0,
            "S03_D04_Cond02_HbO": 4.0,
            "S03_D04_Cond03_HbO": 6.0,
            "S03_D04_Cond04_HbO": 7.0,
        },
        {
            "subject_id": "0003",
            "S04_D02_Cond01_HbR": -3.0,
            "S04_D02_Cond02_HbR": -4.0,
            "S04_D02_Cond03_HbR": -2.0,
            "S04_D02_Cond04_HbR": "",
            "S04_D02_Cond01_HbO": 3.0,
            "S04_D02_Cond02_HbO": 3.5,
            "S04_D02_Cond03_HbO": 4.5,
            "S04_D02_Cond04_HbO": 6.0,
            "S03_D02_Cond01_HbO": 3.0,
            "S03_D02_Cond02_HbO": 4.0,
            "S03_D02_Cond03_HbO": 5.0,
            "S03_D02_Cond04_HbO": 6.0,
            "S03_D04_Cond01_HbO": 5.0,
            "S03_D04_Cond02_HbO": 6.0,
            "S03_D04_Cond03_HbO": 7.0,
            "S03_D04_Cond04_HbO": 8.0,
        },
        {
            "subject_id": "0004",
            "S04_D02_Cond01_HbR": -9.0,
            "S04_D02_Cond02_HbR": -9.0,
            "S04_D02_Cond03_HbR": -9.0,
            "S04_D02_Cond04_HbR": -9.0,
            "S04_D02_Cond01_HbO": 9.0,
            "S04_D02_Cond02_HbO": 9.0,
            "S04_D02_Cond03_HbO": 9.0,
            "S04_D02_Cond04_HbO": 9.0,
            "S03_D02_Cond01_HbO": 9.0,
            "S03_D02_Cond02_HbO": 9.0,
            "S03_D02_Cond03_HbO": 9.0,
            "S03_D02_Cond04_HbO": 9.0,
            "S03_D04_Cond01_HbO": 9.0,
            "S03_D04_Cond02_HbO": 9.0,
            "S03_D04_Cond03_HbO": 9.0,
            "S03_D04_Cond04_HbO": 9.0,
        },
    ]
    pd.DataFrame(rows).to_csv(path, index=False)


def write_toy_roi_json(path: Path) -> None:
    payload = {
        "L_DMPFC": ["S03_D02", "S03_D04", "S04_D02"],
        "Other": ["S01_D01"],
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def write_toy_exclusions(path: Path) -> None:
    path.write_text(json.dumps(["sub_0004"]), encoding="utf-8")


def test_zero_is_missing_and_roi_mean_uses_available_channels(
    csv_path: Path, roi_json: Path
) -> None:
    df = beta_plot.load_merged_input(str(csv_path))
    roi_map = beta_plot.load_roi_definition(str(roi_json))
    members = beta_plot.build_roi_member_long(
        df,
        roi_name="L_DMPFC",
        chrom="HbO",
        roi_map=roi_map,
        zero_is_missing=True,
    )
    roi_mean = beta_plot.build_roi_mean_long(members, roi_name="L_DMPFC", chrom="HbO")
    sub = roi_mean.loc[
        roi_mean["subject_id"].eq(2) & roi_mean["condition_code"].eq("02")
    ]
    assert len(sub) == 1
    assert abs(sub["beta"].iloc[0] - 3.0) < 1e-9


def test_exclusions_and_complete_case_counts(
    csv_path: Path, roi_json: Path, exclusions_json: Path
) -> None:
    df = beta_plot.load_merged_input(str(csv_path))
    excluded_ids = beta_plot.load_excluded_ids(str(exclusions_json))
    df = beta_plot.apply_subject_exclusions(df, excluded_ids, enabled=True)
    roi_map = beta_plot.load_roi_definition(str(roi_json))

    raw_df, summary_df, n_map = beta_plot.build_plot_dataset(
        df,
        roi_map,
        channel="S04_D02",
        channel_chrom="HbR",
        roi="L_DMPFC",
        roi_chrom="HbO",
        ci_level=0.95,
        zero_is_missing=True,
        decompose_roi=True,
    )

    assert n_map["channel"] == 2, "subject 0003 should drop from channel panel complete-case set"
    assert n_map["roi"] == 3, "ROI panel should retain subjects 0001-0003 after exclusion"
    assert 4 not in raw_df["subject_id"].dropna().astype(int).unique()
    assert set(summary_df["panel"].unique()) == {"channel", "roi", "roi_decomposition"}


def test_run_plotting_creates_outputs(
    csv_path: Path, roi_json: Path, exclusions_json: Path, out_dir: Path
) -> None:
    outputs = beta_plot.run_plotting(
        input_csv=str(csv_path),
        roi_json=str(roi_json),
        exclude_subjects_json=str(exclusions_json),
        channel="S04_D02",
        channel_chrom="HbR",
        roi="L_DMPFC",
        roi_chrom="HbO",
        style="raw_means",
        out_dir=str(out_dir),
        dpi=120,
        ci_level=0.95,
        apply_exclusions_flag=True,
        decompose_roi=True,
        zero_is_missing=True,
    )
    assert outputs["figure_path"].exists()
    assert outputs["data_csv_path"].exists()

    exported = pd.read_csv(outputs["data_csv_path"])
    assert {"raw", "summary"} == set(exported["row_type"].unique())
    assert "roi_decomposition" in set(exported["panel"].unique())


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        csv_path = tmp_path / "toy_merged.csv"
        roi_json = tmp_path / "roi_definition.json"
        exclusions_json = tmp_path / "excluded_subjects.json"
        out_dir = tmp_path / "outputs"

        write_toy_merged_csv(csv_path)
        write_toy_roi_json(roi_json)
        write_toy_exclusions(exclusions_json)

        test_zero_is_missing_and_roi_mean_uses_available_channels(csv_path, roi_json)
        test_exclusions_and_complete_case_counts(csv_path, roi_json, exclusions_json)
        test_run_plotting_creates_outputs(csv_path, roi_json, exclusions_json, out_dir)

    print("[PASS] validate_beta_discrepancy_plot_py")


if __name__ == "__main__":
    main()
