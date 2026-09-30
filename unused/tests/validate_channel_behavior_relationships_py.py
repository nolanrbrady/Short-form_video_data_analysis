"""Validation harness for channel-behavior exploratory screening.

Run:
  python tests/validate_channel_behavior_relationships_py.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analyze_channel_behavior_relationships import (  # noqa: E402
    build_condition_matched_results,
    build_behavior_summary,
    classify_behavior_variables,
    compute_pairwise_associations,
    identify_channel_columns,
    parse_channel_metadata,
    sanitize_channel_values,
    summarize_channel_missingness,
)


def make_fixture_dataframe() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "subject_id": [1, 2, 3, 4, 5, 6],
            "homer_subject": [f"sub_{i:04d}" for i in range(1, 7)],
            "S01_D01_Cond01_HbO": [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            "S01_D01_Cond01_HbR": [0.0, -0.2, -0.4, -0.6, -0.8, -1.0],
            "continuous_behavior": [1, 2, 3, 4, 5, 6],
            "binary_behavior": [0, 0, 0, 1, 1, 1],
        }
    )


def test_identify_channel_columns_and_parse_metadata() -> None:
    df = make_fixture_dataframe()
    channels = identify_channel_columns(df.columns)
    assert channels == ["S01_D01_Cond01_HbO", "S01_D01_Cond01_HbR"]
    metadata = parse_channel_metadata("S01_D01_Cond01_HbO")
    assert metadata["source"] == "01"
    assert metadata["detector"] == "01"
    assert metadata["condition"] == "01"
    assert metadata["chromophore"] == "HbO"
    assert metadata["condition_label"] == "short_form_education"


def test_sanitize_channel_values_only_replaces_channel_zeroes() -> None:
    df = make_fixture_dataframe()
    channels = identify_channel_columns(df.columns)
    sanitized = sanitize_channel_values(df, channels)
    assert math.isnan(sanitized.loc[0, "S01_D01_Cond01_HbO"])
    assert math.isnan(sanitized.loc[0, "S01_D01_Cond01_HbR"])
    assert sanitized.loc[0, "binary_behavior"] == 0


def test_behavior_classification_uses_spearman_for_continuous_and_point_biserial_for_binary() -> None:
    df = sanitize_channel_values(make_fixture_dataframe(), ["S01_D01_Cond01_HbO", "S01_D01_Cond01_HbR"])
    behavior_variables = classify_behavior_variables(
        df,
        channel_columns=["S01_D01_Cond01_HbO", "S01_D01_Cond01_HbR"],
        identifier_columns=["subject_id", "homer_subject"],
    )
    method_map = {item.name: item.method for item in behavior_variables}
    assert method_map == {
        "continuous_behavior": "spearman",
        "binary_behavior": "point_biserial",
    }


def test_pairwise_results_apply_fdr_and_preserve_expected_effect_directions() -> None:
    raw = make_fixture_dataframe()
    channels = identify_channel_columns(raw.columns)
    df = sanitize_channel_values(raw, channels)
    behavior_variables = classify_behavior_variables(
        df,
        channel_columns=channels,
        identifier_columns=["subject_id", "homer_subject"],
    )
    results = compute_pairwise_associations(df, channels, behavior_variables)
    assert set(results["behavior"]) == {"continuous_behavior", "binary_behavior"}
    assert results["p_fdr_bh"].notna().sum() == 4

    continuous_hbo = results.loc[
        (results["behavior"] == "continuous_behavior") & (results["channel"] == "S01_D01_Cond01_HbO")
    ].iloc[0]
    assert continuous_hbo["method"] == "spearman"
    assert continuous_hbo["effect_size"] > 0.99
    assert bool(continuous_hbo["fdr_reject_alpha_0_05"])

    binary_hbr = results.loc[
        (results["behavior"] == "binary_behavior") & (results["channel"] == "S01_D01_Cond01_HbR")
    ].iloc[0]
    assert binary_hbr["method"] == "point_biserial"
    assert binary_hbr["effect_size"] < 0


def test_missingness_summary_counts_zeroes_as_pruned_after_sanitization() -> None:
    raw = make_fixture_dataframe()
    channels = identify_channel_columns(raw.columns)
    sanitized = sanitize_channel_values(raw, channels)
    summary = summarize_channel_missingness(raw, sanitized, channels)
    row = summary.loc[summary["channel"] == "S01_D01_Cond01_HbO"].iloc[0]
    assert row["raw_zero_count"] == 1
    assert row["raw_nan_count"] == 0
    assert row["pruned_or_missing_n"] == 1
    assert row["usable_n"] == 5


def test_behavior_summary_picks_best_channel_per_behavior() -> None:
    raw = make_fixture_dataframe()
    channels = identify_channel_columns(raw.columns)
    df = sanitize_channel_values(raw, channels)
    behavior_variables = classify_behavior_variables(
        df,
        channel_columns=channels,
        identifier_columns=["subject_id", "homer_subject"],
    )
    results = compute_pairwise_associations(df, channels, behavior_variables)
    summary = build_behavior_summary(results)
    continuous_row = summary.loc[summary["behavior"] == "continuous_behavior"].iloc[0]
    assert continuous_row["best_channel"] == "S01_D01_Cond01_HbO"
    assert continuous_row["n_fdr_significant"] >= 1


def test_condition_matched_results_filter_to_defined_behavior_condition_pairs() -> None:
    raw = make_fixture_dataframe()
    channels = identify_channel_columns(raw.columns)
    df = sanitize_channel_values(raw, channels)
    renamed = df.rename(
        columns={
            "continuous_behavior": "diff_short_form_education",
            "binary_behavior": "pd_status",
            "S01_D01_Cond01_HbR": "S01_D01_Cond02_HbR",
        }
    )
    renamed_channels = identify_channel_columns(renamed.columns)
    behavior_variables = classify_behavior_variables(
        renamed,
        channel_columns=renamed_channels,
        identifier_columns=["subject_id", "homer_subject"],
    )
    results = compute_pairwise_associations(renamed, renamed_channels, behavior_variables)
    matched = build_condition_matched_results(results)
    assert set(matched["behavior"]) == {"diff_short_form_education"}
    assert set(matched["condition_label"]) == {"short_form_education"}


def test_run_like_output_files_can_be_written() -> None:
    from analyze_channel_behavior_relationships import run_analysis  # noqa: E402

    raw = make_fixture_dataframe()
    with TemporaryDirectory() as tmpdir:
        input_csv = Path(tmpdir) / "fixture.csv"
        out_dir = Path(tmpdir) / "outputs"
        raw.to_csv(input_csv, index=False)
        outputs = run_analysis(input_csv=input_csv, out_dir=out_dir, alpha=0.05, excluded_columns=[])
        for path in outputs.values():
            assert path.exists(), f"Expected output file to exist: {path}"


def main() -> None:
    test_identify_channel_columns_and_parse_metadata()
    test_sanitize_channel_values_only_replaces_channel_zeroes()
    test_behavior_classification_uses_spearman_for_continuous_and_point_biserial_for_binary()
    test_pairwise_results_apply_fdr_and_preserve_expected_effect_directions()
    test_missingness_summary_counts_zeroes_as_pruned_after_sanitization()
    test_behavior_summary_picks_best_channel_per_behavior()
    test_condition_matched_results_filter_to_defined_behavior_condition_pairs()
    test_run_like_output_files_can_be_written()
    print("[PASS] validate_channel_behavior_relationships_py")


if __name__ == "__main__":
    main()
