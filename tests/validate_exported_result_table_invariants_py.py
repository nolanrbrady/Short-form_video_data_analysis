"""Invariant checks for exported result CSVs.

This validator reads `data/results` only. It does not rerun analyses and does
not write any files. The goal is to catch stale or internally inconsistent
publication-facing tables.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "data/results"
ALPHA = 0.05


def assert_true(condition: bool, message: str) -> None:
    if not bool(condition):
        raise AssertionError(message)


def read_csv(relative_path: str) -> pd.DataFrame:
    path = RESULTS / relative_path
    assert_true(path.exists(), f"Missing exported result table: {path}")
    return pd.read_csv(path)


def bh_qvalues(p_values: pd.Series) -> pd.Series:
    p = pd.to_numeric(p_values, errors="coerce")
    q = pd.Series(np.nan, index=p.index, dtype=float)
    finite = p[np.isfinite(p)]
    if finite.empty:
        return q
    ordered = finite.sort_values()
    m = len(ordered)
    adjusted = []
    prev = 1.0
    for rank, (_, value) in reversed(list(enumerate(ordered.items(), start=1))):
        prev = min(prev, (m / rank) * float(value))
        adjusted.append(prev)
    adjusted = list(reversed(adjusted))
    q.loc[ordered.index] = np.clip(adjusted, 0.0, 1.0)
    return q


def holm_qvalues(p_values: pd.Series) -> pd.Series:
    p = pd.to_numeric(p_values, errors="coerce")
    q = pd.Series(np.nan, index=p.index, dtype=float)
    finite = p[np.isfinite(p)]
    if finite.empty:
        return q
    ordered = finite.sort_values()
    m = len(ordered)
    raw = pd.Series(index=ordered.index, dtype=float)
    for rank, (idx, value) in enumerate(ordered.items(), start=1):
        raw.loc[idx] = (m - rank + 1) * float(value)
    running = 0.0
    for idx in ordered.index:
        running = max(running, raw.loc[idx])
        q.loc[idx] = min(running, 1.0)
    return q


def assert_probability_columns(df: pd.DataFrame, table_name: str, columns: list[str]) -> None:
    for col in columns:
        if col not in df.columns:
            continue
        values = pd.to_numeric(df[col], errors="coerce").dropna()
        assert_true(((values >= 0.0) & (values <= 1.0)).all(), f"{table_name}.{col} has values outside [0, 1]")


def assert_ci_contains_estimate(df: pd.DataFrame, table_name: str) -> None:
    rows = df[["estimate", "ci95_low", "ci95_high"]].dropna()
    assert_true(
        ((rows["ci95_low"] <= rows["estimate"]) & (rows["estimate"] <= rows["ci95_high"])).all(),
        f"{table_name} has CI bounds that do not contain estimates",
    )


def assert_estimate_over_se_reconstructs_t(df: pd.DataFrame, table_name: str) -> None:
    rows = df[["estimate", "se", "t"]].dropna()
    rows = rows[rows["se"].abs() > 0]
    assert_true(not rows.empty, f"{table_name} has no finite rows for estimate/se/t audit")
    max_abs = ((rows["estimate"] / rows["se"]) - rows["t"]).abs().max()
    assert_true(max_abs < 1e-10, f"{table_name} estimate/se does not reconstruct t; max_abs={max_abs}")


def assert_neural_tidy_table(df: pd.DataFrame, table_name: str, unit_col: str) -> None:
    required = {
        unit_col,
        "chrom",
        "effect",
        "n_subjects",
        "n_obs",
        "converged",
        "singular_fit",
        "estimate",
        "se",
        "df",
        "t",
        "ci95_low",
        "ci95_high",
        "p_unc",
        "p_fdr",
    }
    assert_true(required.issubset(df.columns), f"{table_name} missing required columns")
    assert_true(set(df["effect"]) == {"format", "content", "interaction"}, f"{table_name} has unexpected effects")
    assert_true((df["n_obs"] == 4 * df["n_subjects"]).all(), f"{table_name} violates n_obs == 4*n_subjects")
    assert_true(df["converged"].isin([True, False]).all(), f"{table_name}.converged is not boolean-like")
    assert_true((df["se"] > 0).all(), f"{table_name} has non-positive SE")
    assert_probability_columns(df, table_name, ["p_unc", "p_fdr"])
    assert_ci_contains_estimate(df, table_name)
    assert_estimate_over_se_reconstructs_t(df, table_name)
    for (chrom, effect), sub in df.groupby(["chrom", "effect"]):
        expected = bh_qvalues(sub["p_unc"])
        observed = pd.to_numeric(sub["p_fdr"], errors="coerce")
        max_abs = (observed - expected).abs().max()
        assert_true(max_abs < 1e-12, f"{table_name} BH-FDR mismatch for {chrom}/{effect}: {max_abs}")


def assert_neural_posthoc_gating(
    tidy: pd.DataFrame,
    posthoc: pd.DataFrame,
    *,
    unit_col: str,
    table_name: str,
) -> None:
    sig_units = tidy.loc[
        tidy["effect"].eq("interaction") & tidy["p_fdr"].lt(ALPHA),
        [unit_col, "chrom"],
    ].drop_duplicates()
    expected_rows = 6 * len(sig_units)
    assert_true(len(posthoc) == expected_rows, f"{table_name} posthoc row count does not match interaction gate")
    if expected_rows == 0:
        return
    observed_units = posthoc[[unit_col, "chrom"]].drop_duplicates()
    merged = observed_units.merge(sig_units, on=[unit_col, "chrom"], how="outer", indicator=True)
    assert_true((merged["_merge"] == "both").all(), f"{table_name} posthoc units do not match FDR interaction gate")
    assert_probability_columns(posthoc, f"{table_name} posthoc", ["p_unc"])
    assert_true((posthoc["se"] > 0).all(), f"{table_name} posthoc has non-positive SE")


def assert_behavior_lmm_table(df: pd.DataFrame, table_name: str) -> None:
    assert_true(set(df["effect"]) == {"length", "content", "interaction"}, f"{table_name} has unexpected effects")
    assert_true((df["n_obs"] == 4 * df["n_subjects"]).all(), f"{table_name} violates n_obs == 4*n_subjects")
    assert_true((df["se"] > 0).all(), f"{table_name} has non-positive SE")
    assert_probability_columns(df, table_name, ["p_unc", "p_fdr", "eta2_p"])
    assert_ci_contains_estimate(df, table_name)
    assert_estimate_over_se_reconstructs_t(df, table_name)
    expected = holm_qvalues(df["p_unc"])
    max_abs = (pd.to_numeric(df["p_fdr"], errors="coerce") - expected).abs().max()
    assert_true(max_abs < 1e-12, f"{table_name} Holm correction mismatch: {max_abs}")


def assert_behavior_posthoc_gating(main_df: pd.DataFrame, posthoc: pd.DataFrame, table_name: str) -> None:
    interaction_p = float(main_df.loc[main_df["effect"].eq("interaction"), "p_fdr"].iloc[0])
    expected_rows = 6 if interaction_p < ALPHA else 0
    assert_true(len(posthoc) == expected_rows, f"{table_name} posthoc row count does not match interaction gate")
    if expected_rows > 0:
        assert_probability_columns(posthoc, f"{table_name} posthoc", ["p_unc", "eta2_p"])
        assert_true((posthoc["se"] > 0).all(), f"{table_name} posthoc has non-positive SE")


def assert_family_adjustment(df: pd.DataFrame, table_name: str) -> None:
    tested = df[df["analysis_status"].eq("tested") & np.isfinite(pd.to_numeric(df["p_unc"], errors="coerce"))]
    assert_probability_columns(df, table_name, ["p_unc", "p_fdr"])
    if tested.empty:
        return
    for family_id, sub in tested.groupby("family_id"):
        method = str(sub["family_adjust_method"].iloc[0]).lower()
        if method in {"bh", "fdr"}:
            expected = bh_qvalues(sub["p_unc"])
        else:
            raise AssertionError(f"{table_name} has unsupported exported family method for audit: {method}")
        observed = pd.to_numeric(sub["p_fdr"], errors="coerce")
        max_abs = (observed - expected).abs().max()
        assert_true(max_abs < 1e-12, f"{table_name} family {family_id} FDR mismatch: {max_abs}")
        assert_true((sub["family_n_tested"] == len(sub)).all(), f"{table_name} family_n_tested mismatch for {family_id}")


def assert_correlation_table(df: pd.DataFrame, table_name: str, estimate_col: str) -> None:
    tested = df[df["analysis_status"].eq("tested")]
    assert_true(not tested.empty, f"{table_name} has no tested rows")
    assert_probability_columns(df, table_name, ["p_unc", "p_fdr"])
    est = pd.to_numeric(tested[estimate_col], errors="coerce")
    assert_true(((est >= -1.0) & (est <= 1.0)).all(), f"{table_name}.{estimate_col} outside [-1, 1]")
    if "r_squared" in tested.columns:
        r2 = pd.to_numeric(tested["r_squared"], errors="coerce").dropna()
        assert_true(((r2 >= 0.0) & (r2 <= 1.0)).all(), f"{table_name}.r_squared outside [0, 1]")
    if "ci95_low" in tested.columns and "ci95_high" in tested.columns:
        ci = tested[["ci95_low", "ci95_high"]].dropna()
        assert_true((ci["ci95_low"] <= ci["ci95_high"]).all(), f"{table_name} has inverted CI bounds")
    if "family_id" in df.columns:
        assert_family_adjustment(df, table_name)


def main() -> None:
    channel_tidy = read_csv("format_content_lmm_main_effects_tidy_r.csv")
    channel_posthoc = read_csv("format_content_lmm_posthoc_pairwise_r.csv")
    assert_neural_tidy_table(channel_tidy, "channelwise neural tidy", "channel")
    assert_neural_posthoc_gating(channel_tidy, channel_posthoc, unit_col="channel", table_name="channelwise neural")

    roi_tidy = read_csv("format_content_lmm_roi_main_effects_tidy_r.csv")
    roi_posthoc = read_csv("format_content_lmm_roi_posthoc_pairwise_r.csv")
    assert_neural_tidy_table(roi_tidy, "ROI neural tidy", "roi")
    assert_neural_posthoc_gating(roi_tidy, roi_posthoc, unit_col="roi", table_name="ROI neural")

    retention = read_csv("retention_format_content_lmm_main_effects_r.csv")
    retention_posthoc = read_csv("retention_format_content_lmm_posthoc_pairwise_r.csv")
    assert_behavior_lmm_table(retention, "retention LMM")
    assert_behavior_posthoc_gating(retention, retention_posthoc, "retention LMM")

    engagement = read_csv("engagement_format_content_lmm_main_effects_r.csv")
    engagement_posthoc = read_csv("engagement_format_content_lmm_posthoc_pairwise_r.csv")
    assert_behavior_lmm_table(engagement, "engagement LMM")
    assert_behavior_posthoc_gating(engagement, engagement_posthoc, "engagement LMM")

    assert_correlation_table(
        read_csv("correlational_relationships/pairwise_correlations_r.csv"),
        "condition-specific correlation follow-up",
        "pearson_r",
    )
    assert_correlation_table(
        read_csv("correlational_relationships_roi_means/pairwise_correlations_r.csv"),
        "ROI-mean correlation follow-up",
        "association_estimate",
    )
    assert_correlation_table(
        read_csv("pooled_mean_correlations/pooled_mean_correlations_r.csv"),
        "pooled-mean correlation follow-up",
        "association_estimate",
    )
    assert_correlation_table(
        read_csv("behavior_pairwise_correlations/behavior_pairwise_correlations_r.csv"),
        "behavior pairwise correlations",
        "pearson_r",
    )
    print("[PASS] validate_exported_result_table_invariants_py")


if __name__ == "__main__":
    main()
