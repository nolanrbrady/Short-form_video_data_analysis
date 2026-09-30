"""Real-data reproducibility harness for exported result CSVs.

This validator reruns result-producing scripts into a temporary directory and
compares the generated CSVs against `data/results`. It intentionally never
writes into `data/results`. Its scope is the four primary LMMs, selected pooled
means, and behavioral pairwise correlations; archived exploratory alternatives
are not invoked. Audit-trail rationale: Sandve et al. (2013), see CITATIONS.md.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from pandas.testing import assert_frame_equal


ROOT = Path(__file__).resolve().parents[1]
INPUT_CSV = "data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv"
EXCLUSIONS_JSON = "data/config/excluded_subjects.json"
ROI_JSON = "data/config/roi_definition.json"
BEHAVIOR_PLAN_JSON = "data/config/behavior_pairwise_correlation_plan.json"


def run_command(args: list[str], *, cwd: Path = ROOT) -> None:
    result = subprocess.run(args, cwd=cwd, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        tail = "\n".join((result.stdout + "\n" + result.stderr).splitlines()[-40:])
        raise AssertionError(f"Command failed: {' '.join(args)}\n{tail}")


def assert_csv_equal(
    exported: Path,
    rerun: Path,
    *,
    ignore_columns: tuple[str, ...] = (),
    exact: bool = False,
) -> None:
    if not exported.exists():
        raise AssertionError(f"Missing exported CSV: {exported}")
    if not rerun.exists():
        raise AssertionError(f"Missing rerun CSV: {rerun}")

    left = pd.read_csv(exported)
    right = pd.read_csv(rerun)
    drop = [col for col in ignore_columns if col in left.columns and col in right.columns]
    try:
        assert_frame_equal(
            left.drop(columns=drop),
            right.drop(columns=drop),
            check_exact=exact,
            rtol=1e-12,
            atol=1e-12,
        )
    except AssertionError as exc:
        raise AssertionError(f"CSV mismatch for {exported} vs {rerun}: {exc}") from exc


def run_primary_lmm_reruns(tmp: Path) -> None:
    run_command(
        [
            "Rscript",
            "analyze_format_content_lmm_channelwise.R",
            "--input_csv",
            INPUT_CSV,
            "--exclude_subjects_json",
            EXCLUSIONS_JSON,
            "--out_main_csv",
            str(tmp / "format_content_lmm_main_effects_r.csv"),
            "--out_main_tidy_csv",
            str(tmp / "format_content_lmm_main_effects_tidy_r.csv"),
            "--out_posthoc_csv",
            str(tmp / "format_content_lmm_posthoc_pairwise_r.csv"),
        ]
    )
    run_command(
        [
            "Rscript",
            "analyze_format_content_lmm_roi.R",
            "--input_csv",
            INPUT_CSV,
            "--roi_json",
            ROI_JSON,
            "--exclude_subjects_json",
            EXCLUSIONS_JSON,
            "--out_main_csv",
            str(tmp / "format_content_lmm_roi_main_effects_r.csv"),
            "--out_main_tidy_csv",
            str(tmp / "format_content_lmm_roi_main_effects_tidy_r.csv"),
            "--out_posthoc_csv",
            str(tmp / "format_content_lmm_roi_posthoc_pairwise_r.csv"),
        ]
    )
    run_command(
        [
            "Rscript",
            "analyze_retention_format_content_lmm.R",
            "--input_csv",
            INPUT_CSV,
            "--exclude_subjects_json",
            EXCLUSIONS_JSON,
            "--out_main_csv",
            str(tmp / "retention_format_content_lmm_main_effects_r.csv"),
            "--out_posthoc_csv",
            str(tmp / "retention_format_content_lmm_posthoc_pairwise_r.csv"),
        ]
    )
    run_command(
        [
            "Rscript",
            "analyze_engagement_format_content_lmm.R",
            "--input_csv",
            INPUT_CSV,
            "--exclude_subjects_json",
            EXCLUSIONS_JSON,
            "--out_main_csv",
            str(tmp / "engagement_format_content_lmm_main_effects_r.csv"),
            "--out_posthoc_csv",
            str(tmp / "engagement_format_content_lmm_posthoc_pairwise_r.csv"),
        ]
    )


def run_correlation_reruns(tmp: Path) -> None:
    pooled = tmp / "pooled_mean_correlations"
    behavior = tmp / "behavior_pairwise_correlations"

    run_command(
        [
            "Rscript",
            "analyze_pooled_mean_correlations.R",
            "--input_csv",
            INPUT_CSV,
            "--roi_json",
            ROI_JSON,
            "--channel_results_csv",
            "data/results/format_content_lmm_main_effects_tidy_r.csv",
            "--roi_results_csv",
            "data/results/format_content_lmm_roi_main_effects_tidy_r.csv",
            "--exclude_subjects_json",
            EXCLUSIONS_JSON,
            "--out_dir",
            str(pooled),
        ]
    )
    run_command(
        [
            "Rscript",
            "analyze_behavior_pairwise_correlations.R",
            "--input_csv",
            INPUT_CSV,
            "--analysis_plan_json",
            BEHAVIOR_PLAN_JSON,
            "--exclude_subjects_json",
            EXCLUSIONS_JSON,
            "--out_dir",
            str(behavior),
        ]
    )


def compare_primary_lmm_outputs(tmp: Path) -> None:
    for name in [
        "format_content_lmm_main_effects_r.csv",
        "format_content_lmm_main_effects_tidy_r.csv",
        "format_content_lmm_posthoc_pairwise_r.csv",
        "format_content_lmm_roi_main_effects_r.csv",
        "format_content_lmm_roi_main_effects_tidy_r.csv",
        "format_content_lmm_roi_posthoc_pairwise_r.csv",
        "retention_format_content_lmm_main_effects_r.csv",
        "retention_format_content_lmm_posthoc_pairwise_r.csv",
        "engagement_format_content_lmm_main_effects_r.csv",
        "engagement_format_content_lmm_posthoc_pairwise_r.csv",
    ]:
        assert_csv_equal(ROOT / "data/results" / name, tmp / name)


def compare_correlation_outputs(tmp: Path) -> None:
    comparisons = [
        (
            "pooled_mean_correlations/selected_pooled_mean_targets_r.csv",
            "pooled_mean_correlations/selected_pooled_mean_targets_r.csv",
            (),
        ),
        (
            "pooled_mean_correlations/subject_level_pooled_mean_pairs_r.csv",
            "pooled_mean_correlations/subject_level_pooled_mean_pairs_r.csv",
            (),
        ),
        (
            "pooled_mean_correlations/pooled_mean_correlations_r.csv",
            "pooled_mean_correlations/pooled_mean_correlations_r.csv",
            ("plot_file",),
        ),
        (
            "behavior_pairwise_correlations/behavior_pairwise_correlations_r.csv",
            "behavior_pairwise_correlations/behavior_pairwise_correlations_r.csv",
            (),
        ),
        (
            "behavior_pairwise_correlations/behavior_pairwise_correlations_fdr_r.csv",
            "behavior_pairwise_correlations/behavior_pairwise_correlations_fdr_r.csv",
            (),
        ),
    ]
    for exported_rel, rerun_rel, ignored in comparisons:
        assert_csv_equal(ROOT / "data/results" / exported_rel, tmp / rerun_rel, ignore_columns=ignored)


def assert_neural_kr_tidy_invariants() -> None:
    for path in [
        ROOT / "data/results/format_content_lmm_main_effects_tidy_r.csv",
        ROOT / "data/results/format_content_lmm_roi_main_effects_tidy_r.csv",
    ]:
        df = pd.read_csv(path)
        audited = df[df[["estimate", "se", "t"]].notna().all(axis=1) & df["se"].abs().gt(0)]
        if audited.empty:
            raise AssertionError(f"No finite rows to audit in {path}")
        max_abs = ((audited["estimate"] / audited["se"]) - audited["t"]).abs().max()
        if max_abs >= 1e-10:
            raise AssertionError(f"KR SE/t invariant failed for {path}: max_abs={max_abs}")


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="rise-real-results-") as tmp_name:
        tmp = Path(tmp_name)
        run_primary_lmm_reruns(tmp)
        run_correlation_reruns(tmp)
        compare_primary_lmm_outputs(tmp)
        compare_correlation_outputs(tmp)
        assert_neural_kr_tidy_invariants()
    print("[PASS] validate_real_result_reproducibility_py")


if __name__ == "__main__":
    main()
