"""Research-correctness tests for Pearson covariate correlation diagnostics."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import covariate_correlation_analysis as covariate_analysis


@pytest.fixture
def no_heatmap(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Record heatmap requests without testing matplotlib rendering."""

    paths: list[str] = []

    def record_heatmap(corr: pd.DataFrame, outpath: Path, title: str) -> None:
        paths.append(Path(outpath).name)

    monkeypatch.setattr(covariate_analysis, "save_heatmap", record_heatmap)
    return paths


def write_csv(path: Path, data: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(path, index=False)


def write_exclusions(path: Path, values: list[str]) -> None:
    path.write_text(json.dumps(values), encoding="utf-8")


def assert_columns(frame: pd.DataFrame, expected: list[str]) -> None:
    assert list(frame.columns) == expected


def test_scipy_is_required_at_module_import() -> None:
    block_scipy_import = """
import builtins
real_import = builtins.__import__
def blocked_import(name, globals=None, locals=None, fromlist=(), level=0):
    if name == "scipy" or name.startswith("scipy."):
        raise ModuleNotFoundError("blocked scipy for test")
    return real_import(name, globals, locals, fromlist, level)
builtins.__import__ = blocked_import
import covariate_correlation_analysis
"""

    result = subprocess.run(
        [sys.executable, "-c", block_scipy_import],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode != 0
    assert "ModuleNotFoundError" in result.stderr


def test_subject_exclusion_normalizes_ids_filters_rows_preserves_order_reports_audit() -> None:
    source = pd.DataFrame(
        {
            "subject_id": ["sub_0001", "2", "sub_0003"],
            "age": [10, 20, 30],
        }
    )
    source_before = source.copy(deep=True)

    filtered, report = covariate_analysis.apply_subject_exclusions_from_ids(
        source,
        subject_column="subject_id",
        excluded_ids=[2, 4],
        context_label="unit",
    )

    assert filtered["subject_id"].tolist() == ["sub_0001", "sub_0003"]
    assert filtered["age"].tolist() == [10, 30]
    assert "_subject_id_norm" not in filtered.columns
    assert report == {
        "listed": 2,
        "matched": 1,
        "removed": 1,
        "remaining": 2,
        "matched_ids": [2],
        "missing_ids": [4],
    }
    assert_frame_equal(source, source_before)


def test_subject_exclusion_requires_subject_ids_when_manifest_is_nonempty() -> None:
    source = pd.DataFrame({"age": [10, 20]})

    unchanged, report = covariate_analysis.apply_subject_exclusions_from_ids(
        source,
        subject_column="subject_id",
        excluded_ids=[],
        context_label="unit",
    )
    assert_frame_equal(unchanged, source)
    assert report["listed"] == 0
    assert report["remaining"] == 2

    with pytest.raises(ValueError, match="cannot apply nonempty subject exclusions"):
        covariate_analysis.apply_subject_exclusions_from_ids(
            source,
            subject_column="subject_id",
            excluded_ids=[1],
            context_label="unit",
        )


def test_subject_count_limit_accepts_boundary_and_rejects_overage() -> None:
    covariate_analysis.assert_subject_count_limit(
        report={"remaining": 48},
        context_label="unit",
        max_subjects=48,
    )

    with pytest.raises(ValueError, match="exceeds configured maximum of 48"):
        covariate_analysis.assert_subject_count_limit(
            report={"remaining": 49},
            context_label="unit",
            max_subjects=48,
        )


def test_recruitment_order_proxy_normalizes_ids_and_does_not_mutate_input() -> None:
    source = pd.DataFrame({"subject_id": ["sub_0007", "8"], "age": [20, 21]})
    source_before = source.copy(deep=True)

    with_proxy = covariate_analysis.add_recruitment_order_proxy(
        source,
        subject_column="subject_id",
    )

    assert with_proxy["recruitment_order_proxy"].tolist() == [7.0, 8.0]
    assert_frame_equal(source, source_before)


def test_recruitment_order_proxy_fails_when_subject_column_absent() -> None:
    with pytest.raises(ValueError, match="subject column 'subject_id' is absent"):
        covariate_analysis.add_recruitment_order_proxy(
            pd.DataFrame({"age": [20, 21]}),
            subject_column="subject_id",
        )


def test_correlation_input_columns_coerces_numeric_and_removes_identifiers_without_imputing() -> None:
    source = pd.DataFrame(
        {
            "subject_id": [1, 2],
            "homer_subject": ["sub_0001", "sub_0002"],
            "age": ["20", "21"],
            "notes": ["alpha", "beta"],
            "true_missing": [np.nan, 4.0],
        }
    )

    corr_inputs = covariate_analysis.correlation_input_columns(source, "subject_id")

    assert_columns(corr_inputs, ["age", "notes", "true_missing"])
    assert corr_inputs["age"].tolist() == [20, 21]
    assert corr_inputs["notes"].isna().all()
    assert pd.isna(corr_inputs.loc[0, "true_missing"])
    assert corr_inputs.loc[1, "true_missing"] == 4.0


def test_prepare_correlation_inputs_applies_exclusion_then_subject_limit_then_optional_proxy() -> None:
    source = pd.DataFrame(
        {
            "subject_id": [1, 2, 3],
            "homer_subject": ["sub_0001", "sub_0002", "sub_0003"],
            "age": [10, 20, 30],
            "behavior_score": [2, 4, 6],
            "notes": ["a", "b", "c"],
        }
    )
    source_before = source.copy(deep=True)

    prepared = covariate_analysis.prepare_correlation_inputs(
        source,
        subject_column="subject_id",
        excluded_ids=[3],
        context_label="unit",
        include_subject_id_correlation=True,
        max_subjects=2,
    )

    assert_columns(
        prepared.correlation_inputs,
        ["age", "behavior_score", "notes", "recruitment_order_proxy"],
    )
    assert prepared.correlation_inputs["age"].tolist() == [10, 20]
    assert prepared.correlation_inputs["behavior_score"].tolist() == [2, 4]
    assert prepared.correlation_inputs["notes"].isna().all()
    assert prepared.correlation_inputs["recruitment_order_proxy"].tolist() == [1.0, 2.0]
    assert prepared.exclusion_report["matched_ids"] == [3]
    assert prepared.exclusion_report["remaining"] == 2
    assert_frame_equal(source, source_before)


def test_prepare_correlation_inputs_fails_when_post_exclusion_subject_count_exceeds_limit() -> None:
    with pytest.raises(ValueError, match="exceeds configured maximum of 2"):
        covariate_analysis.prepare_correlation_inputs(
            pd.DataFrame({"subject_id": [1, 2, 3], "age": [10, 20, 30]}),
            subject_column="subject_id",
            excluded_ids=[],
            context_label="unit",
            include_subject_id_correlation=False,
            max_subjects=2,
        )


def test_build_pearson_outputs_matches_complete_expected_matrix() -> None:
    inputs = pd.DataFrame(
        {
            "x": [1.0, 2.0, 3.0, 4.0],
            "same": [1.0, 2.0, 3.0, 4.0],
            "opposite": [4.0, 3.0, 2.0, 1.0],
            "constant": [5.0, 5.0, 5.0, 5.0],
        }
    )

    outputs = covariate_analysis.build_pearson_outputs(inputs)
    expected_correlations = pd.DataFrame(
        [
            [1.0, 1.0, -1.0, np.nan],
            [1.0, 1.0, -1.0, np.nan],
            [-1.0, -1.0, 1.0, np.nan],
            [np.nan, np.nan, np.nan, np.nan],
        ],
        index=["x", "same", "opposite", "constant"],
        columns=["x", "same", "opposite", "constant"],
    )

    assert_frame_equal(outputs.correlations, expected_correlations)
    assert outputs.pvalues.index.tolist() == expected_correlations.index.tolist()
    assert outputs.pvalues.columns.tolist() == expected_correlations.columns.tolist()
    assert outputs.pvalues.loc["x", "same"] == pytest.approx(0.0)
    assert outputs.pvalues.loc["x", "opposite"] == pytest.approx(0.0)
    assert pd.isna(outputs.pvalues.loc["x", "constant"])


def test_pearson_pvalues_are_symmetric_pairwise_complete_and_do_not_listwise_delete() -> None:
    inputs = pd.DataFrame(
        {
            "x": [1.0, 2.0, 3.0, 4.0],
            "y": [2.0, 4.0, None, 8.0],
            "constant": [1.0, 1.0, 1.0, 1.0],
        }
    )

    pvalues = covariate_analysis.pearson_pvalues(inputs)

    assert pvalues.index.tolist() == ["x", "y", "constant"]
    assert pvalues.columns.tolist() == ["x", "y", "constant"]
    assert pvalues.loc["x", "y"] == pvalues.loc["y", "x"]
    assert pvalues.loc["x", "y"] < 1e-7
    assert pd.isna(pvalues.loc["x", "constant"])


def test_pearson_pvalues_return_nan_when_pair_has_fewer_than_three_complete_values() -> None:
    inputs = pd.DataFrame(
        {
            "x": [1.0, 2.0, np.nan, np.nan],
            "y": [2.0, 4.0, 6.0, np.nan],
        }
    )

    pvalues = covariate_analysis.pearson_pvalues(inputs)

    assert pd.isna(pvalues.loc["x", "y"])
    assert pd.isna(pvalues.loc["y", "x"])


def test_obsolete_output_paths_are_scoped_to_known_covariate_files() -> None:
    out = Path("/tmp/out")
    paths = covariate_analysis.obsolete_output_paths(out)

    assert all(path.parent == out for path in paths)
    assert out / "covariate_correlation_analysis.csv" in paths
    assert out / "covariate_heatmap_spearman.png" in paths
    assert out / "unrelated.csv" not in paths


def test_combined_runner_applies_exclusions_before_pearson(tmp_path: Path, no_heatmap: list[str]) -> None:
    input_csv = tmp_path / "combined.csv"
    output_dir = tmp_path / "out"
    exclusions_json = tmp_path / "excluded_subjects.json"
    write_csv(
        input_csv,
        pd.DataFrame(
            {
                "subject_id": [1, 2, 3, 4],
                "age": [10, 20, 30, 40],
                "behavior_score": [10, 20, 30, -999],
            }
        ),
    )
    write_exclusions(exclusions_json, ["sub_0004"])

    covariate_analysis.run_combined_preset(
        tmp_path,
        input_path=input_csv,
        out_dir=output_dir,
        excluded_subjects_json=exclusions_json,
        subject_column="subject_id",
    )

    correlations = pd.read_csv(output_dir / "covariate_correlation_analysis_pearson.csv")
    assert correlations.loc[0, "behavior_score"] == pytest.approx(1.0)
    assert correlations.loc[1, "behavior_score"] == pytest.approx(1.0)


def test_combined_runner_writes_only_pearson_outputs(tmp_path: Path, no_heatmap: list[str]) -> None:
    input_csv = tmp_path / "combined.csv"
    output_dir = tmp_path / "out"
    exclusions_json = tmp_path / "excluded_subjects.json"
    write_csv(
        input_csv,
        pd.DataFrame(
            {
                "subject_id": [1, 2, 3, 4],
                "age": [1, 2, 3, 4],
                "behavior_score": [1, 2, 3, 100],
            }
        ),
    )
    write_exclusions(exclusions_json, [])

    covariate_analysis.run_combined_preset(
        tmp_path,
        input_path=input_csv,
        out_dir=output_dir,
        excluded_subjects_json=exclusions_json,
        subject_column="subject_id",
    )

    pearson = pd.read_csv(output_dir / "covariate_correlation_analysis_pearson.csv")
    assert pearson.loc[0, "behavior_score"] < 1.0
    assert (output_dir / "covariate_correlation_analysis_pearson_pvalues.csv").exists()
    assert not (output_dir / "covariate_correlation_analysis.csv").exists()
    assert not (output_dir / "covariate_correlation_analysis_spearman.csv").exists()
    assert not (output_dir / "covariate_correlation_heatmap.png").exists()
    assert not (output_dir / "covariate_correlation_heatmap_spearman.png").exists()
    assert no_heatmap == ["covariate_correlation_heatmap_pearson.png"]


def test_covariate_runner_writes_only_pearson_outputs(tmp_path: Path, no_heatmap: list[str]) -> None:
    input_csv = tmp_path / "covariates.csv"
    output_dir = tmp_path / "out"
    exclusions_json = tmp_path / "excluded_subjects.json"
    write_csv(
        input_csv,
        pd.DataFrame(
            {
                "subject_id": [1, 2, 3, 4],
                "age": [1, 2, 3, 4],
                "behavior_score": [1, 2, 3, 100],
            }
        ),
    )
    write_exclusions(exclusions_json, [])

    covariate_analysis.run_covariate_preset(
        tmp_path,
        input_path=input_csv,
        out_dir=output_dir,
        excluded_subjects_json=exclusions_json,
        subject_column="subject_id",
    )

    pearson = pd.read_csv(output_dir / "covariate_correlations_pearson.csv", index_col=0)
    assert pearson.loc["age", "behavior_score"] < 1.0
    assert (output_dir / "covariate_correlations_pearson_pvalues.csv").exists()
    assert not (output_dir / "covariate_correlations_pvalues.csv").exists()
    assert not (output_dir / "covariate_correlations_spearman.csv").exists()
    assert not (output_dir / "covariate_correlations_spearman_pvalues.csv").exists()
    assert not (output_dir / "covariate_heatmap.png").exists()
    assert not (output_dir / "covariate_heatmap_spearman.png").exists()
    assert no_heatmap == ["covariate_heatmap_pearson.png"]


def test_combined_runner_can_include_recruitment_order_proxy(
    tmp_path: Path,
    no_heatmap: list[str],
) -> None:
    input_csv = tmp_path / "combined.csv"
    output_dir = tmp_path / "out"
    exclusions_json = tmp_path / "excluded_subjects.json"
    write_csv(
        input_csv,
        pd.DataFrame(
            {
                "subject_id": [1, 2, 3],
                "age": [10, 20, 30],
                "behavior_score": [3, 2, 1],
            }
        ),
    )
    write_exclusions(exclusions_json, [])

    covariate_analysis.run_combined_preset(
        tmp_path,
        input_path=input_csv,
        out_dir=output_dir,
        excluded_subjects_json=exclusions_json,
        subject_column="subject_id",
        include_subject_id_correlation=True,
    )

    correlations = pd.read_csv(output_dir / "covariate_correlation_analysis_pearson.csv")
    assert "recruitment_order_proxy" in correlations.columns
    assert correlations.loc[0, "recruitment_order_proxy"] == pytest.approx(1.0)
    assert correlations.loc[1, "recruitment_order_proxy"] == pytest.approx(-1.0)


def test_combined_runner_fails_when_more_than_48_subjects_remain(
    tmp_path: Path,
    no_heatmap: list[str],
) -> None:
    input_csv = tmp_path / "combined.csv"
    output_dir = tmp_path / "out"
    exclusions_json = tmp_path / "excluded_subjects.json"
    write_csv(
        input_csv,
        pd.DataFrame(
            {
                "subject_id": list(range(1, 52)),
                "age": list(range(18, 69)),
                "behavior_score": list(range(51)),
            }
        ),
    )
    write_exclusions(exclusions_json, [])

    with pytest.raises(ValueError, match="exceeds configured maximum of 48"):
        covariate_analysis.run_combined_preset(
            tmp_path,
            input_path=input_csv,
            out_dir=output_dir,
            excluded_subjects_json=exclusions_json,
            subject_column="subject_id",
            max_subjects=48,
        )


def test_removed_spearman_outputs_are_deleted_from_existing_output_dir(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    stale_names = [path.name for path in covariate_analysis.obsolete_output_paths(tmp_path / "out")]
    input_csv = tmp_path / "combined.csv"
    output_dir = tmp_path / "out"
    exclusions_json = tmp_path / "excluded_subjects.json"
    write_csv(
        input_csv,
        pd.DataFrame(
            {
                "subject_id": [1, 2, 3],
                "age": [1, 2, 3],
                "behavior_score": [1, 2, 10],
            }
        ),
    )
    write_exclusions(exclusions_json, [])
    output_dir.mkdir()
    for name in stale_names:
        (output_dir / name).write_text("stale", encoding="utf-8")

    monkeypatch.setattr(
        covariate_analysis,
        "save_heatmap",
        lambda corr, outpath, title: Path(outpath).write_text("fresh", encoding="utf-8"),
    )

    covariate_analysis.run_combined_preset(
        tmp_path,
        input_path=input_csv,
        out_dir=output_dir,
        excluded_subjects_json=exclusions_json,
        subject_column="subject_id",
    )

    assert [name for name in stale_names if (output_dir / name).exists()] == []
