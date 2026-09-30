"""Scientific-integrity guards for the active/archived workflow boundary.

Traceable workflow separation follows Sandve et al. (2013),
https://doi.org/10.1371/journal.pcbi.1003285; see CITATIONS.md. These tests check
file layout and verification scope, not change any statistical analysis rules.
"""

from __future__ import annotations

import ast
import csv
import hashlib
import importlib.util
import re
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
ACTIVE_RESULT_SCRIPTS = {
    "pipeline_preprocess_merge.sh",
    "analyze_format_content_lmm_channelwise.R",
    "analyze_format_content_lmm_roi.R",
    "analyze_retention_format_content_lmm.R",
    "analyze_engagement_format_content_lmm.R",
    "analyze_pooled_mean_correlations.R",
    "analyze_behavior_pairwise_correlations.R",
    "create_demographics_table.py",
    "plot_behavior_score_distributions.R",
    "plot_significant_beta_value_distribution.R",
}
ARCHIVED_SCRIPTS = {
    "analyze_sfv_demographics.py",
    "analyze_channel_behavior_relationships.py",
    "analyze_correlational_relationships.R",
    "covariate_correlation_analysis.py",
}


def load_reproducibility_validator():
    """Load verification helpers without invoking costly real-data reruns."""
    spec = importlib.util.spec_from_file_location(
        "current_reproducibility_validator",
        ROOT / "tests/validate_real_result_reproducibility_py.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_archived_files_preserve_bytes_and_original_locations_are_vacant():
    """Every archived input/script retains its pre-move checksum.

    This guards against accidental edits, overwritten destinations, and duplicate
    root entry points that could make the intended analysis version ambiguous.
    """
    with (ROOT / "unused/move_manifest.csv").open(newline="") as stream:
        moves = list(csv.DictReader(stream))
    assert moves
    assert len({row["original_path"] for row in moves}) == len(moves)
    for row in moves:
        original = ROOT / row["original_path"]
        archived = ROOT / row["archived_path"]
        assert not original.exists(), row["original_path"]
        assert archived.is_file(), row["archived_path"]
        assert hashlib.sha256(archived.read_bytes()).hexdigest() == row["sha256_before"]


def test_current_entry_points_helpers_and_optional_qc_remain_available():
    """Preserve production entry points and user-retained QC/robustness tools.

    Helpers are explicitly checked because apparent root clutter can hide a
    dependency used during FIR reconstruction, exclusions, or model inference.
    """
    required = ACTIVE_RESULT_SCRIPTS | {
        "process_engagement.py", "process_sociodemographic.py", "generate_combined_data.py",
        "collapse_homer_fir_to_auc.py", "homer_fir.py", "validate_homer_fir_auc_conversion.py",
        "mask_homer_auc_between_subject_outliers.py", "certify_preprocess_merge_integrity.py",
        "merge_homer3_betas_with_combined_data.R", "demographic/process_recall_assessment.py",
        "r_subject_exclusions.R", "r_lmm_convergence_helpers.R",
        "r_emmeans_posthoc_helpers.R", "r_figure_style.R",
        "audit_check.py", "plot_fir_betas_subjects.py", "analyze_retention_sensitivity.R",
        "check_citations_links.py", "PREPROCESS_NOTES.md", "sfv_data_description.md",
    }
    for relative in required:
        assert (ROOT / relative).is_file(), relative
    assert not any((ROOT / name).exists() for name in ARCHIVED_SCRIPTS)


def test_root_python_imports_do_not_depend_on_archived_modules():
    """Import dependencies of current Python scripts must not point to archives.

    This guards against moving an apparent standalone script that is actually
    imported by preprocessing or another current producer. Comments are ignored.
    """
    archived_modules = {Path(name).stem for name in ARCHIVED_SCRIPTS if name.endswith(".py")}
    for path in ROOT.glob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            modules = []
            if isinstance(node, ast.Import):
                modules = [alias.name.split(".")[0] for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                modules = [node.module.split(".")[0]]
            assert not archived_modules.intersection(modules), path.name


def test_real_reproduction_covers_exactly_the_current_six_analyses(tmp_path):
    """Current verification must rerun all primary analyses and both correlations.

    Expected commands are recorded rather than executed. This proves archived
    screens cannot be invoked and current pooled/behavioral results were not
    accidentally removed from the no-overwrite reproducibility check.
    """
    module = load_reproducibility_validator()
    with patch.object(module, "run_command") as calls:
        module.run_primary_lmm_reruns(tmp_path)
        module.run_correlation_reruns(tmp_path)
    names = [call.args[0][1] for call in calls.call_args_list]
    assert names == [
        "analyze_format_content_lmm_channelwise.R", "analyze_format_content_lmm_roi.R",
        "analyze_retention_format_content_lmm.R", "analyze_engagement_format_content_lmm.R",
        "analyze_pooled_mean_correlations.R", "analyze_behavior_pairwise_correlations.R",
    ]
    for call in calls.call_args_list:
        args = call.args[0]
        for index, value in enumerate(args):
            if value in {"--out_main_csv", "--out_main_tidy_csv", "--out_posthoc_csv", "--out_dir"}:
                assert Path(args[index + 1]).is_relative_to(tmp_path)


def test_current_result_comparisons_do_not_require_archived_outputs(tmp_path):
    """Keep the complete current CSV comparison set without historical screens.

    Ten primary LMM tables plus five pooled/behavioral tables must remain in the
    reproducibility audit. Ignoring an output is not a valid way to fix cleanup.
    """
    module = load_reproducibility_validator()
    with patch.object(module, "assert_csv_equal") as calls:
        module.compare_primary_lmm_outputs(tmp_path)
        module.compare_correlation_outputs(tmp_path)
    relative = [str(call.args[0].relative_to(ROOT / "data/results")) for call in calls.call_args_list]
    assert len(relative) == 15
    assert len(set(relative)) == 15
    assert sum(name.startswith("pooled_mean_correlations/") for name in relative) == 3
    assert sum(name.startswith("behavior_pairwise_correlations/") for name in relative) == 2
    assert not any("correlational_relationships/" in name or "channel_behavior_relationships/" in name for name in relative)


def test_pytest_discovery_excludes_archived_tests():
    """Archived alternatives must not be collected as current scientific QA.

    Explicit discovery configuration prevents archival tests from silently
    reintroducing obsolete dependencies when users invoke pytest at the root.
    """
    config = (ROOT / "pytest.ini").read_text()
    assert re.search(r"^testpaths\s*=\s*tests\s*$", config, re.MULTILINE)
    assert re.search(r"^norecursedirs\s*=.*\bunused\b", config, re.MULTILINE)
