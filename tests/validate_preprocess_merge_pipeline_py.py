"""Synthetic end-to-end validation for collapse -> merge -> certification.

Run:
  python tests/validate_preprocess_merge_pipeline_py.py
"""

from __future__ import annotations

import csv
import json
import shutil
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
from mask_homer_auc_between_subject_outliers import mask_between_subject_outliers
from validate_homer_fir_auc_conversion import validate_homer_fir_auc_conversion

ENGAGEMENT_CONDITION_MAP_JSON = ROOT / "data" / "config" / "engagement_condition_map.json"


def load_engagement_condition_map() -> list[dict[str, object]]:
    payload = json.loads(ENGAGEMENT_CONDITION_MAP_JSON.read_text(encoding="utf-8"))
    conditions = payload.get("conditions")
    if not isinstance(conditions, list) or not conditions:
        raise AssertionError("Engagement condition map fixture must define conditions.")
    return conditions


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
    path.parent.mkdir(parents=True, exist_ok=True)
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
    desired_auc: dict[str, dict[str, float]] = {}
    for subject_num in range(1, 12):
        subject = f"sub_{subject_num:04d}"
        desired_auc[subject] = {
            "S01_D01_Cond01_HbO": 1.0 if subject_num == 11 else 0.0,
            "S01_D01_Cond01_HbR": -0.5 if subject_num % 2 == 0 else 0.5,
            "S01_D01_Cond02_HbO": 0.0,
            "S01_D01_Cond02_HbR": float(subject_num) / 10.0,
            "S01_D01_Cond03_HbO": np.nan if subject_num == 2 else float(subject_num) / 5.0,
            "S01_D01_Cond03_HbR": -float(subject_num) / 8.0,
            "S01_D01_Cond04_HbO": float(subject_num) / 4.0,
            "S01_D01_Cond04_HbR": -float(subject_num) / 6.0,
        }

    path.parent.mkdir(parents=True, exist_ok=True)
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
            "subject_id": [f"{subject_num:04d}" for subject_num in range(1, 12)],
            "covariate_alpha": [subject_num * 10 for subject_num in range(1, 12)],
        }
    )
    df.to_csv(path, index=False)


def write_engagement_input_csv(path: Path, n_subjects: int) -> None:
    condition_map = load_engagement_condition_map()
    base_rating_by_output = {
        "lf_education_engagement": 1.0,
        "lf_entertainment_engagement": 2.0,
        "sf_education_engagement": 3.0,
        "sf_entertainment_engagement": 3.8,
    }
    rows: list[dict[str, object]] = []
    for subject_num in range(1, n_subjects + 1):
        subject_id = f"{subject_num:04d}"
        for condition in condition_map:
            rows.append(
                {
                    "subject_id": subject_id,
                    "Trigger": condition["trigger"],
                    "Category": condition["category"],
                    "Rating": base_rating_by_output[condition["output_column"]] + subject_num / 10.0,
                }
            )
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def write_recall_csv(path: Path, n_subjects: int) -> None:
    df = pd.DataFrame(
        {
            "subject_id": [f"{subject_num:04d}" for subject_num in range(1, n_subjects + 1)],
            "diff_short_form_education": [subject_num * 0.1 for subject_num in range(1, n_subjects + 1)],
            "diff_short_form_entertainment": [subject_num * 0.2 for subject_num in range(1, n_subjects + 1)],
            "diff_long_form_education": [subject_num * 0.3 for subject_num in range(1, n_subjects + 1)],
            "diff_long_form_entertainment": [subject_num * 0.4 for subject_num in range(1, n_subjects + 1)],
        }
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def write_qualtrics_csv(path: Path, n_subjects: int) -> None:
    single_cols = [
        ("Q71", "Study ID", '{"ImportId":"QID71"}'),
        ("Q9", "Age", '{"ImportId":"QID9"}'),
        ("Q16", "PD Status", '{"ImportId":"QID16"}'),
        ("Q80", "SFV Frequency", '{"ImportId":"QID80"}'),
        ("Q81", "SFV Daily Duration", '{"ImportId":"QID81"}'),
        ("Q81_5_TEXT", "SFV Daily Duration - Other Text", '{"ImportId":"QID81_5_TEXT"}'),
    ]
    phq_cols = [(f"QPHQ{i}", "PHQ-9", f'{{"ImportId":"QIDPHQ{i}"}}') for i in range(1, 9)]
    gad_cols = [(f"QGAD{i}", "GAD", f'{{"ImportId":"QIDGAD{i}"}}') for i in range(1, 8)]
    asrs_cols = [(f"QASRS{i}", "ASRS", f'{{"ImportId":"QIDASRS{i}"}}') for i in range(1, 19)]
    yang_pu_cols = [
        (f"QYPU{i}", "Yang Problematic Use", f'{{"ImportId":"QIDYPU{i}"}}') for i in range(1, 22)
    ]
    yang_mot_cols = [
        (f"QYMOT{i}", "Yang Motivation", f'{{"ImportId":"QIDYMOT{i}"}}') for i in range(1, 23)
    ]
    all_cols = single_cols + phq_cols + gad_cols + asrs_cols + yang_pu_cols + yang_mot_cols

    records: list[list[object]] = []
    for subject_num in range(1, n_subjects + 1):
        row: list[object] = [
            f"{subject_num:04d}",
            18 + subject_num,
            "No",
            "Daily",
            "1 - 2 hours",
            "",
        ]
        row.extend(["Several days"] * len(phq_cols))
        row.extend(["Several days"] * len(gad_cols))
        row.extend(["Sometimes"] * len(asrs_cols))
        row.extend(["Neutral"] * len(yang_pu_cols))
        row.extend(["Neutral"] * len(yang_mot_cols))
        records.append(row)

    df = pd.DataFrame(records, columns=pd.MultiIndex.from_tuples(all_cols))
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def build_pipeline_workspace(tmp_path: Path) -> Path:
    workspace = tmp_path / "pipeline_workspace"
    workspace.mkdir(parents=True, exist_ok=True)

    for relative_path in [
        "pipeline_preprocess_merge.sh",
        "process_engagement.py",
        "process_sociodemographic.py",
        "generate_combined_data.py",
        "data/config/engagement_condition_map.json",
        "collapse_homer_fir_to_auc.py",
        "validate_homer_fir_auc_conversion.py",
        "mask_homer_auc_between_subject_outliers.py",
        "merge_homer3_betas_with_combined_data.R",
        "certify_preprocess_merge_integrity.py",
        "homer_fir.py",
    ]:
        src = ROOT / relative_path
        dst = workspace / relative_path
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)

    return workspace


def run_command(args: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(args, cwd=cwd, text=True, capture_output=True, check=False)


def run_process_engagement_fixture(tmp_path: Path, engagement_df: pd.DataFrame) -> subprocess.CompletedProcess[str]:
    workspace = build_pipeline_workspace(tmp_path)
    engagement_csv = workspace / "demographic" / "combined_engagement_data.csv"
    engagement_csv.parent.mkdir(parents=True, exist_ok=True)
    engagement_df.to_csv(engagement_csv, index=False)
    return run_command(["python", "process_engagement.py"], cwd=workspace)


def test_process_engagement_fails_on_missing_subject_condition(tmp_path: Path) -> None:
    base_csv = tmp_path / "engagement_missing_condition.csv"
    write_engagement_input_csv(base_csv, n_subjects=3)
    df = pd.read_csv(base_csv)
    condition_map = load_engagement_condition_map()
    target_category = next(
        condition["category"] for condition in condition_map if condition["output_column"] == "lf_education_engagement"
    )
    df = df[~((df["subject_id"] == 1) & (df["Category"] == target_category))]
    run = run_process_engagement_fixture(tmp_path, df)
    assert run.returncode != 0
    assert "Each subject must have at least one engagement rating for each condition" in (
        run.stderr + run.stdout
    )


def test_process_engagement_fails_on_unknown_category(tmp_path: Path) -> None:
    base_csv = tmp_path / "engagement_unknown_category.csv"
    write_engagement_input_csv(base_csv, n_subjects=3)
    df = pd.read_csv(base_csv)
    df.loc[0, "Category"] = "Short-Form Edu"
    run = run_process_engagement_fixture(tmp_path, df)
    assert run.returncode != 0
    assert "unsupported Category labels" in (run.stderr + run.stdout)


def test_process_engagement_fails_on_trigger_category_mismatch(tmp_path: Path) -> None:
    base_csv = tmp_path / "engagement_bad_trigger.csv"
    write_engagement_input_csv(base_csv, n_subjects=3)
    df = pd.read_csv(base_csv)
    df.loc[0, "Trigger"] = 2
    run = run_process_engagement_fixture(tmp_path, df)
    assert run.returncode != 0
    assert "Trigger/Category mismatches" in (run.stderr + run.stdout)


def test_process_engagement_fails_on_impossible_repeat_count(tmp_path: Path) -> None:
    base_csv = tmp_path / "engagement_too_many_repeats.csv"
    write_engagement_input_csv(base_csv, n_subjects=3)
    df = pd.read_csv(base_csv)
    condition_map = load_engagement_condition_map()
    target_condition = next(
        condition for condition in condition_map if condition["output_column"] == "lf_education_engagement"
    )
    extra_rows = pd.DataFrame(
        [
            {
                "subject_id": 1,
                "Trigger": target_condition["trigger"],
                "Category": target_condition["category"],
                "Rating": 3.0,
            }
            for _ in range(4)
        ]
    )
    df = pd.concat([df, extra_rows], ignore_index=True)
    run = run_process_engagement_fixture(tmp_path, df)
    assert run.returncode != 0
    assert "more repeats than the planned design allows" in (run.stderr + run.stdout)


def test_pipeline_shell_entrypoint_end_to_end(tmp_path: Path) -> None:
    n_subjects = 11
    workspace = build_pipeline_workspace(tmp_path)

    settings_path = workspace / "data" / "config" / "preprocessing_settings.json"
    raw_fir_csv = workspace / "data" / "tabular" / "homer3_glm_betas_wide_fir_pca.csv"
    engagement_csv = workspace / "demographic" / "combined_engagement_data.csv"
    qualtrics_csv = workspace / "qualtrics" / "final_SF_demographic_data.csv"
    recall_csv = workspace / "data" / "tabular" / "generated_data" / "recall_assessment_score_diffs.csv"

    write_settings(settings_path)
    write_raw_fir_csv(raw_fir_csv, settings_path)
    write_engagement_input_csv(engagement_csv, n_subjects=n_subjects)
    write_qualtrics_csv(qualtrics_csv, n_subjects=n_subjects)
    write_recall_csv(recall_csv, n_subjects=n_subjects)

    run = run_command(["bash", "pipeline_preprocess_merge.sh"], cwd=workspace)
    assert run.returncode == 0, run.stderr or run.stdout

    engagement_out = pd.read_csv(workspace / "data" / "tabular" / "generated_data" / "engagement_data_processed.csv")
    sociodem_out = pd.read_csv(workspace / "data" / "tabular" / "generated_data" / "socio_demographic_data_processed.csv")
    combined_out = pd.read_csv(workspace / "data" / "tabular" / "generated_data" / "combined_sfv_data.csv")
    merged_out = pd.read_csv(workspace / "data" / "tabular" / "generated_data" / "homer3_betas_plus_combined_sfv_data_inner_join.csv")
    cert_payload = json.loads((workspace / "data" / "results" / "preprocess_merge_certification.json").read_text(encoding="utf-8"))
    auc_provenance = json.loads((workspace / "data" / "tabular" / "generated_data" / "homer3_glm_betas_wide_auc.provenance.json").read_text(encoding="utf-8"))

    assert engagement_out.shape[0] == n_subjects
    assert sociodem_out.shape[0] == n_subjects
    assert combined_out.shape[0] == n_subjects
    assert merged_out.shape[0] == n_subjects
    assert np.isclose(
        engagement_out.loc[engagement_out["subject_id"] == 1, "lf_education_engagement"].iloc[0],
        1.1,
    )
    assert sociodem_out.loc[sociodem_out["subject_id"] == 1, "phq_total"].iloc[0] == 8
    assert sociodem_out.loc[sociodem_out["subject_id"] == 1, "gad_total"].iloc[0] == 7
    assert sociodem_out.loc[sociodem_out["subject_id"] == 1, "asrs_total"].iloc[0] == 36
    assert cert_payload["passed"] is True
    assert cert_payload["counts"]["merged_rows"] == n_subjects
    assert auc_provenance["input_csv"] == "data/tabular/homer3_glm_betas_wide_fir_pca.csv"


def test_end_to_end_merge_and_certification(tmp_path: Path) -> None:
    settings_path = tmp_path / "settings.json"
    raw_fir_csv = tmp_path / "raw_fir.csv"
    auc_csv = tmp_path / "auc.csv"
    masked_auc_csv = tmp_path / "auc_masked.csv"
    outlier_audit_csv = tmp_path / "auc_outlier_audit.csv"
    outlier_summary_json = tmp_path / "auc_outlier_summary.json"
    combined_csv = tmp_path / "combined.csv"
    merged_csv = tmp_path / "merged.csv"
    cert_json = tmp_path / "cert.json"
    audit_csv = tmp_path / "audit.csv"
    dropped_csv = tmp_path / "dropped.csv"

    write_settings(settings_path)
    write_raw_fir_csv(raw_fir_csv, settings_path)
    write_combined_csv(combined_csv)

    collapse = run_command(
        [
            "python",
            "collapse_homer_fir_to_auc.py",
            "--input-csv",
            str(raw_fir_csv),
            "--output-csv",
            str(auc_csv),
            "--settings-json",
            str(settings_path),
        ],
        cwd=ROOT,
    )
    assert collapse.returncode == 0, collapse.stderr or collapse.stdout

    validate = run_command(
        [
            "python",
            "validate_homer_fir_auc_conversion.py",
            "--raw-fir-csv",
            str(raw_fir_csv),
            "--auc-csv",
            str(auc_csv),
            "--settings-json",
            str(settings_path),
        ],
        cwd=ROOT,
    )
    assert validate.returncode == 0, validate.stderr or validate.stdout
    auc_df = pd.read_csv(auc_csv)
    row1 = auc_df.loc[auc_df["Subject"] == "sub_0001"].iloc[0]
    row2 = auc_df.loc[auc_df["Subject"] == "sub_0002"].iloc[0]
    assert np.isclose(row1["S01_D01_Cond02_HbO"], 0.0)
    assert np.isfinite(row1["S01_D01_Cond02_HbO"])
    assert np.isnan(row2["S01_D01_Cond03_HbO"])
    assert np.isclose(
        auc_df.loc[auc_df["Subject"] == "sub_0011", "S01_D01_Cond01_HbO"].iloc[0],
        1.0,
    )

    mask_between_subject_outliers(
        input_csv=str(auc_csv),
        output_csv=str(masked_auc_csv),
        out_audit_csv=str(outlier_audit_csv),
        out_summary_json=str(outlier_summary_json),
    )
    masked_auc_df = pd.read_csv(masked_auc_csv)
    assert np.isnan(
        masked_auc_df.loc[
            masked_auc_df["Subject"] == "sub_0011", "S01_D01_Cond01_HbO"
        ].iloc[0]
    )
    assert np.isclose(
        masked_auc_df.loc[
            masked_auc_df["Subject"] == "sub_0001", "S01_D01_Cond02_HbO"
        ].iloc[0],
        0.0,
    )
    outlier_audit_df = pd.read_csv(outlier_audit_csv)
    assert outlier_audit_df.shape[0] == 1
    assert outlier_audit_df.loc[0, "Subject"] == "sub_0011"
    assert outlier_audit_df.loc[0, "channel_column"] == "S01_D01_Cond01_HbO"
    outlier_summary = json.loads(outlier_summary_json.read_text(encoding="utf-8"))
    assert outlier_summary["masked_beta_cells"] == 1

    merge = run_command(
        [
            "Rscript",
            "merge_homer3_betas_with_combined_data.R",
            "--homer_csv",
            str(masked_auc_csv),
            "--combined_csv",
            str(combined_csv),
            "--out_csv",
            str(merged_csv),
        ],
        cwd=ROOT,
    )
    assert merge.returncode == 0, merge.stderr or merge.stdout

    merged_df = pd.read_csv(merged_csv)
    assert merged_df.shape[0] == 11
    assert "S01_D01_Cond02_HbO" in merged_df.columns
    assert np.isclose(
        merged_df.loc[merged_df["subject_id"] == 1, "S01_D01_Cond02_HbO"].iloc[0], 0.0
    )
    assert np.isnan(
        merged_df.loc[merged_df["subject_id"] == 11, "S01_D01_Cond01_HbO"].iloc[0]
    )

    certify = run_command(
        [
            "python",
            "certify_preprocess_merge_integrity.py",
            "--combined_csv",
            str(combined_csv),
            "--homer_csv",
            str(masked_auc_csv),
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
        test_process_engagement_fails_on_missing_subject_condition(tmp_path)
        test_process_engagement_fails_on_unknown_category(tmp_path)
        test_process_engagement_fails_on_trigger_category_mismatch(tmp_path)
        test_process_engagement_fails_on_impossible_repeat_count(tmp_path)
        test_pipeline_shell_entrypoint_end_to_end(tmp_path)
        test_end_to_end_merge_and_certification(tmp_path)
        test_auc_lint_fails_if_excluded_channel_is_not_nan(tmp_path)
        test_auc_lint_fails_if_provenance_is_stale(tmp_path)
        test_duplicate_combined_ids_fail_merge(tmp_path)

    print("[PASS] validate_preprocess_merge_pipeline_py")


if __name__ == "__main__":
    main()
