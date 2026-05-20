"""Validation harness for combined tabular dataset generation.

Run:
  python tests/validate_combined_data_generation_py.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from generate_combined_data import build_combined_dataset, validate_subject_id_column
from process_sociodemographic import (
    EDUCATION_CATEGORIES,
    EDUCATION_OUTPUT_NAMES,
    RACE_CATEGORIES,
    RACE_OUTPUT_NAMES,
    col_by_qid_and_label_contains,
    encode_multi_select_indicators,
    encode_single_select_indicators,
    validate_subject_ids,
)


def test_validate_subject_id_column_rejects_duplicates() -> None:
    df = pd.DataFrame({"subject_id": [1, 1], "value": [10, 11]})
    try:
        validate_subject_id_column(df, dataset_name="dup.csv")
    except ValueError as exc:
        assert "duplicate subject_id values" in str(exc)
        assert "1=2" in str(exc)
        return
    raise AssertionError("Expected duplicate subject_id values to fail.")


def test_validate_subject_id_column_rejects_missing_values() -> None:
    df = pd.DataFrame({"subject_id": [1, None], "value": [10, 11]})
    try:
        validate_subject_id_column(df, dataset_name="missing.csv")
    except ValueError as exc:
        assert "missing/non-numeric subject_id values" in str(exc)
        return
    raise AssertionError("Expected missing subject_id values to fail.")


def test_build_combined_dataset_requires_one_row_per_subject() -> None:
    engagement = pd.DataFrame({"subject_id": [1, 2], "engagement": [0.1, 0.2]})
    socio = pd.DataFrame({"subject_id": [1, 1], "age": [20, 21]})
    recall = pd.DataFrame({"subject_id": [1, 2], "recall": [0.5, 0.6]})
    try:
        build_combined_dataset(engagement, socio, recall)
    except ValueError as exc:
        assert "socio_demographic_data_processed.csv contains invalid subject identifiers" in str(exc)
        assert "duplicate subject_id values" in str(exc)
        return
    raise AssertionError("Expected duplicate socio-demographic IDs to fail before merge.")


def test_build_combined_dataset_preserves_expected_subjects() -> None:
    engagement = pd.DataFrame({"subject_id": [1, 2], "engagement": [0.1, 0.2]})
    socio = pd.DataFrame({"subject_id": [1, 2], "age": [20, 21]})
    recall = pd.DataFrame({"subject_id": [1, 2], "recall": [0.5, 0.6]})

    combined = build_combined_dataset(engagement, socio, recall)
    assert list(combined["subject_id"]) == [1, 2]
    assert list(combined.columns) == ["subject_id", "engagement", "age", "recall"]


def test_process_sociodemographic_subject_validation_rejects_duplicates_and_missing() -> None:
    subject_id = pd.Series([1, 1, pd.NA], dtype="Int64", name="subject_id")
    try:
        validate_subject_ids(subject_id, dataset_name="final_SF_demographic_data.csv")
    except ValueError as exc:
        assert "missing/non-numeric subject_id values" in str(exc)
        return
    raise AssertionError("Expected missing subject IDs to fail before duplicate handling.")


def test_sociodemographic_multiselect_race_preserves_all_selected_categories() -> None:
    race = pd.Series(
        [
            "Asian,White/Caucasian",
            "White/Caucasian,Other",
            pd.NA,
        ]
    )
    race_before = race.copy(deep=True)
    encoded = encode_multi_select_indicators(
        race,
        RACE_CATEGORIES,
        RACE_OUTPUT_NAMES,
        name="race (Q12)",
    )

    assert encoded.loc[0, "race_asian"] == 1.0
    assert encoded.loc[0, "race_white_caucasian"] == 1.0
    assert encoded.loc[0, "race_other"] == 0.0
    assert encoded.loc[1, "race_white_caucasian"] == 1.0
    assert encoded.loc[1, "race_other"] == 1.0
    assert encoded.loc[2].isna().all()
    pd.testing.assert_series_equal(race, race_before)


def test_sociodemographic_race_rejects_unmapped_categories() -> None:
    try:
        encode_multi_select_indicators(
            pd.Series(["Asian,Not in codebook"]),
            RACE_CATEGORIES,
            RACE_OUTPUT_NAMES,
            name="race (Q12)",
        )
    except ValueError as exc:
        assert "Unmapped values for race (Q12)" in str(exc)
        assert "Not in codebook" in str(exc)
        return
    raise AssertionError("Expected unmapped race category to fail hard.")


def test_sociodemographic_single_select_indicators_reject_unmapped_values() -> None:
    try:
        encode_single_select_indicators(
            pd.Series(["High school", "Doctorate"]),
            EDUCATION_CATEGORIES,
            EDUCATION_OUTPUT_NAMES,
            name="highest degree completed (Q2)",
        )
    except ValueError as exc:
        assert "Unmapped values for highest degree completed (Q2)" in str(exc)
        assert "Doctorate" in str(exc)
        return
    raise AssertionError("Expected unmapped education category to fail hard.")


def test_sociodemographic_single_select_encoding_does_not_mutate_input() -> None:
    education = pd.Series(["High school", "Associates", "Bachelor's", "Master's", pd.NA])
    education_before = education.copy(deep=True)

    encoded = encode_single_select_indicators(
        education,
        EDUCATION_CATEGORIES,
        EDUCATION_OUTPUT_NAMES,
        name="highest degree completed (Q2)",
    )

    assert encoded.loc[0, "education_high_school"] == 1.0
    assert encoded.loc[0, "education_masters"] == 0.0
    assert encoded.loc[4].isna().all()
    pd.testing.assert_series_equal(education, education_before)


def test_sociodemographic_education_q2_disambiguates_duplicate_qids_by_label() -> None:
    columns = pd.MultiIndex.from_tuples(
        [
            ("Q2", "What is your highest degree completed? - Selected Choice", '{"ImportId":"QID7"}'),
            ("Q2", "PHQ-9", '{"ImportId":"QID2"}'),
        ]
    )
    df = pd.DataFrame([["Bachelor's", "Not at all"]], columns=columns)

    col = col_by_qid_and_label_contains(df, "Q2", "highest degree completed")
    assert col == ("Q2", "What is your highest degree completed? - Selected Choice", '{"ImportId":"QID7"}')


def main() -> None:
    test_validate_subject_id_column_rejects_duplicates()
    test_validate_subject_id_column_rejects_missing_values()
    test_build_combined_dataset_requires_one_row_per_subject()
    test_build_combined_dataset_preserves_expected_subjects()
    test_process_sociodemographic_subject_validation_rejects_duplicates_and_missing()
    test_sociodemographic_multiselect_race_preserves_all_selected_categories()
    test_sociodemographic_race_rejects_unmapped_categories()
    test_sociodemographic_single_select_indicators_reject_unmapped_values()
    test_sociodemographic_single_select_encoding_does_not_mutate_input()
    test_sociodemographic_education_q2_disambiguates_duplicate_qids_by_label()
    print("[PASS] validate_combined_data_generation_py")


if __name__ == "__main__":
    main()
