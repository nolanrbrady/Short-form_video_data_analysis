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
    encode_education_years,
    encode_multi_select_indicators,
    encode_single_select_indicators,
    validate_education_years_encoding,
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


def test_sociodemographic_education_years_maps_current_categories_without_mutating() -> None:
    """Education years must be an explicit codebook proxy, not inferred or imputed.

    Scientific risk guarded against: silently treating degree labels as arbitrary
    integers would make downstream Pearson correlations depend on implementation
    order rather than the study's auditable years-of-education codebook. The
    expected values are valid because they come directly from the planned
    `education_years_encoding.json` mapping for the current Qualtrics labels.
    """
    education = pd.Series(["High school", "Associates", "Bachelor's", "Master's", pd.NA])
    education_before = education.copy(deep=True)
    encoding = {
        "high_school": 12,
        "associates": 14,
        "bachelor": 16,
        "master": 18,
        "phd": 21,
    }

    years = encode_education_years(education, encoding)

    assert years.name == "education_years"
    assert years.tolist()[:4] == [12.0, 14.0, 16.0, 18.0]
    assert pd.isna(years.iloc[4])
    pd.testing.assert_series_equal(education, education_before)


def test_sociodemographic_education_years_rejects_incomplete_or_non_numeric_config() -> None:
    """The education-years codebook must fail hard when it cannot support Q2.

    Scientific risk guarded against: a missing or text-valued codebook entry
    could silently turn a degree category into missing data or an object column,
    changing pairwise complete cases and correlations. The expected failure is
    valid because all current Q2 categories have required config keys.
    """
    incomplete = {"high_school": 12, "associates": 14, "bachelor": 16}
    try:
        validate_education_years_encoding(incomplete)
    except ValueError as exc:
        assert "missing required education-years config keys" in str(exc)
        assert "master" in str(exc)
    else:
        raise AssertionError("Expected incomplete education-years config to fail hard.")

    non_numeric = {"high_school": 12, "associates": 14, "bachelor": "sixteen", "master": 18}
    try:
        validate_education_years_encoding(non_numeric)
    except ValueError as exc:
        assert "non-numeric education-years config values" in str(exc)
        assert "bachelor" in str(exc)
    else:
        raise AssertionError("Expected non-numeric education-years config to fail hard.")


def test_sociodemographic_education_indicators_are_unchanged_by_years_proxy() -> None:
    """Adding the proxy column must not alter existing one-hot education values.

    Scientific risk guarded against: adding a continuous proxy for exploratory
    correlations should be purely additive and must not change the nominal
    indicator encoding used for demographics reporting. The expected matrix is
    the established one-hot representation for the current Q2 categories.
    """
    education = pd.Series(["High school", "Associates", "Bachelor's", "Master's"])
    indicators = encode_single_select_indicators(
        education,
        EDUCATION_CATEGORIES,
        EDUCATION_OUTPUT_NAMES,
        name="highest degree completed (Q2)",
    )
    _ = encode_education_years(
        education,
        {"high_school": 12, "associates": 14, "bachelor": 16, "master": 18},
    )

    expected = pd.DataFrame(
        {
            "education_high_school": [1.0, 0.0, 0.0, 0.0],
            "education_associates": [0.0, 1.0, 0.0, 0.0],
            "education_bachelors": [0.0, 0.0, 1.0, 0.0],
            "education_masters": [0.0, 0.0, 0.0, 1.0],
        }
    )
    pd.testing.assert_frame_equal(indicators, expected)


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
    test_sociodemographic_education_years_maps_current_categories_without_mutating()
    test_sociodemographic_education_years_rejects_incomplete_or_non_numeric_config()
    test_sociodemographic_education_indicators_are_unchanged_by_years_proxy()
    test_sociodemographic_education_q2_disambiguates_duplicate_qids_by_label()
    print("[PASS] validate_combined_data_generation_py")


if __name__ == "__main__":
    main()
