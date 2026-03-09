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
from process_sociodemographic import validate_subject_ids


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


def main() -> None:
    test_validate_subject_id_column_rejects_duplicates()
    test_validate_subject_id_column_rejects_missing_values()
    test_build_combined_dataset_requires_one_row_per_subject()
    test_build_combined_dataset_preserves_expected_subjects()
    test_process_sociodemographic_subject_validation_rejects_duplicates_and_missing()
    print("[PASS] validate_combined_data_generation_py")


if __name__ == "__main__":
    main()
