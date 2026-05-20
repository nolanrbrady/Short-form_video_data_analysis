"""Validation harness for the demographics table exclusion logic.

Run:
  python tests/validate_demographics_table_py.py
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

from create_demographics_table import apply_subject_exclusions, create_demographics_table


def write_exclusions(path: Path, values: list[str]) -> None:
    path.write_text(json.dumps(values), encoding="utf-8")


def test_homer_style_exclusions_remove_numeric_subject_ids() -> None:
    with tempfile.TemporaryDirectory() as tmp_dir:
        exclusions_path = Path(tmp_dir) / "excluded_subjects.json"
        write_exclusions(exclusions_path, ["sub_0002", "sub_0004"])
        data = pd.DataFrame(
            {
                "subject_id": [1, 2, 3, 4],
                "homer_subject": ["sub_0001", "sub_0002", "sub_0003", "sub_0004"],
                "age": [20, 21, 22, 23],
                "S01_D01_Cond01_HbO": [0.1, 0.2, 0.3, 0.4],
            }
        )

        filtered, report = apply_subject_exclusions(
            data,
            subject_column="subject_id",
            excluded_subjects_json=exclusions_path,
        )

    assert filtered["subject_id"].tolist() == [1, 3]
    assert report["listed"] == 2
    assert report["matched"] == 2
    assert report["removed"] == 2
    assert report["remaining"] == 2


def test_demographics_table_excludes_id_and_beta_columns() -> None:
    data = pd.DataFrame(
        {
            "subject_id": [1, 3],
            "homer_subject": ["sub_0001", "sub_0003"],
            "age": [20, 22],
            "S01_D01_Cond01_HbO": [0.1, 0.3],
        }
    )

    table = create_demographics_table(data)

    assert list(table.index) == ["age"]
    assert table.loc["age", "n"] == 2
    assert table.loc["age", "mean"] == 21


def test_duplicate_exclusions_fail_after_normalization() -> None:
    with tempfile.TemporaryDirectory() as tmp_dir:
        exclusions_path = Path(tmp_dir) / "excluded_subjects.json"
        write_exclusions(exclusions_path, ["sub_0002", "2"])
        data = pd.DataFrame(
            {
                "subject_id": [1, 2],
                "homer_subject": ["sub_0001", "sub_0002"],
                "age": [20, 21],
            }
        )

        try:
            apply_subject_exclusions(
                data,
                subject_column="homer_subject",
                excluded_subjects_json=exclusions_path,
            )
        except ValueError as exc:
            assert "Duplicate exclusion IDs detected after normalization" in str(exc)
            return
    raise AssertionError("Expected duplicate normalized exclusion IDs to fail.")


def main() -> None:
    test_homer_style_exclusions_remove_numeric_subject_ids()
    test_demographics_table_excludes_id_and_beta_columns()
    test_duplicate_exclusions_fail_after_normalization()
    print("[PASS] validate_demographics_table_py")


if __name__ == "__main__":
    main()
