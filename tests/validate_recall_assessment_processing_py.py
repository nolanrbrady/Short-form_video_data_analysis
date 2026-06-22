"""Validation harness for recall assessment scoring.

Run:
  python tests/validate_recall_assessment_processing_py.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from demographic.process_recall_assessment import (
    apply_question_aliases,
    build_audit_table,
    compute_subject_condition_means,
    load_invalid_question_config,
    load_question_alias_config,
    unmatched_key_rows_for_assessment,
    validate_invalid_question_config,
)


def make_key_df() -> pd.DataFrame:
    rows = [
        ("Q1", "Short-form Education", "sfe valid"),
        ("Q10", "Short-form Education", "sfe invalid 1"),
        ("Q26", "Short-form Education", "sfe invalid 2"),
        ("Q2", "Short-form Entertainment", "sft valid"),
        ("Q28", "Short-form Entertainment", "sft invalid"),
        ("Q3", "Long-form Education", "lfe valid"),
        ("Q22", "Long-form Education", "lfe alias canonical"),
        ("Q5", "Long-form Education", "lfe invalid 1"),
        ("Q6", "Long-form Education", "lfe invalid 2"),
        ("Q35", "Long-form Education", "lfe invalid 2"),
        ("Q4", "Long-form Entertainment", "lft valid"),
        ("Q7", "Long-form Entertainment", "lft invalid 1"),
        ("Q36", "Long-form Entertainment", "lft invalid 1"),
        ("Q8", "Long-form Entertainment", "lft invalid 2"),
    ]
    return pd.DataFrame(rows, columns=["Question ID", "Condition", "Answer"])


def make_assessment_df() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Q34": "001",
                "Q1": "sfe valid",
                "Q10": "wrong",
                "Q26": "wrong",
                "Q2": "sft valid",
                "Q28": "wrong",
                "Q3": "lfe valid",
                "Q39": "lfe alias canonical",
                "Q5": "wrong",
                "Q6": "wrong",
                "Q35": "wrong",
                "Q4": "lft valid",
                "Q7": "wrong",
                "Q36": "wrong",
                "Q8": "wrong",
            }
        ]
    )


def test_invalid_recall_questions_are_excluded_from_condition_denominators() -> None:
    key_df = make_key_df()
    assessment_df = make_assessment_df()
    invalid_questions = load_invalid_question_config()
    aliases = load_question_alias_config()
    key_with_aliases = apply_question_aliases(key_df, assessment_df, aliases, assessment_label="pre")
    validate_invalid_question_config(key_with_aliases, invalid_questions)

    means = compute_subject_condition_means(
        assessment_df,
        key_with_aliases,
        invalid_question_ids=set(invalid_questions["question_id"]),
    )

    row = means.loc["001"]
    assert row["short_form_education"] == 1.0
    assert row["short_form_entertainment"] == 1.0
    assert row["long_form_education"] == 1.0
    assert row["long_form_entertainment"] == 1.0

    unfiltered = compute_subject_condition_means(assessment_df, key_with_aliases)
    unfiltered_row = unfiltered.loc["001"]
    assert unfiltered_row["short_form_education"] == 1.0 / 3.0
    assert unfiltered_row["short_form_entertainment"] == 0.5
    assert unfiltered_row["long_form_education"] == 0.4
    assert unfiltered_row["long_form_entertainment"] == 0.25


def test_invalid_recall_questions_are_visible_but_not_scored_in_audit() -> None:
    key_df = make_key_df()
    assessment_df = make_assessment_df()
    invalid_questions = load_invalid_question_config()
    aliases = load_question_alias_config()
    key_with_aliases = apply_question_aliases(key_df, assessment_df, aliases, assessment_label="pre")

    audit = build_audit_table(assessment_df, key_with_aliases, invalid_questions=invalid_questions)
    excluded = audit[audit["excluded_from_retention_score"]]

    assert set(excluded["question_id"]) == {"Q5", "Q6", "Q7", "Q8", "Q10", "Q26", "Q28", "Q35", "Q36"}
    assert set(excluded["method"]) == {"excluded_invalid_question"}
    assert excluded["score"].isna().all()
    assert excluded["exclusion_reason"].str.len().gt(0).all()
    assert "Q39" in set(audit["question_id"])
    assert audit.loc[audit["question_id"] == "Q39", "key_answer"].iloc[0] == "lfe alias canonical"


def test_invalid_recall_question_config_rejects_condition_mismatch() -> None:
    key_df = make_key_df()
    invalid_questions = load_invalid_question_config().copy()
    invalid_questions.loc[invalid_questions["question_id"] == "Q5", "condition"] = "Short-form Education"

    try:
        validate_invalid_question_config(key_df, invalid_questions)
    except ValueError as exc:
        assert "Q5" in str(exc)
        assert "does not match key condition" in str(exc)
        return
    raise AssertionError("Expected invalid-question condition mismatch to fail.")


def test_alias_presence_suppresses_false_missing_key_warning() -> None:
    key_df = make_key_df()
    assessment_df = make_assessment_df()
    aliases = load_question_alias_config()
    key_with_aliases = apply_question_aliases(key_df, assessment_df, aliases, assessment_label="pre")

    unmatched = unmatched_key_rows_for_assessment(key_with_aliases, assessment_df)
    assert "Q22" not in set(unmatched["Question ID"])
    assert "Q39" not in set(unmatched["Question ID"])


def run_all_tests() -> None:
    test_invalid_recall_questions_are_excluded_from_condition_denominators()
    test_invalid_recall_questions_are_visible_but_not_scored_in_audit()
    test_invalid_recall_question_config_rejects_condition_mismatch()
    test_alias_presence_suppresses_false_missing_key_warning()


if __name__ == "__main__":
    run_all_tests()
    print("All recall assessment processing validation checks passed.")
