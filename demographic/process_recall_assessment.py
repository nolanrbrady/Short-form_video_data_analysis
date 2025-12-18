import os
import re
import unicodedata
from difflib import SequenceMatcher
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

ASSESSMENT_DIR = "../../Assessment"

drop_columns = ["StartDate", "EndDate", "Status", "Progress", "Duration (in seconds)", "Finished", "RecordedDate", "ResponseId", "DistributionChannel", "UserLanguage"]

_ARROW_RE = re.compile(r"\s*-{1,4}\s*>\s*")
_WS_RE = re.compile(r"\s+")

# These are the four categories you said you want at the end.
# They must match the spelling in Recall_Assessment_Key.csv's `Condition` column.
EXPECTED_CONDITIONS = [
    "Short-form Education",
    "Short-form Entertainment",
    "Long-form Education",
    "Long-form Entertainment",
]

def _condition_to_col(condition: str) -> str:
    """
    Convert key_df's Condition into a stable snake_case column name.
    Example: "Short-form Education" -> "short_form_education"
    """
    s = str(condition).strip().casefold()
    s = s.replace("-", " ")
    s = _WS_RE.sub(" ", s).strip()
    return s.replace(" ", "_")

def _strip_diacritics(text: str) -> str:
    """
    Remove combining marks (accents) from Latin text while keeping non-Latin scripts intact.
    """
    decomposed = unicodedata.normalize("NFKD", text)
    return "".join(ch for ch in decomposed if not unicodedata.combining(ch))

def _normalize_answer(text: str) -> str:
    """
    Normalize free-text answers so comparisons are robust to:
    - smart quotes / curly apostrophes
    - variable arrow formatting ('-->' vs '->')
    - extra whitespace and surrounding quotes
    - diacritics on Latin characters
    """
    if text is None or (isinstance(text, float) and np.isnan(text)):  # handles NaN from pandas
        return ""

    s = str(text).strip()
    if not s:
        return ""

    # Normalize common unicode punctuation variants
    s = s.replace("’", "'").replace("‘", "'").replace("`", "'")
    s = s.replace("“", '"').replace("”", '"')

    # Standardize arrow formatting: "кто --> ktaw" etc.
    s = _ARROW_RE.sub(" -> ", s)

    # Normalize unicode + strip diacritics (helps: "Priv’et" vs "Priv'et", etc.)
    s = _strip_diacritics(s)

    # Collapse whitespace and casefold (unicode-safe lowercasing)
    s = _WS_RE.sub(" ", s).strip().casefold()

    # Drop surrounding quotes if present (common in CSV exports)
    if len(s) >= 2 and ((s[0] == s[-1] == '"') or (s[0] == s[-1] == "'")):
        s = s[1:-1].strip()

    return s

def _split_acceptable_answers(key: str) -> list[str]:
    """
    Support multiple acceptable keys if you ever encode them, e.g.:
    "neil degrasse tyson|neil tyson"
    "bar;bars"
    """
    s = str(key) if key is not None else ""
    parts = re.split(r"[|;]", s)
    return [p.strip() for p in parts if p and p.strip()]

def grade_response(response: str, key: str) -> int:
    """
    Return 1 if the response matches the key (robustly), else 0.

    Matching rules:
    - First try exact match after normalization.
    - If not exact, allow a conservative fuzzy match to tolerate tiny typos/quote differences.
    """
    resp_n = _normalize_answer(response)
    if not resp_n:
        return 0

    key_options = _split_acceptable_answers(key)
    if not key_options:
        return 0

    normalized_keys = [_normalize_answer(k) for k in key_options]

    # Exact match after normalization
    if resp_n in normalized_keys:
        return 1

    # Conservative fuzzy match (keeps false-positives low for short answers like colors)
    for k in normalized_keys:
        if not k:
            continue
        if min(len(resp_n), len(k)) <= 3:
            continue  # require exact match for very short answers
        if SequenceMatcher(None, resp_n, k).ratio() >= 0.90:
            return 1

    return 0

def grade_response_detailed(response: str, key: str) -> dict:
    """
    Detailed version of `grade_response` for auditing.

    Returns a dict including:
    - score: 0/1
    - method: 'empty' | 'exact' | 'fuzzy' | 'no_match' | 'no_key'
    - response_norm, key_norm: normalized strings used for matching
    - fuzzy_ratio: float or NaN (only meaningful when method == 'fuzzy')
    """
    resp_norm = _normalize_answer(response)
    if not resp_norm:
        return {
            "score": 0,
            "method": "empty",
            "response_norm": resp_norm,
            "key_norm": "",
            "fuzzy_ratio": np.nan,
        }

    key_options = _split_acceptable_answers(key)
    if not key_options:
        return {
            "score": 0,
            "method": "no_key",
            "response_norm": resp_norm,
            "key_norm": "",
            "fuzzy_ratio": np.nan,
        }

    normalized_keys = [_normalize_answer(k) for k in key_options]

    # Exact match after normalization
    if resp_norm in normalized_keys:
        return {
            "score": 1,
            "method": "exact",
            "response_norm": resp_norm,
            "key_norm": next((k for k in normalized_keys if k == resp_norm), normalized_keys[0]),
            "fuzzy_ratio": np.nan,
        }

    # Conservative fuzzy match (disabled for very short answers)
    best_ratio = -1.0
    best_key = ""
    for k in normalized_keys:
        if not k:
            continue
        if min(len(resp_norm), len(k)) <= 3:
            continue
        ratio = SequenceMatcher(None, resp_norm, k).ratio()
        if ratio > best_ratio:
            best_ratio = ratio
            best_key = k
        if ratio >= 0.90:
            return {
                "score": 1,
                "method": "fuzzy",
                "response_norm": resp_norm,
                "key_norm": k,
                "fuzzy_ratio": ratio,
            }

    return {
        "score": 0,
        "method": "no_match",
        "response_norm": resp_norm,
        "key_norm": best_key,
        "fuzzy_ratio": (best_ratio if best_ratio >= 0 else np.nan),
    }

def build_key_map(key_df: pd.DataFrame) -> Dict[str, str]:
    """
    Build a mapping from Qualtrics column name (e.g., 'Q1', 'QID14') to correct answer text.
    """
    required = {"Question ID", "Answer"}
    if not required.issubset(set(key_df.columns)):
        raise ValueError("key_df must contain columns: 'Question ID' and 'Answer'")

    # Keep the first occurrence if duplicates exist.
    return (
        key_df[["Question ID", "Answer"]]
        .dropna(subset=["Question ID"])
        .drop_duplicates(subset=["Question ID"])
        .set_index("Question ID")["Answer"]
        .to_dict()
    )

def build_question_condition_map(key_df: pd.DataFrame) -> Dict[str, str]:
    """
    Build a mapping from question id (e.g., 'Q1', 'QID14') -> Condition string
    (e.g., 'Short-form Education').
    """
    required = {"Question ID", "Condition"}
    if not required.issubset(set(key_df.columns)):
        raise ValueError("key_df must contain columns: 'Question ID' and 'Condition'")

    return (
        key_df[["Question ID", "Condition"]]
        .dropna(subset=["Question ID", "Condition"])
        .drop_duplicates(subset=["Question ID"])
        .set_index("Question ID")["Condition"]
        .to_dict()
    )

def compute_subject_condition_means(
    assessment_df: pd.DataFrame,
    key_df: pd.DataFrame,
    *,
    subject_id_col: str = "Q34",
    conditions: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    For each subject, compute mean correctness (0/1) within each Condition category.

    Returns a wide DataFrame indexed by subject_id with one column per condition (snake_case).
    """
    key_map = build_key_map(key_df)
    q_to_cond = build_question_condition_map(key_df)

    if conditions is None:
        conditions = EXPECTED_CONDITIONS
    conditions_set = set(conditions)

    question_cols = [
        c
        for c in assessment_df.columns
        if c in key_map and c in q_to_cond and q_to_cond[c] in conditions_set
    ]

    long_df = assessment_df[[subject_id_col] + question_cols].melt(
        id_vars=[subject_id_col],
        var_name="question_id",
        value_name="response",
    )
    long_df["key_answer"] = long_df["question_id"].map(key_map)
    long_df["condition"] = long_df["question_id"].map(q_to_cond)

    # Grade row-by-row (dataset is small enough that clarity > micro-optimizations)
    long_df["score"] = long_df.apply(
        lambda r: grade_response(r["response"], r["key_answer"]),
        axis=1,
    )

    grouped = (
        long_df.groupby([subject_id_col, "condition"], dropna=False)["score"]
        .mean()
        .reset_index()
    )
    grouped["condition_col"] = grouped["condition"].map(_condition_to_col)

    wide = grouped.pivot(index=subject_id_col, columns="condition_col", values="score")

    # Ensure all expected columns exist (even if empty)
    for cond in conditions:
        col = _condition_to_col(cond)
        if col not in wide.columns:
            wide[col] = np.nan

    # Stable column ordering
    wide = wide[[_condition_to_col(c) for c in conditions]]
    return wide

def build_audit_table(
    assessment_df: pd.DataFrame,
    key_df: pd.DataFrame,
    *,
    subject_id_col: str = "Q34",
    conditions: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Build a full audit table: one row per (subject, question) with raw/normalized text and scoring details.
    """
    key_map = build_key_map(key_df)
    q_to_cond = build_question_condition_map(key_df)

    if conditions is None:
        conditions = EXPECTED_CONDITIONS
    conditions_set = set(conditions)

    question_cols = [
        c
        for c in assessment_df.columns
        if c in key_map and c in q_to_cond and q_to_cond[c] in conditions_set
    ]

    long_df = assessment_df[[subject_id_col] + question_cols].melt(
        id_vars=[subject_id_col],
        var_name="question_id",
        value_name="response_raw",
    )
    long_df["key_answer"] = long_df["question_id"].map(key_map)
    long_df["condition"] = long_df["question_id"].map(q_to_cond)

    details = long_df.apply(
        lambda r: grade_response_detailed(r["response_raw"], r["key_answer"]),
        axis=1,
        result_type="expand",
    )
    out = pd.concat([long_df, details], axis=1)
    return out

def grade_subject_row(
    subject_row: pd.Series,
    key_map: Dict[str, str],
    *,
    question_cols: Optional[List[str]] = None,
    include_per_item_scores: bool = True,
) -> pd.Series:
    """
    Grade all recall questions for a single participant row.

    Returns a Series with:
    - n_items_graded
    - n_correct
    - pct_correct
    - (optional) per-item columns like 'Q1_score'
    """
    if question_cols is None:
        # Grade only columns that have keys, ignoring obvious non-question columns.
        non_question_cols = set(drop_columns + ["Q34", "SC0"])
        question_cols = [c for c in subject_row.index if c in key_map and c not in non_question_cols]

    per_item: Dict[str, int] = {}
    n_graded = 0
    n_correct = 0
    for col in question_cols:
        key = key_map.get(col)
        if key is None:
            continue
        score = grade_response(subject_row.get(col), key)
        n_graded += 1
        n_correct += score
        if include_per_item_scores:
            per_item[f"{col}_score"] = score

    pct = (n_correct / n_graded) if n_graded else np.nan
    out: Dict[str, object] = {
        "n_items_graded": n_graded,
        "n_correct": n_correct,
        "pct_correct": pct,
    }
    out.update(per_item)
    return pd.Series(out)

def _load_assessment_csvs(assessment_dir: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pre_df = pd.read_csv(f"{assessment_dir}/pretask_assessment.csv").iloc[2:]
    post_df = pd.read_csv(f"{assessment_dir}/posttask_assessment.csv").iloc[2:]
    key_df = pd.read_csv(f"{assessment_dir}/Recall_Assessment_Key.csv")

    pre_df = pre_df.drop(columns=drop_columns)
    post_df = post_df.drop(columns=drop_columns)
    return pre_df, post_df, key_df


if __name__ == "__main__":
    pre_df, post_df, key_df = _load_assessment_csvs(ASSESSMENT_DIR)

    print("=== Recall grading: inputs ===")
    print(f"pre_df shape:  {pre_df.shape}")
    print(f"post_df shape: {post_df.shape}")
    print(f"key_df shape:  {key_df.shape}")
    print(f"key_df conditions: {sorted(set(key_df['Condition'].dropna().unique()))}")

    # Coverage checks (critical for catching mismatches)
    key_questions = set(key_df["Question ID"].dropna().unique())
    pre_questions = set(pre_df.columns) & key_questions
    post_questions = set(post_df.columns) & key_questions
    print("=== Recall grading: coverage checks ===")
    print(f"questions in key: {len(key_questions)}")
    print(f"questions present in pre:  {len(pre_questions)}")
    print(f"questions present in post: {len(post_questions)}")
    missing_in_pre = sorted(key_questions - set(pre_df.columns))
    missing_in_post = sorted(key_questions - set(post_df.columns))
    if missing_in_pre:
        print(f"WARNING: {len(missing_in_pre)} key questions missing from pre_df columns (showing up to 15): {missing_in_pre[:15]}")
    if missing_in_post:
        print(f"WARNING: {len(missing_in_post)} key questions missing from post_df columns (showing up to 15): {missing_in_post[:15]}")

    # Compute per-subject mean score in each condition, separately for pre and post.
    pre_means = compute_subject_condition_means(pre_df, key_df)
    post_means = compute_subject_condition_means(post_df, key_df)

    print("=== Recall grading: per-condition question counts used ===")
    q_to_cond = build_question_condition_map(key_df)
    for cond in EXPECTED_CONDITIONS:
        qids = [q for q, c in q_to_cond.items() if c == cond]
        used_pre = [q for q in qids if q in pre_df.columns]
        used_post = [q for q in qids if q in post_df.columns]
        print(f"{cond}: key={len(qids)} used_pre={len(used_pre)} used_post={len(used_post)}")

    # Align subjects and compute post - pre difference per condition.
    diffs = post_means.sub(pre_means, fill_value=np.nan)
    diffs.columns = [f"diff_{c}" for c in diffs.columns]

    out_df = diffs.reset_index().rename(columns={"Q34": "subject_id"})

    # Build audit tables (one row per subject-question) and save alongside outputs.
    print("=== Recall grading: building audit tables (pre/post) ===")
    pre_audit = build_audit_table(pre_df, key_df)
    post_audit = build_audit_table(post_df, key_df)
    print(f"pre_audit rows:  {len(pre_audit)}")
    print(f"post_audit rows: {len(post_audit)}")
    print("Audit method counts (pre):")
    print(pre_audit["method"].value_counts(dropna=False))
    print("Audit method counts (post):")
    print(post_audit["method"].value_counts(dropna=False))

    # Save to CSV next to this script by default.
    out_path = os.path.join("../data/tabular/", "recall_assessment_score_diffs.csv")
    pre_audit_path = os.path.join(os.path.dirname(__file__), "recall_assessment_audit_pre.csv")
    post_audit_path = os.path.join(os.path.dirname(__file__), "recall_assessment_audit_post.csv")
    out_df.to_csv(out_path, index=False)
    pre_audit.to_csv(pre_audit_path, index=False)
    post_audit.to_csv(post_audit_path, index=False)

    print("=== Recall grading: outputs ===")
    print(out_df.head())
    print(f"Wrote: {out_path}")
    print(f"Wrote: {pre_audit_path}")
    print(f"Wrote: {post_audit_path}")
