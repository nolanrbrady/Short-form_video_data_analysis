"""Create a descriptive demographics table from the merged SFV analysis file.

The demographics table is intentionally generated after applying the shared
subject-exclusion manifest. Subject IDs are normalized before filtering so the
same participant is excluded whether represented as `sub_0050`, `0050`, or
`50`. Keeping this as a single auditable exclusion source follows reproducible
workflow guidance (Sandve et al., 2013; see CITATIONS.md).
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Iterable

import pandas as pd

DEFAULT_INPUT_CSV = Path("data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv")
DEFAULT_EXCLUDED_SUBJECTS_JSON = Path("data/config/excluded_subjects.json")


def normalize_subject_id(value: object, *, column_name: str) -> int:
    """Return the numeric participant ID embedded in a subject identifier.

    Accepts both Homer-style IDs (`sub_0041`) and tabular IDs (`0041`, `41`).
    Raises instead of guessing when no numeric ID can be parsed.
    """

    if pd.isna(value):
        raise ValueError(f"Column '{column_name}' contains a missing subject ID.")

    match = re.search(r"\d+", str(value).strip())
    if match is None:
        raise ValueError(
            f"Failed to parse a numeric subject ID from column '{column_name}'. "
            f"Unparseable value: {value!r}"
        )
    return int(match.group(0))


def normalize_subject_ids(values: Iterable[object], *, column_name: str) -> pd.Series:
    """Normalize a sequence of subject identifiers to integer participant IDs."""

    return pd.Series(
        [normalize_subject_id(value, column_name=column_name) for value in values],
        dtype="Int64",
    )


def load_excluded_subject_ids(excluded_subjects_json: Path) -> pd.Series:
    """Load and validate the shared excluded-subject manifest."""

    try:
        with excluded_subjects_json.open(encoding="utf-8") as f:
            payload = json.load(f)
    except FileNotFoundError as exc:
        raise ValueError(f"Exclusion file not found: {excluded_subjects_json}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Failed to parse exclusion JSON at {excluded_subjects_json}: {exc}"
        ) from exc

    if payload is None:
        payload = []
    if not isinstance(payload, list):
        raise ValueError(
            "Exclusion file must be a top-level JSON array of subject IDs. "
            f"File: {excluded_subjects_json}"
        )
    if any(subject is None or str(subject).strip() == "" for subject in payload):
        raise ValueError(f"Exclusion file contains empty subject IDs: {excluded_subjects_json}")

    excluded_ids = normalize_subject_ids(
        payload,
        column_name=f"exclude_subjects_json({excluded_subjects_json})",
    )
    duplicates = sorted(excluded_ids[excluded_ids.duplicated()].astype(int).unique().tolist())
    if duplicates:
        raise ValueError(
            "Duplicate exclusion IDs detected after normalization in "
            f"{excluded_subjects_json}: {duplicates}"
        )
    return excluded_ids


def apply_subject_exclusions(
    data: pd.DataFrame,
    *,
    subject_column: str,
    excluded_subjects_json: Path,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Remove excluded subjects and return the filtered data plus an audit summary."""

    if subject_column not in data.columns:
        raise ValueError(f"Expected subject column '{subject_column}' in merged input.")

    filtered = data.copy()
    filtered["_subject_id_norm"] = normalize_subject_ids(
        filtered[subject_column],
        column_name=subject_column,
    ).to_numpy()
    excluded_ids = load_excluded_subject_ids(excluded_subjects_json)

    present_ids = set(filtered["_subject_id_norm"].astype(int).tolist())
    excluded_set = set(excluded_ids.astype(int).tolist())
    matched_ids = sorted(present_ids & excluded_set)
    missing_ids = sorted(excluded_set - present_ids)

    n_subjects_before = filtered["_subject_id_norm"].nunique()
    filtered = filtered.loc[~filtered["_subject_id_norm"].isin(excluded_set)].copy()
    n_subjects_after = filtered["_subject_id_norm"].nunique()
    filtered = filtered.drop(columns=["_subject_id_norm"])

    report = {
        "listed": len(excluded_set),
        "matched": len(matched_ids),
        "removed": int(n_subjects_before - n_subjects_after),
        "remaining": int(n_subjects_after),
        "matched_ids": matched_ids,
        "missing_ids": missing_ids,
    }
    return filtered, report


def demographic_columns(data: pd.DataFrame) -> list[str]:
    """Return non-identifier, non-beta columns for the descriptive table."""

    return [
        col
        for col in data.columns
        if not col.startswith("S0") and col not in {"subject_id", "homer_subject"}
    ]


def create_demographics_table(data: pd.DataFrame) -> pd.DataFrame:
    """Create the transposed `describe()` table for demographic/covariate columns."""

    columns = demographic_columns(data)
    if not columns:
        raise ValueError("No demographic columns remain after excluding beta and ID columns.")

    demographics_table = data[columns].describe().transpose()
    return demographics_table.rename(columns={"count": "n"})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a descriptive demographics table after subject exclusions."
    )
    parser.add_argument("--input-csv", type=Path, default=DEFAULT_INPUT_CSV)
    parser.add_argument(
        "--excluded-subjects-json",
        type=Path,
        default=DEFAULT_EXCLUDED_SUBJECTS_JSON,
    )
    parser.add_argument("--subject-column", default="homer_subject")
    parser.add_argument("--output-csv", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = pd.read_csv(args.input_csv)
    data, exclusion_report = apply_subject_exclusions(
        data,
        subject_column=args.subject_column,
        excluded_subjects_json=args.excluded_subjects_json,
    )
    table = create_demographics_table(data)

    print(
        "[exclude] demographics: "
        f"listed={exclusion_report['listed']}, "
        f"matched={exclusion_report['matched']}, "
        f"removed={exclusion_report['removed']}, "
        f"remaining={exclusion_report['remaining']}"
    )
    if exclusion_report["missing_ids"]:
        print(
            "[warn] demographics: exclusion IDs not present in dataset: "
            + ", ".join(str(subject_id) for subject_id in exclusion_report["missing_ids"])
        )
    print(table)

    if args.output_csv is not None:
        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(args.output_csv)


if __name__ == "__main__":
    main()
