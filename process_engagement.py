"""Collapse long-format engagement ratings into one mean per subject x condition."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

INPUT_CSV = Path("./demographic/combined_engagement_data.csv")
OUTPUT_CSV = Path("./data/tabular/generated_data/engagement_data_processed.csv")
CONDITION_MAP_JSON = Path("./data/config/engagement_condition_map.json")

REQUIRED_COLUMNS = ("subject_id", "Trigger", "Category", "Rating")
MAX_EXPECTED_REPEATS_PER_CATEGORY = 4
RATING_MIN = 0.0
RATING_MAX = 5.0


def load_condition_map(config_path: Path) -> list[dict[str, object]]:
    payload = json.loads(config_path.read_text(encoding="utf-8"))
    conditions = payload.get("conditions")
    if not isinstance(conditions, list) or not conditions:
        raise ValueError("Engagement condition map must define a non-empty 'conditions' array.")

    required_keys = {"trigger", "category", "output_column"}
    seen_triggers: set[int] = set()
    seen_categories: set[str] = set()
    seen_output_columns: set[str] = set()
    normalized: list[dict[str, object]] = []

    for entry in conditions:
        if not isinstance(entry, dict):
            raise ValueError("Each engagement condition map entry must be a JSON object.")
        missing = sorted(required_keys - set(entry))
        if missing:
            raise ValueError(
                "Engagement condition map entry is missing required keys: "
                + ", ".join(missing)
            )
        trigger = int(entry["trigger"])
        category = str(entry["category"]).strip()
        output_column = str(entry["output_column"]).strip()
        condition_code = str(entry.get("condition_code", "")).strip()
        if trigger in seen_triggers:
            raise ValueError(f"Duplicate trigger in engagement condition map: {trigger}")
        if category in seen_categories:
            raise ValueError(f"Duplicate category in engagement condition map: {category}")
        if output_column in seen_output_columns:
            raise ValueError(
                f"Duplicate output_column in engagement condition map: {output_column}"
            )
        if not category or not output_column:
            raise ValueError("Engagement condition map contains blank category/output_column values.")

        seen_triggers.add(trigger)
        seen_categories.add(category)
        seen_output_columns.add(output_column)
        normalized.append(
            {
                "trigger": trigger,
                "category": category,
                "output_column": output_column,
                "condition_code": condition_code,
            }
        )

    return normalized


def assert_required_columns(df: pd.DataFrame, required_columns: tuple[str, ...]) -> None:
    missing = sorted(set(required_columns) - set(df.columns))
    if missing:
        raise ValueError(
            "Engagement input is missing required columns: " + ", ".join(missing)
        )


def coerce_numeric_strict(series: pd.Series, column_name: str) -> pd.Series:
    raw = series.astype("string").str.strip()
    missing_mask = raw.isna() | (raw == "")
    numeric = pd.to_numeric(raw, errors="coerce")
    bad_mask = ~missing_mask & numeric.isna()
    if bad_mask.any():
        examples = sorted(raw[bad_mask].dropna().unique().tolist())[:10]
        raise ValueError(
            f"Column '{column_name}' contains non-numeric values. Examples: {examples}"
        )
    return numeric


def validate_and_normalize_input(
    df: pd.DataFrame,
    expected_trigger_by_category: dict[str, int],
) -> pd.DataFrame:
    assert_required_columns(df, REQUIRED_COLUMNS)
    normalized = df.loc[:, list(REQUIRED_COLUMNS)].copy()

    normalized["Category"] = normalized["Category"].astype("string").str.strip()
    if normalized["Category"].isna().any() or (normalized["Category"] == "").any():
        raise ValueError("Engagement input contains blank Category values.")

    normalized["subject_id"] = normalized["subject_id"].astype("string").str.strip()
    if normalized["subject_id"].isna().any() or (normalized["subject_id"] == "").any():
        raise ValueError("Engagement input contains blank subject_id values.")

    normalized["Trigger"] = coerce_numeric_strict(normalized["Trigger"], "Trigger")
    normalized["Rating"] = coerce_numeric_strict(normalized["Rating"], "Rating")

    if normalized["Trigger"].isna().any():
        raise ValueError("Engagement input contains missing Trigger values.")
    if normalized["Rating"].isna().any():
        raise ValueError("Engagement input contains missing Rating values.")
    if (~normalized["Trigger"].isin([1, 2, 3, 4])).any():
        bad = sorted(normalized.loc[~normalized["Trigger"].isin([1, 2, 3, 4]), "Trigger"].unique().tolist())
        raise ValueError(f"Engagement input contains unsupported Trigger values: {bad}")
    if (~normalized["Rating"].between(RATING_MIN, RATING_MAX)).any():
        bad = (
            normalized.loc[~normalized["Rating"].between(RATING_MIN, RATING_MAX), "Rating"]
            .drop_duplicates()
            .sort_values()
            .tolist()
        )
        raise ValueError(
            f"Engagement ratings must fall within [{RATING_MIN}, {RATING_MAX}]. Bad values: {bad}"
        )

    unexpected_categories = sorted(
        set(normalized["Category"].dropna().unique()) - set(expected_trigger_by_category)
    )
    if unexpected_categories:
        raise ValueError(
            "Engagement input contains unsupported Category labels: "
            + ", ".join(unexpected_categories)
        )

    expected_trigger = normalized["Category"].map(expected_trigger_by_category)
    mismatch = normalized["Trigger"] != expected_trigger
    if mismatch.any():
        offenders = (
            normalized.loc[mismatch, ["subject_id", "Category", "Trigger"]]
            .assign(expected_trigger=expected_trigger[mismatch].to_numpy())
            .head(10)
            .to_dict("records")
        )
        raise ValueError(
            "Engagement input contains Trigger/Category mismatches. "
            f"Examples: {offenders}"
        )

    counts = (
        normalized.groupby(["subject_id", "Category"], as_index=False)
        .size()
        .rename(columns={"size": "n_rows"})
    )
    expected_index = pd.MultiIndex.from_product(
        [
            sorted(normalized["subject_id"].unique().tolist()),
            list(expected_trigger_by_category.keys()),
        ],
        names=["subject_id", "Category"],
    )
    counts_full = (
        counts.set_index(["subject_id", "Category"])
        .reindex(expected_index, fill_value=0)
        .reset_index()
    )
    missing = counts_full[counts_full["n_rows"] == 0]
    if not missing.empty:
        offenders = missing.head(10).to_dict("records")
        raise ValueError(
            "Each subject must have at least one engagement rating for each condition. "
            f"Missing subject/category cells: {offenders}"
        )

    too_many = counts_full[counts_full["n_rows"] > MAX_EXPECTED_REPEATS_PER_CATEGORY]
    if not too_many.empty:
        offenders = too_many.head(10).to_dict("records")
        raise ValueError(
            "Engagement input contains more repeats than the planned design allows "
            f"({MAX_EXPECTED_REPEATS_PER_CATEGORY} per subject x condition). "
            f"Examples: {offenders}"
        )

    return normalized


def build_processed_output(df: pd.DataFrame, output_column_by_category: dict[str, str]) -> pd.DataFrame:
    processed = (
        df.groupby(["subject_id", "Category"], as_index=False)["Rating"]
        .mean()
        .pivot(index="subject_id", columns="Category", values="Rating")
        .rename(columns=output_column_by_category)
        .reset_index()
    )

    expected_output_columns = ["subject_id", *output_column_by_category.values()]
    missing_output_columns = [
        column for column in expected_output_columns if column not in processed.columns
    ]
    if missing_output_columns:
        raise ValueError(
            "Processed engagement output is missing expected columns: "
            + ", ".join(missing_output_columns)
        )
    if processed["subject_id"].duplicated().any():
        dupes = processed.loc[processed["subject_id"].duplicated(), "subject_id"].tolist()
        raise ValueError(
            "Processed engagement output contains duplicate subject_id values: "
            + ", ".join(dupes[:10])
        )
    value_columns = list(output_column_by_category.values())
    if processed[value_columns].isna().any().any():
        offenders = processed.loc[processed[value_columns].isna().any(axis=1), ["subject_id", *value_columns]]
        raise ValueError(
            "Processed engagement output contains missing subject-condition means. "
            f"Examples: {offenders.head(10).to_dict('records')}"
        )

    return processed.loc[:, expected_output_columns].sort_values("subject_id").reset_index(
        drop=True
    )


def main() -> None:
    condition_map = load_condition_map(CONDITION_MAP_JSON)
    expected_trigger_by_category = {
        entry["category"]: entry["trigger"] for entry in condition_map
    }
    output_column_by_category = {
        entry["category"]: entry["output_column"] for entry in condition_map
    }
    engagement_data = pd.read_csv(INPUT_CSV)
    validated = validate_and_normalize_input(
        engagement_data,
        expected_trigger_by_category=expected_trigger_by_category,
    )
    processed_data = build_processed_output(
        validated,
        output_column_by_category=output_column_by_category,
    )
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    processed_data.to_csv(OUTPUT_CSV, index=False)


if __name__ == "__main__":
    main()
