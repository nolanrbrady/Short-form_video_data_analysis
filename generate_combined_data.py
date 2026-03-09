"""Build the combined tabular dataset used before merging in Homer betas.

This step is intentionally strict because downstream fNIRS analyses assume a
single row per participant. Missing or duplicated `subject_id` values are
treated as hard data-integrity errors rather than being silently dropped or
collapsed.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


GENERATED_TABULAR_DIR = Path("./data/tabular/generated_data")
ENGAGEMENT_CSV = GENERATED_TABULAR_DIR / "engagement_data_processed.csv"
SOCIODEMOGRAPHIC_CSV = GENERATED_TABULAR_DIR / "socio_demographic_data_processed.csv"
RECALL_CSV = GENERATED_TABULAR_DIR / "recall_assessment_score_diffs.csv"
OUTPUT_CSV = GENERATED_TABULAR_DIR / "combined_sfv_data.csv"


def _load_csv(path: Path, *, dataset_name: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"{dataset_name} not found: {path}")
    return pd.read_csv(path)


def validate_subject_id_column(df: pd.DataFrame, *, dataset_name: str) -> None:
    """Fail fast if a dataset is missing a clean one-row-per-subject key."""
    if "subject_id" not in df.columns:
        raise ValueError(f"{dataset_name} is missing required column 'subject_id'.")

    subject_id = pd.to_numeric(df["subject_id"], errors="coerce")
    issues: list[str] = []
    missing_mask = subject_id.isna()
    if missing_mask.any():
        missing_rows = [int(idx) + 2 for idx in df.index[missing_mask][:10]]
        issues.append(
            "missing/non-numeric subject_id values at 1-based CSV row numbers "
            f"{missing_rows}"
        )

    duplicate_counts = subject_id.astype("Int64").value_counts()
    duplicate_counts = duplicate_counts[duplicate_counts > 1]
    if not duplicate_counts.empty:
        examples = ", ".join(f"{int(idx)}={int(count)}" for idx, count in duplicate_counts.head(10).items())
        issues.append(f"duplicate subject_id values {examples}")

    if issues:
        raise ValueError(
            f"{dataset_name} contains invalid subject identifiers: {'; '.join(issues)}. "
            "Resolve subject_id problems before combining datasets."
        )


def build_combined_dataset(
    engagement_data: pd.DataFrame,
    socio_demographic_data: pd.DataFrame,
    recall_assessment_data: pd.DataFrame,
) -> pd.DataFrame:
    """Inner-join the three preprocessed tabular datasets on `subject_id`."""
    validate_subject_id_column(engagement_data, dataset_name="engagement_data_processed.csv")
    validate_subject_id_column(
        socio_demographic_data,
        dataset_name="socio_demographic_data_processed.csv",
    )
    validate_subject_id_column(
        recall_assessment_data,
        dataset_name="recall_assessment_score_diffs.csv",
    )

    combined_data = pd.merge(
        engagement_data,
        socio_demographic_data,
        on="subject_id",
        how="inner",
        validate="one_to_one",
    )
    combined_data = pd.merge(
        combined_data,
        recall_assessment_data,
        on="subject_id",
        how="inner",
        validate="one_to_one",
    )

    if combined_data.empty:
        raise ValueError(
            "Combined dataset is empty after inner joins. "
            "Check subject_id overlap across engagement, socio-demographic, and recall inputs."
        )

    validate_subject_id_column(combined_data, dataset_name="combined_sfv_data.csv")
    return combined_data


def main() -> None:
    engagement_data = _load_csv(ENGAGEMENT_CSV, dataset_name="engagement_data_processed.csv")
    socio_demographic_data = _load_csv(
        SOCIODEMOGRAPHIC_CSV,
        dataset_name="socio_demographic_data_processed.csv",
    )
    recall_assessment_data = _load_csv(RECALL_CSV, dataset_name="recall_assessment_score_diffs.csv")

    combined_data = build_combined_dataset(
        engagement_data,
        socio_demographic_data,
        recall_assessment_data,
    )

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    combined_data.to_csv(OUTPUT_CSV, index=False)

    print(
        f"[combined] wrote {OUTPUT_CSV} with {combined_data.shape[0]} rows and "
        f"{combined_data.shape[1]} columns."
    )


if __name__ == "__main__":
    main()
