"""
Covariate preprocessing for the Short Form Video (SFV) study.

This script is intentionally conservative and paper-friendly:
- Loads a Qualtrics export that has 3 header rows (QID, label, ImportId) into a DataFrame with MultiIndex columns.
- Encodes text responses into numeric (ordinal) values for correlations.
- Computes per-participant composite totals for multi-item scales (PHQ-9, GAD, ASRS, Yang PU, Yang Motivation).
- Outputs a clean covariate dataset (correlation diagnostics are handled in `covariate_correlation_analysis.py`).

Input (currently):
- qualtrics/final_SF_demographic_data.csv

Outputs (to covariate_outputs/):
- covariates_clean.csv
- covariate_missingness.csv
- covariate_column_audit.csv
- sfv_duration_other_audit.csv

Outputs (to data/tabular/):
- socio_demographic_data_processed.csv
"""

from __future__ import annotations

from math import ceil
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

try:
    import pandas as pd
except ModuleNotFoundError as e:  # pragma: no cover
    raise SystemExit(
        "Missing dependency 'pandas'.\n\n"
        "This script expects a scientific Python environment (e.g., conda/venv) with:\n"
        "- pandas\n"
        "- numpy (required by pandas)\n\n"
        "Activate your project environment and re-run."
    ) from e


# =============================================================================
# Configuration
# =============================================================================

# Minimum proportion of scale items required to compute a participant score.
MIN_ITEM_PROPORTION = 0.75

# Ordinal encodings (higher = more / higher severity).
LIKERT_0_TO_3 = {
    "Not at all": 0,
    "Several days": 1,
    "More than half the days": 2,
    "Nearly every day": 3,
}

ASRS_0_TO_4 = {
    "Never": 0,
    "Rarely": 1,
    "Sometimes": 2,
    "Often": 3,
    "Very often": 4,
}

YANG_1_TO_5 = {
    "Completely Disagree": 1,
    "Somewhat Disagree": 2,
    "Neutral": 3,
    "Somewhat Agree": 4,
    "Completely Agree": 5,
}

SFV_FREQUENCY_0_TO_3 = {
    "Never": 0,
    "A few times per week": 1,
    "Daily": 2,
    "Multiple times per day": 3,
}

# Daily SFV duration: higher = more time.
# Note: "Other (please specify)" cannot be safely ordered without parsing the free text.
# We encode it as NaN by default (see `ENCODE_SFV_DURATION_OTHER_AS_NAN`).
SFV_DURATION_0_TO_3 = {
    "Less than 30 minutes": 0,
    "30 - 60 minutes": 1,
    "30-60 minutes": 1,
    "1 - 2 hours": 2,
    "1-2 hours": 2,
    "2 - 3 hours": 3,
    "2-3 hours": 3,
    "Other": float("nan"),
    "Other (please specify)": float("nan"),
}

PD_YES_NO = {"No": 0, "Yes": 1}

# Behavior switches (kept explicit for paper reproducibility).
ENCODE_SFV_DURATION_OTHER_AS_NAN = True  # If False, "Other" / "Other (please specify)" will be encoded as 4.


# =============================================================================
# Qualtrics CSV helpers (MultiIndex columns)
# =============================================================================


def load_qualtrics_multilevel_csv(path: Path) -> pd.DataFrame:
    """
    Load a Qualtrics export with 3 header rows:
      level 0 = QID (e.g., Q80)
      level 1 = label (e.g., "GAD")
      level 2 = ImportId JSON blob
    """
    df = pd.read_csv(path, header=[0, 1, 2])
    if not isinstance(df.columns, pd.MultiIndex) or df.columns.nlevels != 3:
        raise ValueError(
            f"Expected a 3-level MultiIndex column structure from {path.name}, "
            f"got nlevels={getattr(df.columns, 'nlevels', None)}."
        )
    return df


def col_by_qid(df: pd.DataFrame, qid: str) -> Tuple[str, str, str]:
    """Return the unique column tuple that matches a QID (level 0)."""
    matches = [c for c in df.columns if c[0] == qid]
    if len(matches) == 0:
        raise KeyError(f"QID not found: {qid}")
    if len(matches) > 1:
        raise KeyError(
            f"QID matched multiple columns: {qid}. "
            f"Use label-based selection or refine selection. Matches: {matches[:5]}{'...' if len(matches) > 5 else ''}"
        )
    return matches[0]


def cols_by_label(df: pd.DataFrame, label: str) -> List[Tuple[str, str, str]]:
    """Return all columns whose label (level 1) matches `label`."""
    cols = [c for c in df.columns if c[1] == label]
    if len(cols) == 0:
        raise KeyError(f"Label not found: {label}")
    return cols


def merge_other_responses(primary: pd.Series, other_text: pd.Series, other_token: str) -> pd.Series:
    """Replace `other_token` with the free-text response where available."""
    result = primary.copy()
    return result.where(result != other_token, other_text.where(other_text.notna(), result))


# =============================================================================
# Encoding + scoring
# =============================================================================


def encode_ordinal(series: pd.Series, mapping: Dict[str, float]) -> pd.Series:
    """Encode a string response series using an explicit mapping."""
    return series.replace(mapping)


def report_unmapped(series: pd.Series, mapping: Dict[str, float], name: str) -> List[str]:
    """Return a list of unique non-missing values not present in mapping keys."""
    observed = pd.Series(series.dropna().unique())
    unmapped = sorted([v for v in observed.tolist() if v not in mapping])
    if unmapped:
        print(f"[WARN] Unmapped values for {name}: {unmapped}")
    return unmapped


def compute_scale_scores(
    df: pd.DataFrame,
    item_cols: Sequence[Tuple[str, str, str]],
    mapping: Dict[str, float],
    *,
    min_item_proportion: float = MIN_ITEM_PROPORTION,
) -> pd.DataFrame:
    """
    Compute per-participant composite scores for a multi-item scale.
    These are items like PHQ-9, GAD, ASRS, Yang Problematic Use, and Yang Motivation.

    Returns a DataFrame with:
    - total: sum across items (requires >= ceil(p * n_items) non-missing)
    - mean_answered: mean across answered items (requires same threshold)
    - n_answered: number of answered items
    """
    raw = df.loc[:, item_cols]
    numeric = raw.replace(mapping)

    n_items = numeric.shape[1]
    min_items = ceil(n_items * min_item_proportion)
    n_answered = numeric.notna().sum(axis=1)

    total = numeric.sum(axis=1, min_count=min_items)
    mean_answered = numeric.mean(axis=1)
    mean_answered = mean_answered.where(n_answered >= min_items, float("nan"))

    return pd.DataFrame(
        {
            "total": total,
            "mean_answered": mean_answered,
            "n_answered": n_answered,
            "n_items": n_items,
        }
    )


def write_column_audit(
    outpath: Path,
    selections: Dict[str, Sequence[Tuple[str, str, str]]],
) -> None:
    """Write an explicit audit table of which Qualtrics columns were used for each construct."""
    rows = []
    for construct, cols in selections.items():
        for c in cols:
            rows.append({"construct": construct, "qid": c[0], "label": c[1], "import_id": c[2]})
    pd.DataFrame(rows).to_csv(outpath, index=False)


def main() -> None:
    root = Path(__file__).resolve().parent
    in_path = root / "qualtrics" / "final_SF_demographic_data.csv"
    out_dir = root / "covariate_outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    tabular_dir = root / "data" / "tabular"
    tabular_dir.mkdir(parents=True, exist_ok=True)

    df = load_qualtrics_multilevel_csv(in_path)
    print(f"Loaded {in_path.name}: rows={len(df)}, cols={df.shape[1]}")

    # -------------------------
    # Single-item covariates
    # -------------------------
    subject_id_col = col_by_qid(df, "Q71")
    age_col = col_by_qid(df, "Q9")
    pd_col = col_by_qid(df, "Q16")
    sfv_freq_col = col_by_qid(df, "Q80")
    sfv_dur_col = col_by_qid(df, "Q81")
    sfv_dur_other_col = col_by_qid(df, "Q81_5_TEXT")

    # Study ID (Q71): keep as numeric for downstream merges; exclude from correlations.
    # Use pandas' nullable integer dtype so missing / non-numeric values become <NA>.
    subject_id = (
        pd.to_numeric(df[subject_id_col].astype(str).str.strip(), errors="coerce")
        .astype("Int64")
        .rename("subject_id")
    )

    age = pd.to_numeric(df[age_col], errors="coerce").rename("age")

    report_unmapped(df[pd_col], PD_YES_NO, "pd_status (Q16)")
    pd_status = encode_ordinal(df[pd_col], PD_YES_NO).rename("pd_status")

    report_unmapped(df[sfv_freq_col], SFV_FREQUENCY_0_TO_3, "sfv_frequency (Q80)")
    sfv_frequency = encode_ordinal(df[sfv_freq_col], SFV_FREQUENCY_0_TO_3).rename("sfv_frequency")

    # SFV duration: "Other" is kept for audit; the ordinal encoding is applied only to standard choices.
    other_tokens = {"Other", "Other (please specify)"}
    other_mask = df[sfv_dur_col].isin(other_tokens)
    if other_mask.any():
        audit = pd.DataFrame(
            {
                "sfv_duration_choice": df.loc[other_mask, sfv_dur_col],
                "sfv_duration_other_text": df.loc[other_mask, sfv_dur_other_col],
            }
        )
        audit.to_csv(out_dir / "sfv_duration_other_audit.csv")

    duration_map = dict(SFV_DURATION_0_TO_3)
    if not ENCODE_SFV_DURATION_OTHER_AS_NAN:
        for tok in other_tokens:
            duration_map[tok] = 4

    report_unmapped(df[sfv_dur_col], duration_map, "sfv_daily_duration (Q81)")
    sfv_daily_duration = encode_ordinal(df[sfv_dur_col], duration_map).rename("sfv_daily_duration")

    # -------------------------
    # Multi-item scale scores
    # -------------------------
    phq_items = cols_by_label(df, "PHQ-9")
    gad_items = cols_by_label(df, "GAD")
    asrs_items = cols_by_label(df, "ASRS")
    yang_pu_items = cols_by_label(df, "Yang Problematic Use")
    yang_mot_items = cols_by_label(df, "Yang Motivation")

    phq = compute_scale_scores(df, phq_items, LIKERT_0_TO_3).add_prefix("phq_")
    gad = compute_scale_scores(df, gad_items, LIKERT_0_TO_3).add_prefix("gad_")
    asrs = compute_scale_scores(df, asrs_items, ASRS_0_TO_4).add_prefix("asrs_")
    yang_pu = compute_scale_scores(df, yang_pu_items, YANG_1_TO_5).add_prefix("yang_pu_")
    yang_mot = compute_scale_scores(df, yang_mot_items, YANG_1_TO_5).add_prefix("yang_mot_")

    # Persist a paper-friendly audit of exactly which columns fed each construct.
    write_column_audit(
        out_dir / "covariate_column_audit.csv",
        selections={
            "subject_id": [subject_id_col],
            "age": [age_col],
            "pd_status": [pd_col],
            "sfv_frequency": [sfv_freq_col],
            "sfv_daily_duration": [sfv_dur_col],
            "sfv_daily_duration_other_text": [sfv_dur_other_col],
            "PHQ-9 items": phq_items,
            "GAD items": gad_items,
            "ASRS items": asrs_items,
            "Yang Problematic Use items": yang_pu_items,
            "Yang Motivation items": yang_mot_items,
        },
    )

    covariates = pd.concat(
        [
            age,
            pd_status,
            sfv_frequency,
            sfv_daily_duration,
            phq[["phq_total"]],
            gad[["gad_total"]],
            asrs[["asrs_total"]],
            yang_pu[["yang_pu_total"]],
            yang_mot[["yang_mot_total"]],
        ],
        axis=1,
    )

    # Ensure the correlation inputs are purely numeric (prevents silent object-dtype artifacts).
    covariates_numeric = covariates.apply(pd.to_numeric, errors="coerce")

    # -------------------------
    # Save clean dataset + missingness
    # -------------------------
    covariates_numeric.to_csv(out_dir / "covariates_clean.csv", index=False)

    # Also save a copy with the Qualtrics study ID included (kept out of correlations).
    covariates_with_subject_id = pd.concat([subject_id, covariates_numeric], axis=1)
    covariates_with_subject_id.to_csv(tabular_dir / "socio_demographic_data_processed.csv", index=False)

    missingness = (
        covariates_numeric.isna()
        .mean()
        .rename("missing_prop")
        .to_frame()
        .assign(n_missing=covariates_numeric.isna().sum(), n_total=len(covariates_numeric))
        .sort_values("missing_prop", ascending=False)
    )
    missingness.to_csv(out_dir / "covariate_missingness.csv")

    # -------------------------
    # Lightweight console audit
    # -------------------------
    print("\nScale item counts:")
    print(
        pd.Series(
            {
                "PHQ-9 items": len(phq_items),
                "GAD items": len(gad_items),
                "ASRS items": len(asrs_items),
                "Yang PU items": len(yang_pu_items),
                "Yang Motivation items": len(yang_mot_items),
            }
        ).to_string()
    )

    print("\nMissingness (top):")
    print(missingness.head(10).to_string())

    print(f"\nWrote outputs to: {out_dir}")
    print("Tip: run `covariate_correlation_analysis.py` to generate Spearman correlation diagnostics + heatmap.")


if __name__ == "__main__":
    main()