"""Exploratory channel-behavior association screen for the SFV study.

This script evaluates associations between channel beta columns in a merged
subject-level table and all non-identifier behavioral variables.

Method choices for this script are intentionally explicit and citation-backed.
See `CITATIONS.md` for the full reasoning trail and source registry.

- Yucel et al. (2021): transparent handling of fNIRS quality/missingness.
- Spearman (1904): rank correlation for monotonic associations in bounded or
  ordinal-like variables.
- Bonett (2020): point-biserial correlation for dichotomous-vs-continuous
  associations.
- Benjamini and Hochberg (1995): false-discovery-rate control across the full
  screened family of tests.
- Kriegeskorte et al. (2009): exploratory interpretation caution.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import pointbiserialr, spearmanr


CHANNEL_PATTERN = re.compile(
    r"^S(?P<source>\d{2})_D(?P<detector>\d{2})_Cond(?P<condition>\d{2})_(?P<chromophore>HbO|HbR)$"
)

CONDITION_LABELS = {
    "01": "short_form_education",
    "02": "short_form_entertainment",
    "03": "long_form_entertainment",
    "04": "long_form_education",
}

BEHAVIOR_TO_CONDITIONS = {
    "diff_short_form_education": {"short_form_education"},
    "diff_short_form_entertainment": {"short_form_entertainment"},
    "diff_long_form_education": {"long_form_education"},
    "diff_long_form_entertainment": {"long_form_entertainment"},
    "sf_education_engagement": {"short_form_education"},
    "sf_entertainment_engagement": {"short_form_entertainment"},
    "lf_education_engagement": {"long_form_education"},
    "lf_entertainment_engagement": {"long_form_entertainment"},
    "short_form_engagement": {"short_form_education", "short_form_entertainment"},
    "long_form_engagement": {"long_form_education", "long_form_entertainment"},
    "education_engagement": {"short_form_education", "long_form_education"},
    "entertainment_engagement": {"short_form_entertainment", "long_form_entertainment"},
}

DEFAULT_IDENTIFIER_COLUMNS = ("subject_id", "homer_subject")


@dataclass(frozen=True)
class BehaviorVariable:
    """Metadata for a single behavioral variable."""

    name: str
    method: str
    n_missing: int
    n_unique_non_missing: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Inspect exploratory relationships between channel beta columns "
            "(`S##_D##_Cond##_HbO/HbR`) and behavioral variables in a merged CSV."
        )
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=Path("data/tabular/homer3_betas_plus_combined_sfv_data_inner_join.csv"),
        help="Merged CSV containing channel columns and behavioral columns.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("data/results/channel_behavior_relationships"),
        help="Directory for tables, metadata, and the markdown summary.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="FDR threshold used for the global Benjamini-Hochberg correction.",
    )
    parser.add_argument(
        "--exclude-column",
        action="append",
        default=[],
        help="Additional non-behavior columns to exclude from the screen.",
    )
    return parser.parse_args()


def load_dataset(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input CSV does not exist: {path}")
    return pd.read_csv(path)


def identify_channel_columns(columns: Iterable[str]) -> list[str]:
    channel_columns = [column for column in columns if CHANNEL_PATTERN.match(column)]
    if not channel_columns:
        raise ValueError("No channel columns matched the expected S##_D##_Cond##_HbO/HbR schema.")
    return channel_columns


def parse_channel_metadata(channel_name: str) -> dict[str, str]:
    match = CHANNEL_PATTERN.match(channel_name)
    if match is None:
        raise ValueError(f"Invalid channel column name: {channel_name}")
    metadata = match.groupdict()
    metadata["condition_label"] = CONDITION_LABELS.get(metadata["condition"], "unknown_condition")
    return metadata


def sanitize_channel_values(df: pd.DataFrame, channel_columns: list[str]) -> pd.DataFrame:
    """Treat channel 0/NaN placeholders as missing/pruned observations only."""

    sanitized = df.copy()
    sanitized.loc[:, channel_columns] = sanitized.loc[:, channel_columns].replace(0, np.nan)
    return sanitized


def classify_behavior_variables(
    df: pd.DataFrame,
    channel_columns: list[str],
    identifier_columns: Iterable[str],
) -> list[BehaviorVariable]:
    identifier_set = set(identifier_columns)
    behavior_candidates = [column for column in df.columns if column not in channel_columns and column not in identifier_set]
    behavior_variables: list[BehaviorVariable] = []

    for column in behavior_candidates:
        series = df[column]
        if not pd.api.types.is_numeric_dtype(series):
            raise ValueError(
                f"Non-numeric non-identifier behavioral column encountered: {column}. "
                "Exclude it explicitly or encode it numerically before analysis."
            )
        non_missing = series.dropna()
        n_unique = int(non_missing.nunique())
        if n_unique < 2:
            continue
        unique_values = set(non_missing.astype(float).unique().tolist())
        method = "point_biserial" if unique_values.issubset({0.0, 1.0}) else "spearman"
        behavior_variables.append(
            BehaviorVariable(
                name=column,
                method=method,
                n_missing=int(series.isna().sum()),
                n_unique_non_missing=n_unique,
            )
        )

    if not behavior_variables:
        raise ValueError("No analyzable behavioral variables remained after exclusions.")

    return behavior_variables


def summarize_channel_missingness(raw_df: pd.DataFrame, sanitized_df: pd.DataFrame, channel_columns: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for column in channel_columns:
        metadata = parse_channel_metadata(column)
        raw = raw_df[column]
        sanitized = sanitized_df[column]
        rows.append(
            {
                "channel": column,
                **metadata,
                "raw_zero_count": int(raw.eq(0).sum()),
                "raw_nan_count": int(raw.isna().sum()),
                "usable_n": int(sanitized.notna().sum()),
                "pruned_or_missing_n": int(sanitized.isna().sum()),
            }
        )
    summary = pd.DataFrame(rows).sort_values(["condition", "source", "detector", "chromophore"]).reset_index(drop=True)
    summary["pruned_or_missing_fraction"] = summary["pruned_or_missing_n"] / (summary["usable_n"] + summary["pruned_or_missing_n"])
    return summary


def summarize_behavior_variables(df: pd.DataFrame, behavior_variables: list[BehaviorVariable]) -> pd.DataFrame:
    rows = []
    for variable in behavior_variables:
        series = df[variable.name]
        rows.append(
            {
                "behavior": variable.name,
                "method": variable.method,
                "n_missing": variable.n_missing,
                "usable_n": int(series.notna().sum()),
                "n_unique_non_missing": variable.n_unique_non_missing,
            }
        )
    return pd.DataFrame(rows).sort_values("behavior").reset_index(drop=True)


def compute_association(channel: pd.Series, behavior: pd.Series, method: str) -> dict[str, object]:
    aligned = pd.concat([channel, behavior], axis=1).dropna()
    aligned.columns = ["channel", "behavior"]
    n_pairwise = int(len(aligned))
    if n_pairwise < 3:
        return {"effect_size": np.nan, "p_uncorrected": np.nan, "n_pairwise": n_pairwise, "skip_reason": "fewer_than_three_pairwise_observations"}
    if aligned["channel"].nunique() < 2 or aligned["behavior"].nunique() < 2:
        return {"effect_size": np.nan, "p_uncorrected": np.nan, "n_pairwise": n_pairwise, "skip_reason": "insufficient_variation"}

    if method == "spearman":
        effect_size, p_value = spearmanr(aligned["channel"], aligned["behavior"])
        return {"effect_size": float(effect_size), "p_uncorrected": float(p_value), "n_pairwise": n_pairwise, "skip_reason": ""}

    if method == "point_biserial":
        binary_values = set(aligned["behavior"].astype(float).unique().tolist())
        if not binary_values.issubset({0.0, 1.0}):
            return {"effect_size": np.nan, "p_uncorrected": np.nan, "n_pairwise": n_pairwise, "skip_reason": "behavior_not_binary"}
        effect_size, p_value = pointbiserialr(aligned["behavior"], aligned["channel"])
        return {"effect_size": float(effect_size), "p_uncorrected": float(p_value), "n_pairwise": n_pairwise, "skip_reason": ""}

    raise ValueError(f"Unsupported method: {method}")


def compute_pairwise_associations(
    df: pd.DataFrame,
    channel_columns: list[str],
    behavior_variables: list[BehaviorVariable],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for channel_name in channel_columns:
        channel_metadata = parse_channel_metadata(channel_name)
        channel_series = df[channel_name]
        for behavior in behavior_variables:
            stats = compute_association(channel_series, df[behavior.name], method=behavior.method)
            row = {
                "channel": channel_name,
                **channel_metadata,
                "behavior": behavior.name,
                "method": behavior.method,
                "effect_size": stats["effect_size"],
                "abs_effect_size": abs(stats["effect_size"]) if pd.notna(stats["effect_size"]) else np.nan,
                "p_uncorrected": stats["p_uncorrected"],
                "n_pairwise": stats["n_pairwise"],
                "skip_reason": stats["skip_reason"],
            }
            rows.append(row)
    results = pd.DataFrame(rows)
    if results.empty:
        raise ValueError("No pairwise association rows were produced.")
    return apply_bh_fdr(results)


def benjamini_hochberg_adjust(p_values: pd.Series) -> np.ndarray:
    """Return BH-adjusted p-values in the original row order."""

    p = p_values.to_numpy(dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    n_tests = len(ranked)
    adjusted_ranked = ranked * n_tests / np.arange(1, n_tests + 1)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted_ranked = np.clip(adjusted_ranked, 0.0, 1.0)
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = adjusted_ranked
    return adjusted


def apply_bh_fdr(results: pd.DataFrame) -> pd.DataFrame:
    adjusted = results.copy()
    valid_mask = adjusted["p_uncorrected"].notna()
    adjusted["p_fdr_bh"] = np.nan
    adjusted["fdr_reject_alpha_0_05"] = False
    if valid_mask.any():
        adjusted_p = benjamini_hochberg_adjust(adjusted.loc[valid_mask, "p_uncorrected"])
        adjusted.loc[valid_mask, "p_fdr_bh"] = adjusted_p
        adjusted.loc[valid_mask, "fdr_reject_alpha_0_05"] = adjusted_p <= 0.05
    return adjusted


def build_behavior_summary(results: pd.DataFrame) -> pd.DataFrame:
    valid = results.loc[results["p_uncorrected"].notna()].copy()
    if valid.empty:
        return pd.DataFrame(
            columns=[
                "behavior",
                "method",
                "tested_channel_count",
                "best_channel",
                "best_effect_size",
                "best_abs_effect_size",
                "best_p_uncorrected",
                "best_p_fdr_bh",
                "n_fdr_significant",
            ]
        )

    best_hits = (
        valid.sort_values(["behavior", "p_fdr_bh", "p_uncorrected", "abs_effect_size"], ascending=[True, True, True, False])
        .groupby("behavior", as_index=False)
        .first()
    )

    counts = (
        valid.groupby(["behavior", "method"], as_index=False)
        .agg(
            tested_channel_count=("channel", "count"),
            n_fdr_significant=("fdr_reject_alpha_0_05", "sum"),
        )
    )

    merged = counts.merge(
        best_hits[
            [
                "behavior",
                "channel",
                "effect_size",
                "abs_effect_size",
                "p_uncorrected",
                "p_fdr_bh",
            ]
        ],
        on="behavior",
        how="left",
    )

    return merged.rename(
        columns={
            "channel": "best_channel",
            "effect_size": "best_effect_size",
            "abs_effect_size": "best_abs_effect_size",
            "p_uncorrected": "best_p_uncorrected",
            "p_fdr_bh": "best_p_fdr_bh",
        }
    ).sort_values(["best_p_fdr_bh", "best_p_uncorrected", "behavior"]).reset_index(drop=True)


def build_condition_matched_results(results: pd.DataFrame) -> pd.DataFrame:
    def is_condition_matched(row: pd.Series) -> bool:
        allowed_conditions = BEHAVIOR_TO_CONDITIONS.get(row["behavior"])
        return allowed_conditions is not None and row["condition_label"] in allowed_conditions

    matched = results.loc[results["p_uncorrected"].notna()].copy()
    matched = matched.loc[matched.apply(is_condition_matched, axis=1)]
    return matched.sort_values(["p_uncorrected", "abs_effect_size"], ascending=[True, False]).reset_index(drop=True)


def render_markdown_table(df: pd.DataFrame, columns: list[str]) -> list[str]:
    subset = df.loc[:, columns].copy()
    formatted_rows = []
    for row in subset.itertuples(index=False, name=None):
        formatted = []
        for value in row:
            if isinstance(value, float):
                if math.isnan(value):
                    formatted.append("NA")
                else:
                    formatted.append(f"{value:.4f}")
            elif isinstance(value, (np.floating,)):
                if np.isnan(value):
                    formatted.append("NA")
                else:
                    formatted.append(f"{float(value):.4f}")
            else:
                formatted.append(str(value))
        formatted_rows.append(formatted)

    header = "| " + " | ".join(columns) + " |"
    divider = "| " + " | ".join(["---"] * len(columns)) + " |"
    body = ["| " + " | ".join(row) + " |" for row in formatted_rows]
    return [header, divider, *body]


def write_summary_markdown(
    out_path: Path,
    *,
    input_csv: Path,
    raw_df: pd.DataFrame,
    channel_columns: list[str],
    behavior_variables: list[BehaviorVariable],
    results: pd.DataFrame,
    channel_missingness: pd.DataFrame,
    behavior_summary: pd.DataFrame,
    condition_matched_results: pd.DataFrame,
) -> None:
    tested = results.loc[results["p_uncorrected"].notna()].copy()
    significant = tested.loc[tested["fdr_reject_alpha_0_05"]].copy().sort_values(["p_fdr_bh", "p_uncorrected"])
    top_nominal = tested.sort_values(["p_uncorrected", "abs_effect_size"], ascending=[True, False]).head(10)
    most_pruned = channel_missingness.sort_values(["pruned_or_missing_fraction", "pruned_or_missing_n"], ascending=[False, False]).head(10)

    lines = [
        "# Channel-behavior relationship screen",
        "",
        f"- Input CSV: `{input_csv}`",
        f"- Subjects (rows): `{len(raw_df)}`",
        f"- Channel columns analyzed: `{len(channel_columns)}`",
        f"- Behavioral variables analyzed: `{len(behavior_variables)}`",
        f"- Valid pairwise tests: `{len(tested)}`",
        f"- FDR-significant tests at q <= 0.05: `{len(significant)}`",
        "",
        "## Method summary",
        "",
        "- Channel values equal to `0` were re-coded to missing only for channel columns because this dataset encodes pruned channels as `0`/`NaN`; no imputation was performed.",
        "- Continuous and ordinal-like behavioral variables were screened with Spearman rank correlation.",
        "- Binary behavioral variables were screened with point-biserial correlation.",
        "- Benjamini-Hochberg FDR was applied across the full family of valid tests.",
        "- Results are exploratory and should not be interpreted as confirmatory causal evidence.",
        "",
        "## Top nominal associations",
        "",
    ]
    if top_nominal.empty:
        lines.append("No valid channel-behavior tests were available.")
    else:
        lines.extend(
            render_markdown_table(
                top_nominal,
                ["behavior", "channel", "method", "effect_size", "p_uncorrected", "p_fdr_bh", "n_pairwise"],
            )
        )

    lines.extend(["", "## FDR-significant associations", ""])
    if significant.empty:
        lines.append("No associations survived global Benjamini-Hochberg correction at q <= 0.05.")
    else:
        lines.extend(
            render_markdown_table(
                significant.head(25),
                ["behavior", "channel", "method", "effect_size", "p_uncorrected", "p_fdr_bh", "n_pairwise"],
            )
        )

    lines.extend(["", "## Top condition-matched nominal associations", ""])
    if condition_matched_results.empty:
        lines.append("No condition-matched behavior-to-channel rows were defined for this dataset.")
    else:
        lines.extend(
            render_markdown_table(
                condition_matched_results.head(15),
                ["behavior", "channel", "condition_label", "effect_size", "p_uncorrected", "p_fdr_bh", "n_pairwise"],
            )
        )

    lines.extend(["", "## Behavioral summary", ""])
    if behavior_summary.empty:
        lines.append("No behavior-level summary was available.")
    else:
        lines.extend(
            render_markdown_table(
                behavior_summary.head(25),
                [
                    "behavior",
                    "method",
                    "tested_channel_count",
                    "best_channel",
                    "best_effect_size",
                    "best_p_uncorrected",
                    "best_p_fdr_bh",
                    "n_fdr_significant",
                ],
            )
        )

    lines.extend(["", "## Most-pruned channels", ""])
    lines.extend(
        render_markdown_table(
            most_pruned,
            ["channel", "condition_label", "chromophore", "usable_n", "pruned_or_missing_n", "pruned_or_missing_fraction"],
        )
    )

    out_path.write_text("\n".join(lines) + "\n")


def run_analysis(input_csv: Path, out_dir: Path, alpha: float, excluded_columns: Iterable[str]) -> dict[str, Path]:
    raw_df = load_dataset(input_csv)
    channel_columns = identify_channel_columns(raw_df.columns)
    sanitized_df = sanitize_channel_values(raw_df, channel_columns)
    identifier_columns = list(DEFAULT_IDENTIFIER_COLUMNS) + list(excluded_columns)
    behavior_variables = classify_behavior_variables(sanitized_df, channel_columns, identifier_columns)
    channel_missingness = summarize_channel_missingness(raw_df, sanitized_df, channel_columns)
    behavior_variable_summary = summarize_behavior_variables(sanitized_df, behavior_variables)
    results = compute_pairwise_associations(sanitized_df, channel_columns, behavior_variables)
    results["fdr_reject_alpha_0_05"] = results["p_fdr_bh"] <= alpha
    results.loc[results["p_fdr_bh"].isna(), "fdr_reject_alpha_0_05"] = False
    behavior_summary = build_behavior_summary(results)
    condition_matched_results = build_condition_matched_results(results)
    top_hits = results.loc[results["p_uncorrected"].notna()].sort_values(
        ["p_fdr_bh", "p_uncorrected", "abs_effect_size"], ascending=[True, True, False]
    ).head(100)

    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "pairwise_results": out_dir / "channel_behavior_pairwise_results.csv",
        "top_hits": out_dir / "channel_behavior_top_hits.csv",
        "condition_matched_hits": out_dir / "channel_behavior_condition_matched_top_hits.csv",
        "behavior_summary": out_dir / "channel_behavior_behavior_summary.csv",
        "behavior_variables": out_dir / "behavior_variable_profile.csv",
        "channel_missingness": out_dir / "channel_missingness_summary.csv",
        "metadata": out_dir / "analysis_metadata.json",
        "summary_markdown": out_dir / "channel_behavior_summary.md",
    }

    results.to_csv(outputs["pairwise_results"], index=False)
    top_hits.to_csv(outputs["top_hits"], index=False)
    condition_matched_results.to_csv(outputs["condition_matched_hits"], index=False)
    behavior_summary.to_csv(outputs["behavior_summary"], index=False)
    behavior_variable_summary.to_csv(outputs["behavior_variables"], index=False)
    channel_missingness.to_csv(outputs["channel_missingness"], index=False)

    metadata_payload = {
        "input_csv": str(input_csv.resolve()),
        "n_subjects": int(len(raw_df)),
        "n_channel_columns": len(channel_columns),
        "n_behavior_variables": len(behavior_variables),
        "n_valid_tests": int(results["p_uncorrected"].notna().sum()),
        "n_fdr_significant": int(results["fdr_reject_alpha_0_05"].sum()),
        "n_condition_matched_tests": int(len(condition_matched_results)),
        "identifier_columns_excluded": identifier_columns,
        "behavior_variables": [
            {
                "name": variable.name,
                "method": variable.method,
                "n_missing": variable.n_missing,
                "n_unique_non_missing": variable.n_unique_non_missing,
            }
            for variable in behavior_variables
        ],
    }
    outputs["metadata"].write_text(json.dumps(metadata_payload, indent=2))

    write_summary_markdown(
        outputs["summary_markdown"],
        input_csv=input_csv,
        raw_df=raw_df,
        channel_columns=channel_columns,
        behavior_variables=behavior_variables,
        results=results,
        channel_missingness=channel_missingness,
        behavior_summary=behavior_summary,
        condition_matched_results=condition_matched_results,
    )

    return outputs


def main() -> None:
    args = parse_args()
    outputs = run_analysis(
        input_csv=args.input_csv,
        out_dir=args.out_dir,
        alpha=args.alpha,
        excluded_columns=args.exclude_column,
    )
    for name, path in outputs.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
