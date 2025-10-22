"""
Group-level fNIRS analysis using a two-stage 2x2 repeated-measures pipeline.

Stage 1 (global effects)
- Collapse run-level betas to subject-by-channel condition means.
- Average condition means across channels for each subject.
- Test the three within-subject contrasts (Format, Content, Interaction) with
  one-sample t-tests across subjects.
- Adjust the three p-values using the Benjamini-Hochberg method (FDR).

Stage 2 (channel localisation)
- Only for Stage-1 contrasts that survive FDR correction.
- Re-test the same contrasts per channel (paired within-subject; one-sample t).
- Apply Benjamini-Hochberg FDR across channels at alpha' = alpha * (R / m),
  following Benjamini-Bogomolov selective FDR control.

Outputs
- `group_<chroma>_global_effects.csv`: Stage-1 statistics per contrast.
- `group_<chroma>_channel_effects.csv`: Stage-2 channel statistics for
  selected contrasts (empty if no Stage-1 effects survive).

Assumptions
- Balanced 2x2 within-subject design with four conditions present for subjects
  contributing to each channel (handled by complete-case filtering).
- Subject-level betas originate from the GLM in `fnirs_analysis.py`.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import numpy as np
import pandas as pd
from scipy import stats as scistats
from statsmodels.stats.multitest import fdrcorrection, multipletests


# -----------------------
# Configuration defaults
# -----------------------

DEFAULT_COMBINED_PATH = Path("./glm_results/combined_glm_long_runs.csv")
DEFAULT_CHROMA = "hbo"  # Can be "hbo" or "hbr"
DEFAULT_ALPHA = 0.05
MIN_SUBJECTS_PER_CHANNEL = 6  # drop channels with fewer than this many complete subjects
BETA_VALUE_SCALER = 1

EFFECT_KEYS = ("format", "content", "interaction")


@dataclass
class StageOneResult:
    effect_key: str
    effect: str
    contrast_definition: str
    n_subjects: int
    df: float
    t_value: float
    mean_difference: float
    p_value: float
    q_value: float | None = None
    reject: bool = False
    cohen_dz: float | None = None


@dataclass
class StageTwoResult:
    effect_key: str
    effect: str
    contrast_definition: str
    channel: str
    chroma: str
    n_subjects: int
    df: float
    t_value: float
    mean_difference: float
    p_value: float
    q_value: float | None = None
    reject: bool = False
    cohen_dz: float | None = None
    alpha_prime: float | None = None


def _map_condition_to_factors(condition: str) -> tuple[str, str]:
    """Infer Format (Short/Long) and Content (Education/Entertainment) from condition label."""
    cond_lower = condition.lower()
    if "short" in cond_lower:
        format_level = "Short-Form"
    elif "long" in cond_lower:
        format_level = "Long-Form"
    else:
        raise ValueError(f"Could not determine format factor from condition '{condition}'.")

    if "education" in cond_lower:
        content_level = "Education"
    elif "entertainment" in cond_lower:
        content_level = "Entertainment"
    else:
        raise ValueError(f"Could not determine content factor from condition '{condition}'.")
    return format_level, content_level


def _infer_factor_levels(conditions: Sequence[str]) -> tuple[List[str], List[str]]:
    """Determine the unique levels for each within-subject factor."""
    format_levels: List[str] = []
    content_levels: List[str] = []
    for cond in conditions:
        fmt, cont = _map_condition_to_factors(cond)
        if fmt not in format_levels:
            format_levels.append(fmt)
        if cont not in content_levels:
            content_levels.append(cont)
    if len(format_levels) != 2 or len(content_levels) != 2:
        raise ValueError(
            "Two-stage 2x2 design requires exactly two levels per factor. "
            f"Found format levels={format_levels} and content levels={content_levels}."
        )
    return format_levels, content_levels


def _aggregate_subject_channel(
    df_runs: pd.DataFrame,
    conditions: Sequence[str],
    format_levels: Sequence[str],
    content_levels: Sequence[str],
) -> pd.DataFrame:
    """Aggregate run-level betas to subject-channel condition means, keeping complete cases."""
    if df_runs.empty:
        return pd.DataFrame(
            columns=["subject_id", "ch_name", "condition", "theta", "format_factor", "content_factor"]
        )

    aggregated = (
        df_runs.groupby(["subject_id", "ch_name", "condition"], as_index=False)
        .agg(theta=("theta", "mean"))
    )
    aggregated = aggregated[aggregated["condition"].isin(list(conditions))].copy()

    if aggregated.empty:
        return pd.DataFrame(
            columns=["subject_id", "ch_name", "condition", "theta", "format_factor", "content_factor"]
        )

    format_vals, content_vals = zip(*aggregated["condition"].map(_map_condition_to_factors))
    aggregated["format_factor"] = pd.Categorical(format_vals, categories=list(format_levels), ordered=True)
    aggregated["content_factor"] = pd.Categorical(content_vals, categories=list(content_levels), ordered=True)

    counts = (
        aggregated.groupby(["subject_id", "ch_name"])["condition"]
        .nunique()
        .reset_index(name="n_conditions")
    )
    complete_pairs = counts[counts["n_conditions"] == len(conditions)][["subject_id", "ch_name"]]
    if complete_pairs.empty:
        return pd.DataFrame(
            columns=["subject_id", "ch_name", "condition", "theta", "format_factor", "content_factor"]
        )

    aggregated = aggregated.merge(complete_pairs, on=["subject_id", "ch_name"], how="inner")
    aggregated.sort_values(["ch_name", "subject_id", "condition"], inplace=True)
    aggregated.reset_index(drop=True, inplace=True)
    return aggregated


def _pivot_subject_channel(aggregated: pd.DataFrame, conditions: Sequence[str]) -> pd.DataFrame:
    """Pivot subject-channel condition means to wide format."""
    if aggregated.empty:
        return pd.DataFrame(columns=["subject_id", "ch_name", *conditions])
    wide = aggregated.pivot_table(
        index=["subject_id", "ch_name"],
        columns="condition",
        values="theta",
        aggfunc="mean",
    )
    # Ensure column order and drop incomplete rows
    wide = wide.reindex(columns=list(conditions))
    wide = wide.dropna()
    if wide.empty:
        return pd.DataFrame(columns=["subject_id", "ch_name", *conditions])
    wide = wide.reset_index()
    return wide


def _build_contrast_weights(
    conditions: Sequence[str],
    format_levels: Sequence[str],
    content_levels: Sequence[str],
    condition_to_format: Dict[str, str],
    condition_to_content: Dict[str, str],
) -> dict[str, np.ndarray]:
    """Create linear contrast weights for each effect."""
    format_codes = {format_levels[0]: -1.0, format_levels[1]: 1.0}
    content_codes = {content_levels[0]: 1.0, content_levels[1]: -1.0}

    weights: dict[str, np.ndarray] = {}
    for effect in EFFECT_KEYS:
        vec: List[float] = []
        for cond in conditions:
            fmt = condition_to_format[cond]
            cont = condition_to_content[cond]
            if effect == "format":
                weight = format_codes[fmt] / 2.0
            elif effect == "content":
                weight = content_codes[cont] / 2.0
            else:  # interaction
                weight = (format_codes[fmt] * content_codes[cont]) / 2.0
            vec.append(float(weight))
        weights[effect] = np.array(vec, dtype=float)
    return weights


def _make_contrast_definitions(
    format_levels: Sequence[str], content_levels: Sequence[str]
) -> dict[str, str]:
    """Generate human-readable descriptions of each contrast."""
    return {
        "format": f"{format_levels[1]} - {format_levels[0]} (averaged over {', '.join(content_levels)})",
        "content": f"{content_levels[0]} - {content_levels[1]} (averaged over {', '.join(format_levels)})",
        "interaction": (
            f"[{format_levels[1]} {content_levels[0]} - {format_levels[1]} {content_levels[1]}] "
            f"- [{format_levels[0]} {content_levels[0]} - {format_levels[0]} {content_levels[1]}]"
        ),
    }


def _run_stage_one(
    subject_condition_means: pd.DataFrame,
    conditions: Sequence[str],
    weights: dict[str, np.ndarray],
    effect_labels: dict[str, str],
    contrast_definitions: dict[str, str],
    alpha: float,
    min_subjects: int,
) -> List[StageOneResult]:
    """Compute Stage-1 global one-sample t-tests across subjects."""
    records: List[StageOneResult] = []
    if subject_condition_means.empty:
        for effect_key in EFFECT_KEYS:
            records.append(
                StageOneResult(
                    effect_key=effect_key,
                    effect=effect_labels[effect_key],
                    contrast_definition=contrast_definitions[effect_key],
                    n_subjects=0,
                    df=float("nan"),
                    t_value=float("nan"),
                    mean_difference=float("nan"),
                    p_value=float("nan"),
                )
            )
        return records

    condition_values = subject_condition_means[list(conditions)].to_numpy(dtype=float)

    for effect_key in EFFECT_KEYS:
        weight_vec = weights[effect_key]
        contrasts = condition_values @ weight_vec
        series = pd.Series(contrasts, index=subject_condition_means.index)
        n = int(series.notna().sum())

        if n >= min_subjects:
            values = series.dropna().to_numpy(dtype=float)
            t_stat, p_val = scistats.ttest_1samp(values, popmean=0.0)
            mean_diff = float(np.mean(values))
            std = float(np.std(values, ddof=1)) if n > 1 else float("nan")
            cohen_dz = float(mean_diff / std) if n > 1 and std > 0 else float("nan")
            df_val = float(n - 1)
        else:
            t_stat = p_val = mean_diff = cohen_dz = float("nan")
            df_val = float("nan")

        records.append(
            StageOneResult(
                effect_key=effect_key,
                effect=effect_labels[effect_key],
                contrast_definition=contrast_definitions[effect_key],
                n_subjects=n,
                df=df_val,
                t_value=float(t_stat) if np.isfinite(t_stat) else float("nan"),
                mean_difference=mean_diff,
                p_value=float(p_val) if np.isfinite(p_val) else float("nan"),
                cohen_dz=cohen_dz,
            )
        )

    finite_idx = [idx for idx, rec in enumerate(records) if np.isfinite(rec.p_value)]
    if finite_idx:
        raw_p = [records[idx].p_value for idx in finite_idx]
        reject, q_vals, _, _ = multipletests(raw_p, alpha=alpha, method="fdr_bh")
        for dest, q_val, rej in zip(finite_idx, q_vals, reject):
            records[dest].q_value = float(q_val)
            records[dest].reject = bool(rej)

    return records


def _run_stage_two(
    wide_df: pd.DataFrame,
    conditions: Sequence[str],
    weights: dict[str, np.ndarray],
    effect_labels: dict[str, str],
    contrast_definitions: dict[str, str],
    chroma: str,
    selected_effects: Iterable[str],
    alpha_prime: float,
) -> List[StageTwoResult]:
    """Compute Stage-2 channel-level contrasts with FDR control."""
    records: List[StageTwoResult] = []
    if wide_df.empty or not selected_effects:
        return records

    condition_values = wide_df[list(conditions)].to_numpy(dtype=float)

    for effect_key in selected_effects:
        weight_vec = weights[effect_key]
        contrasts = condition_values @ weight_vec
        effect_df = pd.DataFrame(
            {
                "subject_id": wide_df["subject_id"],
                "ch_name": wide_df["ch_name"],
                "contrast": contrasts,
            }
        )

        channel_records: List[StageTwoResult] = []
        for channel, sub_df in effect_df.groupby("ch_name"):
            values = sub_df["contrast"].dropna().to_numpy(dtype=float)
            n = len(values)
            if n == 0:
                t_stat = p_val = mean_diff = cohen_dz = float("nan")
                df_val = float("nan")
            else:
                t_stat, p_val = scistats.ttest_1samp(values, popmean=0.0)
                mean_diff = float(np.mean(values))
                std = float(np.std(values, ddof=1)) if n > 1 else float("nan")
                cohen_dz = float(mean_diff / std) if n > 1 and std > 0 else float("nan")
                df_val = float(n - 1)

            channel_records.append(
                StageTwoResult(
                    effect_key=effect_key,
                    effect=effect_labels[effect_key],
                    contrast_definition=contrast_definitions[effect_key],
                    channel=channel,
                    chroma=chroma,
                    n_subjects=n,
                    df=df_val,
                    t_value=float(t_stat) if np.isfinite(t_stat) else float("nan"),
                    mean_difference=mean_diff,
                    p_value=float(p_val) if np.isfinite(p_val) else float("nan"),
                    cohen_dz=cohen_dz,
                    alpha_prime=alpha_prime if alpha_prime > 0 else None,
                )
            )

        finite_idx = [idx for idx, rec in enumerate(channel_records) if np.isfinite(rec.p_value)]
        if finite_idx and alpha_prime > 0:
            raw_p = [channel_records[idx].p_value for idx in finite_idx]
            reject, q_vals = fdrcorrection(raw_p, alpha=alpha_prime)
            for dest, q_val, rej in zip(finite_idx, q_vals, reject):
                channel_records[dest].q_value = float(q_val)
                channel_records[dest].reject = bool(rej)

        records.extend(channel_records)

    return records


def run_group_analysis(
    combined_runs_path: Path,
    chroma: str,
    alpha: float,
    conditions: Sequence[str] | None = None,
    min_subjects: int = MIN_SUBJECTS_PER_CHANNEL,
) -> tuple[pd.DataFrame, pd.DataFrame, float]:
    """Run the two-stage 2x2 repeated-measures pipeline."""
    df = pd.read_csv(combined_runs_path)
    df["theta"] = df["theta"] * BETA_VALUE_SCALER

    required_cols = ("chroma", "condition", "ch_name", "subject_id", "run")
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        missing_fmt = ", ".join(missing_cols)
        raise ValueError(f"Expected column(s) {missing_fmt} in {combined_runs_path}")

    df["chroma"] = df["chroma"].str.lower()
    df = df[df["chroma"] == chroma.lower()].copy()
    if df.empty:
        raise ValueError(f"No rows found for chroma '{chroma}' in {combined_runs_path}.")

    if conditions is None:
        condition_list = df["condition"].dropna().unique().tolist()
    else:
        condition_list = list(dict.fromkeys(conditions))
        df = df[df["condition"].isin(condition_list)]

    if len(condition_list) != 4:
        raise ValueError(
            "Two-stage repeated-measures analysis expects exactly four unique conditions for the 2x2 design. "
            f"Received {len(condition_list)}: {condition_list}"
        )

    format_levels, content_levels = _infer_factor_levels(condition_list)
    condition_to_format = {cond: _map_condition_to_factors(cond)[0] for cond in condition_list}
    condition_to_content = {cond: _map_condition_to_factors(cond)[1] for cond in condition_list}

    aggregated = _aggregate_subject_channel(
        df_runs=df,
        conditions=condition_list,
        format_levels=format_levels,
        content_levels=content_levels,
    )
    if aggregated.empty:
        raise ValueError(
            "No subject-channel pairs contained all four conditions after preprocessing; "
            "cannot run the repeated-measures analysis."
        )

    wide = _pivot_subject_channel(aggregated, condition_list)
    if wide.empty:
        raise ValueError(
            "All subject-channel combinations were incomplete after pivoting; "
            "cannot compute condition means."
        )

    channel_counts = wide.groupby("ch_name")["subject_id"].nunique()
    valid_channels = channel_counts[channel_counts >= min_subjects].index.tolist()
    if not valid_channels:
        raise ValueError(
            f"No channels retained at least {min_subjects} subjects with complete data; "
            "analysis aborted."
        )
    wide = wide[wide["ch_name"].isin(valid_channels)].reset_index(drop=True)

    subject_condition_means = wide.groupby("subject_id")[condition_list].mean()
    if subject_condition_means.shape[0] < min_subjects:
        raise ValueError(
            f"Fewer than {min_subjects} subjects retained after channel filtering; "
            "cannot run Stage-1 tests."
        )

    weights = _build_contrast_weights(
        conditions=condition_list,
        format_levels=format_levels,
        content_levels=content_levels,
        condition_to_format=condition_to_format,
        condition_to_content=condition_to_content,
    )
    effect_labels = {
        "format": "Length",
        "content": "Content",
        "interaction": "Length:Content",
    }
    contrast_definitions = _make_contrast_definitions(format_levels, content_levels)

    stage1_results = _run_stage_one(
        subject_condition_means=subject_condition_means,
        conditions=condition_list,
        weights=weights,
        effect_labels=effect_labels,
        contrast_definitions=contrast_definitions,
        alpha=alpha,
        min_subjects=min_subjects,
    )
    stage1_df = pd.DataFrame([r.__dict__ for r in stage1_results])

    valid_stage1 = [res for res in stage1_results if np.isfinite(res.p_value)]
    m_total = len(valid_stage1)
    selected_effects = [res.effect_key for res in valid_stage1 if res.reject]

    if selected_effects and m_total > 0:
        alpha_prime = alpha * (len(selected_effects) / m_total)
        stage2_records = _run_stage_two(
            wide_df=wide,
            conditions=condition_list,
            weights=weights,
            effect_labels=effect_labels,
            contrast_definitions=contrast_definitions,
            chroma=chroma,
            selected_effects=selected_effects,
            alpha_prime=alpha_prime,
        )
    else:
        alpha_prime = 0.0
        stage2_records = []

    stage2_df = pd.DataFrame([r.__dict__ for r in stage2_records])

    out_dir = combined_runs_path.parent
    stage1_path = out_dir / f"group_{chroma}_global_effects.csv"
    stage2_path = out_dir / f"group_{chroma}_channel_effects.csv"
    stage1_df.to_csv(stage1_path, index=False)
    stage2_df.to_csv(stage2_path, index=False)

    return stage1_df, stage2_df, alpha_prime


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Two-stage group-level fNIRS analysis (global 2x2 contrasts with selective channel follow-up)."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_COMBINED_PATH,
        help="Path to runs-level combined long-format CSV (from combine_glm_output.py).",
    )
    parser.add_argument(
        "--chroma",
        choices=["hbo", "hbr"],
        default=DEFAULT_CHROMA,
        help="Chromophore to analyze at group level.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=DEFAULT_ALPHA,
        help="Alpha for Benjamini-Hochberg FDR correction (Stage 1). Stage-2 alpha' follows Benjamini-Bogomolov.",
    )
    parser.add_argument(
        "--conditions",
        nargs="*",
        default=[
            "Short-Form Education",
            "Short-Form Entertainment",
            "Long-Form Education",
            "Long-Form Entertainment",
        ],
        help="Condition names (order preserved) defining the 2x2 within-subject design.",
    )
    parser.add_argument(
        "--min-subjects",
        type=int,
        default=MIN_SUBJECTS_PER_CHANNEL,
        help="Minimum number of complete-case subjects required per channel (and for Stage 1).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    stage1_df, stage2_df, alpha_prime = run_group_analysis(
        combined_runs_path=args.input,
        chroma=args.chroma,
        alpha=args.alpha,
        conditions=args.conditions,
        min_subjects=args.min_subjects,
    )

    if stage1_df.empty:
        print("Stage 1: No effects were evaluated.")
    else:
        print("Stage 1: Global within-subject contrasts (Benjamini-Hochberg FDR)")
        for row in stage1_df.itertuples():
            t_val = f"{row.t_value:.3f}" if pd.notna(row.t_value) else "nan"
            df_val = f"{row.df:.1f}" if pd.notna(row.df) else "nan"
            mean_diff = f"{row.mean_difference:.3f}" if pd.notna(row.mean_difference) else "nan"
            p_val = f"{row.p_value:.3g}" if pd.notna(row.p_value) else "nan"
            q_val = f"{row.q_value:.3g}" if pd.notna(row.q_value) else "nan"
            dz = f"{row.cohen_dz:.3f}" if pd.notna(row.cohen_dz) else "nan"
            mark = "*" if bool(row.reject) else ""
            print(
                f"  {row.effect}{mark}: n={row.n_subjects}, t={t_val} (df={df_val}), "
                f"mean_diff={mean_diff}, p={p_val}, q={q_val}, dz={dz}"
            )

    selected_effects = stage1_df[stage1_df["reject"] == True]["effect_key"].tolist() if not stage1_df.empty else []  # noqa: E712
    if not selected_effects:
        print("Stage 2: Not run (no Stage-1 contrasts survived Holm correction).")
    elif stage2_df.empty:
        print("Stage 2: No channel results available after filtering.")
    else:
        print(
            f"Stage 2: Channel localisation for selected contrasts "
            f"(Benjamini-Hochberg at alpha' = {alpha_prime:.4f})"
        )
        for effect_key in selected_effects:
            subset = stage2_df[
                (stage2_df["effect_key"] == effect_key) & (stage2_df["reject"] == True)  # noqa: E712
            ]
            label = stage1_df.loc[stage1_df["effect_key"] == effect_key, "effect"].iloc[0]
            if subset.empty:
                print(f"  {label}: no channels passed FDR.")
                continue
            print(f"  {label}:")
            for row in subset.itertuples():
                t_val = f"{row.t_value:.3f}" if pd.notna(row.t_value) else "nan"
                p_val = f"{row.p_value:.3g}" if pd.notna(row.p_value) else "nan"
                q_val = f"{row.q_value:.3g}" if pd.notna(row.q_value) else "nan"
                dz = f"{row.cohen_dz:.3f}" if pd.notna(row.cohen_dz) else "nan"
                mean_diff = f"{row.mean_difference:.3f}" if pd.notna(row.mean_difference) else "nan"
                print(
                    f"    {row.channel}: n={row.n_subjects}, t={t_val}, mean_diff={mean_diff}, "
                    f"p={p_val}, q={q_val}, dz={dz}"
                )

    print(f"Stage-1 results saved: {args.input.parent / f'group_{args.chroma}_global_effects.csv'}")
    if stage2_df.empty:
        print("Stage-2 results file is empty (no contrasts selected or no channels retained).")
    else:
        print(f"Stage-2 results saved: {args.input.parent / f'group_{args.chroma}_channel_effects.csv'}")


if __name__ == "__main__":
    main()
