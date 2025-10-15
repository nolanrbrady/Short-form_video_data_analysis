"""
Combines subject-level GLM outputs into analysis-ready tables.

What it does
- Reads each subject's `*_ALL_runs_glm_results.csv` under `glm_results/<subject_id>/`.
- Filters to task regressors (drops derivatives, constant, and drift terms).
- Optionally filters by chromophore (HbO/HbR).
- Aggregates across runs per subject using inverse-variance weights (1 / se^2).
- Writes:
  - Long-format combined table for group analysis: one row per subject×channel×condition×chroma.
  - Wide matrices for convenience: one CSV per condition with subjects as rows and channels as columns.

Notes
- The subject-level GLM beta (theta) is already an estimate across repeated blocks per condition
  within a run; this script only aggregates across multiple runs for the same subject.
- Inverse-variance weighting improves efficiency when run-level standard errors differ.

Author: Nolan Brady
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List, Sequence

import numpy as np
import pandas as pd


# -----------------------
# Configuration defaults
# -----------------------

GLM_RESULTS_ROOT = Path("./glm_results")

# Default chromophores to include: both HbO and HbR
DEFAULT_CHROMAS: Sequence[str] = ("hbo", "hbr")

# Condition name filters. If None, infer from data by excluding nuisance terms.
DEFAULT_CONDITIONS: Sequence[str] | None = (
    "Short-Form Education",
    "Short-Form Entertainment",
    "Long-Form Education",
    "Long-Form Entertainment",
)


def _discover_subject_all_run_csvs(root: Path) -> List[Path]:
    """Return paths to each subject's ALL_runs GLM CSV under `root`.

    Expected layout: `<root>/<subject_id>/*_ALL_runs_glm_results.csv`
    """
    csvs: List[Path] = []
    if not root.exists():
        return csvs
    for subj_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        for csv in subj_dir.glob("*_ALL_runs_glm_results.csv"):
            csvs.append(csv)
    return csvs


def _is_task_regressor(row: pd.Series) -> bool:
    """True for task regressors (conditions) and False for nuisances.

    Drops:
    - Derivative terms (suffix `_derivative`)
    - Constant term (`constant`)
    - Cosine drifts (prefix `drift_`)
    """
    cond: str = str(row["Condition"]) if "Condition" in row else ""
    if not cond:
        return False
    if cond == "constant":
        return False
    if cond.startswith("drift_"):
        return False
    if cond.endswith("_derivative"):
        return False
    return True


def _weighted_mean_theta(df: pd.DataFrame) -> pd.Series:
    """Compute inverse-variance weighted mean of theta and count of runs.

    Assumes `theta` and `se` columns exist. If any `se` is 0 or NaN, falls
    back to an unweighted mean for numerical stability.
    """
    thetas = df["theta"].astype(float).to_numpy()
    ses = df["se"].astype(float).replace(0.0, np.nan).to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        weights = 1.0 / np.square(ses)
    if not np.all(np.isfinite(weights)) or np.nansum(weights) == 0.0:
        # Fallback: simple mean
        theta_hat = float(np.nanmean(thetas))
        n_runs = int(df["run"].nunique()) if "run" in df else int(len(df))
        return pd.Series({"theta": theta_hat, "n_runs": n_runs})
    theta_hat = float(np.nansum(weights * thetas) / np.nansum(weights))
    n_runs = int(df["run"].nunique()) if "run" in df else int(len(df))
    return pd.Series({"theta": theta_hat, "n_runs": n_runs})


def combine_subject_glm_runs(
    csv_paths: Iterable[Path],
    include_chromas: Sequence[str],
    conditions: Sequence[str] | None,
) -> pd.DataFrame:
    """Combine subject-level GLM CSVs without aggregating across runs.

    Output columns:
    - subject_id, run, ch_name, chroma, condition, theta, se (if available)
    """
    frames: List[pd.DataFrame] = []
    for csv in csv_paths:
        df = pd.read_csv(csv)
        df = df[df.apply(_is_task_regressor, axis=1)].copy()

        if "ch_name" not in df.columns:
            if {"Source", "Detector", "Chroma"}.issubset(df.columns):
                df["ch_name"] = df["Source"].astype(str).radd("S").str.cat(
                    df["Detector"].astype(str).radd("_D"), sep=" ").str.cat(
                    df["Chroma"].astype(str), sep=" ")
            else:
                raise ValueError(f"Cannot infer ch_name for rows in {csv}")

        # Filter by chromophore
        if "Chroma" in df.columns:
            df["Chroma"] = df["Chroma"].str.lower()
            df = df[df["Chroma"].isin([c.lower() for c in include_chromas])]
        else:
            df["Chroma"] = df["ch_name"].str.split().str[-1].str.lower()
            df = df[df["Chroma"].isin([c.lower() for c in include_chromas])]

        if conditions is not None:
            df = df[df["Condition"].isin(list(conditions))]

        if "subject_id" not in df.columns:
            df["subject_id"] = csv.parent.name

        # Keep only required columns and rename
        needed = [
            "subject_id",
            "run" if "run" in df.columns else None,
            "ch_name",
            "Chroma",
            "Condition",
            "theta",
            "se" if "se" in df.columns else None,
        ]
        needed = [c for c in needed if c is not None]
        df = df[needed].copy()
        df = df.rename(columns={"Chroma": "chroma", "Condition": "condition"})
        frames.append(df)

    if not frames:
        raise FileNotFoundError("No subject GLM CSVs found to combine (runs-level).")

    all_runs = pd.concat(frames, axis=0, ignore_index=True)
    # Sort for reproducibility
    all_runs = all_runs.sort_values(["subject_id", "run", "ch_name", "chroma", "condition"]).reset_index(drop=True)
    return all_runs


def combine_subject_glm(
    csv_paths: Iterable[Path],
    include_chromas: Sequence[str],
    conditions: Sequence[str] | None,
) -> pd.DataFrame:
    """Combine subject-level GLM CSVs into one long-format DataFrame.

    Output columns:
    - subject_id, ch_name, chroma, condition, theta, n_runs
    """
    frames: List[pd.DataFrame] = []
    for csv in csv_paths:
        # Explicit fail on read errors: surface problems early rather than proceed silently
        df = pd.read_csv(csv)

        # Keep only task regressors
        df = df[df.apply(_is_task_regressor, axis=1)].copy()

        # Align/rename expected columns
        # ch_name: e.g., "S1_D1 hbo"; chroma stored in `Chroma`
        if "ch_name" not in df.columns:
            # Derive from Source/Detector/Chroma if needed
            if {"Source", "Detector", "Chroma"}.issubset(df.columns):
                df["ch_name"] = df["Source"].astype(str).radd("S").str.cat(
                    df["Detector"].astype(str).radd("_D"), sep=" ").str.cat(
                    df["Chroma"].astype(str), sep=" ")
            else:
                raise ValueError(f"Cannot infer ch_name for rows in {csv}")

        # Filter by chromophore
        if "Chroma" in df.columns:
            df["Chroma"] = df["Chroma"].str.lower()
            df = df[df["Chroma"].isin([c.lower() for c in include_chromas])]
        else:
            # Try to extract chroma from ch_name suffix
            df["Chroma"] = df["ch_name"].str.split().str[-1].str.lower()
            df = df[df["Chroma"].isin([c.lower() for c in include_chromas])]

        # Optionally filter to specified conditions
        if conditions is not None:
            df = df[df["Condition"].isin(list(conditions))]

        # Ensure subject_id present
        if "subject_id" not in df.columns:
            # Fallback: parent folder name
            df["subject_id"] = csv.parent.name

        frames.append(df)

    if not frames:
        raise FileNotFoundError("No subject GLM CSVs found to combine.")

    all_df = pd.concat(frames, axis=0, ignore_index=True)

    # Aggregate across runs using inverse-variance weights
    grouped = (
        all_df.groupby(["subject_id", "ch_name", "Chroma", "Condition"], dropna=False)
        .apply(_weighted_mean_theta)
        .reset_index()
        .rename(columns={"Chroma": "chroma", "Condition": "condition"})
    )

    # Sort for pleasant, reproducible output
    grouped = grouped.sort_values(["subject_id", "ch_name", "chroma", "condition"]).reset_index(drop=True)
    return grouped


def write_outputs(grouped_long: pd.DataFrame, out_root: Path, all_runs_long: pd.DataFrame) -> None:
    """Write combined long table and convenience wide matrices per condition.

    Files written under `<out_root>`:
    - `combined_glm_long.csv` (subject×channel×condition×chroma)
    - `combined_matrices/<chroma>_<condition>.csv` (subjects rows, channels columns)
    - `combined_matrices/<chroma>_all_conditions.csv` (subjects rows, MultiIndex columns)
    """
    out_root.mkdir(parents=True, exist_ok=True)
    long_path = out_root / "combined_glm_long.csv"
    grouped_long.to_csv(long_path, index=False)
    # Save the non-aggregated runs-level long table for mixed-effects analysis
    long_runs_path = out_root / "combined_glm_long_runs.csv"
    all_runs_long.to_csv(long_runs_path, index=False)

    # Create wide matrices per condition and an all-conditions matrix (MultiIndex columns)
    matrices_dir = out_root / "combined_matrices"
    matrices_dir.mkdir(parents=True, exist_ok=True)

    # One file per chroma × condition
    chromas = grouped_long["chroma"].unique().tolist()
    conditions = grouped_long["condition"].unique().tolist()

    # Also collect wide with MultiIndex columns (condition, channel)
    for chroma in chromas:
        wide_all = (
            grouped_long[grouped_long["chroma"] == chroma]
            .pivot_table(index="subject_id", columns=["condition", "ch_name"], values="theta")
            .sort_index(axis=1)
        )
        wide_all.to_csv(matrices_dir / f"{chroma}_all_conditions.csv")

        for cond in conditions:
            wide = (
                grouped_long[(grouped_long["chroma"] == chroma) & (grouped_long["condition"] == cond)]
                .pivot(index="subject_id", columns="ch_name", values="theta")
                .sort_index(axis=1)
            )
            # Only write if non-empty
            if not wide.empty:
                # Replace spaces in condition name for file naming
                safe_cond = cond.replace(" ", "_")
                wide.to_csv(matrices_dir / f"{chroma}_{safe_cond}.csv")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Combine subject-level GLM outputs for group analysis.")
    parser.add_argument(
        "--root",
        type=Path,
        default=GLM_RESULTS_ROOT,
        help="Root directory containing <subject_id>/*_ALL_runs_glm_results.csv",
    )
    parser.add_argument(
        "--chroma",
        dest="chromas",
        choices=["hbo", "hbr", "both"],
        default="both",
        help="Chromophore(s) to include in outputs.",
    )
    parser.add_argument(
        "--conditions",
        nargs="*",
        default=list(DEFAULT_CONDITIONS) if DEFAULT_CONDITIONS is not None else None,
        help="Explicit list of condition names to include. Defaults to the four task conditions.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=GLM_RESULTS_ROOT,
        help="Output directory for combined files.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    chroma_choice = args.chromas
    if chroma_choice == "both":
        chromas = ("hbo", "hbr")
    else:
        chromas = (chroma_choice,)

    csvs = _discover_subject_all_run_csvs(args.root)
    if not csvs:
        raise SystemExit(f"No *_ALL_runs_glm_results.csv files found under {args.root}")

    combined_runs = combine_subject_glm_runs(csvs, include_chromas=chromas, conditions=args.conditions)
    combined_long = combine_subject_glm(csvs, include_chromas=chromas, conditions=args.conditions)
    write_outputs(combined_long, args.out, combined_runs)

    print(f"Combined long table written to: {args.out / 'combined_glm_long.csv'}")
    print(f"Per-condition matrices written under: {args.out / 'combined_matrices'}")


if __name__ == "__main__":
    main()
