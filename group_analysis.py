"""
Group-level analysis pipeline using linear mixed-effects models (LMM).

Workflow
- Input: combined long-format table from `combine_glm_output.py` with columns
  `subject_id, ch_name, chroma, condition, theta, n_runs`.
- For each channel (and specified chroma):
  1) Fit an LMM on run-level betas: `theta ~ C(condition)` with random intercept for subject.
  2) Test the main effect of condition via a joint Wald test of all condition coefficients.
  3) Correct omnibus p-values across channels using FDR.
  4) For channels passing omnibus FDR, run uncorrected LMM pairwise contrasts for condition differences.

Stats choices (following docstring guidance)
- Omnibus: LMM joint Wald test of condition fixed effects using `statsmodels.formula.api.mixedlm`.
- Multiple comparisons: Benjamini–Hochberg FDR across channels for the omnibus tests.
- Post-hoc: LMM-based pairwise Wald contrasts (uncorrected, gated by omnibus).

Assumptions and notes
- Uses aggregated subject-level betas (weighted across runs) as the dependent variable.
- Subjects with any missing condition for a given channel are dropped from that channel's test.
- AnovaRM assumes sphericity; for maximum rigor, consider complementing with a
  mixed-effects model as a sensitivity analysis (not implemented here to honor the docstring).

Author: Nolan Brady
"""

from __future__ import annotations

import argparse
import itertools
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats as scistats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection


# -----------------------
# Configuration defaults
# -----------------------

DEFAULT_COMBINED_PATH = Path("./glm_results/combined_glm_long_runs.csv")
DEFAULT_CHROMA = "hbo"  # typical primary analyte
DEFAULT_ALPHA = 0.05
MIN_SUBJECTS_PER_CHANNEL = 6  # drop channels with fewer than this many complete subjects


@dataclass
class OmnibusResult:
    channel: str
    chroma: str
    n_subjects: int
    df: float
    chi2: float
    p_value: float
    q_value: float | None = None


@dataclass
class PosthocResult:
    channel: str
    chroma: str
    n_subjects: int
    condition_a: str
    condition_b: str
    t_value: float
    p_value: float
    q_value: float | None
    cohen_dz: float


def _prepare_runs_for_channel(
    df_runs: pd.DataFrame, channel: str, chroma: str, conditions: Sequence[str]
) -> pd.DataFrame:
    """Filter to a single channel/chroma and subset to specified conditions.

    Returns columns: subject_id, run, condition, theta.
    """
    sub = df_runs[(df_runs["ch_name"] == channel) & (df_runs["chroma"] == chroma)].copy()
    sub = sub[sub["condition"].isin(list(conditions))]
    needed_cols = ["subject_id", "run", "condition", "theta"]
    for c in needed_cols:
        if c not in sub.columns:
            raise ValueError(f"Expected column '{c}' in runs-level data")
    return sub[needed_cols].copy()


def _fit_lmm_and_wald_omnibus(sub_runs: pd.DataFrame) -> Tuple[float, float, float, int, object]:
    """Fit LMM: theta ~ C(condition) + (1|subject_id). Return (chi2, p, df, n_subj, result).

    Uses a joint Wald test that all condition coefficients equal zero.
    """
    if sub_runs.empty:
        return (np.nan, np.nan, np.nan, 0, None)
    n_subj = int(sub_runs["subject_id"].nunique())
    model = smf.mixedlm("theta ~ C(condition)", data=sub_runs, groups=sub_runs["subject_id"], re_formula="1")
    result = model.fit(reml=True, method="lbfgs", maxiter=200, disp=False)
    # Identify fixed effect columns for condition terms
    fe_names = list(result.fe_params.index)
    cond_cols = [i for i, name in enumerate(fe_names) if name.startswith("C(condition)[T.")]
    if len(cond_cols) == 0:
        # No variation in condition or only baseline present
        return (np.nan, np.nan, np.nan, n_subj, result)
    import numpy as _np
    R = _np.zeros((len(cond_cols), len(fe_names)))
    for r, idx in enumerate(cond_cols):
        R[r, idx] = 1.0
    w = result.wald_test(R)
    chi2 = float(w.statistic)
    df = float(w.df_num)
    pval = float(w.pvalue)
    return (chi2, pval, df, n_subj, result)


def _lmm_pairwise_contrasts(result, conditions: Sequence[str]) -> List[Tuple[str, str, float, float, float]]:
    """Compute LMM pairwise contrasts for condition means using Wald tests.

    Returns list of (cond_a, cond_b, estimate, se, p_value).
    """
    fe = result.fe_params
    cov = result.cov_params()
    names = list(fe.index)

    # Helper to build fixed-effects row vector for a given condition value under treatment coding
    def row_for_condition(cond: str) -> np.ndarray:
        x = np.zeros(len(names), dtype=float)
        # Intercept term
        if "Intercept" in names:
            x[names.index("Intercept")] = 1.0
        # Treatment coding: columns like C(condition)[T.<level>]
        for i, nm in enumerate(names):
            if nm.startswith("C(condition)[T.") and nm.endswith("]"):
                level = nm[len("C(condition)[T.") : -1]
                if cond == level:
                    x[i] = 1.0
        return x

    results: List[Tuple[str, str, float, float, float]] = []
    for a, b in itertools.combinations(conditions, 2):
        xa = row_for_condition(a)
        xb = row_for_condition(b)
        L = xb - xa  # tests mean(b) - mean(a) = 0
        est = float(L @ fe)
        se = float(np.sqrt(L @ cov.values @ L))
        # Wald chi-square with df=1; convert to two-sided p-value via normal approx
        z = est / se if se > 0 else np.nan
        pval = float(2 * (1 - scistats.norm.cdf(abs(z)))) if np.isfinite(z) else np.nan
        results.append((a, b, est, se, pval))
    return results


def run_group_analysis(
    combined_runs_path: Path,
    chroma: str,
    alpha: float,
    conditions: Sequence[str] | None = None,
    min_subjects: int = MIN_SUBJECTS_PER_CHANNEL,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Run the LMM-based group pipeline and return (omnibus_df, posthoc_df).

    Saves CSV outputs alongside the input file under `combined_runs_path.parent`.
    """
    df = pd.read_csv(combined_runs_path)
    # Normalize column names
    for col in ("chroma", "condition", "ch_name", "subject_id"):
        if col not in df.columns:
            raise ValueError(f"Expected column '{col}' in {combined_long_path}")
    df["chroma"] = df["chroma"].str.lower()

    # Filter chroma and conditions
    df = df[df["chroma"] == chroma.lower()].copy()
    if conditions is None:
        # Infer from data by excluding typical nuisances (already done upstream)
        conditions = (
            df["condition"].dropna().unique().tolist()
        )
    else:
        df = df[df["condition"].isin(list(conditions))]

    channels = sorted(df["ch_name"].unique().tolist())
    omni_rows: List[OmnibusResult] = []

    # Omnibus tests per channel
    for ch in channels:
        sub_runs = _prepare_runs_for_channel(df, channel=ch, chroma=chroma, conditions=conditions)
        n_subj = int(sub_runs["subject_id"].nunique()) if not sub_runs.empty else 0
        if n_subj < min_subjects:
            omni_rows.append(
                OmnibusResult(channel=ch, chroma=chroma, n_subjects=n_subj, df=np.nan, chi2=np.nan, p_value=np.nan)
            )
            continue
        chi2, p_val, df_dof, _, _res = _fit_lmm_and_wald_omnibus(sub_runs)
        omni_rows.append(
            OmnibusResult(channel=ch, chroma=chroma, n_subjects=n_subj, df=df_dof, chi2=chi2, p_value=p_val)
        )

    anova_df = pd.DataFrame([r.__dict__ for r in omni_rows])

    # FDR across channels (omnibus)
    valid_mask = anova_df["p_value"].notna()
    if valid_mask.any():
        reject, qvals = fdrcorrection(anova_df.loc[valid_mask, "p_value"].to_numpy(), alpha=alpha)
        anova_df.loc[valid_mask, "q_value"] = qvals
        anova_df.loc[valid_mask, "reject_fdr"] = reject
    else:
        anova_df["q_value"] = np.nan
        anova_df["reject_fdr"] = False

    # Post-hoc within significant channels
    posthoc_records: List[PosthocResult] = []
    sig_channels = anova_df[(anova_df["reject_fdr"] == True) & (anova_df["p_value"].notna())]["channel"].tolist()  # noqa: E712
    for ch in sig_channels:
        sub_runs = _prepare_runs_for_channel(df, channel=ch, chroma=chroma, conditions=conditions)
        _, _, _, n_subj, res = _fit_lmm_and_wald_omnibus(sub_runs)
        if res is None:
            continue
        pairs = _lmm_pairwise_contrasts(res, conditions=conditions)
        for (a, b, est, se, pval) in pairs:
            # For LMM contrasts, report z-stat implied by estimate/SE and leave q_value as None
            z = est / se if se > 0 else np.nan
            posthoc_records.append(
                PosthocResult(
                    channel=ch,
                    chroma=chroma,
                    n_subjects=n_subj,
                    condition_a=a,
                    condition_b=b,
                    t_value=z,
                    p_value=pval,
                    q_value=None,
                    cohen_dz=np.nan,
                )
            )

    posthoc_df = pd.DataFrame([r.__dict__ for r in posthoc_records])

    # Write outputs next to input
    out_dir = combined_runs_path.parent
    main_path = out_dir / f"group_{chroma}_main_effects.csv"
    posthoc_path = out_dir / f"group_{chroma}_posthoc_pairs.csv"
    anova_df.to_csv(main_path, index=False)
    posthoc_df.to_csv(posthoc_path, index=False)

    return anova_df, posthoc_df


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Group-level fNIRS analysis: LMM omnibus with FDR and gated LMM contrasts.")
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
        help="Family-wise alpha for FDR procedures.",
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
        help="Explicit condition names and order to analyze.",
    )
    parser.add_argument(
        "--min-subjects",
        type=int,
        default=MIN_SUBJECTS_PER_CHANNEL,
        help="Minimum number of complete-case subjects required per channel.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    anova_df, posthoc_df = run_group_analysis(
        combined_runs_path=args.input,
        chroma=args.chroma,
        alpha=args.alpha,
        conditions=args.conditions,
        min_subjects=args.min_subjects,
    )

    reject_col = anova_df["reject_fdr"] if "reject_fdr" in anova_df.columns else pd.Series(False, index=anova_df.index)
    sig_main = anova_df[(reject_col == True) & (anova_df["p_value"].notna())]  # noqa: E712
    if not sig_main.empty:
        print("Significant omnibus effects (FDR-corrected):")
        for row in sig_main.itertuples():
            chi2 = f"{row.chi2:.3f}" if pd.notna(row.chi2) else "nan"
            p_val = f"{row.p_value:.3g}" if pd.notna(row.p_value) else "nan"
            q_val = f"{row.q_value:.3g}" if pd.notna(row.q_value) else "nan"
            print(f"  {row.channel} (n={row.n_subjects}, chi2={chi2}, p={p_val}, q={q_val})")
    else:
        print("No significant omnibus effects after FDR correction.")

    sig_posthoc = posthoc_df[(posthoc_df["p_value"].notna()) & (posthoc_df["p_value"] < args.alpha)]
    if not sig_posthoc.empty:
        print(f"Significant post-hoc contrasts (uncorrected p < {args.alpha:g}):")
        for row in sig_posthoc.itertuples():
            t_val = f"{row.t_value:.3f}" if pd.notna(row.t_value) else "nan"
            p_val = f"{row.p_value:.3g}" if pd.notna(row.p_value) else "nan"
            print(
                f"  {row.channel}: {row.condition_b} vs {row.condition_a} "
                f"(n={row.n_subjects}, z={t_val}, p={p_val})"
            )
    else:
        print(f"No post-hoc contrasts met p < {args.alpha:g}.")

    print(f"Main-effects saved: {args.input.parent / f'group_{args.chroma}_main_effects.csv'}")
    if not posthoc_df.empty:
        print(f"Post-hocs saved: {args.input.parent / f'group_{args.chroma}_posthoc_pairs.csv'}")
    else:
        print("No channels passed omnibus FDR; no post-hocs run.")


if __name__ == "__main__":
    main()
