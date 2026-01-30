"""
Channelwise within-subject inference for Format x Content effects on prefrontal activation.

Inputs
  - Homer3 GLM betas (wide): `data/tabular/homer3_glm_betas_wide.csv`
      - ID column: `Subject` (e.g., `sub_0001`)
      - Beta columns: `S##_D##_Cond##_HbO` / `S##_D##_Cond##_HbR`
  - Combined tabular data (wide): `data/tabular/combined_sfv_data.csv`
      - ID column: `subject_id` (may be zero-padded in some exports)

Design
  - Within-subjects 2x2 factorial:
      Format:  Short vs Long
      Content: Education vs Entertainment
  - Condition mapping (trigger codes):
      Cond01 = Short-Form Education
      Cond02 = Short-Form Entertainment
      Cond03 = Long-Form Entertainment
      Cond04 = Long-Form Education

Model (per channel x chromophore)
  - Linear mixed model (random intercept for subject):
      beta ~ format_c * content_c + (1 | subject_id)
    where:
      format_c  = -0.5 (Short), +0.5 (Long)
      content_c = -0.5 (Entertainment), +0.5 (Education)

Missingness / pruned channels
  - Per repo policy, BOTH 0 and NaN in `homer3_glm_betas_wide.csv` can indicate pruned channels.
    This script treats:
      beta == 0  -> missing
      beta == NaN -> missing
  - Default behavior is complete-case within-channel: subjects missing any of the 4 conditions for a
    given (channel, chromophore) are dropped for that model.

Multiple testing correction
  - Benjamini–Hochberg FDR (BH) performed:
      separately per chromophore (HbO vs HbR),
      separately per effect (Format, Content, Interaction),
      across channels.

Post-hoc (only if interaction is FDR-significant for that channel/chromophore)
  - All pairwise within-subject comparisons among the 4 conditions (6 contrasts),
    using paired t-tests (two-sided), uncorrected (per study instruction).

Citations (see CITATIONS.md)
  - Benjamini & Hochberg (1995): BH-FDR correction.
  - Laird & Ware (1982): mixed-effects models.

Note: This Python script requires `statsmodels` to fit MixedLM. If not installed, run with
`--dry-run` to validate parsing/reshaping only, or install dependencies before full analysis.
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal

import numpy as np
import pandas as pd
from scipy import stats as scistats


CHANNEL_COL_RE = re.compile(
    r"^(?P<channel>S\d+_D\d+)_Cond(?P<cond>\d{2})_(?P<chrom>HbO|HbR)$"
)
_DIGITS_RE = re.compile(r"(\d+)")


ConditionCode = Literal["01", "02", "03", "04"]


COND_MAP: dict[ConditionCode, dict[str, object]] = {
    "01": {"condition": "SF_Edu", "format": "Short", "content": "Education", "format_c": -0.5, "content_c": +0.5},
    "02": {"condition": "SF_Ent", "format": "Short", "content": "Entertainment", "format_c": -0.5, "content_c": -0.5},
    "03": {"condition": "LF_Ent", "format": "Long", "content": "Entertainment", "format_c": +0.5, "content_c": -0.5},
    "04": {"condition": "LF_Edu", "format": "Long", "content": "Education", "format_c": +0.5, "content_c": +0.5},
}


def _bh_fdr(p_values: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR (BH) q-values for a 1D array; preserves NaNs."""
    p = np.asarray(p_values, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    finite_mask = np.isfinite(p)
    if finite_mask.sum() == 0:
        return q
    p_f = p[finite_mask]
    order = np.argsort(p_f)
    ranked = p_f[order]
    m = ranked.size
    # BH: q_i = min_{j>=i} (m/j) * p_j
    q_ranked = np.empty_like(ranked)
    prev = 1.0
    for i in range(m - 1, -1, -1):
        rank = i + 1
        val = (m / rank) * ranked[i]
        prev = min(prev, val)
        q_ranked[i] = prev
    q_ranked = np.clip(q_ranked, 0.0, 1.0)
    # Undo sorting
    q_f = np.empty_like(q_ranked)
    q_f[order] = q_ranked
    q[finite_mask] = q_f
    return q


def _normalize_subject_id(series: pd.Series, *, column_name: str) -> pd.Series:
    """Extract numeric subject id digits; returns pandas nullable Int64."""
    as_str = series.astype("string")
    extracted = as_str.str.extract(_DIGITS_RE, expand=False)
    if extracted.isna().any():
        bad = as_str[extracted.isna()].dropna().unique().tolist()[:10]
        raise ValueError(
            f"Failed to parse numeric IDs from column '{column_name}'. "
            f"Examples of unparseable values: {bad}"
        )
    return extracted.astype("Int64")


def _load_and_merge(homer_csv: Path, combined_csv: Path) -> pd.DataFrame:
    homer = pd.read_csv(homer_csv)
    combined = pd.read_csv(combined_csv)

    if "Subject" not in homer.columns:
        raise ValueError(f"Expected column 'Subject' in {homer_csv}")
    if "subject_id" not in combined.columns:
        raise ValueError(f"Expected column 'subject_id' in {combined_csv}")

    homer = homer.copy()
    combined = combined.copy()
    homer["subject_id"] = _normalize_subject_id(homer["Subject"], column_name="Subject")
    combined["subject_id"] = _normalize_subject_id(combined["subject_id"], column_name="subject_id")

    # Fail hard if either input contains duplicate subject IDs. The analysis spec assumes one row per subject.
    for df_name, df in [("combined", combined), ("homer", homer)]:
        dup_mask = df["subject_id"].duplicated(keep=False)
        if dup_mask.any():
            counts = df.loc[dup_mask, "subject_id"].astype(int).value_counts().sort_index()
            examples = counts.head(10).to_dict()
            raise ValueError(
                f"Duplicate subject_id values detected in {df_name} dataset after normalization. "
                f"Expected exactly one row per subject. "
                f"Example duplicate counts (subject_id: n_rows): {examples}. "
                f"Fix the upstream CSV before running the channelwise analysis."
            )

    homer = homer.rename(columns={"Subject": "homer_subject"})
    merged = combined.merge(homer, how="inner", on="subject_id", validate="one_to_one")
    return merged


def _extract_beta_columns(df: pd.DataFrame) -> list[str]:
    beta_cols: list[str] = []
    for c in df.columns:
        if CHANNEL_COL_RE.match(c):
            beta_cols.append(c)
    if not beta_cols:
        raise ValueError("No beta columns matched expected pattern like 'S01_D01_Cond01_HbO'.")
    return beta_cols


def _reshape_to_long(df_merged: pd.DataFrame) -> pd.DataFrame:
    if "subject_id" not in df_merged.columns:
        raise ValueError("Expected column 'subject_id' in merged dataset.")

    beta_cols = _extract_beta_columns(df_merged)
    long = df_merged[["subject_id"] + beta_cols].melt(
        id_vars=["subject_id"], var_name="beta_col", value_name="beta"
    )

    extracted = long["beta_col"].str.extract(CHANNEL_COL_RE)
    if extracted.isna().any().any():
        bad = long.loc[extracted.isna().any(axis=1), "beta_col"].unique().tolist()[:10]
        raise ValueError(f"Failed to parse channel/cond/chrom from some beta columns. Examples: {bad}")

    long = pd.concat([long.drop(columns=["beta_col"]), extracted], axis=1)
    long["beta"] = pd.to_numeric(long["beta"], errors="coerce")

    # Treat 0 and NaN as pruned/missing per repo policy
    long.loc[long["beta"] == 0, "beta"] = np.nan

    if not set(long["cond"].unique()).issubset(set(COND_MAP.keys())):
        bad_conds = sorted(set(long["cond"].unique()) - set(COND_MAP.keys()))
        raise ValueError(f"Unexpected condition codes in betas: {bad_conds}")

    long["condition"] = long["cond"].map(lambda c: COND_MAP[c]["condition"])
    long["format"] = long["cond"].map(lambda c: COND_MAP[c]["format"])
    long["content"] = long["cond"].map(lambda c: COND_MAP[c]["content"])
    long["format_c"] = long["cond"].map(lambda c: float(COND_MAP[c]["format_c"]))
    long["content_c"] = long["cond"].map(lambda c: float(COND_MAP[c]["content_c"]))
    long = long.drop(columns=["cond"])

    return long


def _complete_case_subjects(sub: pd.DataFrame) -> pd.DataFrame:
    """Keep only subjects with all 4 conditions present (non-missing beta)."""
    non_missing = sub.dropna(subset=["beta"]).copy()
    counts = non_missing.groupby("subject_id")["condition"].nunique()
    keep_ids = counts[counts == 4].index
    out = non_missing[non_missing["subject_id"].isin(keep_ids)].copy()
    return out


@dataclass(frozen=True)
class EffectRow:
    channel: str
    chrom: str
    n_subjects: int
    n_obs: int
    estimate: float
    se: float
    stat: float
    p_unc: float
    ci95_low: float
    ci95_high: float


def _fit_mixedlm(sub: pd.DataFrame):
    try:
        import statsmodels.formula.api as smf
    except Exception as exc:  # pragma: no cover
        raise RuntimeError(
            "Missing dependency: statsmodels. Install it to run the LMM analysis, "
            "or rerun with --dry-run."
        ) from exc

    model = smf.mixedlm("beta ~ format_c * content_c", data=sub, groups=sub["subject_id"], re_formula="1")
    # Use REML; random intercept only
    result = model.fit(reml=True, method=["lbfgs", "bfgs", "powell"], maxiter=200, disp=False)
    return result


def _extract_effect(result, term: str) -> tuple[float, float, float, float, float, float]:
    if term not in result.params.index:
        raise ValueError(f"Expected term '{term}' in fitted model params. Found: {list(result.params.index)}")
    est = float(result.params[term])
    se = float(result.bse[term])
    stat = float(result.tvalues[term])
    p = float(result.pvalues[term])
    zcrit = 1.959963984540054  # ~N(0,1) 97.5%
    ci_low = est - zcrit * se
    ci_high = est + zcrit * se
    return est, se, stat, p, ci_low, ci_high


def _paired_posthoc(sub_complete: pd.DataFrame) -> pd.DataFrame:
    """
    All pairwise within-subject comparisons among the 4 conditions using paired t-tests.

    Returns columns:
      condition_a, condition_b, mean_diff (a-b), ci95_low, ci95_high, t, df, p_unc
    """
    wide = sub_complete.pivot_table(index="subject_id", columns="condition", values="beta", aggfunc="mean")
    needed = ["SF_Edu", "SF_Ent", "LF_Ent", "LF_Edu"]
    missing_cols = [c for c in needed if c not in wide.columns]
    if missing_cols:
        raise ValueError(f"Missing condition columns in posthoc wide table: {missing_cols}")
    wide = wide[needed].dropna(axis=0, how="any")
    if wide.empty:
        raise ValueError("No complete subjects available for posthoc contrasts.")

    rows = []
    pairs = [
        ("SF_Edu", "SF_Ent"),
        ("SF_Edu", "LF_Ent"),
        ("SF_Edu", "LF_Edu"),
        ("SF_Ent", "LF_Ent"),
        ("SF_Ent", "LF_Edu"),
        ("LF_Ent", "LF_Edu"),
    ]
    for a, b in pairs:
        diffs = (wide[a] - wide[b]).to_numpy(dtype=float)
        n = diffs.size
        if n < 2:
            continue
        mean_diff = float(np.mean(diffs))
        se = float(np.std(diffs, ddof=1) / math.sqrt(n))
        t_stat, p_val = scistats.ttest_rel(wide[a], wide[b], nan_policy="raise")
        df = n - 1
        tcrit = float(scistats.t.ppf(0.975, df))
        ci_low = mean_diff - tcrit * se
        ci_high = mean_diff + tcrit * se
        rows.append(
            {
                "condition_a": a,
                "condition_b": b,
                "n_subjects": n,
                "mean_diff": mean_diff,
                "ci95_low": ci_low,
                "ci95_high": ci_high,
                "t": float(t_stat),
                "df": float(df),
                "p_unc": float(p_val),
            }
        )
    return pd.DataFrame(rows)


def run_analysis(
    *,
    homer_csv: Path,
    combined_csv: Path,
    out_main_csv: Path,
    out_posthoc_csv: Path,
    alpha: float,
    min_subjects: int,
    dry_run: bool,
    fail_fast: bool,
) -> None:
    merged = _load_and_merge(homer_csv=homer_csv, combined_csv=combined_csv)
    long = _reshape_to_long(merged)

    channels = sorted(long["channel"].unique().tolist())
    chroms = ["HbO", "HbR"]
    print(f"[data] merged subjects: {merged['subject_id'].nunique()}")
    print(f"[data] channels detected: {len(channels)}")
    print(f"[data] chromophores detected: {sorted(long['chrom'].unique().tolist())}")

    if dry_run:
        # Sanity diagnostics for parsing and missingness
        non_missing = long.dropna(subset=["beta"])
        summary = (
            non_missing.groupby(["chrom", "channel"])["subject_id"].nunique().reset_index(name="n_subjects_nonmissing")
        )
        print(f"[dry-run] non-missing channel/chrom pairs: {len(summary)}")
        print("[dry-run] done (no model fits).")
        return

    # Model per channel x chrom
    main_rows = []
    errors = []
    gated_count = 0
    gated_examples = []

    for chrom in chroms:
        for channel in channels:
            sub = long[(long["chrom"] == chrom) & (long["channel"] == channel)].copy()
            sub = _complete_case_subjects(sub)
            n_subjects = int(sub["subject_id"].nunique())
            if n_subjects < min_subjects:
                gated_count += 1
                if len(gated_examples) < 10:
                    gated_examples.append(f"{chrom} {channel} (n_subjects={n_subjects} < min_subjects={min_subjects})")
                continue
            n_obs = int(len(sub))
            try:
                result = _fit_mixedlm(sub)
                for effect_name, term in [
                    ("format", "format_c"),
                    ("content", "content_c"),
                    ("interaction", "format_c:content_c"),
                ]:
                    est, se, stat, p, ci_low, ci_high = _extract_effect(result, term)
                    main_rows.append(
                        {
                            "channel": channel,
                            "chrom": chrom,
                            "effect": effect_name,
                            "n_subjects": n_subjects,
                            "n_obs": n_obs,
                            "estimate": est,
                            "se": se,
                            "stat": stat,
                            "p_unc": p,
                            "ci95_low": ci_low,
                            "ci95_high": ci_high,
                        }
                    )
            except Exception as exc:
                msg = f"{chrom} {channel}: {exc}"
                if fail_fast:
                    raise RuntimeError(msg) from exc
                errors.append(msg)

    if not main_rows:
        msg = "No models were fit (possibly due to min_subjects threshold or missing data)."
        if gated_count > 0:
            msg += (
                f" {gated_count} channel/chrom pairs were skipped due to min_subjects gating (min_subjects={min_subjects})."
            )
            if gated_examples:
                msg += f" Examples: {gated_examples}"
        if errors:
            msg += f" {len(errors)} channel/chrom fits failed. First few: {errors[:10]}"
        raise ValueError(msg)

    main_df = pd.DataFrame(main_rows)

    # BH-FDR: separate per chromophore, separate per effect, across channels
    main_df["p_fdr"] = np.nan
    for chrom in chroms:
        for effect in ["format", "content", "interaction"]:
            mask = (main_df["chrom"] == chrom) & (main_df["effect"] == effect)
            main_df.loc[mask, "p_fdr"] = _bh_fdr(main_df.loc[mask, "p_unc"].to_numpy(dtype=float))

    # Gate posthoc by interaction q-value
    sig_interactions = main_df[(main_df["effect"] == "interaction") & (main_df["p_fdr"] < alpha)].copy()
    print(f"[results] interaction FDR-significant (q<{alpha}): {len(sig_interactions)} channel/chrom pairs")

    if gated_count > 0:
        print(
            f"[warn] skipped {gated_count} channel/chrom models due to min_subjects gating "
            f"(min_subjects={min_subjects}); showing up to 10 examples: {gated_examples}"
        )

    posthoc_rows = []
    for row in sig_interactions.itertuples(index=False):
        chrom = row.chrom
        channel = row.channel
        sub = long[(long["chrom"] == chrom) & (long["channel"] == channel)].copy()
        sub = _complete_case_subjects(sub)
        ph = _paired_posthoc(sub)
        if ph.empty:
            continue
        ph.insert(0, "channel", channel)
        ph.insert(1, "chrom", chrom)
        posthoc_rows.append(ph)

    out_main_csv.parent.mkdir(parents=True, exist_ok=True)
    main_df = main_df.sort_values(["chrom", "effect", "channel"]).reset_index(drop=True)
    main_df.to_csv(out_main_csv, index=False)
    print(f"[write] main effects: {out_main_csv}")

    if posthoc_rows:
        posthoc_df = pd.concat(posthoc_rows, ignore_index=True)
        posthoc_df = posthoc_df.sort_values(["chrom", "channel", "condition_a", "condition_b"]).reset_index(drop=True)
        posthoc_df.to_csv(out_posthoc_csv, index=False)
        print(f"[write] posthoc:      {out_posthoc_csv}")
    else:
        # Still write an empty CSV with headers for predictable downstream automation
        empty = pd.DataFrame(
            columns=[
                "channel",
                "chrom",
                "condition_a",
                "condition_b",
                "n_subjects",
                "mean_diff",
                "ci95_low",
                "ci95_high",
                "t",
                "df",
                "p_unc",
            ]
        )
        empty.to_csv(out_posthoc_csv, index=False)
        print(f"[write] posthoc:      {out_posthoc_csv} (empty; no significant interactions)")

    if errors:
        print("[warn] some channel/chrom models failed:")
        for e in errors[:25]:
            print("  -", e)
        if len(errors) > 25:
            print(f"  ... {len(errors) - 25} more")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--homer-csv", type=Path, default=Path("data/tabular/homer3_glm_betas_wide.csv"))
    parser.add_argument("--combined-csv", type=Path, default=Path("data/tabular/combined_sfv_data.csv"))
    parser.add_argument("--out-main-csv", type=Path, default=Path("data/results/format_content_lmm_main_effects_python.csv"))
    parser.add_argument(
        "--out-posthoc-csv", type=Path, default=Path("data/results/format_content_lmm_posthoc_pairwise_python.csv")
    )
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--min-subjects", type=int, default=6)
    parser.add_argument("--dry-run", action="store_true", help="Parse/reshape only; do not fit models.")
    parser.add_argument("--fail-fast", action="store_true", help="Stop on first model failure.")
    args = parser.parse_args()

    run_analysis(
        homer_csv=args.homer_csv,
        combined_csv=args.combined_csv,
        out_main_csv=args.out_main_csv,
        out_posthoc_csv=args.out_posthoc_csv,
        alpha=args.alpha,
        min_subjects=args.min_subjects,
        dry_run=args.dry_run,
        fail_fast=args.fail_fast,
    )


if __name__ == "__main__":
    main()
