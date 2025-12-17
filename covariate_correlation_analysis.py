"""
Correlation diagnostics (Spearman) for the Short Form Video (SFV) study.

This script intentionally lives outside `process_sociodemographic.py` so preprocessing stays focused
and correlation/plotting stays modular.

Presets
-------
- covariates (default):
  - Input:  ./covariate_outputs/covariates_clean.csv
  - Output: ./covariate_outputs/covariate_correlations_spearman.csv
            ./covariate_outputs/covariate_correlations_pvalues.csv (if SciPy available; else NaN)
            ./covariate_outputs/covariate_heatmap.png

- combined:
  - Input:  ./data/tabular/combined_sfv_data.csv
  - Output: ./data/tabular/covariate_correlation_analysis.csv
            ./data/tabular/covariate_correlation_heatmap.png
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def spearman_pvalues(df_numeric: pd.DataFrame) -> pd.DataFrame:
    """
    Pairwise Spearman p-values.
    Uses SciPy if available; otherwise returns NaNs (with a warning).
    """
    try:
        from scipy.stats import spearmanr  # type: ignore
    except Exception:
        print("[WARN] SciPy not available; p-values will be written as NaN.")
        return pd.DataFrame(float("nan"), index=df_numeric.columns, columns=df_numeric.columns)

    cols = list(df_numeric.columns)
    pvals = pd.DataFrame(float("nan"), index=cols, columns=cols)

    for i, c1 in enumerate(cols):
        for j, c2 in enumerate(cols):
            if j < i:
                continue
            x = df_numeric[c1]
            y = df_numeric[c2]
            mask = x.notna() & y.notna()
            if mask.sum() < 3:
                p = float("nan")
            else:
                _, p = spearmanr(x[mask], y[mask])
            pvals.loc[c1, c2] = p
            pvals.loc[c2, c1] = p

    return pvals


def save_heatmap(corr: pd.DataFrame, outpath: Path, title: str) -> None:
    """
    Save a readable correlation heatmap.

    Layout improvements vs a naive heatmap:
    - Plot only the lower triangle (avoids duplicating information).
    - Disable cell annotations automatically for larger matrices.
    - Scale figure size + tick label styling based on matrix size.
    """
    n = int(corr.shape[0])
    # Heuristic sizing that stays readable without becoming enormous.
    # Give a bit more horizontal space so we can use larger tick-label fonts.
    fig_w = min(28, max(12, 0.80 * n + 7))
    fig_h = min(22, max(9, 0.65 * n + 5))

    # Only show one triangle to reduce clutter.
    mask = pd.DataFrame(
        [[j > i for j in range(n)] for i in range(n)],
        index=corr.index,
        columns=corr.columns,
    )

    # Annotating every cell gets unreadable quickly; keep it for small matrices only.
    annot = n <= 12

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        corr,
        mask=mask,
        annot=annot,
        fmt=".2f" if annot else "",
        annot_kws={"size": 8} if annot else None,
        cmap="vlag",
        center=0,
        vmin=-1,
        vmax=1,
        square=True,
        linewidths=0.5,
        cbar_kws={"shrink": 0.8},
        ax=ax,
    )
    ax.set_title(title, pad=12)

    # Tick label styling: rotate x labels and reduce font size for readability.
    # Slightly larger fonts help readability in papers/slides; still scale down for big matrices.
    tick_fs = 14 if n <= 22 else 12 if n <= 30 else 10
    ax.tick_params(axis="x", labelrotation=50, labelsize=tick_fs)
    ax.tick_params(axis="y", labelrotation=0, labelsize=tick_fs)
    for lbl in ax.get_xticklabels():
        lbl.set_horizontalalignment("right")

    # Extra bottom margin helps keep rotated x-labels readable with larger font sizes.
    fig.subplots_adjust(bottom=0.22)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def run_covariate_preset(root: Path, *, input_path: Optional[Path], out_dir: Optional[Path]) -> None:
    in_path = input_path or (root / "covariate_outputs" / "covariates_clean.csv")
    out = out_dir or (root / "covariate_outputs")
    out.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(in_path)
    df = df.apply(pd.to_numeric, errors="coerce")

    # Defensive: ensure subject_id (if ever present) is not included in correlations.
    corr_inputs = df.drop(columns=["subject_id"], errors="ignore")

    corr = corr_inputs.corr(method="spearman")
    corr.to_csv(out / "covariate_correlations_spearman.csv")

    pvals = spearman_pvalues(corr_inputs)
    pvals.to_csv(out / "covariate_correlations_pvalues.csv")

    save_heatmap(
        corr,
        out / "covariate_heatmap.png",
        title="Covariate Spearman correlations (encoded ordinal / totals)",
    )

    print(f"Loaded: {in_path}")
    print(f"Wrote:  {out / 'covariate_correlations_spearman.csv'}")
    print(f"Wrote:  {out / 'covariate_correlations_pvalues.csv'}")
    print(f"Wrote:  {out / 'covariate_heatmap.png'}")


def run_combined_preset(root: Path, *, input_path: Optional[Path], out_dir: Optional[Path]) -> None:
    in_path = input_path or (root / "data" / "tabular" / "combined_sfv_data.csv")
    out = out_dir or (root / "data" / "tabular")
    out.mkdir(parents=True, exist_ok=True)

    combined_data = pd.read_csv(in_path)
    correlations = combined_data.corr(method="spearman")

    # Preserve existing output format (no index column).
    correlations.to_csv(out / "covariate_correlation_analysis.csv", index=False)
    save_heatmap(
        correlations,
        out / "covariate_correlation_heatmap.png",
        title="Correlation Heatmap (Spearman)",
    )

    print(f"Loaded: {in_path}")
    print(f"Wrote:  {out / 'covariate_correlation_analysis.csv'}")
    print(f"Wrote:  {out / 'covariate_correlation_heatmap.png'}")


def main() -> None:
    root = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(description="Spearman correlation diagnostics + heatmaps.")
    parser.add_argument(
        "--preset",
        choices=["covariates", "combined"],
        default="combined",
        help="Which input/output convention to use.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Override the default input CSV path for the selected preset.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default="./covariate_outputs/",
        help="Override the default output directory for the selected preset.",
    )
    args = parser.parse_args()

    if args.preset == "covariates":
        print("[INFO] Running covariate correlation analysis...")
        run_covariate_preset(root, input_path=args.input, out_dir=args.out_dir)
    else:
        print("[INFO] Running combined correlation analysis...")
        run_combined_preset(root, input_path=args.input, out_dir=args.out_dir)


if __name__ == "__main__":
    main()