"""
Pearson correlation diagnostics for the Short Form Video (SFV) study.

This script intentionally lives outside `process_sociodemographic.py` so preprocessing stays focused
and correlation/plotting stays modular.

Presets
-------
- covariates (default):
  - Input:  ./covariate_outputs/covariates_clean.csv
  - Output: ./covariate_outputs/covariate_correlations_pearson.csv
            ./covariate_outputs/covariate_correlations_pearson_pvalues.csv
            ./covariate_outputs/covariate_heatmap_pearson.png

- combined:
  - Input:  ./data/tabular/generated_data/combined_sfv_data.csv
  - Output: ./data/tabular/generated_data/covariate_correlation_analysis_pearson.csv
            ./data/tabular/generated_data/covariate_correlation_analysis_pearson_pvalues.csv
            ./data/tabular/generated_data/covariate_correlation_heatmap_pearson.png
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

try:
    import seaborn as sns
except ModuleNotFoundError:
    sns = None

from create_demographics_table import load_excluded_subject_ids, normalize_subject_ids


DEFAULT_EXCLUDED_SUBJECTS_JSON = Path("data/config/excluded_subjects.json")
DEFAULT_SUBJECT_COLUMN = "subject_id"
DEFAULT_MAX_SUBJECTS = 48
OBSOLETE_OUTPUT_FILENAMES = {
    "covariate_correlation_analysis.csv",
    "covariate_correlation_analysis_spearman.csv",
    "covariate_correlation_heatmap.png",
    "covariate_correlation_heatmap_spearman.png",
    "covariate_correlations_spearman.csv",
    "covariate_correlations_spearman_pvalues.csv",
    "covariate_correlations_pvalues.csv",
    "covariate_heatmap.png",
    "covariate_heatmap_spearman.png",
}


@dataclass(frozen=True)
class PreparedCorrelationInputs:
    """Pure correlation-input preparation result."""

    correlation_inputs: pd.DataFrame
    exclusion_report: dict[str, object]


@dataclass(frozen=True)
class PearsonOutputs:
    """Pure Pearson correlation result tables."""

    correlations: pd.DataFrame
    pvalues: pd.DataFrame


def pearson_pvalues(df_numeric: pd.DataFrame) -> pd.DataFrame:
    """Return pairwise Pearson p-values."""

    cols = list(df_numeric.columns)
    pvals = pd.DataFrame(float("nan"), index=cols, columns=cols)

    for i, c1 in enumerate(cols):
        for j, c2 in enumerate(cols):
            if j < i:
                continue
            x = df_numeric[c1]
            y = df_numeric[c2]
            mask = x.notna() & y.notna()
            if mask.sum() < 3 or x[mask].nunique() < 2 or y[mask].nunique() < 2:
                p = float("nan")
            else:
                _, p = pearsonr(x[mask], y[mask])
            pvals.loc[c1, c2] = p
            pvals.loc[c2, c1] = p

    return pvals


def apply_subject_exclusions_from_ids(
    data: pd.DataFrame,
    *,
    subject_column: str,
    excluded_ids: Iterable[object],
    context_label: str,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Remove manifest-listed subjects before correlation diagnostics.

    Sandve et al. (2013; see CITATIONS.md) support keeping analysis decisions
    centralized and auditable. This pure function accepts preloaded exclusion
    IDs and fails rather than filtering by row order when subject IDs are
    unavailable.
    """

    excluded_series = normalize_subject_ids(
        excluded_ids,
        column_name="excluded_ids",
    )
    excluded_set = set(excluded_series.astype(int).tolist())

    if subject_column not in data.columns:
        if excluded_set:
            raise ValueError(
                f"{context_label}: cannot apply nonempty subject exclusions because "
                f"subject column '{subject_column}' is absent from the input. "
                "Use an input CSV with subject IDs or pass an empty exclusion manifest "
                "only for a deliberately unfiltered diagnostic run."
            )
        return data.copy(), {
            "listed": 0,
            "matched": 0,
            "removed": 0,
            "remaining": int(len(data)),
            "matched_ids": [],
            "missing_ids": [],
        }

    filtered = data.copy()
    filtered["_subject_id_norm"] = normalize_subject_ids(
        filtered[subject_column],
        column_name=subject_column,
    ).to_numpy()

    present_ids = set(filtered["_subject_id_norm"].astype(int).tolist())
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


def print_exclusion_report(context_label: str, report: dict[str, object]) -> None:
    """Print an auditable subject-exclusion summary."""

    print(
        f"[exclude] {context_label}: "
        f"listed={report['listed']}, "
        f"matched={report['matched']}, "
        f"removed={report['removed']}, "
        f"remaining={report['remaining']}"
    )
    if report["missing_ids"]:
        print(
            f"[warn] {context_label}: exclusion IDs not present in dataset: "
            + ", ".join(str(subject_id) for subject_id in report["missing_ids"])
        )


def assert_subject_count_limit(
    *,
    report: dict[str, object],
    context_label: str,
    max_subjects: int,
) -> None:
    """Fail if too many participants remain after exclusions.

    Sandve et al. (2013; see CITATIONS.md) support fail-fast analysis
    safeguards for reproducible, auditable workflows.
    """

    remaining = int(report["remaining"])
    if remaining > max_subjects:
        raise ValueError(
            f"{context_label}: {remaining} subjects remain after exclusions, which "
            f"exceeds configured maximum of {max_subjects}. Confirm that "
            "data/config/excluded_subjects.json was applied before correlations."
        )


def correlation_input_columns(data: pd.DataFrame, subject_column: str) -> pd.DataFrame:
    """Return numeric correlation inputs with identifier columns removed."""

    df_numeric = data.apply(pd.to_numeric, errors="coerce")
    return df_numeric.drop(
        columns=[subject_column, "subject_id", "homer_subject"],
        errors="ignore",
    )


def add_recruitment_order_proxy(data: pd.DataFrame, *, subject_column: str) -> pd.DataFrame:
    """Add a normalized numeric subject ID for exploratory recruitment-order checks.

    Simmons et al. (2011; see CITATIONS.md) motivate explicit labeling of this
    analysis as exploratory: subject ID is only interpretable if it encodes
    recruitment/order, and correlations should not be read as primary evidence.
    """

    if subject_column not in data.columns:
        raise ValueError(
            "Cannot include subject ID in correlations because subject column "
            f"'{subject_column}' is absent from the input."
        )
    out = data.copy()
    out["recruitment_order_proxy"] = normalize_subject_ids(
        out[subject_column],
        column_name=subject_column,
    ).astype(float).to_numpy()
    return out


def prepare_correlation_inputs(
    data: pd.DataFrame,
    *,
    subject_column: str,
    excluded_ids: Iterable[object],
    context_label: str,
    include_subject_id_correlation: bool,
    max_subjects: int,
) -> PreparedCorrelationInputs:
    """Purely prepare numeric Pearson inputs after exclusions and guard checks."""

    filtered, exclusion_report = apply_subject_exclusions_from_ids(
        data,
        subject_column=subject_column,
        excluded_ids=excluded_ids,
        context_label=context_label,
    )
    assert_subject_count_limit(
        report=exclusion_report,
        context_label=context_label,
        max_subjects=max_subjects,
    )
    if include_subject_id_correlation:
        filtered = add_recruitment_order_proxy(filtered, subject_column=subject_column)
    return PreparedCorrelationInputs(
        correlation_inputs=correlation_input_columns(filtered, subject_column),
        exclusion_report=exclusion_report,
    )


def build_pearson_outputs(corr_inputs: pd.DataFrame) -> PearsonOutputs:
    """Purely build Pearson correlation and p-value tables from numeric inputs."""

    return PearsonOutputs(
        correlations=corr_inputs.corr(method="pearson"),
        pvalues=pearson_pvalues(corr_inputs),
    )


def obsolete_output_paths(out: Path) -> list[Path]:
    """Return known stale Spearman/generic covariate-correlation output paths."""

    return [out / filename for filename in sorted(OBSOLETE_OUTPUT_FILENAMES)]


def remove_obsolete_outputs(out: Path) -> None:
    """Remove known stale Spearman/generic covariate-correlation outputs."""

    for path in obsolete_output_paths(out):
        if path.exists():
            path.unlink()


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
    if sns is not None:
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
    else:
        corr_masked = corr.mask(mask)
        image = ax.imshow(corr_masked.to_numpy(dtype=float), cmap="coolwarm", vmin=-1, vmax=1)
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(corr.columns)
        ax.set_yticklabels(corr.index)
        fig.colorbar(image, ax=ax, shrink=0.8)
        if annot:
            for i in range(n):
                for j in range(n):
                    value = corr_masked.iat[i, j]
                    if pd.notna(value):
                        ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=8)
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


def write_covariate_correlation_outputs(corr_inputs: pd.DataFrame, out: Path) -> None:
    """Write covariates preset Pearson table, p-values, and heatmap."""

    remove_obsolete_outputs(out)
    outputs = build_pearson_outputs(corr_inputs)
    outputs.correlations.to_csv(out / "covariate_correlations_pearson.csv")
    outputs.pvalues.to_csv(out / "covariate_correlations_pearson_pvalues.csv")
    save_heatmap(
        outputs.correlations,
        out / "covariate_heatmap_pearson.png",
        title="Covariate Pearson correlations (encoded ordinal / totals)",
    )


def write_combined_correlation_outputs(corr_inputs: pd.DataFrame, out: Path) -> None:
    """Write combined preset Pearson table, p-values, and heatmap."""

    remove_obsolete_outputs(out)
    outputs = build_pearson_outputs(corr_inputs)
    outputs.correlations.to_csv(out / "covariate_correlation_analysis_pearson.csv", index=False)
    outputs.pvalues.to_csv(out / "covariate_correlation_analysis_pearson_pvalues.csv")
    save_heatmap(
        outputs.correlations,
        out / "covariate_correlation_heatmap_pearson.png",
        title="Correlation Heatmap (Pearson)",
    )


def run_covariate_preset(
    root: Path,
    *,
    input_path: Optional[Path],
    out_dir: Optional[Path],
    excluded_subjects_json: Path,
    subject_column: str,
    include_subject_id_correlation: bool = False,
    max_subjects: int = DEFAULT_MAX_SUBJECTS,
) -> None:
    in_path = input_path or (root / "covariate_outputs" / "covariates_clean.csv")
    out = out_dir or (root / "covariate_outputs")
    out.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(in_path)
    excluded_ids = load_excluded_subject_ids(excluded_subjects_json)
    prepared = prepare_correlation_inputs(
        df,
        subject_column=subject_column,
        excluded_ids=excluded_ids,
        context_label="covariates",
        include_subject_id_correlation=include_subject_id_correlation,
        max_subjects=max_subjects,
    )
    write_covariate_correlation_outputs(prepared.correlation_inputs, out)

    print(f"Loaded: {in_path}")
    print_exclusion_report("covariates", prepared.exclusion_report)
    print(f"Wrote:  {out / 'covariate_correlations_pearson.csv'}")
    print(f"Wrote:  {out / 'covariate_correlations_pearson_pvalues.csv'}")
    print(f"Wrote:  {out / 'covariate_heatmap_pearson.png'}")


def run_combined_preset(
    root: Path,
    *,
    input_path: Optional[Path],
    out_dir: Optional[Path],
    excluded_subjects_json: Path,
    subject_column: str,
    include_subject_id_correlation: bool = False,
    max_subjects: int = DEFAULT_MAX_SUBJECTS,
) -> None:
    in_path = input_path or (root / "data" / "tabular" / "generated_data" / "combined_sfv_data.csv")
    out = out_dir or (root / "data" / "tabular" / "generated_data")
    out.mkdir(parents=True, exist_ok=True)

    combined_data = pd.read_csv(in_path)
    excluded_ids = load_excluded_subject_ids(excluded_subjects_json)
    prepared = prepare_correlation_inputs(
        combined_data,
        subject_column=subject_column,
        excluded_ids=excluded_ids,
        context_label="combined",
        include_subject_id_correlation=include_subject_id_correlation,
        max_subjects=max_subjects,
    )
    write_combined_correlation_outputs(prepared.correlation_inputs, out)

    print(f"Loaded: {in_path}")
    print_exclusion_report("combined", prepared.exclusion_report)
    print(f"Wrote:  {out / 'covariate_correlation_analysis_pearson.csv'}")
    print(f"Wrote:  {out / 'covariate_correlation_analysis_pearson_pvalues.csv'}")
    print(f"Wrote:  {out / 'covariate_correlation_heatmap_pearson.png'}")


def main() -> None:
    root = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(description="Pearson correlation diagnostics + heatmaps.")
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
    parser.add_argument(
        "--excluded-subjects-json",
        type=Path,
        default=DEFAULT_EXCLUDED_SUBJECTS_JSON,
        help="Shared subject-exclusion manifest applied before correlations.",
    )
    parser.add_argument(
        "--subject-column",
        default=DEFAULT_SUBJECT_COLUMN,
        help="Column containing participant IDs for subject exclusions.",
    )
    parser.add_argument(
        "--include-subject-id-correlation",
        action="store_true",
        help=(
            "Include normalized numeric subject ID as recruitment_order_proxy "
            "for exploratory recruitment-order diagnostics."
        ),
    )
    parser.add_argument(
        "--max-subjects",
        type=int,
        default=DEFAULT_MAX_SUBJECTS,
        help="Fail if more than this many subjects remain after exclusions.",
    )
    args = parser.parse_args()

    if args.max_subjects < 1:
        raise ValueError("--max-subjects must be a positive integer.")

    if args.preset == "covariates":
        print("[INFO] Running covariate correlation analysis...")
        run_covariate_preset(
            root,
            input_path=args.input,
            out_dir=args.out_dir,
            excluded_subjects_json=args.excluded_subjects_json,
            subject_column=args.subject_column,
            include_subject_id_correlation=args.include_subject_id_correlation,
            max_subjects=args.max_subjects,
        )
    else:
        print("[INFO] Running combined correlation analysis...")
        run_combined_preset(
            root,
            input_path=args.input,
            out_dir=args.out_dir,
            excluded_subjects_json=args.excluded_subjects_json,
            subject_column=args.subject_column,
            include_subject_id_correlation=args.include_subject_id_correlation,
            max_subjects=args.max_subjects,
        )


if __name__ == "__main__":
    main()
