"""Analyze demographic and survey data from the Short Form Video study."""

from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import pandas as pd

# =============================================================================
# Likert Scale Mappings
# =============================================================================

LIKERT_0_TO_3 = {
    "Not at all": 0,
    "Several days": 1,
    "More than half the days": 2,
    "Nearly every day": 3,
}

ASRS_SCALE = {
    "Never": 0,
    "Rarely": 1,
    "Sometimes": 2,
    "Often": 3,
    "Very often": 4,
}

YANG_SCALE = {
    "Completely Disagree": 1,
    "Somewhat Disagree": 2,
    "Neutral": 3,
    "Somewhat Agree": 4,
    "Completely Agree": 5,
}

# Minimum proportion of items required for valid scale scores
MIN_ITEM_PROPORTION = 0.75


# =============================================================================
# Column Lookup Helpers
# =============================================================================


def col_by_label(df: pd.DataFrame, label: str) -> Tuple:
    """Find column tuple by matching the second level (label)."""
    for col in df.columns:
        if col[1] == label:
            return col
    raise KeyError(label)


def col_by_qid(df: pd.DataFrame, qid: str) -> Tuple:
    """Find column tuple by matching the first level (question ID)."""
    for col in df.columns:
        if col[0] == qid:
            return col
    raise KeyError(qid)


# =============================================================================
# Summary Formatting Helpers
# =============================================================================


def categorical_summary(series: pd.Series, label: str, n_total: int) -> str:
    """Generate a formatted summary for a categorical variable."""
    counts = series.value_counts(dropna=False).sort_index()
    lines = [f"{label} (N={series.notna().sum()}):"]
    for level, count in counts.items():
        level_str = "Missing / No response" if pd.isna(level) else str(level)
        pct = (count / n_total) * 100
        lines.append(f"  {level_str}: {count} ({pct:.1f}%)")
    return "\n".join(lines)


def continuous_summary(series: pd.Series, label: str) -> str:
    """Generate a formatted summary for a continuous variable."""
    s = series.dropna()
    if s.empty:
        return f"{label}: no non-missing data"
    return (
        f"{label}: N={s.shape[0]}, mean={s.mean():.2f}, SD={s.std():.2f}, "
        f"min={s.min():.2f}, max={s.max():.2f}"
    )


# =============================================================================
# Plotting Helpers
# =============================================================================


def setup_plot_style() -> None:
    """Apply consistent plot styling."""
    try:
        plt.style.use("seaborn-v0_8")
    except OSError:
        pass


def save_histogram(
    data: pd.Series,
    xlabel: str,
    title: str,
    filepath: Path,
    bins: int = 10,
    color: str = "steelblue",
    figsize: Tuple[int, int] = (6, 4),
) -> None:
    """Create and save a histogram figure."""
    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(data, bins=bins, color=color, edgecolor="black")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Number of participants")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(filepath, dpi=300)
    plt.close(fig)


def save_bar_plot(
    series: pd.Series,
    title: str,
    xlabel: str,
    filepath: Path,
) -> None:
    """Create and save a bar plot for categorical data."""
    counts = series.value_counts(dropna=False)
    labels = [
        "Missing / No response" if pd.isna(level) else str(level)
        for level in counts.index
    ]
    values = counts.values

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(range(len(values)), values, color="steelblue", edgecolor="black")
    ax.set_xticks(range(len(values)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Number of participants")
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(filepath, dpi=300)
    plt.close(fig)


# =============================================================================
# Data Loading and Column Extraction
# =============================================================================


def load_data(data_path: Path) -> pd.DataFrame:
    """Load the survey data with multi-level headers."""
    return pd.read_csv(data_path, header=[0, 1, 2])


def extract_demographic_columns(df: pd.DataFrame) -> dict:
    """Extract demographic variable columns from the dataframe."""
    return {
        "age": col_by_label(df, "What is your current age in years?"),
        "degree": col_by_label(df, "What is your highest degree completed? - Selected Choice"),
        "income": col_by_label(
            df, "What was your total household income before taxes during the past year?"
        ),
        "ethnicity": col_by_label(df, "Are you of Hispanic or Latino descent?"),
        "race": col_by_label(
            df,
            "Regardless of your answer to the prior question, please indicate how you identify yourself - Selected Choice",
        ),
        "race_other": col_by_label(
            df,
            "Regardless of your answer to the prior question, please indicate how you identify yourself - Other - Text",
        ),
        "sex_birth": col_by_label(df, "What is your sex at birth? - Selected Choice"),
        "gender": col_by_label(df, "What is your gender? - Selected Choice"),
        "relationship": col_by_label(df, "What is your relationship status? - Selected Choice"),
        "psych_dx": col_by_label(
            df,
            "Have you ever been diagnosed with ADHD, Autism spectrum disorder, anxiety, depression or any other psychiatric diagnosis? - Selected Choice",
        ),
        "adhd_med": col_by_label(
            df, "Are you currently taking any medications for ADHD? - Selected Choice"
        ),
    }


def extract_sfv_usage_columns(df: pd.DataFrame) -> dict:
    """Extract short-form video usage columns from the dataframe."""
    return {
        "time_of_day": col_by_qid(df, "Q79"),
        "freq": col_by_qid(df, "Q80"),
        "duration": col_by_qid(df, "Q81"),
        "duration_other": col_by_qid(df, "Q81_5_TEXT"),
        "weekday_weekend": col_by_qid(df, "Q82"),
        "platform": col_by_qid(df, "Q83"),
        "platform_other": col_by_qid(df, "Q83_5_TEXT"),
    }


def extract_language_columns(df: pd.DataFrame) -> dict:
    """Extract language exposure columns from the dataframe."""
    return {
        "arabic": col_by_qid(df, "Q73"),
        "hawaiian": col_by_qid(df, "Q74"),
        "russian": col_by_qid(df, "Q75"),
        "swedish": col_by_qid(df, "Q76"),
    }


def extract_scale_columns(df: pd.DataFrame) -> dict:
    """Extract psychological scale columns from the dataframe."""
    gad_cols = [c for c in df.columns if c[1] == "GAD"]
    gad_text_cols = [c for c in df.columns if "Over the last two weeks" in c[1]]

    return {
        "phq": [c for c in df.columns if c[1] == "PHQ-9"],
        "asrs": [c for c in df.columns if c[1] == "ASRS"],
        "gad": gad_cols + gad_text_cols,
        "yang_prob": [c for c in df.columns if c[1] == "Yang Problematic Use"],
        "yang_mot": [c for c in df.columns if c[1] == "Yang Motivation"],
    }


# =============================================================================
# Data Processing
# =============================================================================


def merge_other_responses(primary: pd.Series, other: pd.Series) -> pd.Series:
    """Replace 'Other' responses with text from the corresponding other column."""
    result = primary.copy()
    return result.where(result != "Other", other.where(other.notna(), result))


def compute_scale_scores(
    df: pd.DataFrame,
    columns: list,
    scale_map: dict,
) -> Tuple[pd.Series, pd.Series]:
    """
    Compute total and mean item scores for a psychological scale.

    Returns tuple of (total_score, mean_item_score).
    """
    numeric_data = df[columns].replace(scale_map)
    min_items = int(len(columns) * MIN_ITEM_PROPORTION)
    total_score = numeric_data.sum(axis=1, min_count=min_items)
    mean_item_score = total_score / len(columns)
    return total_score, mean_item_score


# =============================================================================
# Report Generation
# =============================================================================


def generate_summary_report(
    df: pd.DataFrame,
    demographics: dict,
    sfv_usage: dict,
    language_vars: dict,
    scale_scores: dict,
) -> str:
    """Generate the complete demographic summary report."""
    summary = []
    n = len(df)

    # Total participants
    summary.append(f"Total participants with survey data: {n}")

    # Age statistics
    age_valid = demographics["age"].dropna()
    summary.append(
        "Age (years): N={n_age}, mean={mean:.2f}, SD={sd:.2f}, min={min_val:.1f}, max={max_val:.1f}".format(
            n_age=age_valid.shape[0],
            mean=age_valid.mean(),
            sd=age_valid.std(),
            min_val=age_valid.min(),
            max_val=age_valid.max(),
        )
    )

    # Demographic variables
    demographic_vars = [
        (demographics["degree"], "Highest degree completed"),
        (demographics["income"], "Household income (past year)"),
        (demographics["ethnicity"], "Hispanic / Latino ethnicity"),
        (demographics["race"], "Racial identification"),
        (demographics["sex_birth"], "Sex at birth"),
        (demographics["gender"], "Gender identity"),
        (demographics["relationship"], "Relationship status"),
        (demographics["psych_dx"], "Any prior psychiatric diagnosis"),
        (demographics["adhd_med"], "Currently taking ADHD medication"),
    ]
    for series, label in demographic_vars:
        summary.append(categorical_summary(series, label, n))

    # SFV usage variables
    sfv_vars = [
        (sfv_usage["time_of_day"], "Typical time of day for SFV use"),
        (sfv_usage["freq"], "Frequency of SFV use"),
        (sfv_usage["duration"], "Daily duration of SFV use"),
        (sfv_usage["weekday_weekend"], "SFV use: weekdays vs. weekends"),
        (sfv_usage["platform"], "Primary SFV platform"),
    ]
    for series, label in sfv_vars:
        summary.append(categorical_summary(series, label, n))

    # Language exposure variables
    language_vars_list = [
        (language_vars["arabic"], "Arabic language exposure"),
        (language_vars["hawaiian"], "Hawaiian language exposure"),
        (language_vars["russian"], "Russian language exposure"),
        (language_vars["swedish"], "Swedish language exposure"),
    ]
    for series, label in language_vars_list:
        summary.append(categorical_summary(series, label, n))

    # Psychological scale scores
    scale_summaries = [
        (scale_scores["phq_total"], "PHQ total score (sum of 8 items, 0-24 range)"),
        (scale_scores["phq_mean"], "PHQ mean item score (0-3 scale)"),
        (scale_scores["gad_total"], "GAD total score (sum of 7 items, 0-21 range)"),
        (scale_scores["gad_mean"], "GAD mean item score (0-3 scale)"),
        (scale_scores["asrs_total"], "ASRS total score (sum of 18 items, 0-72 range)"),
        (scale_scores["asrs_mean"], "ASRS mean item score (0-4 scale)"),
        (
            scale_scores["yang_prob_total"],
            "Yang Short-Form Video Problematic Use total (21 items, 21-105 range)",
        ),
        (
            scale_scores["yang_prob_mean"],
            "Yang Short-Form Video Problematic Use mean item score (1-5 scale)",
        ),
        (
            scale_scores["yang_mot_total"],
            "Yang Short-Form Video Motivation total (22 items, 22-110 range)",
        ),
        (
            scale_scores["yang_mot_mean"],
            "Yang Short-Form Video Motivation mean item score (1-5 scale)",
        ),
    ]
    for series, label in scale_summaries:
        summary.append(continuous_summary(series, label))

    return "\n\n".join(summary)


# =============================================================================
# Figure Generation
# =============================================================================


def generate_figures(
    root: Path,
    demographics: dict,
    sfv_usage: dict,
    scale_scores: dict,
) -> None:
    """Generate all summary figures."""
    setup_plot_style()

    # Figure 1: Age distribution
    age_valid = demographics["age"].dropna()
    save_histogram(
        age_valid,
        xlabel="Age (years)",
        title="Distribution of participant ages",
        filepath=root / "fig1_age_distribution.png",
    )

    # Figures 2-5: Categorical distributions
    bar_plots = [
        (demographics["gender"], "Gender identity", "Gender category", "fig2_gender_distribution.png"),
        (demographics["race"], "Racial identification", "Race / ethnicity category", "fig3_race_ethnicity_distribution.png"),
        (demographics["degree"], "Highest degree completed", "Education level", "fig4_education_distribution.png"),
        (sfv_usage["duration"], "Daily time spent on short-form video", "Daily SFV duration", "fig5_sfv_daily_duration.png"),
    ]
    for series, title, xlabel, filename in bar_plots:
        save_bar_plot(series, title, xlabel, root / filename)

    # Figure 6: PHQ and GAD distributions
    phq_valid = scale_scores["phq_total"].dropna()
    gad_valid = scale_scores["gad_total"].dropna()
    if not phq_valid.empty and not gad_valid.empty:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)

        axes[0].hist(
            phq_valid,
            bins=range(int(phq_valid.min()), int(phq_valid.max()) + 2),
            color="steelblue",
            edgecolor="black",
        )
        axes[0].set_xlabel("PHQ total score")
        axes[0].set_ylabel("Number of participants")
        axes[0].set_title("Distribution of PHQ total scores")

        axes[1].hist(
            gad_valid,
            bins=range(int(gad_valid.min()), int(gad_valid.max()) + 2),
            color="indianred",
            edgecolor="black",
        )
        axes[1].set_xlabel("GAD total score")
        axes[1].set_title("Distribution of GAD total scores")

        fig.tight_layout()
        fig.savefig(root / "fig6_phq_gad_distributions.png", dpi=300)
        plt.close(fig)

    # Figure 7: Yang problematic use distribution
    yang_prob_valid = scale_scores["yang_prob_total"].dropna()
    save_histogram(
        yang_prob_valid,
        xlabel="Yang problematic-use total score",
        title="Distribution of short-form video problematic-use scores",
        filepath=root / "fig7_yang_problematic_distribution.png",
        color="mediumpurple",
    )


# =============================================================================
# Main Entry Point
# =============================================================================


def main() -> None:
    """Run the demographic analysis pipeline."""
    root = Path(__file__).resolve().parent
    data_path = root / "SFV_demo_data.csv"

    # Load data
    df = load_data(data_path)

    # Extract column references
    demo_cols = extract_demographic_columns(df)
    sfv_cols = extract_sfv_usage_columns(df)
    lang_cols = extract_language_columns(df)
    scale_cols = extract_scale_columns(df)

    # Process demographic variables
    demographics = {
        "age": pd.to_numeric(df[demo_cols["age"]], errors="coerce"),
        "degree": df[demo_cols["degree"]],
        "income": df[demo_cols["income"]],
        "ethnicity": df[demo_cols["ethnicity"]],
        "race": merge_other_responses(df[demo_cols["race"]], df[demo_cols["race_other"]]),
        "sex_birth": df[demo_cols["sex_birth"]],
        "gender": df[demo_cols["gender"]],
        "relationship": df[demo_cols["relationship"]],
        "psych_dx": df[demo_cols["psych_dx"]],
        "adhd_med": df[demo_cols["adhd_med"]],
    }

    # Process SFV usage variables
    sfv_usage = {
        "time_of_day": df[sfv_cols["time_of_day"]],
        "freq": df[sfv_cols["freq"]],
        "duration": merge_other_responses(df[sfv_cols["duration"]], df[sfv_cols["duration_other"]]),
        "weekday_weekend": df[sfv_cols["weekday_weekend"]],
        "platform": merge_other_responses(df[sfv_cols["platform"]], df[sfv_cols["platform_other"]]),
    }

    # Process language exposure variables
    language_vars = {
        "arabic": df[lang_cols["arabic"]],
        "hawaiian": df[lang_cols["hawaiian"]],
        "russian": df[lang_cols["russian"]],
        "swedish": df[lang_cols["swedish"]],
    }

    # Compute psychological scale scores
    phq_total, phq_mean = compute_scale_scores(df, scale_cols["phq"], LIKERT_0_TO_3)
    gad_total, gad_mean = compute_scale_scores(df, scale_cols["gad"], LIKERT_0_TO_3)
    asrs_total, asrs_mean = compute_scale_scores(df, scale_cols["asrs"], ASRS_SCALE)
    yang_prob_total, yang_prob_mean = compute_scale_scores(df, scale_cols["yang_prob"], YANG_SCALE)
    yang_mot_total, yang_mot_mean = compute_scale_scores(df, scale_cols["yang_mot"], YANG_SCALE)

    scale_scores = {
        "phq_total": phq_total,
        "phq_mean": phq_mean,
        "gad_total": gad_total,
        "gad_mean": gad_mean,
        "asrs_total": asrs_total,
        "asrs_mean": asrs_mean,
        "yang_prob_total": yang_prob_total,
        "yang_prob_mean": yang_prob_mean,
        "yang_mot_total": yang_mot_total,
        "yang_mot_mean": yang_mot_mean,
    }

    # Generate summary report
    summary_text = generate_summary_report(
        df, demographics, sfv_usage, language_vars, scale_scores
    )
    print(summary_text)

    summary_path = root / "demographics_summary.txt"
    summary_path.write_text(summary_text, encoding="utf-8")

    # Generate figures
    generate_figures(root, demographics, sfv_usage, scale_scores)


if __name__ == "__main__":
    main()
