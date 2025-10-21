
import pandas as pd
import os
import statsmodels.formula.api as smf
import statsmodels.api as sm
import numpy as np
from statsmodels.stats.anova import AnovaRM

def run_two_by_two_anova(df):
    """
    2×2 Repeated-Measures ANOVA
    Use per-subject cell means for AnovaRM
    """

    agg = (
        df.groupby(["subject_id", "Content", "Length"], as_index=False)
        .agg(Rating_mean=("Rating", "mean"),
            n=("Rating", "size"))
    )

    rm = AnovaRM(
        data=agg,
        depvar="Rating_mean",
        subject="subject_id",
        within=["Content", "Length"]
    ).fit()

    anova_rm_table = rm.anova_table.copy()
    return anova_rm_table


def run_fixed_effect_ols(df):
    """
    Runs a fixed-effects model to get the ANOVA table and cell means.
    """
    # Step 1: Run a regression that includes subject-level controls.
    model = smf.ols("Rating ~ C(Content)*C(Length) + C(subject_id)", data=df).fit()

    # Step 2: Convert the regression results into a classic ANOVA table.
    anova_results = sm.stats.anova_lm(model, typ=3)

    # Step 3: Calculate the average rating for each of the four conditions.
    mean_ratings = df.groupby(["Content", "Length"])["Rating"].mean().reset_index()

    return anova_results, mean_ratings

def run_mixed_linear_model(df):
    # The formula defines the fixed effects
    formula = "Rating ~ C(Category)"

    # The 'groups' argument defines the random effect (a random intercept for each subject)
    model = smf.mixedlm(formula, data=df, groups=df["subject_id"]).fit()
    return model

