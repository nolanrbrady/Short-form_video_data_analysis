"""
This file is processing the engagement data collected during participants fNIRs sessions.

Author: Nolan Brady
"""

import pandas as pd
import os
import statsmodels.formula.api as smf
import statsmodels.api as sm
import numpy as np
from statsmodels.stats.anova import AnovaRM
from engagement_stats import run_two_by_two_anova, run_fixed_effect_ols, run_mixed_linear_model

DATA_DIR = "../../Engagement"

all_data = pd.DataFrame()
for file in os.listdir(DATA_DIR):
    sub_id = file.split("_engagement")[0].split("_")[-1]
    df = pd.read_csv(f"{DATA_DIR}/{file}")
    df = df[['Rating', 'Trigger', 'Category']]
    df['subject_id'] = int(sub_id)
    all_data = pd.concat([all_data, df], ignore_index=True)

# Write the aggregated enagement data to a CSV
all_data.to_csv("combined_engagement_data.csv")

#================================
# Statistical Analysis
#================================

df = all_data

# Derive factors from Category
def parse_content(cat):
    return "Education" if "Educat" in cat else "Entertainment"

def parse_length(cat: str) -> str:
    return "Short" if "Short-Form" in cat else "Long"

df["Content"] = df["Category"].apply(parse_content).astype("category")
df["Length"] = df["Category"].apply(parse_length).astype("category")
df["Condition"] = df["Category"].astype("category")
df["subject_id"] = df["subject_id"].astype("category")
df['Rating'] = pd.to_numeric(df['Rating'])
print(df.head())
# 2 x 2 ANOVA
results_table = run_two_by_two_anova(df)
print(results_table)

# # Fixed Effects Ordinary Least Squares Model
anova_result, mean_rating = run_fixed_effect_ols(df)
print("ANOVA Results: ")
print(anova_result)

# Mixed linear model
model = run_mixed_linear_model(df)
print("Mixed Effects Model: ")
print(model.summary())

