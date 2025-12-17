import numpy as np
import pandas as pd


"""
Process the engagement data from the Short Form Video study. We want to consolidate the into a few categories per subject.

Categories:
- Long-form Education
- Long-form Entertainment
- Short-form Education
- Short-form Entertainment

We also want to get the data isolating the specifics of the stimulus.
- Education
- Entertainment
- Long-Form
- Short-Form

The result should be that each subject has a row with the columns being the above categories

Input:
- ./demographic/combined_engagement_data.csv

Output:
- ./data/tabular/engagement_data_processed.csv
"""

# combined_engagement_data.csv is the aggregated engagement data from the Short Form Video study
# It relies on combine_engagement.py to be run first.
engagement_data = pd.read_csv("./demographic/combined_engagement_data.csv")
processed_data = pd.DataFrame()

def get_average_rating(subject_data, category):
    return subject_data[subject_data["Category"] == category]["Rating"].mean()

def get_condition_averages(subject_data, condition):
    """
    condition: one of {"Long-Form", "Short-Form", "Education", "Entertainment"}
    Returns the mean Rating for that subject restricted to the condition,
    marginalizing over the other dimension.
    """
    cond = condition.strip().lower()
    cat = subject_data["Category"].astype(str).str.lower()

    if cond == "long-form":
        mask = cat.str.contains("long-form", na=False)
    elif cond == "short-form":
        mask = cat.str.contains("short-form", na=False)
    elif cond == "education":
        mask = cat.str.contains("educ", na=False)          # education/educational
    elif cond == "entertainment":
        mask = cat.str.contains("entertain", na=False)
    else:
        raise ValueError(f"Unknown condition: {condition}")

    return subject_data.loc[mask, "Rating"].mean()
    

for subject_id in engagement_data["subject_id"].unique():
    subject_data = engagement_data[engagement_data["subject_id"] == subject_id]
    print(subject_data)

    # Get the class averages for the subject
    lf_education = get_average_rating(subject_data, "Long-Form Educational")
    lf_entertainment = get_average_rating(subject_data, "Long-Form Entertainment")
    sf_education = get_average_rating(subject_data, "Short-Form Education")
    sf_entertainment = get_average_rating(subject_data, "Short-Form Entertainment")
    print(lf_education, lf_entertainment, sf_education, sf_entertainment)

    # Get the condition averages for the subject
    long_form = get_condition_averages(subject_data, "Long-Form")
    short_form = get_condition_averages(subject_data, "Short-Form")
    education = get_condition_averages(subject_data, "Education")
    entertainment = get_condition_averages(subject_data, "Entertainment")
    print(long_form, short_form, education, entertainment)

    # Add the data to the processed data
    new_row = pd.DataFrame([{
        "subject_id": subject_id,
        "lf_education_engagement": lf_education,
        "lf_entertainment_engagement": lf_entertainment,
        "sf_education_engagement": sf_education,
        "sf_entertainment_engagement": sf_entertainment,
        "long_form_engagement": long_form,
        "short_form_engagement": short_form,
        "education_engagement": education,
        "entertainment_engagement": entertainment,
    }])
    processed_data = pd.concat([processed_data, new_row], ignore_index=True)

processed_data.to_csv("./data/tabular/engagement_data_processed.csv", index=False)