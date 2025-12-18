"""

Combines the engagement and socio-demographic data into a single file.
Future work will add in the beta values from the GLM analysis.

Input:
- ./data/tabular/engagement_data_processed.csv
- ./data/tabular/socio_demographic_data_processed.csv

Output:
- ./data/tabular/combined_sfv_data.csv
"""

import pandas as pd


engagement_data = pd.read_csv("./data/tabular/engagement_data_processed.csv")
socio_demographic_data = pd.read_csv("./data/tabular/socio_demographic_data_processed.csv")
recall_assessment_data = pd.read_csv("./data/tabular/recall_assessment_score_diffs.csv")

# Merge the data on the subject_id column
combined_data = pd.merge(engagement_data, socio_demographic_data, on="subject_id", how="inner")
combined_data = pd.merge(combined_data, recall_assessment_data, on="subject_id", how="inner")

# Save the combined data
combined_data.to_csv("./data/tabular/combined_sfv_data.csv", index=False)

print(combined_data)