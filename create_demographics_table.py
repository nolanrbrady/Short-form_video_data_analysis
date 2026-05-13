"""
Using the table in ./data/tabular/homer3_betas_plus_combined_sfv_data_inner_join.csv we want to create a demographics table.
The table should report n, mean, std, min, 25%, 50%, 75%, and max for each of the columns that are not specific to the beta values.

Ensure that excluded subjects (./data/config/excluded_subjects.json) are removed from the dataset.
"""

import json
import pandas as pd


# Load the data
data = pd.read_csv('./data/tabular/homer3_betas_plus_combined_sfv_data_inner_join.csv')

# Remove excluded subjects
with open('./data/config/excluded_subjects.json') as f:
    excluded_subjects = json.load(f)

# Remove rows where 'homer_subject' is in the list of excluded subjects
data = data[~data['homer_subject'].isin(excluded_subjects)]

# Select only the columns that are not specific to the beta values
demographics_columns = [col for col in data.columns if not col.startswith('S0') and col not in ['subject_id', 'homer_subject']]
print(demographics_columns)

# Create the demographics table
demographics_table = data[demographics_columns].describe().transpose()
demographics_table.rename(columns={'count': 'n'}, inplace=True)
print(demographics_table)