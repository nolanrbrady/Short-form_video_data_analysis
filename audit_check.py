
import pandas as pd
import numpy as np

def check_audit_file(filepath):
    print(f"Checking {filepath}...")
    df = pd.read_csv(filepath)
    
    issues = []

    for idx, row in df.iterrows():
        score = row['score']
        method = row['method']
        fuzzy_ratio = row['fuzzy_ratio']
        response_norm = str(row['response_norm']) if pd.notna(row['response_norm']) else ""
        key_norm = str(row['key_norm']) if pd.notna(row['key_norm']) else ""
        
        # 1. Exact Match Consistency
        if method == 'exact':
            if score != 1:
                issues.append(f"Row {idx}: Method is exact but score is {score}")
            if response_norm != key_norm:
                # Sometimes key_norm might be one of multiple options, but in this script 
                # key_norm in audit seems to be the matched key. 
                # Let's verify if response_norm is truly equal to key_norm.
                issues.append(f"Row {idx}: Method is exact but norms differ: '{response_norm}' vs '{key_norm}'")

        # 2. Fuzzy Match Consistency
        elif method == 'fuzzy':
            if score != 1:
                issues.append(f"Row {idx}: Method is fuzzy but score is {score}")
            if pd.isna(fuzzy_ratio) or fuzzy_ratio < 0.90:
                issues.append(f"Row {idx}: Method is fuzzy but ratio is {fuzzy_ratio} (expected >= 0.90)")
            if response_norm == key_norm:
                 issues.append(f"Row {idx}: Method is fuzzy but norms are identical (should be exact): '{response_norm}'")

        # 3. No Match Consistency
        elif method == 'no_match':
            if score != 0:
                issues.append(f"Row {idx}: Method is no_match but score is {score}")
            if pd.notna(fuzzy_ratio) and fuzzy_ratio >= 0.90:
                 # This is tricky. If length difference is small (<=3), we force exact match.
                 # So high ratio could still be a mismatch if string is short.
                 pass 

        # 4. Empty / No Key Consistency
        elif method in ['empty', 'no_key']:
            if score != 0:
                issues.append(f"Row {idx}: Method is {method} but score is {score}")

    if not issues:
        print("  No logic inconsistencies found in scoring.")
    else:
        print(f"  Found {len(issues)} inconsistencies:")
        for i in issues[:10]:
            print("    " + i)
        if len(issues) > 10:
            print(f"    ... and {len(issues)-10} more.")

    # Summary Stats
    print("\n  Summary:")
    print(df['method'].value_counts())
    print("\n  Fuzzy Matches (Sample):")
    fuzzy_rows = df[df['method'] == 'fuzzy']
    if not fuzzy_rows.empty:
        print(fuzzy_rows[['response_raw', 'key_answer', 'response_norm', 'key_norm', 'fuzzy_ratio']].head())
    else:
        print("  No fuzzy matches found.")
    print("-" * 40)

if __name__ == "__main__":
    check_audit_file("demographic/recall_assessment_audit_post.csv")
    check_audit_file("demographic/recall_assessment_audit_pre.csv")
