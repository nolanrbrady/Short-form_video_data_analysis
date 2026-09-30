# Historical validation instructions

Archived on 2026-09-30. These tests are retained for provenance, not invoked by
current result validation. Restore the original file layout and dependencies in
an isolated workspace before attempting historical test commands.

## Correlation Follow-up Pipeline (R): reverted predictor-by-condition validation

Runs `analyze_correlational_relationships.R` on a synthetic merged dataset and verifies:
- The old predictor-by-condition plan schema still runs cleanly against the reverted script
- One output CSV is written with the expected predictor x target x condition rows
- Pearson summary statistics and BH-FDR family bookkeeping are present
- A known positive synthetic row recovers a perfect Pearson correlation
- Families are defined across the four condition-specific rows for each predictor x neural target
- The correlation output folder is rebuilt at run start, removing stale CSVs and PNGs
- Fail-hard behavior for duplicate IDs and malformed required inputs

Command:

```bash
Rscript tests/validate_correlational_relationships_r.R
```

## Covariate Correlation Diagnostics (Python): shared subject-exclusion validation

Runs `covariate_correlation_analysis.py` helper functions on synthetic data and verifies:
- manifest-listed subjects are removed before Spearman correlations are computed
- nonempty exclusion manifests fail hard when the input lacks the configured subject ID column
- normalized subject ID can be included as `recruitment_order_proxy` only when explicitly requested
- the analysis fails if more than 48 subjects remain after exclusions
Command:

```bash
pytest -q tests/test_covariate_correlation_analysis.py
```

## Beta Discrepancy Plotting (Python): channel-vs-ROI descriptive validation

Runs `plot_beta_discrepancy_dynamics.py` on a synthetic merged beta table and verifies:
- shared subject exclusions are applied before plotting
- exact-zero betas are treated as pruned/missing when configured
- ROI means use the available non-missing member channels rather than filling missing values
- complete-case panel counts match the intended channel-vs-ROI comparison logic
- the composite PNG and audit CSV are both created

Command:

```bash
python tests/validate_beta_discrepancy_plot_py.py
```
