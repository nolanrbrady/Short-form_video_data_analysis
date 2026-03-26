# Tests / Validation

This folder contains lightweight validation harnesses that aim to verify scientific/analysis integrity against `ANALYSIS_SPEC.md`.

## Pipeline C (R): synthetic end-to-end validation

Runs `analyze_format_content_lmm_channelwise.R` on a synthetic dataset that:
- Uses Homer-style subject IDs (`sub_0001`) and combined IDs (`0001`)
- Uses varying subject ages and checks the age-adjusted omnibus model against a direct reference fit
- Includes both HbO/HbR and multiple channels
- Includes a pruned beta (`0`) to verify pruned-channel missingness + complete-case logic
- Verifies coefficient recovery (within tolerance) for known ground-truth fixed effects
- Verifies explicit inferential TP/TN outcomes (known significant interaction channel vs known null channel)
- Verifies post-hoc gating only triggers for interaction-FDR-significant channels
- Verifies duplicate `subject_id` fails hard
- Verifies missing/non-numeric `age` fails hard

Command:

```bash
Rscript tests/validate_pipeline_c_r.R
```

## Pipeline C ROI (R): synthetic end-to-end validation

Runs `analyze_format_content_lmm_roi.R` on a synthetic dataset that:
- Uses Homer-style subject IDs (`sub_0001`) and combined IDs (`0001`)
- Uses varying subject ages and checks the age-adjusted omnibus model against a direct reference fit
- Defines ROIs from a strict JSON ROI map (`ROI -> [channels]`)
- Verifies ROI mean aggregation over available channels (`0` treated as pruned/missing)
- Verifies complete-case behavior when all ROI channels are pruned for a condition
- Verifies coefficient recovery for ROI-level known generating effects
- Verifies explicit inferential TP/TN outcomes (known significant ROI interaction vs known null ROI interaction)
- Verifies BH-FDR correctness across ROIs (per chrom/effect family)
- Verifies interaction-gated post-hoc behavior
- Verifies fail-hard behavior for invalid JSON, missing/non-numeric `age`, missing ROI channels, and overlapping ROI assignments

Command:

```bash
Rscript tests/validate_pipeline_c_roi_r.R
```

## Retention Pipeline (R): deterministic + synthetic validation

Runs `analyze_retention_format_content_lmm.R` on generated retention datasets and verifies:
- Deterministic analytic ground truth for Length/Content/Interaction contrasts
- End-to-end coefficient recovery for known generating parameters
- Direct agreement with an age-adjusted reference omnibus fit
- Explicit TP/TN inferential checks using strong-effect and null-effect synthetic datasets
- Holm-correction correctness across the 3 omnibus effects
- Interaction-gated post-hoc behavior (on/off cases)
- Zero-as-valid retention handling (not treated as missing)
- Complete-case subject drop only for true `NA`
- Fail-hard integrity checks (duplicate IDs, missing required columns, missing/non-numeric `age`, non-numeric retention values)

Command:

```bash
Rscript tests/validate_retention_pipeline_r.R
```

## Engagement Pipeline (R): deterministic + synthetic validation

Runs `analyze_engagement_format_content_lmm.R` on generated engagement datasets and verifies:
- Deterministic analytic ground truth for Length/Content/Interaction contrasts
- End-to-end coefficient recovery for known generating parameters
- Direct agreement with an age-adjusted reference omnibus fit
- Explicit TP/TN inferential checks using strong-effect and null-effect synthetic datasets
- Holm-correction correctness across the 3 omnibus effects
- Interaction-gated post-hoc behavior (on/off cases)
- Zero-as-valid engagement handling (not treated as missing)
- Complete-case subject drop only for true `NA`
- Fail-hard integrity checks (duplicate IDs, missing required columns, missing/non-numeric `age`, non-numeric engagement values)

Command:

```bash
Rscript tests/validate_engagement_pipeline_r.R
```

## Correlation Follow-up Pipeline (R): targeted pairwise validation

Runs `analyze_correlational_relationships.R` on a synthetic merged dataset and verifies:
- Requested predictor discovery (`sfv_daily_duration`, `sfv_frequency`, `diff_*`, `*_engagement`)
- Requested neural target coverage for `S04_D02` and the specified ROI/chrom pairs across all 4 conditions
- Pairwise-complete-case behavior only for the affected correlation, with no imputation
- ROI averaging across available non-missing channels
- Explicit positive/negative/null correlation cases
- BH-FDR correctness within each configured family (default study config: neural-target x chromophore x predictor, i.e. across the 4 conditions)
- Output artifact creation (CSV plus one figure per tested pair)
- Fail-hard behavior for malformed ROI JSON, duplicate IDs, and non-numeric required inputs

Command:

```bash
Rscript tests/validate_correlational_relationships_r.R
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

## Type-I Error Calibration (R): Monte Carlo null simulations across all pipelines

Runs repeated null-effect synthetic datasets through all four inferential scripts:
- `analyze_format_content_lmm_channelwise.R`
- `analyze_format_content_lmm_roi.R`
- `analyze_retention_format_content_lmm.R`
- `analyze_engagement_format_content_lmm.R`

It estimates empirical false-positive rates from adjusted p-values (`p_fdr`) and fails when any
pipeline/effect exceeds a configured upper bound.
The synthetic generators now include varying `age` so the required omnibus covariate path is exercised during calibration.

Command:

```bash
Rscript tests/calibrate_type1_error_r.R
```

Useful overrides:

```bash
# Faster local smoke run
Rscript tests/calibrate_type1_error_r.R --n_reps 20 --type1_upper_bound 0.20

# Stricter run for manuscript QA
Rscript tests/calibrate_type1_error_r.R --n_reps 200 --type1_upper_bound 0.10
```

## Type-II Error Calibration (R): Monte Carlo power simulations across all pipelines

Runs repeated non-null synthetic datasets through all four inferential scripts:
- `analyze_format_content_lmm_channelwise.R`
- `analyze_format_content_lmm_roi.R`
- `analyze_retention_format_content_lmm.R`
- `analyze_engagement_format_content_lmm.R`

It estimates empirical power and type-II error from adjusted p-values (`p_fdr`) and fails when any
pipeline/effect exceeds a configured type-II upper bound.
The synthetic generators now include varying `age` so the required omnibus covariate path is exercised during calibration.

Command:

```bash
Rscript tests/calibrate_type2_error_r.R
```

Useful overrides:

```bash
# Faster local smoke run
Rscript tests/calibrate_type2_error_r.R --n_reps 20 --type2_upper_bound 0.60

# Stricter run for manuscript QA
Rscript tests/calibrate_type2_error_r.R --n_reps 200 --type2_upper_bound 0.25
```
