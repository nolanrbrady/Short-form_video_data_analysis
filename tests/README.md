# Tests / Validation

This folder contains lightweight validation harnesses that aim to verify scientific/analysis integrity against `ANALYSIS_SPEC.md`.

## Pipeline C (R): synthetic end-to-end validation

Runs `analyze_format_content_lmm_channelwise.R` on a synthetic dataset that:
- Uses Homer-style subject IDs (`sub_0001`) and combined IDs (`0001`)
- Uses varying subject ages and checks the age-adjusted omnibus model against a direct reference fit
- Includes both HbO/HbR and multiple channels
- Includes a pruned beta (`NA`) to verify pruned-channel missingness + complete-case logic
- Verifies literal `0` remains a valid observed beta rather than being treated as missing
- Verifies the output `converged` flag is present and TRUE for the clean synthetic fits
- Verifies the tidy output is sorted by ascending `p_unc`
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
- Verifies ROI mean aggregation over available channels with pruned values represented as `NA`
- Verifies literal `0` remains a valid observed beta rather than being treated as missing
- Verifies the output `converged` flag is present and TRUE for the clean synthetic fits
- Verifies the tidy output is sorted by ascending `p_unc`
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

## Mixed-model convergence helper (R): warning classification validation

Runs `r_lmm_convergence_helpers.R` through a focused validation that:
- Verifies canonical `lme4` non-convergence warnings are detected
- Verifies singular-fit and unrelated warnings are not misclassified as non-convergence
- Verifies duplicate warnings are deduplicated and the returned `converged` flag behaves as expected

Command:

```bash
Rscript tests/validate_lmm_convergence_helpers_r.R
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

## Correlation Follow-up ROI Means (R): standalone pooled ROI subset validation

Runs `analyze_correlational_relationships_roi_means.R` and verifies:
- the standalone ROI-focused script exits cleanly on the study inputs
- it writes combined, Pearson-only, and Spearman-only CSVs
- it uses its dedicated ROI-means analysis plan rather than the broader correlation-plan schema
- it contains only pooled behavioral rows and ROI neural rows
- it clears stale ROI-means CSV and figure artifacts before rerun
- it writes figures for uncorrected-significant Pearson or Spearman rows under the default ROI-means plan

Command:

```bash
Rscript tests/validate_correlational_relationships_roi_means_r.R
```

## Behavior Pairwise Correlations (R): standalone behavioral screen validation

Runs `analyze_behavior_pairwise_correlations.R` on a synthetic merged dataset and verifies:
- the standalone behavioral-correlation script exits cleanly
- it writes a full CSV and a significant-only CSV
- it uses pairwise complete cases for each variable pair
- known positive and negative synthetic pairs recover the expected Pearson correlations
- constant-input pairs are skipped with an explicit reason
- binary-vs-continuous pairs remain analyzable under the Pearson-only policy
- it clears stale CSV and figure artifacts before rerun
- uncorrected-significant rows generate figures

Command:

```bash
Rscript tests/validate_behavior_pairwise_correlations_r.R
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
