# Tests / Validation

This folder contains lightweight validation harnesses that aim to verify scientific/analysis integrity against `ANALYSIS_SPEC.md`.

## Real Result Reproducibility (Python): no-overwrite exported-results check

Runs the primary real-data result-producing scripts into a temporary directory and compares their CSV outputs to the exported files under `data/results`.
This validator intentionally never writes to `data/results`; it only reads exported results for comparison.

It verifies:
- channelwise, ROI, retention, and engagement LMM CSVs reproduce from the documented real inputs
- main exploratory/correlation CSV data reproduce, including the behavior-pairwise all-attempted and BH-FDR CSVs, ignoring only expected temp-vs-exported `plot_file` paths where those columns exist
- channel-behavior screening CSVs and core metadata reproduce
- neural tidy outputs satisfy `estimate / se == t` for Kenward-Roger-consistent 1-df rows

Command:

```bash
python tests/validate_real_result_reproducibility_py.py
```

## Exported Result Table Invariants (Python): no-overwrite CSV audit

Reads exported CSVs under `data/results` and checks internal statistical/reporting invariants without rerunning analyses or writing files.

It verifies:
- neural tidy tables have valid required columns, complete-case `n_obs == 4 * n_subjects`, valid p/q values, CIs containing estimates, BH-FDR by chromophore/effect family, and `estimate / se == t`
- neural posthoc rows match the FDR-significant interaction gate
- retention and engagement tables have valid Holm corrections, CIs containing estimates, and posthoc gating
- correlation tables have valid effect-size ranges, p/q ranges, BH family adjustments, and family sizes

Command:

```bash
python tests/validate_exported_result_table_invariants_py.py
```

## Real Neural Model Diagnostics (R): no-overwrite model audit

Refits the real channelwise and ROI neural primary LMMs in memory, compares the
fixed-effect statistics to the exported tidy result tables, and checks that
residuals/fitted values are finite and non-degenerate. This validator reads
`data/results` but intentionally writes no files.

It verifies:
- channelwise and ROI real-data neural models reproduce exported estimates, SEs, KR df, signed t values, p-values, and KR-df confidence intervals
- exported `n_subjects`, `n_obs`, `converged`, and `singular_fit` flags match the refit models
- real model residuals and fitted values are finite and non-degenerate
- real neural primary models do not emit captured non-convergence warnings

The Shapiro-Wilk residual p-value is printed as a diagnostic, not a hard
correctness gate, because mixed-model residual normality tests can be sensitive
and should be interpreted with plots and sensitivity checks.

Command:

```bash
Rscript tests/validate_real_model_diagnostics_r.R
```

## Recall assessment processing (Python)

Runs synthetic checks for `demographic/process_recall_assessment.py` and verifies:
- configured invalid recall questions are excluded from condition mean denominators
- configured pre/post Qualtrics ID aliases are applied before scoring
- excluded questions remain visible in audit output with `method = excluded_invalid_question`
- invalid-question manifest condition mismatches fail hard

Command:

```bash
python tests/validate_recall_assessment_processing_py.py
```

## Pipeline C (R): synthetic end-to-end validation

Runs `analyze_format_content_lmm_channelwise.R` on a synthetic dataset that:
- Uses Homer-style subject IDs (`sub_0001`) and combined IDs (`0001`)
- Uses varying subject ages and a non-collinear education-years pattern, then checks the age- and education-adjusted omnibus model against a direct reference fit
- Includes both HbO/HbR and multiple channels
- Includes a pruned beta (`NA`) to verify pruned-channel missingness + complete-case logic
- Verifies literal `0` remains a valid observed beta rather than being treated as missing
- Verifies the output `converged` flag is present and TRUE for the clean synthetic fits
- Verifies the tidy output is sorted by ascending `p_unc`
- Verifies coefficient recovery (within tolerance) for known ground-truth fixed effects
- Verifies reported SE values reconstruct the reported signed Kenward-Roger t statistics
- Verifies explicit inferential TP/TN outcomes (known significant interaction channel vs known null channel)
- Verifies post-hoc gating only triggers for interaction-FDR-significant channels
- Verifies duplicate `subject_id` fails hard
- Verifies missing/non-numeric `age` or `education_years` fails hard

Command:

```bash
Rscript tests/validate_pipeline_c_r.R
```

## Pipeline C ROI (R): synthetic end-to-end validation

Runs `analyze_format_content_lmm_roi.R` on a synthetic dataset that:
- Uses Homer-style subject IDs (`sub_0001`) and combined IDs (`0001`)
- Uses varying subject ages and a non-collinear education-years pattern, then checks the age- and education-adjusted omnibus model against a direct reference fit
- Defines ROIs from a strict JSON ROI map (`ROI -> [channels]`)
- Verifies ROI mean aggregation requires at least 2 of 3 available channels, with pruned values represented as `NA`
- Verifies literal `0` remains a valid observed beta rather than being treated as missing
- Verifies the output `converged` flag is present and TRUE for the clean synthetic fits
- Verifies the tidy output is sorted by ascending `p_unc`
- Verifies a participant is excluded from an ROI × chromophore when any condition has fewer than 2 of 3 available channels
- Verifies non-three-channel ROI definitions fail hard
- Verifies coefficient recovery for ROI-level known generating effects
- Verifies reported SE values reconstruct the reported signed Kenward-Roger t statistics
- Verifies explicit inferential TP/TN outcomes (known significant ROI interaction vs known null ROI interaction)
- Verifies BH-FDR correctness across ROIs (per chrom/effect family)
- Verifies interaction-gated post-hoc behavior
- Verifies fail-hard behavior for invalid JSON, missing/non-numeric `age` or `education_years`, missing ROI channels, and overlapping ROI assignments

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
- Direct agreement with an age- and education-adjusted reference omnibus fit
- Explicit TP/TN inferential checks using strong-effect and null-effect synthetic datasets
- Holm-correction correctness across the 3 omnibus effects
- Interaction-gated post-hoc behavior (on/off cases)
- Zero-as-valid retention handling (not treated as missing)
- Complete-case subject drop only for true `NA`
- Fail-hard integrity checks (duplicate IDs, missing required columns, missing/non-numeric `age` or `education_years`, non-numeric retention values)

Command:

```bash
Rscript tests/validate_retention_pipeline_r.R
```

## Retention equal-question-count sensitivity (R)

`tests/validate_retention_sensitivity_r.R` uses the repository's scored audits and
merged input, with corrupted fixtures and all generated outputs confined to a
temporary directory. It verifies all seven common-item omissions, independent
participant/item score calculations, known positive/negative/zero item changes,
unchanged covariates/other conditions, and a fixed complete-case cohort. It checks
audit/manifest consistency, paired answer keys, missing/duplicate coverage,
missing/nonbinary valid scores, stale merged outcomes, independent balanced
factorial contrasts, manual three-effect Holm adjustment, Wald intervals, posthoc
gating and combined reporting in both branches, convergence/rank-deficiency
failures, and baseline agreement with saved primary results. Regression tolerances
are `1e-10` for coefficients, `1e-3` for Satterthwaite df and `1e-7` for other
statistics; score reconstruction uses `1e-12`.

```bash
Rscript tests/validate_retention_sensitivity_r.R
Rscript tests/validate_retention_pipeline_r.R
```

The primary suite supplies additional synthetic true-positive/true-negative
inference and both posthoc-gate branches. The sensitivity script reuses that
primary analysis directly. Run from the repository root after generating current
scored audits, merged outcomes and primary results. A changed question-count
design requires review; tests must not silently adopt a different analysis.

## Engagement Pipeline (R): deterministic + synthetic validation

Runs `analyze_engagement_format_content_lmm.R` on generated engagement datasets and verifies:
- Deterministic analytic ground truth for Length/Content/Interaction contrasts
- End-to-end coefficient recovery for known generating parameters
- Direct agreement with an age- and education-adjusted reference omnibus fit
- Explicit TP/TN inferential checks using strong-effect and null-effect synthetic datasets
- Holm-correction correctness across the 3 omnibus effects
- Interaction-gated post-hoc behavior (on/off cases)
- Zero-as-valid engagement handling (not treated as missing)
- Complete-case subject drop only for true `NA`
- Fail-hard integrity checks (duplicate IDs, missing required columns, missing/non-numeric `age` or `education_years`, non-numeric engagement values)

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
- only existing planned ROI/chromophore comparisons present in the supplied ROI JSON are analyzed, and skipped undefined targets are reported
- updated channel membership changes the corresponding ROI mean without adding unplanned ROIs or chromophores
- hand-calculated condition and format means preserve the existing pruned-channel and missing-condition rules
- missing beta columns for retained targets still fail; a configuration with no eligible targets fails before clearing existing outputs
- output counts and manual BH-adjusted p-values follow the eligible target set within the unchanged families
- it writes combined and Pearson-only CSVs, and does not write a Spearman CSV
- it uses its dedicated ROI-means analysis plan rather than the broader correlation-plan schema
- it contains only pooled behavioral rows and ROI neural rows
- it clears stale ROI-means CSV and figure artifacts before rerun
- it writes figures for uncorrected-significant Pearson rows under the default ROI-means plan

Command:

```bash
Rscript tests/validate_correlational_relationships_roi_means_r.R
```

## Pooled-Mean Correlations (R): standalone pooled-target validation

Runs `analyze_pooled_mean_correlations.R` on a synthetic dataset and verifies:
- Target selection keeps only significant `format`/`content` rows from the tidy channel/ROI LMM outputs
- Stale significant `M_DMPFC` and `M_VMPFC` ROI rows are explicitly excluded even when matching beta columns are available
- Pool gating restricts `format` targets to `short/long` and `content` targets to `education/entertainment`
- Channel pooled neural values recover exact known values
- ROI condition means require at least 2 of 3 good channels in every condition
- Participants failing the 2-of-3 rule in any condition are excluded from every pooled row for that ROI
- A selected ROI with anything other than exactly three configured and available channels fails hard
- Non-significant and interaction-only targets are excluded from the exported target set
- The exported Pearson correlations agree with independently known gated reference cases
- Channel `0` placeholders are treated as missing rather than true beta values
- BH-FDR correctness within each configured family
- Figure emission obeys the default `significant_only` policy
- Duplicate normalized `subject_id` values fail hard

Command:

```bash
Rscript tests/validate_pooled_mean_correlations_r.R
```

## Behavior Pairwise Correlations (R): standalone behavioral screen validation

Runs `analyze_behavior_pairwise_correlations.R` on a synthetic merged dataset and on the real merged publication input, then verifies:
- the standalone behavioral-correlation script exits cleanly
- it writes a full all-attempted CSV and an all-tested BH-FDR CSV
- it emits the lower-triangle matrix as PNG and PDF
- incomplete variable-figure-label configs fail hard and identify the omitted label
- the exact requested variable set is used and both `pd_status` and `recruitment_order_proxy` are absent
- plans that attempt to include `recruitment_order_proxy` fail hard
- it uses pairwise complete cases for each variable pair
- known positive and negative synthetic pairs recover the expected Pearson correlations
- the exploratory `education_years` proxy is included in the declared matrix and treated as numeric Pearson input
- ordinal SFV-use variables are treated as numeric Pearson inputs
- an underpowered synthetic pair is skipped with `n_complete<6`, excluded from the FDR CSV, and left without q-value/significance flags
- global BH-FDR q-values match `p.adjust(..., method = "BH")`
- `significant_fdr` is based on `p_fdr < alpha`
- it clears stale CSV and figure artifacts before rerun
- every real-data pair's `n_complete`, Pearson `r`, raw p-value, Fisher CI, BH-FDR q-value, and FDR flag match an independently recomputed reference table

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

## Behavioral Score Distribution Plotting (R): engagement + retention validation

Runs `plot_behavior_score_distributions.R` on a synthetic final merged dataset and verifies:
- shared subject exclusions are applied before plotting
- engagement and recall/retention condition columns map to the intended four conditions
- content main-effect marginal means are computed within subject as Education and Entertainment averages across Short and Long after complete-case filtering
- behavioral zero values remain valid observations
- missing condition cells trigger domain-specific complete-case exclusion without imputation
- non-numeric score tokens and missing required columns fail hard
- PNG figures plus audit and summary CSVs are created, including retention length-marginal recall figures

Command:

```bash
Rscript tests/validate_behavior_score_distribution_plot_r.R
```

## Demographics Table (Python): shared subject-exclusion validation

Runs `create_demographics_table.py` helper functions on synthetic merged data and verifies:
- Homer-style exclusion IDs remove matching numeric `subject_id` rows after normalization
- ID and beta columns are excluded from the descriptive demographics table
- duplicate exclusion IDs fail after normalization instead of silently removing a participant twice

Command:

```bash
python tests/validate_demographics_table_py.py
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

## Type-I Error Calibration (R): Monte Carlo null simulations across all pipelines

Runs repeated null-effect synthetic datasets through all four inferential scripts:
- `analyze_format_content_lmm_channelwise.R`
- `analyze_format_content_lmm_roi.R`
- `analyze_retention_format_content_lmm.R`
- `analyze_engagement_format_content_lmm.R`

It estimates empirical false-positive rates from adjusted p-values (`p_fdr`) and fails when any
pipeline/effect exceeds a configured upper bound.
The synthetic generators include varying `age` and a non-collinear `education_years` pattern so both required omnibus covariate paths are exercised during calibration.

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
The synthetic generators include varying `age` and a non-collinear `education_years` pattern so both required omnibus covariate paths are exercised during calibration.

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
