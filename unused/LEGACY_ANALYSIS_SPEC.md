# Historical analysis specification

Archived on 2026-09-30; this is not the current pooled-mean follow-up.
Its target configuration, missingness rules, and multiplicity/output descriptions
are historical and must not be used to describe current manuscript results.

# ANALYSIS_SPEC — Correlational Follow-up Format Effects (post-hoc)
Last updated: 2026-04-03

## Scope

Goal: Evaluate whether subject-level pooled **long-form** and **short-form** neural means track the matching pooled **long-form** and **short-form** behavioral means in an explicitly exploratory post-hoc analysis, while also retaining supplementary raw behavioral rows that can localize which task cells appear to drive a pattern without reverting the neural side to condition-specific values.

Behavior runs:
- `engagement`
- `retention`
- Supplementary raw-value runs:
  - `sf_education_engagement`
  - `sf_entertainment_engagement`
  - `lf_entertainment_engagement`
  - `lf_education_engagement`
  - `diff_short_form_education`
  - `diff_short_form_entertainment`
  - `diff_long_form_entertainment`
  - `diff_long_form_education`

Neural targets:
- `S04_D02` (`HbO`, `HbR`)
- `R_DLPFC` (`HbR`)
- `L_DLPFC` (`HbO`)
- `M_DMPFC` (`HbO`)
- `L_DMPFC` (`HbO`)

## Inputs

Primary CSV input:
- `data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv`

Required behavior columns:
- `sf_education_engagement`
- `sf_entertainment_engagement`
- `lf_education_engagement`
- `lf_entertainment_engagement`
- `diff_short_form_education`
- `diff_short_form_entertainment`
- `diff_long_form_education`
- `diff_long_form_entertainment`

Required support files:
- `data/config/roi_definition.json`
- `data/config/correlational_analysis_plan.json`
- `data/config/excluded_subjects.json`

## Effect construction

Behavior pooled values:
- `engagement_long = ((lf_education_engagement + lf_entertainment_engagement) / 2)`
- `engagement_short = ((sf_education_engagement + sf_entertainment_engagement) / 2)`
- `retention_long = ((diff_long_form_education + diff_long_form_entertainment) / 2)`
- `retention_short = ((diff_short_form_education + diff_short_form_entertainment) / 2)`
- Supplementary raw-value rows keep one behavior value per task cell.

Neural pooled values:
- First compute condition-level values for each selected target.
- Then compute:
  - `neural_long = mean(LF_Ent, LF_Edu)`
  - `neural_short = mean(SF_Edu, SF_Ent)`
- Supplementary raw-value rows reuse the corresponding pooled neural value:
  - long raw behavioral cells are tested against `neural_long`
  - short raw behavioral cells are tested against `neural_short`

ROI rule:
- ROI condition values are arithmetic means across available non-missing member channels.

## Missingness and data integrity

- `subject_id` is normalized by extracting digits and converting to integer.
- Input must contain exactly one row per normalized `subject_id`; duplicates are a hard error.
- Required behavioral source columns and analyzed beta columns must all exist; missing columns are a hard error.
- Required behavior and beta columns must be numeric/coercible to numeric; non-numeric values are a hard error.
- Channel beta value `0` is treated as a pruned/missing observation for this analysis, consistent with the repo's Homer import note.
- Channel beta value `NA` is treated as missing/pruned.
- A subject contributes a pooled `long` row only when both long-form cells needed for that pooled mean are present.
- A subject contributes a pooled `short` row only when both short-form cells needed for that pooled mean are present.
- A subject contributes a raw-value row only when that behavior cell and its matched pooled neural value are both present.
- After effect construction, association tests use pairwise complete cases only.
- No imputation is allowed.

## Primary statistical outputs

Association methods:
- Pearson correlation only

Reported quantities per tested row:
- `behavior_run`
- `format_pool`
- `behavior_run_type`
- `behavior_condition_code`, `behavior_condition_label`
- `association_method`
- `association_method_tier`
- `analysis_tier`
- `neural_level`, `neural_name`, `chrom`
- `neural_effect`, `neural_condition_label`
- `association_estimate`
- `n_complete`
- `p_unc`
- `p_fdr`
- `family_id`, `family_n_tested`
- Pearson-only fields: `ci95_low`, `ci95_high`, `slope`, `intercept`, `r_squared`

## Multiple testing correction

- Apply **BH-FDR** within each configured `analysis_tier x behavior_run x format_pool x association_method` family.
- Under the default study plan, each family contains the `6` selected neural targets for that behavior run, pooled format, and association method.
- The default plan therefore yields `24` families total:
  - `8` pooled-format families: `2` behavior runs x `2` format pools x `2` association methods
  - `16` raw-value families: `8` raw behavior runs x `2` association methods

## Figures

- Generate figures only for Pearson rows selected by `figures.policy` in the analysis plan.
- Under the default study plan, plots are emitted only for Pearson rows with `p_unc < 0.05`.

## Outputs

- The script clears `data/results/correlational_relationships/` before each run so stale result files and figures cannot persist.
- Results:
  - `data/results/correlational_relationships/pairwise_correlations_r.csv`
  - `data/results/correlational_relationships/pairwise_correlations_r_pearson.csv`
  - `data/results/correlational_relationships/pairwise_correlations_r_spearman.csv`
- Figures:
  - `data/results/correlational_relationships/figures/`

Output ordering:
- all correlation CSVs are sorted in descending `association_estimate` order
- the combined CSV retains both association methods
- the method-specific CSVs contain only the requested association method

## Validation requirements

Validation script:
- `tests/validate_correlational_relationships_r.R`

Must verify:
- exact recovery of known pooled long and pooled short behavioral means
- exact recovery of known pooled long and pooled short channel and ROI neural means
- exact recovery of known raw behavioral task-cell values
- raw behavioral rows use the correct pooled long or pooled short neural values
- zero beta placeholders are treated as missing rather than true activation
- ROI means use available member channels when only a subset is pruned
- all-missing pooled-format ROI requirements invalidate only the affected subject-target row
- Pearson and Spearman rows agree with manual reference calculations
- manual BH agreement within each configured family
- figures are emitted only for Pearson rows allowed by the plan
- fail-hard behavior for malformed ROI/config JSON, duplicate IDs, missing columns, and non-numeric values

---
