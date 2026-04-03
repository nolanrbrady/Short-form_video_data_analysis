# ANALYSIS_SPEC — Homer3 betas + Format×Content (channelwise)
Last updated: 2026-03-30

This document captures the exact specifications agreed **before** implementation of the merge and
channelwise statistical analysis scripts.

If this spec conflicts with any code, the spec should be treated as authoritative until explicitly revised.

---

## Scope

Goal: Evaluate whether **the effect of Format depends on Content (Format×Content interaction)** at the
level of **prefrontal activation**, using **Homer3 subject-level GLM betas**, while properly accounting for the
**within-subjects** design.

Both chromophores are analyzed and reported:
- **HbO**
- **HbR**

Inference is performed **per channel**.

---

## Inputs

Primary CSV inputs (repo-local):
- `data/tabular/generated_data/homer3_glm_betas_wide_auc_outliers_masked.csv`
  - Subject ID column is named `Subject` (e.g., `sub_0001`)
  - Beta columns are wide and follow the pattern:
    - `S##_D##_Cond##_HbO`
    - `S##_D##_Cond##_HbR`
    - Example: `S01_D01_Cond01_HbO`
- `data/tabular/generated_data/homer3_glm_betas_wide_auc.csv`
  - Raw post-collapse AUC table retained for provenance and pre-mask validation
- `data/tabular/generated_data/combined_sfv_data.csv`
  - Subject ID column is named `subject_id`
  - `subject_id` may be **zero-padded** in some sources (e.g., `0017` vs `17`)

Outputs are written under:
- `data/results/`

Toy inputs (for quick validation only):
- `data/toy/toy_homer3_glm_betas_wide.csv`
- `data/toy/toy_combined_sfv_data.csv`

---

## Subject ID normalization + merge specifications

Join type:
- **INNER JOIN** between Homer3 and combined tabular datasets.

ID handling:
- Both `combined_sfv_data.csv:subject_id` and `homer3_glm_betas_wide_auc_outliers_masked.csv:Subject` are treated as **numeric IDs**
  (even if stored as strings).
- IDs are normalized by **extracting digits** and converting to integer (handles `0017` vs `17`, and `sub_0017`).

Expected row structure:
- `homer3_glm_betas_wide_auc_outliers_masked.csv` contains **one row per subject**.
- Result of merge is **one row per subject** containing:
  - all relevant combined tabular columns (demographics/behavior)
  - all Homer beta columns

No imputation during merge:
- The merge step must not silently impute missing values.

Implementation scripts:
- Python: `merge_homer3_betas_with_combined_data.py`
- R: `merge_homer3_betas_with_combined_data.R`

---

## Missingness / pruned-channel policy (critical)

Per repo policy, Homer betas can include values that stand in for **pruned channels**:
- in the raw FIR export (`homer3_glm_betas_wide_fir_pca.csv`), `0` values may indicate a pruned channel
- in the raw FIR export (`homer3_glm_betas_wide_fir_pca.csv`), `NaN` values may indicate a pruned channel

Required handling:
- In the derived single-beta AUC table (`homer3_glm_betas_wide_auc.csv`), carry pruned channels forward as **`NaN`**.
- In the between-subject outlier-masked AUC table (`homer3_glm_betas_wide_auc_outliers_masked.csv`), carry both pruned channels and censored outlier values as **`NaN`**.
- Do **not** treat these values as true zero activation.
- Do **not** silently impute these values.

Downstream modeling must handle this explicitly as missingness.

Upstream derivation note:
- `homer3_glm_betas_wide_auc.csv` is produced from the raw FIR export `data/tabular/homer3_glm_betas_wide_fir_pca.csv`
  by reconstructing the latent HRF from the Gaussian basis weights and computing a baseline-corrected
  task-window AUC. The merged/statistical analysis scripts consume the derived single-beta table, not the raw FIR table.
- `homer3_glm_betas_wide_auc_outliers_masked.csv` is produced from `homer3_glm_betas_wide_auc.csv` by screening each exact
  channel x condition x chromophore column across subjects and masking values outside `mean +/- 3 SD`.
- Because sample mean/SD screening cannot detect a `3 SD` outlier when fewer than 11 observed subjects are available,
  undersized columns are reported as skipped rather than silently treated as screened.
- Production settings are:
  - `idxBasis = 1`
  - reconstruction support `[-10, 130]`
  - basis spacing `0.5 s`
  - Gaussian sigma `0.5 s`
  - baseline window `[-10, 0]`
  - AUC window `[0, 120]`

---

## Condition mapping (experiment design)

Condition codes are embedded in beta column names as `Cond01`, `Cond02`, `Cond03`, `Cond04`.

Mapping to experimental conditions:
- `Cond01` → **Short-Form Education**
- `Cond02` → **Short-Form Entertainment**
- `Cond03` → **Long-Form Entertainment**
- `Cond04` → **Long-Form Education**

Factor definitions:
- **Format**: Short vs Long
- **Content**: Education vs Entertainment

---

## Data reshaping for analysis (per channel × chromophore)

Each beta column represents:
- a single beta value for a given (subject × channel × chromophore × condition).

Analysis requires conversion from wide → long with fields:
- `subject_id`
- `age`
- `channel` (e.g., `S01_D01`)
- `chrom` (`HbO` or `HbR`)
- `condition` (one of: `SF_Edu`, `SF_Ent`, `LF_Ent`, `LF_Edu`)
- `beta`
- derived predictors:
  - `format_c`
  - `content_c`

Complete-case rule (within channel/chromophore):
- For a given (channel, chromophore), **only subjects with all 4 conditions present (non-missing beta)** are included.
- `age` is a required subject-level omnibus covariate: it must exist, be numeric, and be complete after subject exclusions or the script fails hard.

---

## Primary statistical model (main effects + interaction)

Model form:
- One model per **(channel × chromophore)**.
- Linear mixed model with random intercept for subject:
  - `beta ~ format_c * content_c + age + (1 | subject_id)`

Coding (required):
- Use numeric sum/effect coding with ±0.5:
  - `format_c = -0.5` for Short, `+0.5` for Long
  - `content_c = -0.5` for Entertainment, `+0.5` for Education

Reported quantities (per channel × chromophore × effect):
- Fixed-effect estimate (Format, Content, and **Format×Content interaction**)
- 95% CI
- p-value (uncorrected and FDR-corrected)

Significance threshold:
- Gate significance on **FDR-corrected p < 0.05**.

R implementation notes:
- LMM via `lme4::lmer`, with fixed-effect p-values/df from Kenward-Roger Type-III tests via `lmerTest` + `pbkrtest`.
- The current omnibus covariate adjustment includes `age` only; `sfv_daily_duration` is deferred until its missingness is resolved upstream.
- For numerical conditioning, the implemented R script may fit the neural response after multiplying beta by one fixed global constant (`1e6`), but reported estimates/CIs are back-transformed into the original beta units before output.
- Implemented outputs also include a boolean `converged` flag based on captured mixed-model convergence warnings so any numerically suspect fits remain auditable in the result tables.
- Post-hoc via `emmeans`, using the existing condition-only follow-up model.

Python implementation notes:
- LMM intended via `statsmodels` `MixedLM` (if installed).
- If `statsmodels` is not installed, Python script supports a parse/reshape validation via `--dry-run`.

References (see `CITATIONS.md`):
- Mixed models: Laird & Ware (1982); Bates et al. (2015).
- Kenward-Roger inference: Kenward & Roger (1997); Halekoh & Højsgaard (2014); Kuznetsova et al. (2017).

---

## Multiple testing correction (FDR)

Correction method:
- **Benjamini–Hochberg (BH) FDR**.

Correction families (agreed):
- HbO and HbR are treated as **separate hypothesis families**.
- Within each chromophore:
  - Apply BH-FDR **separately per effect** across channels:
    - Format main effect
    - Content main effect
    - Format×Content interaction

Note to revisit later (documented in `README.md`):
- We explicitly deferred a decision on whether a broader family (e.g., across effects and/or chromophores)
  is preferable for the final manuscript reporting plan.

Reference:
- Benjamini & Hochberg (1995) — see `CITATIONS.md`.

---

## Post-hoc analysis (only if interaction significant)

Gate:
- Perform post-hoc tests **only** for (channel × chromophore) where the **interaction** is
  FDR-significant (q < 0.05).

Post-hoc set (Option B):
- All pairwise comparisons among the 4 conditions (6 contrasts):
  - `SF_Edu` vs `SF_Ent`
  - `SF_Edu` vs `LF_Ent`
  - `SF_Edu` vs `LF_Edu`
  - `SF_Ent` vs `LF_Ent`
  - `SF_Ent` vs `LF_Edu`
  - `LF_Ent` vs `LF_Edu`

Multiplicity handling:
- **No multiple-test correction** is applied in the post-hoc analysis (per instruction).

Reference:
- Estimated marginal means / contrasts: Lenth (2016); Searle et al. (1980) — see `CITATIONS.md`.
- Paired t-test framework (Python post-hoc implementation): Student (1908) — see `CITATIONS.md`.

---

## Output requirements

Merge output:
- A merged CSV with one row per subject and columns from both inputs.
- Destination under `data/results/` (filename chosen at run time via script argument defaults).

Analysis outputs:
- Main effects table (CSV) including, at minimum:
  - `channel`, `chrom`, `effect`
  - `estimate`, `ci95_low`, `ci95_high`
  - `p_unc`, `p_fdr`
  - sample sizes (`n_subjects`, `n_obs`)
- Post-hoc table (CSV) containing the 6 contrasts **only for gated channel/chrom pairs**.

Console reporting:
- Scripts should print key counts (subjects, channels) and where outputs were written.

---

## Explicit non-requirements / exclusions

- Do not rely on or extend the existing/older statistical analysis scripts in `fnirs_analysis/`
  for this pipeline.
- Do not impute pruned channels.
- Do not “collapse” channels into an ROI summary (inference is per-channel).

---

# ANALYSIS_SPEC — Homer3 betas + Format×Content (ROI-wise)
Last updated: 2026-03-30

## Scope

Goal: Evaluate whether **the effect of Format depends on Content (Format×Content interaction)** at the
level of **ROI-wise prefrontal activation**, using the same Homer3 subject-level GLM betas and within-subject
design as the channelwise pipeline.

Both chromophores are analyzed and reported:
- **HbO**
- **HbR**

Inference is performed **per ROI**.

---

## Inputs

Primary CSV input:
- `data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv`
  - Must include one row per subject (`subject_id`), a numeric `age` column, and Homer beta columns matching:
    - `S##_D##_Cond##_HbO`
    - `S##_D##_Cond##_HbR`

ROI definition input:
- `data/config/roi_definition.json`
  - Must be strict JSON.
  - Top-level object maps ROI names to arrays of channel IDs:
    - Example: `"VMPFC": ["S01_D01", "S01_D02"]`

---

## ROI definition integrity rules

- ROI JSON must parse without coercion/fallback.
- ROI names must be non-empty.
- Each ROI must contain at least one channel.
- Channel IDs must match Homer naming (`S##_D##` after normalization).
- A channel cannot be assigned to multiple ROIs.
- ROI channels not present in the input beta columns are a hard error.

---

## Missingness / pruned-channel policy (critical)

Per repo policy:
- In the derived FIR-to-AUC beta table, pruned channels are encoded as `NaN`.
- Do not impute.

ROI summary construction:
- For each `subject × ROI × chrom × condition`, ROI beta is the arithmetic mean over
  available (non-missing) channels in that ROI.
- If all channels are missing for that cell, ROI beta is missing.

Complete-case inclusion rule (within ROI/chrom):
- Keep only subjects with non-missing ROI beta in all 4 conditions.
- `age` is a required subject-level omnibus covariate: it must exist, be numeric, and be complete after subject exclusions or the script fails hard.

---

## Condition mapping and model

Condition mapping from beta columns:
- `Cond01` → `SF_Edu`
- `Cond02` → `SF_Ent`
- `Cond03` → `LF_Ent`
- `Cond04` → `LF_Edu`

Effect coding:
- `format_c = -0.5` (Short), `+0.5` (Long)
- `content_c = -0.5` (Entertainment), `+0.5` (Education)

Primary model (per ROI × chrom):
- `beta ~ format_c * content_c + age + (1 | subject_id)`

Inference and post-hoc:
- Main effects reported for Format, Content, and Interaction.
- The current omnibus covariate adjustment includes `age` only; `sfv_daily_duration` is deferred until its missingness is resolved upstream.
- For numerical conditioning, the implemented R script may fit the neural response after multiplying beta by one fixed global constant (`1e6`), but reported estimates/CIs are back-transformed into the original beta units before output.
- Implemented outputs also include a boolean `converged` flag based on captured mixed-model convergence warnings so any numerically suspect fits remain auditable in the result tables.
- Post-hoc pairwise condition contrasts (6 total) run only when ROI/chrom interaction
  is FDR-significant.
- Post-hoc p-values are uncorrected (`adjust = "none"`).

---

## Multiple testing correction

- Use BH-FDR.
- Families are defined separately by chromophore and effect, across ROIs:
  - HbO / Format across ROIs
  - HbO / Content across ROIs
  - HbO / Interaction across ROIs
  - HbR / Format across ROIs
  - HbR / Content across ROIs
  - HbR / Interaction across ROIs

---

## Output requirements

Main effects (wide):
- CSV with one row per ROI × chrom and per-effect estimate/statistics fields.

Main effects (tidy/spec):
- CSV with one row per ROI × chrom × effect including:
  - `roi`, `chrom`, `effect`
  - `estimate`, `ci95_low`, `ci95_high`
  - `p_unc`, `p_fdr`
  - `n_subjects`, `n_obs`
  - `singular_fit`

Post-hoc:
- CSV with pairwise contrasts only for gated ROI/chrom pairs.

---

# ANALYSIS_SPEC — Retention Length×Content (subject-level LMM)
Last updated: 2026-03-16

## Scope

Goal: Evaluate whether **retention improvement** (post - pre) differs by:
- **Length** (Short vs Long),
- **Content** (Education vs Entertainment),
- and their **Length×Content interaction**,
while accounting for repeated measures (4 within-subject conditions).

## Inputs

Primary CSV input:
- `data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv`

Required columns:
- `subject_id`
- `age`
- `diff_short_form_education`
- `diff_short_form_entertainment`
- `diff_long_form_education`
- `diff_long_form_entertainment`

## Subject ID + data integrity

- `subject_id` is normalized by extracting digits and converting to integer.
- Input must contain exactly one row per normalized `subject_id`; duplicates are a hard error.
- Required retention columns and `age` must all exist; missing columns are a hard error.
- Retention columns and `age` must be numeric/coercible to numeric; non-numeric values are a hard error.
- `age` must be complete after subject exclusions; any remaining missing value is a hard error.

## Condition mapping and coding

Condition mapping:
- `diff_short_form_education` -> `SF_Edu`
- `diff_short_form_entertainment` -> `SF_Ent`
- `diff_long_form_entertainment` -> `LF_Ent`
- `diff_long_form_education` -> `LF_Edu`

Effect coding:
- `length_c = -0.5` (Short), `+0.5` (Long)
- `content_c = -0.5` (Entertainment), `+0.5` (Education)

## Missingness policy

- Complete-case by subject across the 4 retention conditions:
  include only subjects with all 4 non-missing retention values.
- Retention value `0` is treated as a valid observed value (not missing).
- No imputation is allowed.

## Primary statistical model

Model:
- One subject-level LMM:
  - `retention_diff ~ length_c * content_c + age + (1 | subject_id)`

Reported quantities (for `length_c`, `content_c`, `length_c:content_c`):
- estimate, SE, df, t, uncorrected p, Holm-adjusted p, Wald 95% CI
- sample size fields: `n_subjects`, `n_obs`
- singular-fit flag

R implementation notes:
- Fit with `lmerTest::lmer` (REML).
- p-values from `lmerTest` (Satterthwaite df approximation).
- The current omnibus covariate adjustment includes `age` only; `sfv_daily_duration` is deferred until its missingness is resolved upstream.

## Multiple testing correction

- Apply **Holm correction** across the **three omnibus effects**:
  - Length
  - Content
  - Length×Content interaction

Significance gate:
- Use adjusted `p < alpha` for inferential gating.

## Post-hoc analysis

Gate:
- Run post-hoc contrasts only when interaction adjusted p-value is significant (`p_adj_interaction < alpha`).

Method:
- Fit condition model: `retention_diff ~ condition + (1 | subject_id)`
- Compute all 6 pairwise condition contrasts via `emmeans`, uncorrected (`adjust = "none"`):
  - `SF_Edu` vs `SF_Ent`
  - `SF_Edu` vs `LF_Ent`
  - `SF_Edu` vs `LF_Edu`
  - `SF_Ent` vs `LF_Ent`
  - `SF_Ent` vs `LF_Edu`
  - `LF_Ent` vs `LF_Edu`

Output includes:
- `condition_a`, `condition_b`, `mean_diff` (`condition_a - condition_b`), `se`, `df`, `t`, `p_unc`, `stat_type`

## Outputs

- Main effects:
  - `data/results/retention_format_content_lmm_main_effects_r.csv`
- Post-hoc pairwise:
  - `data/results/retention_format_content_lmm_posthoc_pairwise_r.csv`

## Validation requirements

Validation script:
- `tests/validate_retention_pipeline_r.R`

Must verify:
- deterministic analytic recovery of known Length/Content/Interaction effects
- end-to-end coefficient recovery in synthetic data with known generating parameters
- direct agreement with an age-adjusted reference omnibus fit
- manual Holm agreement with output adjusted p-values
- post-hoc gating behavior (on/off)
- complete-case behavior for `NA`
- retention `0` handling as valid (not missing)
- fail-hard behavior for duplicates, missing columns, missing `age`, and non-numeric values

---

# ANALYSIS_SPEC — Engagement Length×Content (subject-level LMM)
Last updated: 2026-03-16

## Scope

Goal: Evaluate whether **engagement ratings** differ by:
- **Length** (Short vs Long),
- **Content** (Education vs Entertainment),
- and their **Length×Content interaction**,
while accounting for repeated measures (4 within-subject conditions).

## Inputs

Primary CSV input:
- `data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv`

Required columns:
- `subject_id`
- `age`
- `sf_education_engagement`
- `sf_entertainment_engagement`
- `lf_education_engagement`
- `lf_entertainment_engagement`

## Subject ID + data integrity

- `subject_id` is normalized by extracting digits and converting to integer.
- Input must contain exactly one row per normalized `subject_id`; duplicates are a hard error.
- Required engagement columns and `age` must all exist; missing columns are a hard error.
- Engagement columns and `age` must be numeric/coercible to numeric; non-numeric values are a hard error.
- `age` must be complete after subject exclusions; any remaining missing value is a hard error.

## Condition mapping and coding

Condition mapping:
- `sf_education_engagement` -> `SF_Edu`
- `sf_entertainment_engagement` -> `SF_Ent`
- `lf_entertainment_engagement` -> `LF_Ent`
- `lf_education_engagement` -> `LF_Edu`

Effect coding:
- `length_c = -0.5` (Short), `+0.5` (Long)
- `content_c = -0.5` (Entertainment), `+0.5` (Education)

## Missingness policy

- Complete-case by subject across the 4 engagement conditions:
  include only subjects with all 4 non-missing engagement values.
- Engagement value `0` is treated as a valid observed value (not missing).
- No imputation is allowed.

## Primary statistical model

Model:
- One subject-level LMM:
  - `engagement ~ length_c * content_c + age + (1 | subject_id)`

Reported quantities (for `length_c`, `content_c`, `length_c:content_c`):
- estimate, SE, df, t, uncorrected p, Holm-adjusted p, Wald 95% CI
- sample size fields: `n_subjects`, `n_obs`
- singular-fit flag

R implementation notes:
- Fit with `lmerTest::lmer` (REML).
- p-values from `lmerTest` (Satterthwaite df approximation).
- The current omnibus covariate adjustment includes `age` only; `sfv_daily_duration` is deferred until its missingness is resolved upstream.

## Multiple testing correction

- Apply **Holm correction** across the **three omnibus effects**:
  - Length
  - Content
  - Length×Content interaction

Significance gate:
- Use adjusted `p < alpha` for inferential gating.

## Post-hoc analysis

Gate:
- Run post-hoc contrasts only when interaction adjusted p-value is significant (`p_adj_interaction < alpha`).

Method:
- Fit condition model: `engagement ~ condition + (1 | subject_id)`
- Compute all 6 pairwise condition contrasts via `emmeans`, uncorrected (`adjust = "none"`):
  - `SF_Edu` vs `SF_Ent`
  - `SF_Edu` vs `LF_Ent`
  - `SF_Edu` vs `LF_Edu`
  - `SF_Ent` vs `LF_Ent`
  - `SF_Ent` vs `LF_Edu`
  - `LF_Ent` vs `LF_Edu`

Output includes:
- `condition_a`, `condition_b`, `mean_diff` (`condition_a - condition_b`), `se`, `df`, `t`, `p_unc`, `stat_type`

## Outputs

- Main effects:
  - `data/results/engagement_format_content_lmm_main_effects_r.csv`
- Post-hoc pairwise:
  - `data/results/engagement_format_content_lmm_posthoc_pairwise_r.csv`

## Validation requirements

Validation script:
- `tests/validate_engagement_pipeline_r.R`

Must verify:
- deterministic analytic recovery of known Length/Content/Interaction effects
- end-to-end coefficient recovery in synthetic data with known generating parameters
- direct agreement with an age-adjusted reference omnibus fit
- manual Holm agreement with output adjusted p-values
- post-hoc gating behavior (on/off)
- complete-case behavior for `NA`
- engagement `0` handling as valid (not missing)
- fail-hard behavior for duplicates, missing columns, missing `age`, and non-numeric values

---

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
- `data/tabular/homer3_betas_plus_combined_sfv_data_inner_join.csv`

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
- Primary: Pearson correlation
- Sensitivity: Spearman correlation

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
