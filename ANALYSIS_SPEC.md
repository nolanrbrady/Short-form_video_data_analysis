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
- `education_years`
- `channel` (e.g., `S01_D01`)
- `chrom` (`HbO` or `HbR`)
- `condition` (one of: `SF_Edu`, `SF_Ent`, `LF_Ent`, `LF_Edu`)
- `beta`
- derived predictors:
  - `format_c`
  - `content_c`

Complete-case rule (within channel/chromophore):
- For a given (channel, chromophore), **only subjects with all 4 conditions present (non-missing beta)** are included.
- `age` and `education_years` are required subject-level omnibus covariates: both must exist, be numeric, and be complete after subject exclusions or the script fails hard.

---

## Primary statistical model (main effects + interaction)

Model form:
- One model per **(channel × chromophore)**.
- Linear mixed model with random intercept for subject:
  - `beta ~ format_c * content_c + age + education_years + (1 | subject_id)`

Coding (required):
- Use numeric sum/effect coding with ±0.5:
  - `format_c = -0.5` for Short, `+0.5` for Long
  - `content_c = -0.5` for Entertainment, `+0.5` for Education

Reported quantities (per channel × chromophore × effect):
- Fixed-effect estimate (Format, Content, and **Format×Content interaction**)
- Kenward-Roger-consistent standard error, Kenward-Roger denominator df, and signed t-statistic
- 95% CI from the Kenward-Roger-consistent standard error and denominator df
- p-value (uncorrected and FDR-corrected)

Significance threshold:
- Gate significance on **FDR-corrected p < 0.05**.

R implementation notes:
- LMM via `lme4::lmer`, with fixed-effect p-values/df from Kenward-Roger Type-III tests via `lmerTest` + `pbkrtest`.
- For the 1-df omnibus terms, the R script derives the reported SE from the same Kenward-Roger F statistic used for signed t, so `estimate / se` reconstructs `t` in publication tables.
- The current omnibus covariate adjustment includes `age` and the study-codebook `education_years` proxy; `sfv_daily_duration` remains deferred until its missingness is resolved upstream.
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
  - Must include one row per subject (`subject_id`), numeric `age` and `education_years` columns, and Homer beta columns matching:
    - `S##_D##_Cond##_HbO`
    - `S##_D##_Cond##_HbR`

ROI definition input:
- `data/config/roi_definition.json`
  - Must be strict JSON.
  - Top-level object maps ROI names to arrays of channel IDs:
    - Example: `"VMPFC": ["S01_D01", "S01_D02", "S01_D03"]`

---

## ROI definition integrity rules

- ROI JSON must parse without coercion/fallback.
- ROI names must be non-empty.
- Each inferential ROI must contain exactly three channels.
- Channel IDs must match Homer naming (`S##_D##` after normalization).
- A channel cannot be assigned to multiple ROIs.
- ROI channels not present in the input beta columns are a hard error.

---

## Missingness / pruned-channel policy (critical)

Per repo policy:
- In the derived FIR-to-AUC beta table, pruned channels are encoded as `NaN`.
- Do not impute.

ROI summary construction:
- Every inferential ROI must contain exactly 3 channels; non-three-channel definitions fail hard.
- For each `subject × ROI × chrom × condition`, ROI beta is the arithmetic mean over
  available channels only when at least 2 of 3 channels are non-missing.
- If fewer than 2 channels are available for that cell, ROI beta is missing.
- The exact 2-of-3 cutoff is a predeclared study-specific conservative rule informed by published fNIRS good-channel inclusion precedents; it is not presented as a universal cutoff.

Complete-case inclusion rule (within ROI/chrom):
- Keep only subjects satisfying the 2-of-3 channel rule in all 4 conditions.
- `age` and `education_years` are required subject-level omnibus covariates: both must exist, be numeric, and be complete after subject exclusions or the script fails hard.

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
- `beta ~ format_c * content_c + age + education_years + (1 | subject_id)`

Inference and post-hoc:
- Main effects reported for Format, Content, and Interaction.
- The ROI tidy main-effects output reports estimate, Kenward-Roger-consistent
  standard error, Kenward-Roger denominator df, signed t-statistic, 95% CI,
  uncorrected p-value, and BH-FDR q-value for each ROI × chromophore × effect row.
  For the 1-df omnibus terms, `estimate / se` reconstructs the reported signed t.
- The current omnibus covariate adjustment includes `age` and the study-codebook `education_years` proxy; `sfv_daily_duration` remains deferred until its missingness is resolved upstream.
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
- `education_years`
- `diff_short_form_education`
- `diff_short_form_entertainment`
- `diff_long_form_education`
- `diff_long_form_entertainment`

## Subject ID + data integrity

- `subject_id` is normalized by extracting digits and converting to integer.
- Input must contain exactly one row per normalized `subject_id`; duplicates are a hard error.
- Required retention columns, `age`, and `education_years` must all exist; missing columns are a hard error.
- Retention inputs must come from the recall preprocessing script using the invalid-question manifest (`data/config/recall_invalid_questions.json`) so `Q5`, `Q6`, `Q7`, `Q8`, `Q10`, `Q26`, and `Q28` are excluded from both pre-task and post-task scoring denominators, including pre-task aliases `Q35` for invalid `Q6` and `Q36` for invalid `Q7`.
- Retention preprocessing must also apply `data/config/recall_question_aliases.json`; currently this maps pre-task `Q39` to canonical post-task/key item `Q22`.
- Retention columns, `age`, and `education_years` must be numeric/coercible to numeric; non-numeric values are a hard error.
- `age` and `education_years` must be complete after subject exclusions; any remaining missing value is a hard error.

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
  - `retention_diff ~ length_c * content_c + age + education_years + (1 | subject_id)`

Reported quantities (for `length_c`, `content_c`, `length_c:content_c`):
- estimate, SE, df, t, uncorrected p, Holm-adjusted p, Wald 95% CI
- sample size fields: `n_subjects`, `n_obs`
- singular-fit flag

R implementation notes:
- Fit with `lmerTest::lmer` (REML).
- p-values from `lmerTest` (Satterthwaite df approximation).
- The current omnibus covariate adjustment includes `age` and the study-codebook `education_years` proxy; `sfv_daily_duration` remains deferred until its missingness is resolved upstream.

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
- direct agreement with an age- and education-adjusted reference omnibus fit
- manual Holm agreement with output adjusted p-values
- post-hoc gating behavior (on/off)
- complete-case behavior for `NA`
- retention `0` handling as valid (not missing)
- fail-hard behavior for duplicates, missing columns, missing `age` or `education_years`, and non-numeric values

---

# ANALYSIS_SPEC — Retention equal-question-count sensitivity

## Scope and fixed design

`analyze_retention_sensitivity.R` is a supplementary robustness check of the
retention analysis above, not a replacement primary analysis or null permutation
test. Baseline valid question counts are 6, 7, 6, 6 in Short Education, Short
Entertainment, Long Education, Long Entertainment order. Enumerate all seven
six-of-seven subsets of Short Entertainment (`Q3`, `Q4`, `Q11`, `Q12`, `Q19`,
`Q20`, `Q27` in the current audits). Each omitted item is removed for every
participant and both phases. Compute post mean minus pre mean over six retained
items. Preserve the other three outcomes, covariates, exclusions and primary
complete-case cohort. No RNG, independent participant-specific omissions, or
selection of a favorable subset is used.

## Integrity and inference

- Validate pre/post scored audits against `recall_invalid_questions.json`, with
  consistent item keys/conditions and participant coverage; each phase must have
  eight original items per condition and valid counts 6/7/6/6. The resampled
  condition must have matching pre/post item IDs and answer keys.
- Excluded scores must be missing; all valid scores must be binary and present.
  Upstream valid blank responses already scored zero retain that score. No
  regrading, imputation, or denominator changes on missing scores are performed.
- Reconstruct all four baseline outcomes for participants remaining after shared
  exclusions; require agreement with every observed merged outcome at absolute
  tolerance `1e-12`. Preserve the primary complete-case rule when a merged outcome
  is missing, and hold the resulting cohort fixed for all eight models.
- Invoke the primary R script's analysis function for each input: same REML model,
  coding, Satterthwaite tests, unadjusted Wald CIs, effect sizes, three-effect Holm
  family, and interaction-gated uncorrected six-pair posthoc comparisons.
- Do not pool the 24 p-values into a new family. Summaries are descriptive
  sensitivity ranges, not an additional inferential test. No subset is promoted
  to the primary result. Stop for nonconvergence, nonfinite inference, or dropped
  fixed-effect columns; retain and flag singular fits.

## Reporting and limits

Write all scenario inputs and results to a separate new/empty output directory,
with a 24-row combined effects table, three-row summary, omission plan,
participant scores, diagnostics, gated contrasts, input/code hashes and session
versions. Name the combined adjusted-p column `p_holm`; the original-format
scenario tables retain the primary script's historical `p_fdr` name.

Report coefficient and adjusted-p ranges, plus counts of direction/significance
agreement. Overlapping subsets are not independent replications. Their mean
score equals the baseline arithmetically; this is a validation identity, not
evidence of unbiased inference. Balancing counts does not address unequal item
difficulty, omitted-video confounding, or stimulus generalization (Judd et al.,
2012). No exchangeability-based null permutations occur (Winkler et al., 2014).
See the retention-sensitivity section in `CITATIONS.md`.

Validate using `tests/validate_retention_sensitivity_r.R` and the existing primary
retention validation suite. Cover independent paired-item scoring and balanced
factorial contrasts, exhaustive enumeration, stable cohort/covariates, valid zeros,
manual Holm agreement, saved baseline regression, fail-fast corrupt/stale input
cases, and protection against overwriting a prior sensitivity run.

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
- `education_years`
- `sf_education_engagement`
- `sf_entertainment_engagement`
- `lf_education_engagement`
- `lf_entertainment_engagement`

## Subject ID + data integrity

- `subject_id` is normalized by extracting digits and converting to integer.
- Input must contain exactly one row per normalized `subject_id`; duplicates are a hard error.
- Required engagement columns, `age`, and `education_years` must all exist; missing columns are a hard error.
- Engagement columns, `age`, and `education_years` must be numeric/coercible to numeric; non-numeric values are a hard error.
- `age` and `education_years` must be complete after subject exclusions; any remaining missing value is a hard error.

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
  - `engagement ~ length_c * content_c + age + education_years + (1 | subject_id)`

Reported quantities (for `length_c`, `content_c`, `length_c:content_c`):
- estimate, SE, df, t, uncorrected p, Holm-adjusted p, Wald 95% CI
- sample size fields: `n_subjects`, `n_obs`
- singular-fit flag

R implementation notes:
- Fit with `lmerTest::lmer` (REML).
- p-values from `lmerTest` (Satterthwaite df approximation).
- The current omnibus covariate adjustment includes `age` and the study-codebook `education_years` proxy; `sfv_daily_duration` remains deferred until its missingness is resolved upstream.

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
- direct agreement with an age- and education-adjusted reference omnibus fit
- manual Holm agreement with output adjusted p-values
- post-hoc gating behavior (on/off)
- complete-case behavior for `NA`
- engagement `0` handling as valid (not missing)
- fail-hard behavior for duplicates, missing columns, missing `age` or `education_years`, and non-numeric values

---

# ANALYSIS_SPEC — Pooled-Mean Neural-Behavior Correlations
Last updated: 2026-09-30

## Scope

Goal: Evaluate whether pooled neural activation for `short`, `long`, `education`, and `entertainment` tracks pooled engagement and retention for the same pools across subjects in an explicitly exploratory follow-up analysis.

This is the selected pooled neural-behavior follow-up for the manuscript. The statistical rules below are unchanged; this selection does not alter the tested targets or multiplicity correction.

Selection and gating rule:
- Read `data/results/format_content_lmm_main_effects_tidy_r.csv` and keep rows with `p_fdr < 0.05`.
- Read `data/results/format_content_lmm_roi_main_effects_tidy_r.csv` and keep rows with `p_fdr < 0.05`.
- Keep only `format` and `content` effects for this pooled-main-effect follow-up.
- Gate eligible pools by selected main effect:
  - `format -> short, long`
  - `content -> education, entertainment`
- Exclude pure `interaction` hits from this workflow.
- Treat the resulting target set as exploratory because the same dataset is used for target selection and pooled follow-up testing.

Modeling target:
- For each selected `analysis_level x unit_id x chrom x selected_effect x behavior_domain x pool_name`, compute a pooled neural mean and a pooled behavioral mean, then fit a Pearson correlation across subjects.

Behavior domains:
- `engagement`
- `retention`

## Inputs

Primary CSV input:
- `data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv`

Required support files:
- `data/config/roi_definition.json`
- `data/results/format_content_lmm_main_effects_tidy_r.csv`
- `data/results/format_content_lmm_roi_main_effects_tidy_r.csv`
- `data/config/excluded_subjects.json`

Required behavior source columns:
- `sf_education_engagement`
- `sf_entertainment_engagement`
- `lf_entertainment_engagement`
- `lf_education_engagement`
- `diff_short_form_education`
- `diff_short_form_entertainment`
- `diff_long_form_entertainment`
- `diff_long_form_education`

## Pooled construction

Condition order:
- `SF_Edu`
- `SF_Ent`
- `LF_Ent`
- `LF_Edu`

Pools:
- `short = mean(SF_Edu, SF_Ent)` when both cells are present
- `long = mean(LF_Edu, LF_Ent)` when both cells are present
- `education = mean(SF_Edu, LF_Edu)` when both cells are present
- `entertainment = mean(SF_Ent, LF_Ent)` when both cells are present

ROI rule:
- Each selected ROI must contain exactly three configured channels, and all three must be represented in the beta input.
- A condition-level ROI value is the arithmetic mean when at least 2 of 3 channels are non-missing.
- A participant contributes to a selected ROI only when that 2-of-3 rule is satisfied in all four conditions; otherwise all pooled rows for that participant and ROI are excluded.

Target-scope rule:
- `M_DMPFC` and `M_VMPFC` are excluded from pooled-mean target selection even if a stale ROI LMM result table still contains significant rows for them.

## Missingness and data integrity

- `subject_id` is normalized by extracting digits and converting to integer.
- Input must contain exactly one row per normalized `subject_id`; duplicates are a hard error.
- Required behavior columns and analyzed beta columns must all exist; missing columns are a hard error.
- Required behavior and beta columns must be numeric/coercible to numeric; non-numeric values are a hard error.
- Channel beta value `0` is treated as a pruned/missing observation for this workflow.
- Channel beta value `NA` is treated as missing/pruned.
- A subject contributes to a pooled target/domain row only when both constituent neural cells and both constituent behavioral cells for that pool are present.
- No imputation is allowed.

## Statistical outputs

Reported quantities:
- `analysis_level`, `unit_id`, `chrom`
- `target_origin`, `roi_member_count`
- `selected_effect`
- `selection_estimate`, `selection_p_fdr`, `selection_source`
- `behavior_domain`
- `pool_name`
- `analysis_status`, `skip_reason`
- `n_complete`
- `association_estimate` (Pearson `r`)
- `r_squared`
- `p_unc`
- `ci95_low`, `ci95_high`
- `slope`, `intercept`
- `p_fdr`
- `family_id`, `family_n_tested`
- `plot_file`

## Multiple testing correction

- Apply **BH-FDR** within each `behavior_domain x pool_name` family across all tested neural targets.
- Families may differ in size because available targets can differ by chromophore and analysis level.

## Figures

- The script clears `data/results/pooled_mean_correlations/` before each run so stale result files and figures cannot persist.
- `figure_policy = significant_only` emits one condition-faceted plot per tested target/domain group when any exported row for that group has `p_unc < alpha`.
- `figure_policy = all_tested` emits one condition-faceted plot per tested target/domain group.

## Outputs

- `data/results/pooled_mean_correlations/selected_pooled_mean_targets_r.csv`
- `data/results/pooled_mean_correlations/subject_level_pooled_mean_pairs_r.csv`
- `data/results/pooled_mean_correlations/pooled_mean_correlations_r.csv`
- `data/results/pooled_mean_correlations/figures/`

Compatibility note:
- The script filename is retained for continuity, and the main results CSV retains the legacy filename `pooled_mean_correlations_r.csv` even though the workflow exports pooled Pearson correlations.

## Validation requirements

Validation script:
- `tests/validate_pooled_mean_correlations_r.R`

Must verify:
- target selection keeps only `p_fdr < 0.05` format/content rows from the channel/ROI tidy LMM outputs
- stale significant `M_DMPFC` and `M_VMPFC` ROI rows cannot re-enter the pooled-mean target set
- ROI pooling requires at least 2 of 3 channels in every condition and excludes the participant from that ROI when any condition fails
- pool gating restricts format targets to `short/long` and content targets to `education/entertainment`
- channel pooled neural values recover exact known values
- ROI condition means retain exact values when one of three members is pruned
- behavioral pooled rows join to the correct target/domain rows
- channel `0` placeholders are treated as missing rather than as true beta values
- non-significant and interaction-only targets are excluded from the exported target set
- the exported Pearson correlations agree with independently known reference cases
- manual BH agreement within each configured family
- figures are emitted only for rows permitted by the current figure policy
- fail-hard behavior for duplicate IDs and malformed required inputs

# ANALYSIS_SPEC — Standalone Pairwise Behavioral Correlations
Last updated: 2026-08-06

## Scope

Goal: Screen exploratory pairwise associations among the declared behavioral variables in the same merged subject-level CSV used by the primary LMMs and `analyze_pooled_mean_correlations.R`, with a global Benjamini-Hochberg FDR correction across all tested pairs.

Default behavioral variables:
- `sf_education_engagement`
- `sf_entertainment_engagement`
- `lf_entertainment_engagement`
- `lf_education_engagement`
- `diff_short_form_education`
- `diff_short_form_entertainment`
- `diff_long_form_education`
- `diff_long_form_entertainment`
- `age`
- `education_years`
- `sfv_frequency`
- `sfv_daily_duration`
- `asrs_total`
- `yang_pu_total`
- `yang_mot_total`
- `phq_total`
- `gad_total`

## Inputs

Primary CSV input:
- `data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv`

Required support files:
- `data/config/behavior_pairwise_correlation_plan.json`
- `data/config/variable_figure_names.json`
- `data/config/excluded_subjects.json`

## Missingness and data integrity

- `subject_id` is normalized by extracting digits and converting to integer.
- Recruitment order and subject-ID-derived proxies are excluded from the behavioral correlation family.
- Input must contain exactly one row per normalized `subject_id`; duplicates are a hard error.
- All declared behavioral variables must exist in the merged CSV; missing columns are a hard error.
- All declared behavioral variables must be numeric/coercible to numeric; non-numeric values are a hard error.
- Each tested pair uses pairwise complete cases only.
- No imputation is allowed.
- `education_years` is an exploratory approximate-years proxy from the study codebook, not a directly measured continuous education-duration variable.

## Statistical outputs

- One result row per unique unordered behavioral variable pair.
- Pearson correlation is used for every tested pair.
- Report:
  - `var_x`, `var_y`
  - `analysis_status`, `skip_reason`
  - `n_complete`
  - `pearson_r`
  - `p_unc`
  - `p_fdr`
  - `significant_fdr`
  - `ci95_low`, `ci95_high`
- Apply BH-FDR once across every tested row in this workflow.
- `sfv_frequency` and `sfv_daily_duration` are ordinal 0-3 codes but are intentionally treated as numeric, equally spaced scores in this Pearson diagnostic screen.
- `pd_status` is not part of this workflow.

## Figures

- The script clears `data/results/behavior_pairwise_correlations/` before each run so stale result files and figures cannot persist.
- Emit one lower-triangle matrix as both PNG and PDF.
- Matrix cells show Pearson `r` with global BH-FDR `q` values in parentheses.
- Axis labels are loaded from `data/config/variable_figure_names.json`; missing labels for plotted variables are a hard error.
- The upper triangle and diagonal are omitted.

## Outputs

- `data/results/behavior_pairwise_correlations/behavior_pairwise_correlations_r.csv`
- `data/results/behavior_pairwise_correlations/behavior_pairwise_correlations_fdr_r.csv`
- `data/results/behavior_pairwise_correlations/figures/behavior_pairwise_correlation_lower_triangle.png`
- `data/results/behavior_pairwise_correlations/figures/behavior_pairwise_correlation_lower_triangle.pdf`

Output ordering:
- rows are sorted by ascending `p_fdr`
- ties are broken by ascending `p_unc`, then descending absolute `pearson_r`

## Validation requirements

Validation script:
- `tests/validate_behavior_pairwise_correlations_r.R`

Must verify:
- the exact requested variable set is used and both `pd_status` and `recruitment_order_proxy` are absent
- a plan that attempts to include `recruitment_order_proxy` fails hard
- a known positive synthetic pair recovers a perfect positive Pearson correlation
- a known negative synthetic pair recovers a perfect negative Pearson correlation
- ordinal SFV variables are treated as numeric Pearson inputs
- pairwise complete-case handling drops only the subjects missing one member of a pair
- global `p.adjust(..., method = "BH")` matches exported `p_fdr`
- `significant_fdr` is based on `p_fdr < alpha`
- the FDR CSV contains every tested pair, not only significant rows
- the real merged publication input is audited against an independently recomputed reference table for every pair's `n_complete`, Pearson `r`, raw p-value, Fisher CI, BH-FDR q-value, and FDR flag
- the output directory is rebuilt at run start, removing stale CSVs and PNGs
- the lower-triangle PNG and PDF are emitted
