# ANALYSIS_SPEC — Homer3 betas + Format×Content (channelwise)
Last updated: 2026-01-30

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
- `data/tabular/homer3_glm_betas_wide.csv`
  - Subject ID column is named `Subject` (e.g., `sub_0001`)
  - Beta columns are wide and follow the pattern:
    - `S##_D##_Cond##_HbO`
    - `S##_D##_Cond##_HbR`
    - Example: `S01_D01_Cond01_HbO`
- `data/tabular/combined_sfv_data.csv`
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
- Both `combined_sfv_data.csv:subject_id` and `homer3_glm_betas_wide.csv:Subject` are treated as **numeric IDs**
  (even if stored as strings).
- IDs are normalized by **extracting digits** and converting to integer (handles `0017` vs `17`, and `sub_0017`).

Expected row structure:
- `homer3_glm_betas_wide.csv` contains **one row per subject**.
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
- `0` values may indicate a pruned channel
- `NaN` values may indicate a pruned channel

Required handling:
- Treat **both `0` and `NaN`** in beta columns as **missing/pruned**.
- Do **not** treat these values as true zero activation.
- Do **not** silently impute these values.

Downstream modeling must handle this explicitly as missingness.

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
- `channel` (e.g., `S01_D01`)
- `chrom` (`HbO` or `HbR`)
- `condition` (one of: `SF_Edu`, `SF_Ent`, `LF_Ent`, `LF_Edu`)
- `beta`
- derived predictors:
  - `format_c`
  - `content_c`

Complete-case rule (within channel/chromophore):
- For a given (channel, chromophore), **only subjects with all 4 conditions present (non-missing beta)** are included.

---

## Primary statistical model (main effects + interaction)

Model form:
- One model per **(channel × chromophore)**.
- Linear mixed model with random intercept for subject:
  - `beta ~ format_c * content_c + (1 | subject_id)`

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
- LMM via `lme4::lmer` with p-values via `lmerTest` (Satterthwaite df approximation).
- Post-hoc via `emmeans`.

Python implementation notes:
- LMM intended via `statsmodels` `MixedLM` (if installed).
- If `statsmodels` is not installed, Python script supports a parse/reshape validation via `--dry-run`.

References (see `CITATIONS.md`):
- Mixed models: Laird & Ware (1982); Bates et al. (2015).
- Satterthwaite df: Satterthwaite (1946); Kuznetsova et al. (2017).

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

# ANALYSIS_SPEC — Retention Length×Content (subject-level LMM)
Last updated: 2026-02-09

## Scope

Goal: Evaluate whether **retention improvement** (post - pre) differs by:
- **Length** (Short vs Long),
- **Content** (Education vs Entertainment),
- and their **Length×Content interaction**,
while accounting for repeated measures (4 within-subject conditions).

## Inputs

Primary CSV input:
- `data/results/homer3_betas_plus_combined_sfv_data_inner_join.csv`

Required columns:
- `subject_id`
- `diff_short_form_education`
- `diff_short_form_entertainment`
- `diff_long_form_education`
- `diff_long_form_entertainment`

## Subject ID + data integrity

- `subject_id` is normalized by extracting digits and converting to integer.
- Input must contain exactly one row per normalized `subject_id`; duplicates are a hard error.
- Required retention columns must all exist; missing columns are a hard error.
- Retention columns must be numeric/coercible to numeric; non-numeric values are a hard error.

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
  - `retention_diff ~ length_c * content_c + (1 | subject_id)`

Reported quantities (for `length_c`, `content_c`, `length_c:content_c`):
- estimate, SE, df, t, uncorrected p, Holm-adjusted p, Wald 95% CI
- sample size fields: `n_subjects`, `n_obs`
- singular-fit flag

R implementation notes:
- Fit with `lmerTest::lmer` (REML).
- p-values from `lmerTest` (Satterthwaite df approximation).

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
- manual Holm agreement with output adjusted p-values
- post-hoc gating behavior (on/off)
- complete-case behavior for `NA`
- retention `0` handling as valid (not missing)
- fail-hard behavior for duplicates, missing columns, and non-numeric retention values

---

# ANALYSIS_SPEC — Engagement Length×Content (subject-level LMM)
Last updated: 2026-02-09

## Scope

Goal: Evaluate whether **engagement ratings** differ by:
- **Length** (Short vs Long),
- **Content** (Education vs Entertainment),
- and their **Length×Content interaction**,
while accounting for repeated measures (4 within-subject conditions).

## Inputs

Primary CSV input:
- `data/results/homer3_betas_plus_combined_sfv_data_inner_join.csv`

Required columns:
- `subject_id`
- `sf_education_engagement`
- `sf_entertainment_engagement`
- `lf_education_engagement`
- `lf_entertainment_engagement`

## Subject ID + data integrity

- `subject_id` is normalized by extracting digits and converting to integer.
- Input must contain exactly one row per normalized `subject_id`; duplicates are a hard error.
- Required engagement columns must all exist; missing columns are a hard error.
- Engagement columns must be numeric/coercible to numeric; non-numeric values are a hard error.

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
  - `engagement ~ length_c * content_c + (1 | subject_id)`

Reported quantities (for `length_c`, `content_c`, `length_c:content_c`):
- estimate, SE, df, t, uncorrected p, Holm-adjusted p, Wald 95% CI
- sample size fields: `n_subjects`, `n_obs`
- singular-fit flag

R implementation notes:
- Fit with `lmerTest::lmer` (REML).
- p-values from `lmerTest` (Satterthwaite df approximation).

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
- manual Holm agreement with output adjusted p-values
- post-hoc gating behavior (on/off)
- complete-case behavior for `NA`
- engagement `0` handling as valid (not missing)
- fail-hard behavior for duplicates, missing columns, and non-numeric engagement values
