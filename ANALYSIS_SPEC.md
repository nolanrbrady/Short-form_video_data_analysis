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

