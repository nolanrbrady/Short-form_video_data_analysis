## Short-form Video Study — Analysis Repo
Last updated: 03-09-2026
Updated by: Nolan Brady

This repo contains two primary analysis “tracks”:

- **Tabular / behavioral + survey preprocessing** → derived outputs in `data/tabular/generated_data/`
- **fNIRS preprocessing + subject-level GLM + group stats** → outputs in `glm_results/`

### Trigger markers (task conditions)

These are the stimulus/trigger codes used in the fNIRS recordings and the analysis code.

- **Short-Form Education**: `1`
- **Short-Form Entertainment**: `2`
- **Long-Form Entertainment**: `3`
- **Long-Form Education**: `4`
---

## Repo layout (what lives where)

- **`qualtrics/`**: raw Qualtrics export(s)
  - `qualtrics/final_SF_demographic_data.csv`: Qualtrics export with **3 header rows** (MultiIndex columns)
- **`demographic/`**: scripts + outputs for engagement + recall assessment preprocessing
- **`data/tabular/`**: imported/raw tabular files copied into the repo for analysis
- **`data/tabular/generated_data/`**: preprocessing outputs and merged analysis-ready CSVs (merged by `subject_id`)
- **`fnirs_analysis/`**: fNIRS GLM pipeline + second-level (group) inference scripts (MNE / MNE-NIRS)
- **`glm_results/`**: subject GLM outputs + combined tables + group-level outputs
- **`covariate_outputs/`**: covariate-only outputs (clean covariates, missingness audits, correlation heatmaps)

---

## Environment / dependencies (high level)

There is no packaged Python project here; most scripts are “run as a script”.

Common dependencies across scripts:

- Tabular processing: `pandas`, `numpy`
- Plotting: `matplotlib`, `seaborn` (and optionally `scipy` for p-values)
- fNIRS: `mne`, `mne-nirs`, `numpy`, `pandas`
- Stats models (engagement + LME): `statsmodels`, `scipy`

---

## Pipeline A — Tabular preprocessing (demographics + engagement + recall)

All tabular merges are keyed on **`subject_id`** (the Qualtrics study ID from `qualtrics/final_SF_demographic_data.csv`, Q71).

If you want to merge tabular data with fNIRS betas, ensure you have a consistent ID scheme (or a mapping table), because fNIRS pipelines often use IDs like `sub-XXXX` while Qualtrics uses a numeric `subject_id`.

### A1) Engagement preprocessing

**Goal:** Convert per-trigger engagement ratings into per-subject condition means.

1) **Aggregate raw engagement files → `demographic/combined_engagement_data.csv`**

- Script: `demographic/combine_engagement.py`
- **Important:** this script expects raw engagement exports to live at `../../Engagement/` (outside this repo by default).
- **Important:** it uses `DATA_DIR = "../../Engagement"` as a *relative path*, so run it from inside `demographic/` (or edit `DATA_DIR`).

Command (recommended):

```bash
cd demographic
python combine_engagement.py
cd ..
```

Outputs:

- `demographic/combined_engagement_data.csv`: long-format (`subject_id`, `Trigger`, `Category`, `Rating`)

2) **Compute per-subject engagement features → `data/tabular/generated_data/engagement_data_processed.csv`**

- Script: `process_engagement.py`
- Input: `demographic/combined_engagement_data.csv`
- Output: `data/tabular/generated_data/engagement_data_processed.csv`

Command:

```bash
python process_engagement.py
```

Output columns include:

- Condition means: `lf_education_engagement`, `lf_entertainment_engagement`, `sf_education_engagement`, `sf_entertainment_engagement`
- Fail-fast checks: required columns, exact condition labels, trigger/category consistency, finite 0-5 ratings, and at least one observation per subject x condition

### A2) Recall assessment preprocessing

**Goal:** Grade free-text recall answers, compute per-condition improvement (post − pre).

- Script: `demographic/process_recall_assessment.py`
- Inputs (expected on disk): `../../Assessment/pretask_assessment.csv`, `../../Assessment/posttask_assessment.csv`, `../../Assessment/Recall_Assessment_Key.csv`
  - **Important:** the script defines `ASSESSMENT_DIR = "../../Assessment"` as a *relative path*, so run it from inside `demographic/` (or edit `ASSESSMENT_DIR`).

Command (recommended):

```bash
cd demographic
python process_recall_assessment.py
cd ..
```

Outputs:

- `data/tabular/generated_data/recall_assessment_score_diffs.csv`: one row per `subject_id`, columns like `diff_short_form_education`, etc.
- `demographic/recall_assessment_audit_pre.csv`, `demographic/recall_assessment_audit_post.csv`: detailed per-question audits (raw text, normalized text, exact/fuzzy match method)

Optional audit sanity-check:

```bash
python audit_check.py
```

### A3) Socio-demographic + covariate preprocessing

**Goal:** Turn the Qualtrics export into a clean numeric covariate table (ordinal encodings + scale totals).

- Script: `process_sociodemographic.py`
- Input: `qualtrics/final_SF_demographic_data.csv`
- Outputs:
  - `data/tabular/generated_data/socio_demographic_data_processed.csv` (includes `subject_id` for merges)
  - `covariate_outputs/covariates_clean.csv` (covariates only; excludes `subject_id`)
  - `covariate_outputs/covariate_missingness.csv`, `covariate_outputs/covariate_column_audit.csv`, `covariate_outputs/sfv_duration_other_audit.csv`
- Fails hard if the Qualtrics study ID column contains missing or duplicate `subject_id` values.

Command:

```bash
python process_sociodemographic.py
```

### A4) Combine tabular sources into one dataset

- Script: `generate_combined_data.py`
- Inputs:
  - `data/tabular/generated_data/engagement_data_processed.csv`
  - `data/tabular/generated_data/socio_demographic_data_processed.csv`
  - `data/tabular/generated_data/recall_assessment_score_diffs.csv`
- Output:
  - `data/tabular/generated_data/combined_sfv_data.csv`
- Fails hard if any input lacks a unique, non-missing `subject_id` key or if the inner join would not remain one row per subject.

Command:

```bash
python generate_combined_data.py
```

### A5) Import Homer3 FIR basis weights (optional; produced externally)

This repo does **not** run Homer3. If you run Homer3 elsewhere, copy the exported wide table into:

- `data/tabular/homer3_glm_betas_wide_fir_pca.csv`

Current format on disk:

- ID column: `Subject` (values like `sub_0001`)
- Feature columns (FIR): `S##_D##_Cond##_HbO_Basis###` / `S##_D##_Cond##_HbR_Basis###`
- Upstream Homer settings for this export:
  - `idxBasis = 1`
  - `trange = [-10, 130]`
  - `basis spacing = 0.5 s`
  - `Gaussian sigma = 0.5 s`

**Important missingness note:** this file can contain both `0` and `NaN` values that are *stand-ins for channels that were pruned during preprocessing*. Do **not** interpret these as “true zero activation”; treat them as missing/pruned channels in downstream modeling.

To merge with `data/tabular/generated_data/combined_sfv_data.csv`, you will need a shared key:

- Either export Homer betas with a `subject_id` column that matches Qualtrics, **or**
- Create a mapping between `Subject` (e.g., `sub_0001`) and Qualtrics `subject_id` and merge using that.

### A5b) Collapse FIR basis weights to single AUC betas

- Python: `collapse_homer_fir_to_auc.py`
- Shared settings: `data/config/preprocessing_settings.json`

What it does:
- Reads `data/tabular/homer3_glm_betas_wide_fir_pca.csv`.
- Reconstructs the latent HRF from the exported Gaussian basis weights using the shared settings file.
- Baseline-corrects each HRF by subtracting the configured pre-onset mean.
- Computes task-window AUC with trapezoidal integration over the configured time window.
- Writes `data/tabular/generated_data/homer3_glm_betas_wide_auc.csv` with single-beta columns like `S01_D01_Cond01_HbO`.
- Writes `data/tabular/generated_data/homer3_glm_betas_wide_auc.provenance.json` to lock the AUC table to the raw FIR input, settings file, and exact basis configuration used to generate it.
- Treats both `0` and `NaN` as pruned/missing only when the entire basis vector is all-zero or all-`NaN`.
- Fails explicitly on partial missing basis vectors or malformed FIR schemas.

Configured production settings:
- Reconstruction support: `-10 s` to `130 s`
- Baseline window: `-10 s` to `0 s`
- AUC window: `0 s` to `120 s`

Command:

```bash
python collapse_homer_fir_to_auc.py \
  --input-csv data/tabular/homer3_glm_betas_wide_fir_pca.csv \
  --output-csv data/tabular/generated_data/homer3_glm_betas_wide_auc.csv \
  --settings-json data/config/preprocessing_settings.json
```

### A5c) Mask between-subject AUC outliers within each beta column

- Python: `mask_homer_auc_between_subject_outliers.py`

What it does:
- Reads `data/tabular/generated_data/homer3_glm_betas_wide_auc.csv`.
- Treats each `S##_D##_Cond##_HbO/HbR` column independently across subjects.
- Computes the column mean and sample SD using only observed non-missing subject values.
- Replaces values outside `mean +/- 3 SD` with `NaN`.
- Ignores existing `NaN` values from pruned channels when estimating the mean and SD.
- Writes `data/tabular/generated_data/homer3_glm_betas_wide_auc_outliers_masked.csv` for downstream merge/modeling.
- Writes `data/results/homer_auc_outlier_audit.csv` with one row per censored subject-column cell.
- Writes `data/results/homer_auc_outlier_summary.json` with per-column screening counts, skipped-column reasons, and input/output hashes.

Important limitation:
- For a `3 SD` rule, columns with fewer than `11` observed subjects cannot mathematically yield a detected outlier when the sample mean and sample SD are computed from the same data. Those columns are reported as skipped rather than being silently treated as screened.

Command:

```bash
python mask_homer_auc_between_subject_outliers.py \
  --input-csv data/tabular/generated_data/homer3_glm_betas_wide_auc.csv \
  --output-csv data/tabular/generated_data/homer3_glm_betas_wide_auc_outliers_masked.csv \
  --out-audit-csv data/results/homer_auc_outlier_audit.csv \
  --out-summary-json data/results/homer_auc_outlier_summary.json
```

### A5d) Plot reconstructed FIR HRFs for selected subjects (HbO + HbR on same graph)

- Python: `plot_fir_betas_subjects.py`
- Shared settings: `data/config/preprocessing_settings.json`

What it does:
- Reads `data/tabular/homer3_glm_betas_wide_fir_pca.csv` in a streaming/selective way.
- Uses top-of-file variables (no CLI) to choose:
  - `TARGET_SUBJECTS` (default: `sub_0001`, `sub_0002`)
  - `TARGET_CONDITION` (default: `02`)
- Reconstructs the Homer `idxBasis=1` HRF from the exported Gaussian basis weights using the shared settings file, then plots `HbO` and `HbR` overlaid in one figure.
- Treats both `0` and `NaN` as pruned/missing (not true zero activation; no imputation).
- Fails explicitly if any selected subject/channel/chromophore has partially missing basis weights, because exact HRF reconstruction is not possible without the full coefficient vector.

Output (default directory):
- `data/results/fir_beta_plots/`
  - `sub_0001_cond02_S01_D01_fir_hrf_hbo_hbr.png`
  - `sub_0002_cond02_S01_D01_fir_hrf_hbo_hbr.png`

Command:

```bash
python plot_fir_betas_subjects.py
```

### A6) Single entry-point preprocessing + merge + certification

If you already have the required upstream inputs on disk and want one command to
run preprocessing/merge in the correct order with strict integrity checks:

```bash
bash pipeline_preprocess_merge.sh
```

The key file locations for the FIR-to-AUC, merge, and certification steps are
declared as variables at the top of `pipeline_preprocess_merge.sh`, so those
paths can now be changed from the pipeline entry point without editing the
Python/R scripts.

What this entry-point does (in order):
1. Clears `data/results/` at run start.
2. Runs `process_engagement.py`.
3. Runs `process_sociodemographic.py`.
   Fails hard on missing/duplicate study IDs in the Qualtrics-derived covariate table.
4. Runs `generate_combined_data.py`.
   Fails hard if any tabular input violates the one-row-per-subject merge contract.
5. Runs `collapse_homer_fir_to_auc.py`.
6. Runs `validate_homer_fir_auc_conversion.py` and fails hard if excluded FIR basis vectors are not represented as `NaN` in the derived AUC table or if the AUC provenance sidecar does not match the current raw FIR export + settings JSON.
7. Runs `mask_homer_auc_between_subject_outliers.py` and writes a separate outlier-masked AUC table plus audit artifacts.
8. Runs `merge_homer3_betas_with_combined_data.R` using `data/tabular/generated_data/homer3_glm_betas_wide_auc_outliers_masked.csv`.
9. Runs `certify_preprocess_merge_integrity.py` and fails hard if merge invariants are violated.

Required inputs for this entry-point:
- `demographic/combined_engagement_data.csv`
- `data/tabular/generated_data/recall_assessment_score_diffs.csv`
- `qualtrics/final_SF_demographic_data.csv`
- the raw FIR CSV pointed to by `HOMER_RAW_FIR_CSV` in `pipeline_preprocess_merge.sh`

Certification outputs:
- `data/results/preprocess_merge_certification.json`
- `data/results/preprocess_merge_id_audit.csv`
- `data/results/preprocess_merge_dropped_ids.csv`
- `data/results/homer_auc_outlier_audit.csv`
- `data/results/homer_auc_outlier_summary.json`

---

## Pipeline C — Homer3 betas + Format×Content (channelwise LMM)

**Goal:** Merge the externally-produced Homer3 betas table with the combined tabular dataset and test whether
**Format depends on Content (and vice versa)** at the level of **channelwise prefrontal activation**, for **HbO** and **HbR**.

### C0) Prerequisites / required prior work

Inputs required:
- `data/tabular/homer3_glm_betas_wide_fir_pca.csv` (externally produced; this is the production FIR export used for downstream results)
- `data/tabular/generated_data/homer3_glm_betas_wide_auc.csv` (generated locally by `collapse_homer_fir_to_auc.py`)
- `data/tabular/generated_data/homer3_glm_betas_wide_auc.provenance.json` (generated locally by `collapse_homer_fir_to_auc.py`)
- `data/tabular/generated_data/homer3_glm_betas_wide_auc_outliers_masked.csv` (generated locally by `mask_homer_auc_between_subject_outliers.py`)
- `data/tabular/generated_data/combined_sfv_data.csv` (produced by the tabular preprocessing pipeline)

Requirements / invariants:
- `combined_sfv_data.csv` must contain **exactly one row per subject** (unique `subject_id`).
- `homer3_glm_betas_wide_auc.csv` must contain **exactly one row per subject** (unique `Subject` after normalization).
- `homer3_glm_betas_wide_auc.provenance.json` must match the current raw FIR export and `data/config/preprocessing_settings.json`.
- `homer3_glm_betas_wide_auc_outliers_masked.csv` must contain **exactly one row per subject** and preserve the same beta schema as the raw AUC table.
- In the raw derived AUC table, pruned channels are carried forward as `NaN`.
- In the outlier-masked AUC table, pruned channels and between-subject `mean +/- 3 SD` censored values are both encoded as `NaN`; downstream modeling treats both as missing (do not impute).

Recommended prior step (if you haven’t generated it yet):
- Run `bash pipeline_preprocess_merge.sh` for end-to-end preprocessing + merge + certification.
- Or run `generate_combined_data.py` to (re)build `data/tabular/generated_data/combined_sfv_data.csv` if orchestrating manually.

### C1) Inner-join Homer3 betas with combined tabular data

- R: `merge_homer3_betas_with_combined_data.R`

What it does:
- Normalizes IDs by extracting digits (handles `0017` vs `17`, and `sub_0017`-style IDs).
- Performs an **INNER JOIN** and writes a merged “one row per subject” CSV for inspection / downstream use.

Notes:
- IDs are normalized by extracting digits (handles `0017` vs `17`, and `sub_0017`-style IDs).
- Output is one row per subject containing both demographics/behavior and beta columns.
- Scripts **fail hard** if either input contains duplicate `subject_id` values after normalization (expected exactly one row per subject).

Example (R):

```bash
Rscript merge_homer3_betas_with_combined_data.R \
  --homer_csv data/tabular/generated_data/homer3_glm_betas_wide_auc_outliers_masked.csv \
  --combined_csv data/tabular/generated_data/combined_sfv_data.csv \
  --out_csv data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv
```

### C1b) QC report from Homer3 beta-wide table (subject-level channel exclusion)

- Python: `fnirs_analysis/homer_betas_qc.py`

Goal:
- Quantify subject-level excluded/pruned channel burden directly from the imported Homer3 wide table.

Definitions used:
- A beta entry is **excluded** if it is `0` or `NaN` (per project data-integrity note).
- A channel is **available** for a condition only if **both** `HbO` and `HbR` are non-excluded.
- Primary bad-channel rule: excluded in **>= 2 of 4** conditions (`--bad-channel-min-excluded-conds 2`).
- Task pass rule: condition has **strictly > 50%** channels available.

Sensitivity outputs:
- The script also reports bad-channel counts/lists for:
  - any-condition excluded (>=1/4)
  - all-condition excluded (4/4)

Outputs:
- `data/results/homer3_betas_qc_subject_level.csv` (one row per subject)
- `data/results/homer3_betas_qc_cohort_summary.csv` (single-row cohort summary)

Example:

```bash
python fnirs_analysis/homer_betas_qc.py \
  --input-csv data/tabular/homer3_glm_betas_wide_fir_pca.csv \
  --output-csv data/results/homer3_betas_qc_subject_level.csv \
  --summary-csv data/results/homer3_betas_qc_cohort_summary.csv
```

### C1c) Centralized participant exclusions for inferential scripts

Use a single exclusion manifest so subject filtering is consistent across all inferential endpoints.

- File: `data/config/excluded_subjects.json`
- Format: top-level JSON array of subject IDs (e.g., `["sub_0041", "sub_0044", "sub_0050"]`)
- Matching rule: IDs are normalized by numeric component (e.g., `sub_0041`, `0041`, and `41` are treated as the same participant)
- Missing-ID behavior: if an ID in the exclusion file is not present in a given analysis input, the script prints a warning and continues

This manifest is consumed by:
- `analyze_format_content_lmm_channelwise.R`
- `analyze_format_content_lmm_roi.R`
- `analyze_retention_format_content_lmm.R`
- `analyze_engagement_format_content_lmm.R`

Override path (optional):

```bash
Rscript analyze_format_content_lmm_channelwise.R \
  --exclude_subjects_json data/config/excluded_subjects.json
```

### C2) Channelwise within-subject inference: Format, Content, and Format×Content

- Python: `analyze_format_content_lmm_channelwise.py`
- R: `analyze_format_content_lmm_channelwise.R`

Order / dependencies:
- Assumes Pipeline C0 prerequisites are satisfied.
- Preferred input is the pre-merged table `data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv`
  so covariates and beta columns are available in one file.

Model (per channel × chromophore):
- LMM: `beta ~ format_c * content_c + (1 | subject_id)`
- Coding: `format_c = -0.5 (Short), +0.5 (Long)`; `content_c = -0.5 (Entertainment), +0.5 (Education)`

Pruned channels / missingness policy:
- In the derived FIR-to-AUC beta table, pruned channels are encoded as **`NaN`** (do **not** impute).
- Default behavior is **complete-case within channel**: drop subjects missing any of the 4 conditions for that channel.

Multiple testing:
- BH-FDR is applied **separately** per chromophore (**HbO**, **HbR**) and **separately per effect**
  (Format, Content, Interaction), across channels.
- **Reminder (ask before publication):** consider whether you want a broader correction family
  (e.g., across effects and/or chromophores) depending on the final reporting plan.

Post-hoc (only if the interaction is FDR-significant for that channel/chromophore):
- All 6 pairwise contrasts among the 4 conditions.
- **No multiple-test correction** in post-hoc contrasts (per study instruction).
- Post-hoc `mean_diff` is reported as `condition_a - condition_b`.
- Post-hoc outputs include `stat_type` (`t` vs `z`) to indicate whether emmeans used a t-statistic (finite df) or asymptotic z.

Minimum sample gating:
- Models are only fit for channel/chrom pairs with at least `min_subjects` complete-case subjects (default: 6).
- Scripts print a warning summary (count + examples) for models skipped due to `min_subjects`.

Large-sample df limits (R only):
- `emmeans` may disable some denominator-df adjustments when the number of observations is large (prints a note).
- If you explicitly want to enable those adjustments (may be slow / memory-heavy), pass:
  - `--pbkrtest_limit <N>` and/or `--lmerTest_limit <N>`

Example (R; preferred merged input path):

```bash
Rscript analyze_format_content_lmm_channelwise.R \
  --input_csv data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv \
  --exclude_subjects_json data/config/excluded_subjects.json
```

Outputs:
- `data/results/format_content_lmm_main_effects_*.csv`
- `data/results/format_content_lmm_main_effects_tidy_r.csv` (R: spec-compliant tidy main-effects table)
- `data/results/format_content_lmm_posthoc_pairwise_*.csv`
- Subject exclusions are applied from `data/config/excluded_subjects.json` (or `--exclude_subjects_json` override).

Optional output filtering (R only):
- `analyze_format_content_lmm_channelwise.R` contains `FILTER_MAIN_EFFECTS_TO_FDR_SIG_ONLY` (default `FALSE`) to write only rows with `p_fdr < 0.05` to the main-effects CSVs for quick review.

### C2b) ROI-wise within-subject inference: Format, Content, and Format×Content

- R: `analyze_format_content_lmm_roi.R`

Order / dependencies:
- Assumes Pipeline C0 prerequisites are satisfied.
- Uses the same pre-merged input table as C2:
  `data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv`

ROI definition input:
- `data/config/roi_definition.json` (strict JSON object: `ROI -> [channel_ids]`)
- Channel IDs must match Homer naming (example: `S01_D01`).
- Script fails fast on malformed JSON, overlapping channel assignments, or ROI channels absent from the data.

ROI beta construction:
- For each `subject × ROI × chrom × condition`, ROI beta is the arithmetic mean
  across available (non-missing) channels in that ROI.
- In the derived FIR-to-AUC beta table, pruned channels are encoded as `NaN` and are not imputed.

Model / inference:
- LMM (per ROI × chrom): `beta ~ format_c * content_c + (1 | subject_id)`
- Same coding and interaction-gated post-hoc policy as C2.
- BH-FDR is applied separately per chromophore and per effect across ROIs.

Example:

```bash
Rscript analyze_format_content_lmm_roi.R \
  --input_csv data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv \
  --roi_json data/config/roi_definition.json \
  --exclude_subjects_json data/config/excluded_subjects.json
```

Outputs:
- `data/results/format_content_lmm_roi_main_effects_r.csv`
- `data/results/format_content_lmm_roi_main_effects_tidy_r.csv`
- `data/results/format_content_lmm_roi_posthoc_pairwise_r.csv`

Validation:

```bash
Rscript tests/validate_pipeline_c_roi_r.R
```

### C3) Retention within-subject inference: Length, Content, and Length×Content

- R: `analyze_retention_format_content_lmm.R`

Input:
- `data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv`
  (must contain `subject_id` and:
  `diff_short_form_education`, `diff_short_form_entertainment`,
  `diff_long_form_education`, `diff_long_form_entertainment`)

Model:
- LMM: `retention_diff ~ length_c * content_c + (1 | subject_id)`
- Coding: `length_c = -0.5 (Short), +0.5 (Long)`; `content_c = -0.5 (Entertainment), +0.5 (Education)`

Missingness policy:
- Complete-case by subject across the 4 retention conditions.
- Retention `0` values are treated as valid values (not missing).
- Subject exclusions are applied from `data/config/excluded_subjects.json` (or `--exclude_subjects_json` override).

Multiple testing:
- Holm correction across the 3 planned omnibus effects:
  - Length
  - Content
  - Length×Content interaction

Post-hoc:
- Run all 6 pairwise condition contrasts only if interaction adjusted p `< alpha`.
- Pairwise p-values are uncorrected (`adjust="none"`).

Example:

```bash
Rscript analyze_retention_format_content_lmm.R \
  --input_csv data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv \
  --out_main_csv data/results/retention_format_content_lmm_main_effects_r.csv \
  --out_posthoc_csv data/results/retention_format_content_lmm_posthoc_pairwise_r.csv
```

Outputs:
- `data/results/retention_format_content_lmm_main_effects_r.csv`
- `data/results/retention_format_content_lmm_posthoc_pairwise_r.csv`

Validation:

```bash
Rscript tests/validate_retention_pipeline_r.R
```

### C4) Engagement within-subject inference: Length, Content, and Length×Content

- R: `analyze_engagement_format_content_lmm.R`

Input:
- `data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv`
  (must contain `subject_id` and:
  `sf_education_engagement`, `sf_entertainment_engagement`,
  `lf_education_engagement`, `lf_entertainment_engagement`)

Model:
- LMM: `engagement ~ length_c * content_c + (1 | subject_id)`
- Coding: `length_c = -0.5 (Short), +0.5 (Long)`; `content_c = -0.5 (Entertainment), +0.5 (Education)`

Missingness policy:
- Complete-case by subject across the 4 engagement conditions.
- Engagement `0` values are treated as valid values (not missing).
- Subject exclusions are applied from `data/config/excluded_subjects.json` (or `--exclude_subjects_json` override).

Multiple testing:
- Holm correction across the 3 planned omnibus effects:
  - Length
  - Content
  - Length×Content interaction

Post-hoc:
- Run all 6 pairwise condition contrasts only if interaction adjusted p `< alpha`.
- Pairwise p-values are uncorrected (`adjust="none"`).

Example:

```bash
Rscript analyze_engagement_format_content_lmm.R \
  --input_csv data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv \
  --out_main_csv data/results/engagement_format_content_lmm_main_effects_r.csv \
  --out_posthoc_csv data/results/engagement_format_content_lmm_posthoc_pairwise_r.csv
```

Outputs:
- `data/results/engagement_format_content_lmm_main_effects_r.csv`
- `data/results/engagement_format_content_lmm_posthoc_pairwise_r.csv`

Validation:

```bash
Rscript tests/validate_engagement_pipeline_r.R
```

### C4b) Targeted Pearson correlations for selected channel/ROI follow-up

- R: `analyze_correlational_relationships.R`

Purpose:
- Run targeted follow-up Pearson correlations between selected sociodemographic variables and raw-condition neural values from:
  - `S04_D02` for both `HbO` and `HbR`
  - `R_DLPFC (HbR)`, `L_DLPFC (HbO)`, `M_DMPFC (HbO)`, `L_DMPFC (HbO)`
- Use the same merged input table and subject-exclusion manifest as the other R analyses.

Predictors:
- `age`
- `sfv_daily_duration`
- `sfv_frequency`
- `phq_total`
- `asrs_total`
- `gad_total`
- `yang_pu_total`
- `yang_mot_total`
- all `diff_*` columns
- all `*_engagement` columns

Neural target construction:
- Channel targets are the raw merged beta columns for `S04_D02`, split by chromophore and condition.
- ROI targets are arithmetic means across available non-missing channels in the ROI for each `subject x chrom x condition`.
- ROI channel membership is read from `data/config/roi_definition.json`.

Missingness policy:
- Pairwise complete cases only: for each `predictor x neural target` correlation, drop subjects missing either value for that pair.
- Do not impute pruned channels or missing covariates.
- Subjects excluded in `data/config/excluded_subjects.json` are removed before any pairwise filtering.

Multiple testing:
- BH-FDR is applied within each configured neural-target x chromophore family.
- Family membership and the explicit predictor list are declared in `data/config/correlational_analysis_plan.json`.
- Output includes both raw `p_unc` and adjusted `p_fdr`.

Example:

```bash
Rscript analyze_correlational_relationships.R \
  --input_csv data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv \
  --roi_json data/config/roi_definition.json \
  --analysis_plan_json data/config/correlational_analysis_plan.json \
  --exclude_subjects_json data/config/excluded_subjects.json \
  --out_csv data/results/correlational_relationships/pairwise_correlations_r.csv \
  --out_fig_dir data/results/correlational_relationships/figures
```

Outputs:
- `data/results/correlational_relationships/pairwise_correlations_r.csv`
- `data/results/correlational_relationships/figures/` (one scatterplot with linear fit per tested pair)

Validation:

```bash
Rscript tests/validate_correlational_relationships_r.R
```

### C5) Monte Carlo type-I error calibration (all inferential pipelines)

- Script: `tests/calibrate_type1_error_r.R`
- Purpose:
  - Run repeated **null-effect** synthetic datasets through all four R inferential scripts:
    - `analyze_format_content_lmm_channelwise.R`
    - `analyze_format_content_lmm_roi.R`
    - `analyze_retention_format_content_lmm.R`
    - `analyze_engagement_format_content_lmm.R`
  - Estimate empirical type-I error rates from adjusted p-values (`p_fdr`) per effect.
  - Fail if any observed rate exceeds a configured upper bound.

Example:

```bash
# Default calibration run
Rscript tests/calibrate_type1_error_r.R

# Faster smoke run
Rscript tests/calibrate_type1_error_r.R --n_reps 20 --type1_upper_bound 0.20

# Stricter manuscript QA run
Rscript tests/calibrate_type1_error_r.R --n_reps 200 --type1_upper_bound 0.10
```

### C6) Monte Carlo type-II error calibration (all inferential pipelines)

- Script: `tests/calibrate_type2_error_r.R`
- Purpose:
  - Run repeated **non-null** synthetic datasets through all four R inferential scripts.
  - Estimate empirical power and type-II error (`1 - power`) from adjusted p-values (`p_fdr`) per effect.
  - Fail if any observed type-II error exceeds a configured upper bound.

Example:

```bash
# Default calibration run
Rscript tests/calibrate_type2_error_r.R

# Faster smoke run
Rscript tests/calibrate_type2_error_r.R --n_reps 20 --type2_upper_bound 0.60

# Stricter manuscript QA run
Rscript tests/calibrate_type2_error_r.R --n_reps 200 --type2_upper_bound 0.25
```

### A7) Correlation diagnostics / heatmaps (optional)

- Script: `covariate_correlation_analysis.py`
- Two presets:
  - **`covariates` preset**: correlations on `covariate_outputs/covariates_clean.csv`
  - **`combined` preset**: correlations on `data/tabular/generated_data/combined_sfv_data.csv`

Recommended commands:

```bash
# Covariates-only correlations + p-values (if SciPy installed)
python covariate_correlation_analysis.py --preset covariates --out-dir covariate_outputs

# Combined dataset correlations + heatmap saved alongside the combined dataset
python covariate_correlation_analysis.py --preset combined --out-dir data/tabular/generated_data
```

---

## Pipeline B — fNIRS preprocessing + first-level GLM (MNE / MNE-NIRS)

### B1) Subject-level GLM from SNIRF

- Script: `fnirs_analysis/fnirs_analysis.py`
- Purpose:
  - Find `.snirf` files under a configurable root
  - Rename triggers into the four task condition labels
  - Preprocess intensity → OD → (optional SCI pruning) → (optional TDDR) → Beer–Lambert → (optional filter)
  - Build a first-level design matrix and fit a GLM per run
- **Configuration:** this script is configured via constants at the top (e.g., `DATA_ROOT`, `OUTPUT_ROOT`, `STIMULUS_DURATION_SEC`, `SUBJECT_ID_REGEX`)
  - Tip: if you store SNIRF files inside this repo, a common choice is setting `DATA_ROOT = "./data/fnirs"`.

Command (after configuring paths in the script):

```bash
python fnirs_analysis/fnirs_analysis.py
```

Outputs (per subject under `glm_results/<subject_id>/`):

- `*_glm.h5`: serialized GLM object
- `*_glm_results.csv`: tidy dataframe of GLM estimates per channel/condition/chromophore
- `*_ALL_runs_glm_results.csv`: append-only “all runs” table for that subject

### B2) Quality Control & Exclusion Criteria

The pipeline uses the following criteria for subject-level exclusion (generated via `fnirs_analysis/qc_check.py`):

1.  **Scalp Coupling Index (SCI)**: Subject average SCI must be **≥ 0.8**.
2.  **Bad Channel Count**: Subjects with **> 50% bad channels** (where a channel is bad if its average SCI < 0.8) are excluded.
3.  **Minimum Usable Trials**: At least **50% of trials (2/4)** per condition must be usable.
    - A trial is "usable" if its window-level SCI ≥ 0.8 and it does not exceed the bad channel threshold.

For methodological justifications and citations, see `fnirs_analysis/fnirs_preprocess_justifications.md`.

### B3) Combine first-level outputs across subjects (and across runs within subject)

- Script: `fnirs_analysis/combine_glm_output.py`
- Input: `glm_results/<subject_id>/*_ALL_runs_glm_results.csv`
- Outputs:
  - `glm_results/combined_glm_long.csv`: aggregated across runs within subject using inverse-variance weighting
  - `glm_results/combined_glm_long_runs.csv`: run-level long table (preferred for run-level LME/group methods)
  - `glm_results/combined_matrices/*.csv`: convenience wide matrices per condition/chroma

Command:

```bash
python fnirs_analysis/combine_glm_output.py --root glm_results --out glm_results --chroma both
```

### B4) Group-level inference options

There are two group pipelines provided; both read the combined runs-level table.

1) **Two-stage 2×2 contrasts with selective channel follow-up**

- Script: `fnirs_analysis/group_analysis_anova.py`
- Input default: `glm_results/combined_glm_long_runs.csv`
- Outputs (in `glm_results/`):
  - `group_<chroma>_global_effects.csv`
  - `group_<chroma>_channel_effects.csv`

Example:

```bash
python fnirs_analysis/group_analysis_anova.py --input glm_results/combined_glm_long_runs.csv --chroma hbo
```

2) **Per-channel linear mixed effects (LMM) with FDR across channels + gated post-hocs**

- Script: `fnirs_analysis/group_analysis_lme.py`
- Input default: `glm_results/combined_glm_long_runs.csv`
- Outputs (in `glm_results/`):
  - `group_<chroma>_main_effects.csv`
  - `group_<chroma>_posthoc_pairs.csv`

Example:

```bash
python fnirs_analysis/group_analysis_lme.py --input glm_results/combined_glm_long_runs.csv --chroma hbo
```

Methodology notes / planned improvements live in:

- `fnirs_analysis/FNIRS_TODO.md`

---

## “What does each file do?” (quick reference)

### Tabular / survey / engagement

- `process_sociodemographic.py`: Qualtrics demographics + psych scales → numeric covariates + audits (`data/tabular/generated_data/` + `covariate_outputs/`)
- `demographic/combine_engagement.py`: raw engagement CSVs (external `../../Engagement`) → `demographic/combined_engagement_data.csv` (+ runs basic statsmodels analyses)
- `process_engagement.py`: per-subject engagement condition means → `data/tabular/generated_data/engagement_data_processed.csv`
- `demographic/process_recall_assessment.py`: grade pre/post recall (external `../../Assessment`) → diffs CSV + detailed audit CSVs
- `generate_combined_data.py`: merge engagement + sociodemographics + recall diffs → `data/tabular/generated_data/combined_sfv_data.csv`
- `data/tabular/homer3_glm_betas_wide_fir_pca.csv`: externally produced production Homer3 FIR basis-weight table (wide table; copied into this repo; contains `0` and `NaN` as stand-ins for pruned channels)
- `data/tabular/generated_data/homer3_glm_betas_wide_auc.csv`: locally derived single-beta table produced by `collapse_homer_fir_to_auc.py` from the FIR basis weights
- `data/tabular/generated_data/homer3_glm_betas_wide_auc.provenance.json`: sidecar provenance record for the raw derived AUC table
- `data/tabular/generated_data/homer3_glm_betas_wide_auc_outliers_masked.csv`: between-subject outlier-masked AUC table consumed by merge/modeling
- `data/results/homer_auc_outlier_audit.csv`: row-level audit of censored subject-column AUC cells
- `data/results/homer_auc_outlier_summary.json`: machine-readable summary of between-subject AUC screening
- `validate_homer_fir_auc_conversion.py`: hard-fail lint that checks excluded FIR basis vectors map to `NaN` in the derived AUC CSV and that the provenance sidecar matches the current raw FIR export + settings JSON
- `data/config/excluded_subjects.json`: central participant-exclusion manifest consumed by inferential R analyses
- `covariate_correlation_analysis.py`: Spearman correlation tables + heatmaps (covariates-only or combined dataset)
- `analyze_correlational_relationships.R`: targeted Pearson follow-up correlations for selected raw-condition channel/ROI neural targets, with per-pair figures
- `plot_fir_betas_subjects.py`: plots selected-subject FIR betas for one condition with HbO/HbR overlaid (streaming/selective read; top-of-file config)
- `audit_check.py`: consistency checks for the recall assessment audit CSVs
- `demographic/engagement_stats.py`: helper functions used by `combine_engagement.py` for ANOVA/OLS/mixed model

### fNIRS (MNE/MNE-NIRS)

- `fnirs_analysis/fnirs_analysis.py`: SNIRF preprocessing + first-level GLM per run; writes per-subject outputs under `glm_results/`
- `fnirs_analysis/homer_betas_qc.py`: subject-level QC from imported Homer beta-wide table (`0`/`NaN` treated as excluded channels); writes subject and cohort QC CSVs
- `r_subject_exclusions.R`: shared R helpers that enforce centralized subject exclusions from JSON across inferential scripts
- `fnirs_analysis/combine_glm_output.py`: combine subject GLM CSVs; writes `combined_glm_long*.csv` + `combined_matrices/`
- `fnirs_analysis/group_analysis_anova.py`: group-level two-stage contrast testing + selective channel inference
- `fnirs_analysis/group_analysis_lme.py`: per-channel mixed-effects + FDR + gated post-hocs
- `fnirs_analysis/FNIRS_TODO.md`: methodological notes and a change list for the fNIRS pipeline

### Misc / exploratory

- `analyze_sfv_demographics.py`: standalone demographic summary + figures for an input `SFV_demo_data.csv` colocated with the script (older/one-off analysis)
- `analyze_format_content_lmm_channelwise.R`: within-subject channelwise LMM for Format/Content/Interaction + interaction-gated post-hoc contrasts
- `analyze_format_content_lmm_roi.R`: within-subject ROI-wise LMM using `data/config/roi_definition.json` + interaction-gated post-hoc contrasts
- `analyze_retention_format_content_lmm.R`: within-subject retention LMM for Length/Content/Interaction + interaction-gated post-hoc contrasts
- `analyze_engagement_format_content_lmm.R`: within-subject engagement LMM for Length/Content/Interaction + interaction-gated post-hoc contrasts
- `data_quality_evaluation.py`: placeholder for a future subject-level fNIRS QC report (see `TODO.md`)
- `TODO.md`: high-level analysis TODO list
