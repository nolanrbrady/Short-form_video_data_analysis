## Short-form Video Study — Analysis Repo
Last updated: 2026-09-30
Updated by: Codex

This repo contains two primary analysis “tracks”:

- **Tabular / behavioral + survey preprocessing** → derived outputs in `data/tabular/generated_data/`
- **Imported Homer3 FIR weights → AUC summaries + group statistics** → outputs in `data/results/`

### Trigger markers (task conditions)

These are the stimulus/trigger codes used in the fNIRS recordings and the analysis code.
The single source of truth for the engagement preprocessing mapping is
`data/config/engagement_condition_map.json`.

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
- **`covariate_outputs/`**: current clean covariates and preprocessing/missingness audits
- **`unused/`**: archived analyses, dedicated tests/configuration, old diagnostic outputs, and historical planning material; see `unused/README.md` and `unused/move_manifest.csv`

---

## Current result-producing workflow

The active analysis uses `pipeline_preprocess_merge.sh`, the channel/ROI,
retention, and engagement R mixed-model scripts, `analyze_behavior_pairwise_correlations.R`,
`analyze_pooled_mean_correlations.R`, `create_demographics_table.py`, and the two
R publication-distribution plotting scripts. Their raw inputs, configuration,
shared helpers, provenance checks, and validation tests remain in place.
Optional recall audits, FIR plots, retention sensitivity, and citation checks
remain available as QC/robustness tools, not additional primary results.

Historical alternative analyses and instructions are archived under `unused/`;
they are not dependencies of the current results. The older covariate-sensitivity
decision record is retained at `unused/COVARIATE_SENSITIVITY_ANALYSIS.md` for
provenance, not as a table of current numerical findings.

## Environment / dependencies (high level)

There is no packaged Python project here; most scripts are “run as a script”.

Common dependencies across scripts:

- Tabular processing: `pandas`, `numpy`
- Plotting/correlation diagnostics: `matplotlib`, `seaborn`, `scipy`
- Current neural summary/analysis: Python `numpy`/`pandas`; R `lme4`, `lmerTest`, `emmeans`, `pbkrtest`
- Stats models (engagement + LME): `statsmodels`, `scipy`

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
- Condition map: `data/config/engagement_condition_map.json`
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
- Invalid-item exclusion manifest: `data/config/recall_invalid_questions.json`
  - Excluded from both pre-task and post-task retention scoring: `Q5`, `Q6`, `Q7`, `Q8`, `Q10`, `Q26`, `Q28`, plus pre-task aliases `Q35` for invalid `Q6` and `Q36` for invalid `Q7`.
  - The excluded items remain in the pre/post audit CSVs with `method = excluded_invalid_question`, but they are not included in the per-condition mean-correctness denominator.
- Question-alias manifest: `data/config/recall_question_aliases.json`
  - Maps pre-task `Q39` to canonical/keyed item `Q22` so the same `Wa anta` item is scored in both pre-task and post-task exports.

Command (recommended):

```bash
cd demographic
python process_recall_assessment.py
cd ..
```

Outputs:

- `data/tabular/generated_data/recall_assessment_score_diffs.csv`: one row per `subject_id`, columns like `diff_short_form_education`, etc.
- `demographic/recall_assessment_audit_pre.csv`, `demographic/recall_assessment_audit_post.csv`: detailed per-question audits (raw text, normalized text, exact/fuzzy/excluded match method)

Optional audit sanity-check:

```bash
python audit_check.py
```

### A3) Socio-demographic + covariate preprocessing

**Goal:** Turn the Qualtrics export into a clean numeric covariate table (indicator encodings, ordinal encodings, and scale totals).

- Script: `process_sociodemographic.py`
- Input: `qualtrics/final_SF_demographic_data.csv`
- Outputs:
  - `data/tabular/generated_data/socio_demographic_data_processed.csv` (includes `subject_id` for merges)
  - `covariate_outputs/covariates_clean.csv` (covariates only; excludes `subject_id`)
  - `covariate_outputs/covariate_missingness.csv`, `covariate_outputs/covariate_column_audit.csv`, `covariate_outputs/sfv_duration_other_audit.csv`
- Fails hard if the Qualtrics study ID column contains missing or duplicate `subject_id` values.
- Race/ethnicity, sex, and education are encoded as explicit indicator columns. Q12 race supports
  comma-separated multiple selections, so selected race categories are preserved as separate
  indicators; Q11 Hispanic/Latino identity is retained as its own binary ethnicity indicator.
  `education_years` is also added as an exploratory approximate-years proxy from
  `data/config/education_years_encoding.json`; the original education indicators are retained.
  Unknown category labels fail hard instead of being silently coerced.

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

### A4b) Create descriptive demographics table

- Script: `create_demographics_table.py`
- Input: `data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv`
- Shared exclusion manifest: `data/config/excluded_subjects.json`
- Subject IDs are normalized before exclusions, so Homer-style IDs such as `sub_0050` match numeric tabular IDs such as `50`.
- The script prints an exclusion audit (`listed`, `matched`, `removed`, `remaining`) before printing the table.

Command:

```bash
python create_demographics_table.py
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
- Treats a basis vector as pruned/missing only when it is **exactly all-zero** or all-`NaN`; individual zero coefficients and arbitrarily small finite nonzero coefficients remain valid. Pruning uses no approximate-zero amplitude tolerance.
- Fails explicitly on partial missing basis vectors or malformed FIR schemas.

This sentinel rule follows the imported-data policy in `AGENTS.md`, not a
physiological amplitude cutoff. Transparent missingness/QC reporting follows
Yücel et al. (2021); see `CITATIONS.md`. Regression checks include tiny nonzero
vectors that the former `np.allclose(beta, 0.0)` check incorrectly pruned:

```bash
python tests/validate_fir_auc_adapter_py.py
```

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

### A5f) Plot beta-value distributions for significant channelwise/ROI LMM hits

- R: `plot_significant_beta_value_distribution.R`
- Inputs:
  - `data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv`
  - `data/config/roi_definition.json`
  - `data/config/excluded_subjects.json`
  - `data/results/format_content_lmm_main_effects_tidy_r.csv`
  - `data/results/format_content_lmm_roi_main_effects_tidy_r.csv`

What it does:
- Reads the tidy LMM result tables and selects rows with `p_fdr < 0.05` by default.
- Applies the shared subject-exclusion manifest so the plotted subject set matches the inferential analyses.
- Rebuilds channel-level beta rows from the merged wide beta table and recreates ROI betas using the inferential analysis's 2-of-3 channel-completeness rule.
- Uses the same complete-case rule as the LMM scripts within each plotted `channel x chrom` or `roi x chrom` unit.
- Fails hard if literal `0` values appear in beta columns, because this project treats zero placeholders as invalid stand-ins for pruned channels rather than true zero activation.
- For significant `format` or `content` main effects, plots subject-level marginal means collapsed across the orthogonal factor.
- For significant `interaction` effects, plots the four raw condition beta distributions (`SF_Edu`, `SF_Ent`, `LF_Ent`, `LF_Edu`) without collapsing.
- Uses violin density envelopes with jittered subject-level points plus mean and +/- 1 SD overlays; subject trajectories are not connected by lines.
- Writes one PNG per significant hit plus a combined audit CSV of the raw rows and plotted values used in each figure.
- Clears previously generated beta-distribution PNGs and the audit CSV before rebuilding them so plots for effects that no longer pass FDR cannot remain stale.

Default output directory:
- `data/results/beta_value_distribution/`

Command:

```bash
Rscript plot_significant_beta_value_distribution.R
```

### A5g) Plot behavioral score distributions by condition

- R: `plot_behavior_score_distributions.R`
- Inputs:
  - `data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv`
  - `data/config/excluded_subjects.json`

What it does:
- Applies the shared subject-exclusion manifest before plotting so behavioral figures use the same excluded-participant policy as the inferential analyses.
- Plots engagement scores and recall/retention improvement scores (`post - pre`) as separate publication-ready PNG figures across the four conditions.
- Adds content main-effect PNG figures for engagement and recall/retention by averaging Education and Entertainment scores across Short and Long within each complete-case subject.
- Adds retention length main-effect PNG figures by averaging content conditions within short and long blocks (Short = SF_Edu + SF_Ent; Long = LF_Edu + LF_Ent) for each complete-case subject.
- Uses domain-specific complete-case filtering: subjects need all four engagement values for the engagement figure and all four retention values for the retention figure.
- Treats behavioral `0` values as valid observed scores and does not impute missing values.
- Uses violin density envelopes with jittered subject-level points plus mean and +/- 1 SD overlays.
- Writes PNG figures plus audit and summary CSVs.

Default output directory:
- `data/results/behavior_score_distribution/`

Command:

```bash
Rscript plot_behavior_score_distributions.R
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
2. Runs `demographic/process_recall_assessment.py` from its required `demographic/`
   working directory. Regenerates recall differences and both question-level
   audits using the current invalid-question and alias manifests. A scoring
   failure stops the pipeline before any merge; an existing recall CSV is never
   used as a fallback.
3. Runs `process_engagement.py`.
4. Runs `process_sociodemographic.py`.
   Fails hard on missing/duplicate study IDs in the Qualtrics-derived covariate table.
5. Runs `generate_combined_data.py`.
   Fails hard if any tabular input violates the one-row-per-subject merge contract.
6. Runs `collapse_homer_fir_to_auc.py`.
7. Runs `validate_homer_fir_auc_conversion.py` and fails hard if excluded FIR basis vectors are not represented as `NaN` in the derived AUC table or if the AUC provenance sidecar does not match the current raw FIR export + settings JSON.
8. Runs `mask_homer_auc_between_subject_outliers.py` and writes a separate outlier-masked AUC table plus audit artifacts.
9. Runs `merge_homer3_betas_with_combined_data.R` using `data/tabular/generated_data/homer3_glm_betas_wide_auc_outliers_masked.csv`.
10. Runs `certify_preprocess_merge_integrity.py` and fails hard if merge invariants are violated.

Required inputs for this entry-point:
- `demographic/combined_engagement_data.csv`
- `../Assessment/pretask_assessment.csv`, `../Assessment/posttask_assessment.csv`,
  and `../Assessment/Recall_Assessment_Key.csv` (relative to the repository root)
- `data/config/recall_invalid_questions.json` and `data/config/recall_question_aliases.json`
- `qualtrics/final_SF_demographic_data.csv`
- the raw FIR CSV pointed to by `HOMER_RAW_FIR_CSV` in `pipeline_preprocess_merge.sh`

Required files are checked before results cleanup. A pre-existing
`data/tabular/generated_data/recall_assessment_score_diffs.csv` is not required.
Regeneration preserves the scoring rules; if raw responses, the answer key or the
recall manifests change, the resulting scores and downstream analyses can change.

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

- R: `analyze_format_content_lmm_channelwise.R`

Order / dependencies:
- Assumes Pipeline C0 prerequisites are satisfied.
- Preferred input is the pre-merged table `data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv`
  so covariates and beta columns are available in one file.

Model (per channel × chromophore):
- Omnibus LMM: `beta ~ format_c * content_c + age + education_years + (1 | subject_id)`
- Coding: `format_c = -0.5 (Short), +0.5 (Long)`; `content_c = -0.5 (Entertainment), +0.5 (Education)`
- `age` and `education_years` are required, must be numeric, and must be complete after subject exclusions or the script fails hard.
- The current omnibus covariate adjustment includes `age` and the study-codebook `education_years` proxy; `sfv_daily_duration` remains deferred until its missingness is resolved upstream.
- For numerical conditioning, the R script fits the neural models on a fixed internal response scale (`beta * 1e6`) and back-transforms reported estimates, CIs, and post-hoc mean differences into the original beta units before writing outputs.
- Output tables include a boolean `converged` column based on captured mixed-model convergence warnings so flagged fits remain auditable without being silently dropped.
- The R tidy main-effects table reports estimate, Kenward-Roger-consistent standard error, Kenward-Roger denominator df, signed t-statistic, 95% CI, uncorrected p-value, and BH-FDR q-value for each channel x chromophore x effect row. For the 1-df omnibus terms, the reported SE is derived from the same Kenward-Roger F statistic used for the signed t value, so `estimate / se` reconstructs `t`.

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
- Every inferential ROI must contain exactly 3 channels; other ROI sizes fail hard so the completeness denominator cannot change silently.
- Channel IDs must match Homer naming (example: `S01_D01`).
- Script fails fast on malformed JSON, overlapping channel assignments, or ROI channels absent from the data.

ROI beta construction:
- For each `subject × ROI × chrom × condition`, ROI beta is the arithmetic mean
  across available channels only when at least 2 of the ROI's 3 channels are non-missing.
- If fewer than 2 channels are available for any condition, that condition-level ROI beta is missing; the participant is then excluded from that ROI × chromophore model because all 4 conditions are required.
- The 2-of-3 cutoff is a study-specific conservative completeness decision informed by published fNIRS good-channel inclusion precedents, not a universal threshold.
- In the derived FIR-to-AUC beta table, pruned channels are encoded as `NaN` and are not imputed.

Model / inference:
- Omnibus LMM (per ROI × chrom): `beta ~ format_c * content_c + age + education_years + (1 | subject_id)`
- Same coding and interaction-gated post-hoc policy as C2.
- `age` and `education_years` are required, must be numeric, and must be complete after subject exclusions or the script fails hard.
- The current omnibus covariate adjustment includes `age` and the study-codebook `education_years` proxy; `sfv_daily_duration` remains deferred until its missingness is resolved upstream.
- BH-FDR is applied separately per chromophore and per effect across ROIs.
- For numerical conditioning, the R script fits the neural models on a fixed internal response scale (`beta * 1e6`) and back-transforms reported estimates, CIs, and post-hoc mean differences into the original beta units before writing outputs.
- Output tables include a boolean `converged` column based on captured mixed-model convergence warnings so flagged fits remain auditable without being silently dropped.
- The R tidy main-effects table reports estimate, Kenward-Roger-consistent standard error, Kenward-Roger denominator df, signed t-statistic, 95% CI, uncorrected p-value, and BH-FDR q-value for each ROI x chromophore x effect row. For the 1-df omnibus terms, the reported SE is derived from the same Kenward-Roger F statistic used for the signed t value, so `estimate / se` reconstructs `t`.

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
  (must contain `subject_id`, `age`, `education_years`, and:
  `diff_short_form_education`, `diff_short_form_entertainment`,
  `diff_long_form_education`, `diff_long_form_entertainment`)

Model:
- Omnibus LMM: `retention_diff ~ length_c * content_c + age + education_years + (1 | subject_id)`
- Coding: `length_c = -0.5 (Short), +0.5 (Long)`; `content_c = -0.5 (Entertainment), +0.5 (Education)`
- `age` and `education_years` are required, must be numeric, and must be complete after subject exclusions or the script fails hard.
- The current omnibus covariate adjustment includes `age` and the study-codebook `education_years` proxy; `sfv_daily_duration` remains deferred until its missingness is resolved upstream.

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

### C3b) Retention sensitivity to equal question counts

`analyze_retention_sensitivity.R` reproduces the exhaustive six-question sensitivity
check. The valid counts are 6/7/6/6 for Short Education, Short Entertainment, Long
Education, and Long Entertainment. It fits the original baseline and all seven
ways to retain six Short Entertainment questions, omitting the same question for
every participant in both pre- and post-task assessments. The other conditions,
covariates, participant exclusions, and complete-case cohort stay fixed.

Run from the repository root after recall scoring and the merge pipeline:

```bash
Rscript analyze_retention_sensitivity.R
# To preserve an earlier run, select a new output directory:
Rscript analyze_retention_sensitivity.R --out_dir /tmp/retention_sensitivity_review
Rscript tests/validate_retention_sensitivity_r.R
```

The script reads the scored `demographic/recall_assessment_audit_pre.csv` and
`recall_assessment_audit_post.csv`; it does not regrade raw responses. It validates
the invalid-question manifest, binary scores, shared item coverage, paired target
question IDs/answer keys, and agreement of reconstructed scores with observed
merged outcomes (absolute tolerance `1e-12`). Missing valid-item scores are errors;
valid zeros remain valid. Changed item counts require explicit design review.

The existing retention R script supplies the complete model/reporting path:
REML with age and education covariates, Satterthwaite tests, unadjusted 95% Wald
CIs, Holm over three effects **within each fit**, and the existing interaction-gated
posthoc procedure. Nonconvergence or a rank-deficient design stops the sensitivity
run; singularity is reported. Primary data/results are not overwritten.

Outputs default to a new or empty `data/results/retention_sensitivity/` directory:

- `main_effects.csv`: all 24 effect rows; the adjusted column is explicitly named
  `p_holm` (the primary script's legacy `p_fdr` column also contains Holm values).
- `summary.csv`: baseline results and the seven subsets' coefficient/p-value
  ranges, direction/significance agreement counts, and singular-fit counts.
- `omission_plan.csv`, `subject_scores.csv`, `diagnostics.csv`, and `posthoc.csv`:
  retained questions/counts, participant outcomes/covariates, model diagnostics,
  and any gated contrasts.
- `scenarios/`: each fit's exact input and original-format main/posthoc tables.
- `metadata.json` and `session_info.txt`: settings, cohort counts, input/code hashes,
  interpretation limits, and R/package versions.

This is an **exhaustive subset sensitivity analysis, not a permutation test**.
Significance counts describe seven overlapping analyses, not independent
replications or a probability of robustness. The average of all seven reduced
scores equals the baseline score by construction. Equal item counts cannot remove
missing-video confounding, equalize item difficulty, or establish generalization
over stimuli (Judd et al., 2012; Winkler et al., 2014; see `CITATIONS.md`).

### C4) Engagement within-subject inference: Length, Content, and Length×Content

- R: `analyze_engagement_format_content_lmm.R`

Input:
- `data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv`
  (must contain `subject_id`, `age`, `education_years`, and:
  `sf_education_engagement`, `sf_entertainment_engagement`,
  `lf_education_engagement`, `lf_entertainment_engagement`)

Model:
- Omnibus LMM: `engagement ~ length_c * content_c + age + education_years + (1 | subject_id)`
- Coding: `length_c = -0.5 (Short), +0.5 (Long)`; `content_c = -0.5 (Entertainment), +0.5 (Education)`
- `age` and `education_years` are required, must be numeric, and must be complete after subject exclusions or the script fails hard.
- The current omnibus covariate adjustment includes `age` and the study-codebook `education_years` proxy; `sfv_daily_duration` remains deferred until its missingness is resolved upstream.

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

### C4c) Exploratory pooled-mean neural-behavior correlations

- R: `analyze_pooled_mean_correlations.R`

This is the selected pooled neural-behavior follow-up for the manuscript. Its target-selection, missingness, pooling, and BH-FDR rules remain unchanged.

Purpose:
- Select channel and ROI targets with FDR-significant `format` or `content` main effects from the tidy LMM outputs.
- Explicitly exclude the retired `M_DMPFC` and `M_VMPFC` ROI targets, including when they remain in a stale ROI LMM result table.
- Gate eligible pooled follow-up tests by the selected main effect:
  - `format` targets test `short` and `long`
  - `content` targets test `education` and `entertainment`
- Do not carry pure `interaction` hits into this pooled-main-effect follow-up.
- Reconstruct subject-level pooled neural means for:
  - `short`
  - `long`
  - `education`
  - `entertainment`
- Reconstruct the matching pooled engagement and retention means for the same four pools.
- Correlate matched neural and behavioral pools:
  - neural `short` vs behavioral `short`
  - neural `long` vs behavioral `long`
  - neural `education` vs behavioral `education`
  - neural `entertainment` vs behavioral `entertainment`
- Keep the workflow exploratory because the same dataset is used for target selection and pooled follow-up testing.
- Clear `data/results/pooled_mean_correlations/` before each run so stale CSVs and PNGs do not persist.

Inputs:
- `data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv`
- `data/config/roi_definition.json`
- `data/results/format_content_lmm_main_effects_tidy_r.csv`
- `data/results/format_content_lmm_roi_main_effects_tidy_r.csv`
- `data/config/excluded_subjects.json`

Pool construction:
- `short = mean(SF_Edu, SF_Ent)` when both cells are present
- `long = mean(LF_Edu, LF_Ent)` when both cells are present
- `education = mean(SF_Edu, LF_Edu)` when both cells are present
- `entertainment = mean(SF_Ent, LF_Ent)` when both cells are present
- Each selected ROI must contain exactly three configured channels, all three must be present in the beta input, and a condition-level ROI mean requires at least 2 of those 3 channels to be non-missing.
- A participant contributes pooled rows for an ROI only when the 2-of-3 rule is satisfied in all four conditions.

Missingness and quality policy:
- Channel `0` and `NA` beta values are treated as pruned/missing observations.
- A pooled mean requires both constituent condition cells for that subject.
- ROI condition means may use two or three good channels; fewer than two excludes that participant from the ROI entirely.
- No imputation is performed.

Correlation outputs:
- `association_estimate`: Pearson `r`
- `r_squared`: coefficient of determination from `lm(neural_value ~ behavior_value)`
- `p_unc`, `ci95_low`, `ci95_high`
- `slope`, `intercept`
- BH-FDR is applied:
  - across all tested neural targets within each `behavior_domain x pool_name` family
- Figure emission under `figure_policy = significant_only` is based on uncorrected `p_unc < alpha`.

Example:

```bash
Rscript analyze_pooled_mean_correlations.R \
  --input_csv data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv \
  --roi_json data/config/roi_definition.json \
  --channel_results_csv data/results/format_content_lmm_main_effects_tidy_r.csv \
  --roi_results_csv data/results/format_content_lmm_roi_main_effects_tidy_r.csv \
  --exclude_subjects_json data/config/excluded_subjects.json \
  --out_dir data/results/pooled_mean_correlations
```

Outputs:
- `data/results/pooled_mean_correlations/selected_pooled_mean_targets_r.csv`
- `data/results/pooled_mean_correlations/subject_level_pooled_mean_pairs_r.csv`
- `data/results/pooled_mean_correlations/pooled_mean_correlations_r.csv`
- `data/results/pooled_mean_correlations/figures/`

Notes:
- The script retains only FDR-significant `format` and `content` hits as pooled follow-up targets and gates each target to the matching pool family.
- The results CSV stores pooled Pearson correlations despite the legacy filename `pooled_mean_correlations_r.csv`.

Validation:

```bash
Rscript tests/validate_pooled_mean_correlations_r.R
```

### C4d) Standalone pairwise behavioral correlations

- R: `analyze_behavior_pairwise_correlations.R`

Purpose:
- Screen pairwise associations among the explicitly declared behavioral variables in the same merged CSV used by the primary LMMs and `analyze_pooled_mean_correlations.R`.
- Keep the analysis exploratory and apply one global Benjamini-Hochberg FDR correction across all tested behavioral pairs.
- Exclude recruitment order and subject-ID-derived proxies from the behavioral correlation family.
- Generate one lower-triangle correlation matrix with Pearson `r` and BH-FDR `q` values in each tested cell.
- Clear `data/results/behavior_pairwise_correlations/` before each run so stale CSVs and PNGs do not persist.

Input:
- `data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv`
- `data/config/behavior_pairwise_correlation_plan.json`
- `data/config/variable_figure_names.json` for lower-triangle matrix axis labels

Behavioral variable set under the default plan:
- `sf_education_engagement`
- `sf_entertainment_engagement`
- `lf_entertainment_engagement`
- `lf_education_engagement`
- `diff_short_form_education`
- `diff_short_form_entertainment`
- `diff_long_form_education`
- `diff_long_form_entertainment`
- `age`
- `education_years` (exploratory approximate-years proxy from the study codebook)
- `sfv_frequency`
- `sfv_daily_duration`
- `asrs_total`
- `yang_pu_total`
- `yang_mot_total`
- `phq_total`
- `gad_total`

Method and missingness policy:
- Pearson correlation for every tested pair, with Fisher-z confidence intervals and global BH-FDR q-values.
- `education_years` is treated as an exploratory approximate-years proxy, not a directly measured continuous education-duration variable.
- `sfv_frequency` and `sfv_daily_duration` are ordinal 0-3 codes but are intentionally treated as numeric, equally spaced scores in this Pearson diagnostic screen.
- Pairwise complete cases only for each variable pair.
- No imputation is performed.

Example:

```bash
Rscript analyze_behavior_pairwise_correlations.R \
  --input_csv data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv \
  --analysis_plan_json data/config/behavior_pairwise_correlation_plan.json \
  --variable_figure_names_json data/config/variable_figure_names.json \
  --exclude_subjects_json data/config/excluded_subjects.json \
  --out_dir data/results/behavior_pairwise_correlations
```

Outputs:
- `data/results/behavior_pairwise_correlations/behavior_pairwise_correlations_r.csv`
- `data/results/behavior_pairwise_correlations/behavior_pairwise_correlations_fdr_r.csv`
- `data/results/behavior_pairwise_correlations/figures/behavior_pairwise_correlation_lower_triangle.png`
- `data/results/behavior_pairwise_correlations/figures/behavior_pairwise_correlation_lower_triangle.pdf`

Output ordering:
- rows are sorted by ascending `p_fdr`, then ascending `p_unc`, then descending absolute `pearson_r`

Validation:

```bash
Rscript tests/validate_behavior_pairwise_correlations_r.R
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

## “What does each file do?” (quick reference)

### Tabular / survey / engagement

- `process_sociodemographic.py`: Qualtrics demographics + psych scales → numeric covariates + audits (`data/tabular/generated_data/` + `covariate_outputs/`)
- `demographic/combine_engagement.py`: raw engagement CSVs (external `../../Engagement`) → `demographic/combined_engagement_data.csv` (+ runs basic statsmodels analyses)
- `process_engagement.py`: per-subject engagement condition means → `data/tabular/generated_data/engagement_data_processed.csv`
- `demographic/process_recall_assessment.py`: grade pre/post recall (external `../../Assessment`) → diffs CSV + detailed audit CSVs; uses `data/config/recall_invalid_questions.json` to exclude known invalid assessment items from the normalized retention denominator and `data/config/recall_question_aliases.json` for explicit pre/post Qualtrics ID aliases
- `generate_combined_data.py`: merge engagement + sociodemographics + recall diffs → `data/tabular/generated_data/combined_sfv_data.csv`
- `data/tabular/homer3_glm_betas_wide_fir_pca.csv`: externally produced production Homer3 FIR basis-weight table (wide table; copied into this repo; contains `0` and `NaN` as stand-ins for pruned channels)
- `data/tabular/generated_data/homer3_glm_betas_wide_auc.csv`: locally derived single-beta table produced by `collapse_homer_fir_to_auc.py` from the FIR basis weights
- `data/tabular/generated_data/homer3_glm_betas_wide_auc.provenance.json`: sidecar provenance record for the raw derived AUC table
- `data/tabular/generated_data/homer3_glm_betas_wide_auc_outliers_masked.csv`: between-subject outlier-masked AUC table consumed by merge/modeling
- `data/results/homer_auc_outlier_audit.csv`: row-level audit of censored subject-column AUC cells
- `data/results/homer_auc_outlier_summary.json`: machine-readable summary of between-subject AUC screening
- `validate_homer_fir_auc_conversion.py`: hard-fail lint that checks excluded FIR basis vectors map to `NaN` in the derived AUC CSV and that the provenance sidecar matches the current raw FIR export + settings JSON
- `data/config/excluded_subjects.json`: central participant-exclusion manifest consumed by inferential R analyses
- `analyze_behavior_pairwise_correlations.R`: standalone Pearson screen across declared behavioral-variable pairs in the merged SFV dataset, with global BH-FDR and a lower-triangle matrix figure
- `plot_fir_betas_subjects.py`: plots selected-subject FIR betas for one condition with HbO/HbR overlaid (streaming/selective read; top-of-file config)
- `plot_significant_beta_value_distribution.R`: plots simple beta-value point distributions for the FDR-significant channelwise and ROI LMM hits, with one audit CSV covering every plotted row
- `plot_behavior_score_distributions.R`: plots engagement and recall/retention improvement score distributions by condition, plus content and retention-length marginal means from the final merged Homer3 + SFV dataset, applying the shared participant-exclusion manifest
- `audit_check.py`: consistency checks for the recall assessment audit CSVs

### Misc / exploratory

- `analyze_format_content_lmm_channelwise.R`: within-subject channelwise LMM for Format/Content/Interaction + interaction-gated post-hoc contrasts
- `analyze_format_content_lmm_roi.R`: within-subject ROI-wise LMM using `data/config/roi_definition.json` + interaction-gated post-hoc contrasts
- `analyze_retention_format_content_lmm.R`: within-subject retention LMM for Length/Content/Interaction + interaction-gated post-hoc contrasts
- `analyze_engagement_format_content_lmm.R`: within-subject engagement LMM for Length/Content/Interaction + interaction-gated post-hoc contrasts
