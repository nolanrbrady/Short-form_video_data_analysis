## Short-form Video Study — Analysis Repo
Last updated: 01-30-2026
Updated by: Nolan Brady

This repo contains two primary analysis “tracks”:

- **Tabular / behavioral + survey preprocessing** → outputs in `data/tabular/`
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
- **`data/tabular/`**: “analysis-ready” CSVs produced by preprocessing scripts (merged by `subject_id`)
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

2) **Compute per-subject engagement features → `data/tabular/engagement_data_processed.csv`**

- Script: `process_engagement.py`
- Input: `demographic/combined_engagement_data.csv`
- Output: `data/tabular/engagement_data_processed.csv`

Command:

```bash
python process_engagement.py
```

Output columns include:

- Condition means: `lf_education_engagement`, `lf_entertainment_engagement`, `sf_education_engagement`, `sf_entertainment_engagement`
- Marginal means: `long_form_engagement`, `short_form_engagement`, `education_engagement`, `entertainment_engagement`

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

- `data/tabular/recall_assessment_score_diffs.csv`: one row per `subject_id`, columns like `diff_short_form_education`, etc.
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
  - `data/tabular/socio_demographic_data_processed.csv` (includes `subject_id` for merges)
  - `covariate_outputs/covariates_clean.csv` (covariates only; excludes `subject_id`)
  - `covariate_outputs/covariate_missingness.csv`, `covariate_outputs/covariate_column_audit.csv`, `covariate_outputs/sfv_duration_other_audit.csv`

Command:

```bash
python process_sociodemographic.py
```

### A4) Combine tabular sources into one dataset

- Script: `generate_combined_data.py`
- Inputs:
  - `data/tabular/engagement_data_processed.csv`
  - `data/tabular/socio_demographic_data_processed.csv`
  - `data/tabular/recall_assessment_score_diffs.csv`
- Output:
  - `data/tabular/combined_sfv_data.csv`

Command:

```bash
python generate_combined_data.py
```

### A5) Import Homer3 GLM betas (optional; produced externally)

This repo does **not** run Homer3. If you run Homer3 elsewhere, copy the exported wide table into:

- `data/tabular/homer3_glm_betas_wide.csv`

Current format on disk:

- ID column: `Subject` (values like `sub_0001`)
- Feature columns: `S##_D##_Cond##_HbO` / `S##_D##_Cond##_HbR`

**Important missingness note:** this file currently contains both `0` and `NaN` values that are *stand-ins for channels that were pruned during preprocessing*. Do **not** interpret these as “true zero activation”; treat them as missing/pruned channels in downstream modeling.

To merge with `data/tabular/combined_sfv_data.csv`, you will need a shared key:

- Either export Homer betas with a `subject_id` column that matches Qualtrics, **or**
- Create a mapping between `Subject` (e.g., `sub_0001`) and Qualtrics `subject_id` and merge using that.

---

## Pipeline C — Homer3 betas + Format×Content (channelwise LMM)

**Goal:** Merge the externally-produced Homer3 betas table with the combined tabular dataset and test whether
**Format depends on Content (and vice versa)** at the level of **channelwise prefrontal activation**, for **HbO** and **HbR**.

### C0) Prerequisites / required prior work

Inputs required:
- `data/tabular/homer3_glm_betas_wide.csv` (externally produced; see “Homer3 GLM betas export” above)
- `data/tabular/combined_sfv_data.csv` (produced by the tabular preprocessing pipeline)

Requirements / invariants:
- `combined_sfv_data.csv` must contain **exactly one row per subject** (unique `subject_id`).
- `homer3_glm_betas_wide.csv` must contain **exactly one row per subject** (unique `Subject` after normalization).
- Beta columns can include `0`/`NaN` as pruned-channel stand-ins; downstream modeling treats these as missing (do not impute).

Recommended prior step (if you haven’t generated it yet):
- Run `generate_combined_data.py` to (re)build `data/tabular/combined_sfv_data.csv` from the raw sources.

### C1) Inner-join Homer3 betas with combined tabular data

- Python: `merge_homer3_betas_with_combined_data.py`
- R: `merge_homer3_betas_with_combined_data.R`

What it does:
- Normalizes IDs by extracting digits (handles `0017` vs `17`, and `sub_0017`-style IDs).
- Performs an **INNER JOIN** and writes a merged “one row per subject” CSV for inspection / downstream use.

Notes:
- IDs are normalized by extracting digits (handles `0017` vs `17`, and `sub_0017`-style IDs).
- Output is one row per subject containing both demographics/behavior and beta columns.
- Scripts **fail hard** if either input contains duplicate `subject_id` values after normalization (expected exactly one row per subject).

Example (Python):

```bash
python merge_homer3_betas_with_combined_data.py \
  --homer-csv data/tabular/homer3_glm_betas_wide.csv \
  --combined-csv data/tabular/combined_sfv_data.csv \
  --out-csv data/results/homer3_betas_plus_combined_sfv_data_inner_join.csv
```

### C2) Channelwise within-subject inference: Format, Content, and Format×Content

- Python: `analyze_format_content_lmm_channelwise.py`
- R: `analyze_format_content_lmm_channelwise.R`

Order / dependencies:
- Assumes Pipeline C0 prerequisites are satisfied.
- Can be run directly from the two input CSVs (it merges internally); C1 is still recommended to create an auditable merged artifact.

Model (per channel × chromophore):
- LMM: `beta ~ format_c * content_c + (1 | subject_id)`
- Coding: `format_c = -0.5 (Short), +0.5 (Long)`; `content_c = -0.5 (Entertainment), +0.5 (Education)`

Pruned channels / missingness policy:
- Treat **both `0` and `NaN`** in Homer betas as **pruned/missing** (do **not** impute).
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

Minimum sample gating:
- Models are only fit for channel/chrom pairs with at least `min_subjects` complete-case subjects (default: 6).
- Scripts print a warning summary (count + examples) for models skipped due to `min_subjects`.

Example (Python dry-run to validate parsing without fitting models):

```bash
python analyze_format_content_lmm_channelwise.py --dry-run
```

Outputs:
- `data/results/format_content_lmm_main_effects_*.csv`
- `data/results/format_content_lmm_main_effects_tidy_r.csv` (R: spec-compliant tidy main-effects table)
- `data/results/format_content_lmm_posthoc_pairwise_*.csv`

### A6) Correlation diagnostics / heatmaps (optional)

- Script: `covariate_correlation_analysis.py`
- Two presets:
  - **`covariates` preset**: correlations on `covariate_outputs/covariates_clean.csv`
  - **`combined` preset**: correlations on `data/tabular/combined_sfv_data.csv`

Recommended commands:

```bash
# Covariates-only correlations + p-values (if SciPy installed)
python covariate_correlation_analysis.py --preset covariates --out-dir covariate_outputs

# Combined dataset correlations + heatmap saved alongside the combined dataset
python covariate_correlation_analysis.py --preset combined --out-dir data/tabular
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

- `process_sociodemographic.py`: Qualtrics demographics + psych scales → numeric covariates + audits (`data/tabular/` + `covariate_outputs/`)
- `demographic/combine_engagement.py`: raw engagement CSVs (external `../../Engagement`) → `demographic/combined_engagement_data.csv` (+ runs basic statsmodels analyses)
- `process_engagement.py`: per-subject engagement condition means → `data/tabular/engagement_data_processed.csv`
- `demographic/process_recall_assessment.py`: grade pre/post recall (external `../../Assessment`) → diffs CSV + detailed audit CSVs
- `generate_combined_data.py`: merge engagement + sociodemographics + recall diffs → `data/tabular/combined_sfv_data.csv`
- `data/tabular/homer3_glm_betas_wide.csv`: externally produced Homer3 GLM betas (wide table; copied into this repo for downstream modeling; contains `0` and `NaN` as stand-ins for pruned channels)
- `covariate_correlation_analysis.py`: Spearman correlation tables + heatmaps (covariates-only or combined dataset)
- `audit_check.py`: consistency checks for the recall assessment audit CSVs
- `demographic/engagement_stats.py`: helper functions used by `combine_engagement.py` for ANOVA/OLS/mixed model

### fNIRS (MNE/MNE-NIRS)

- `fnirs_analysis/fnirs_analysis.py`: SNIRF preprocessing + first-level GLM per run; writes per-subject outputs under `glm_results/`
- `fnirs_analysis/combine_glm_output.py`: combine subject GLM CSVs; writes `combined_glm_long*.csv` + `combined_matrices/`
- `fnirs_analysis/group_analysis_anova.py`: group-level two-stage contrast testing + selective channel inference
- `fnirs_analysis/group_analysis_lme.py`: per-channel mixed-effects + FDR + gated post-hocs
- `fnirs_analysis/FNIRS_TODO.md`: methodological notes and a change list for the fNIRS pipeline

### Misc / exploratory

- `analyze_sfv_demographics.py`: standalone demographic summary + figures for an input `SFV_demo_data.csv` colocated with the script (older/one-off analysis)
- `data_quality_evaluation.py`: placeholder for a future subject-level fNIRS QC report (see `TODO.md`)
- `TODO.md`: high-level analysis TODO list
