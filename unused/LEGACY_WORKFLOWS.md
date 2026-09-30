# Historical workflow instructions

These sections were relocated from the root README on 2026-09-30. They are
not current execution instructions or a declaration of current results. Some
referenced files were already absent before this cleanup. Do not run historical
commands against current outputs. Restore original paths from the move manifest
and use an isolated workspace if revisiting an archived workflow.

### C4b) Exploratory pooled long/short ROI/channel correlations with pooled and raw behavioral follow-up

- R: `analyze_correlational_relationships.R`

Purpose:
- Run exploratory post-hoc associations between:
  - pooled long-form and pooled short-form behavioral means and the matching pooled long-form and pooled short-form neural means
  - raw behavioral task-cell values and the matching pooled long-form and pooled short-form neural means
  from:
  - `S04_D02` for both `HbO` and `HbR`
  - `R_DLPFC (HbR)`, `L_DLPFC (HbO)`, `M_DMPFC (HbO)`, `L_DMPFC (HbO)`
- Use the same merged input table and subject-exclusion manifest as the other R analyses.
- Interpret these results as exploratory rather than confirmatory when the channel/ROI target set was selected from this same dataset.
- The pooled long/short rows are the primary outputs.
- The raw condition rows are supplementary localization checks and do not directly test whether the long-form association differs from the short-form association.
- The raw condition rows do not support condition-specific neural statements because the neural side remains pooled within `long` or `short`.

Behavior runs:
- `engagement`
  - `engagement_long = ((lf_education_engagement + lf_entertainment_engagement) / 2)`
  - `engagement_short = ((sf_education_engagement + sf_entertainment_engagement) / 2)`
- `retention`
  - `retention_long = ((diff_long_form_education + diff_long_form_entertainment) / 2)`
  - `retention_short = ((diff_short_form_education + diff_short_form_entertainment) / 2)`
- Supplementary raw-value runs:
  - engagement: `sf_education_engagement`, `sf_entertainment_engagement`, `lf_entertainment_engagement`, `lf_education_engagement`
  - retention: `diff_short_form_education`, `diff_short_form_entertainment`, `diff_long_form_entertainment`, `diff_long_form_education`
- Stored pooled engagement columns such as `long_form_engagement` and `short_form_engagement` are not used because they can differ from the cell-wise arithmetic means implied by the 2x2 design.

Neural target construction:
- Channel targets are collapsed to one pooled `long` mean and one pooled `short` mean per subject and chromophore.
- ROI targets are arithmetic means across available non-missing channels in the ROI for each `subject x chrom x condition`, then collapsed to one pooled `long` mean and one pooled `short` mean per subject and chromophore.
- Supplementary rows reuse those same pooled neural `long` or `short` means rather than reverting to condition-specific neural values.
- ROI channel membership is read from `data/config/roi_definition.json`.

Missingness policy:
- Channel `0` and `NA` beta values are treated as pruned/missing observations for this workflow.
- A subject contributes a pooled `long` row only when both long-form cells for that pooled mean are present.
- A subject contributes a pooled `short` row only when both short-form cells for that pooled mean are present.
- A subject contributes a raw-value row only when that specific behavior cell and its matched pooled neural mean are both present.
- After effect construction, each `behavior_run x neural target` association uses pairwise complete cases only.
- Subjects excluded in `data/config/excluded_subjects.json` are removed before any effect construction.

Association methods:
- Pearson correlation is the primary analysis.
- Spearman correlation is emitted as a sensitivity analysis for bounded behavioral outcomes.
- Output includes both uncorrected `p_unc` and BH-adjusted `p_fdr`.

Multiple testing:
- BH-FDR is applied within each configured `analysis_tier x behavior_run x format_pool x association_method` family.
- Under the current default config, that means:
  - `6` tests per family, because each family contains the same 6 selected neural targets
  - `24` families total:
    - `8` pooled-format families: `2` behavior runs x `2` format pools x `2` association methods
    - `16` raw-value families: `8` raw behavior runs x `2` association methods
- Family membership and behavior-run definitions are declared in `data/config/correlational_analysis_plan.json`.
- Output includes both raw `p_unc` and adjusted `p_fdr` in a single combined CSV.
- The script clears `data/results/correlational_relationships/` before each run so stale CSVs and PNGs do not persist.

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
- `data/results/correlational_relationships/figures/` (scatterplots only for tested rows selected by the current figure policy in `data/config/correlational_analysis_plan.json`; under the current default plan that means FDR-significant Pearson rows)

CSV ordering:
- the combined CSV is ordered by neural target, condition, and predictor

Validation:

```bash
Rscript tests/validate_correlational_relationships_r.R
```

### C4f) Exploratory channel-behavior screening across all channel columns

- Python: `analyze_channel_behavior_relationships.py`

Purpose:
- Screen every channel column matching `S##_D##_Cond##_HbO/HbR` against every non-identifier behavioral variable in the merged CSV.
- Treat `subject_id` and `homer_subject` as identifiers rather than behavioral predictors.
- Keep the analysis exploratory and reproducible, with explicit missingness handling and global FDR control.

Input:
- `data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv`

Missingness and quality policy:
- Channel `0` and `NaN` values are treated as missing/pruned observations, consistent with the repo's Homer beta import note.
- No imputation is performed.
- Tests are computed on pairwise complete cases only for each `behavior x channel` pair.

Association methods:
- Spearman rank correlation for non-binary behavioral variables.
- Point-biserial correlation for binary behavioral variables such as `pd_status`.
- Benjamini-Hochberg FDR is applied across the full family of valid tests in the run.

Example:

```bash
python analyze_channel_behavior_relationships.py \
  --input-csv data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv \
  --out-dir data/results/channel_behavior_relationships
```

Outputs:
- `channel_behavior_pairwise_results.csv`: one row per tested `behavior x channel` pair with method, effect size, `p_uncorrected`, `p_fdr_bh`, and pairwise `n`
- `channel_behavior_top_hits.csv`: top-ranked results after sorting by FDR, then raw p-value
- `channel_behavior_condition_matched_top_hits.csv`: subset of results where the behavioral variable and channel condition refer to the same condition family
- `channel_behavior_behavior_summary.csv`: one summary row per behavioral variable
- `behavior_variable_profile.csv`: method assignment and missingness profile for each included behavioral variable
- `channel_missingness_summary.csv`: per-channel counts for raw zero placeholders, raw NaNs, and usable observations after pruning
- `analysis_metadata.json`: run-level analysis metadata
- `channel_behavior_summary.md`: concise markdown report with top hits, FDR-significant hits, and pruning summaries

Validation:

```bash
python tests/validate_channel_behavior_relationships_py.py
```

### A7) Correlation diagnostics / heatmaps (optional)

- Script: `covariate_correlation_analysis.py`
- Two presets:
  - **`covariates` preset**: correlations on `covariate_outputs/covariates_clean.csv`
  - **`combined` preset**: correlations on `data/tabular/generated_data/combined_sfv_data.csv`
- Each run writes Pearson-only correlation outputs:
  - `covariate_correlations_pearson.csv`, `covariate_correlations_pearson_pvalues.csv`, and `covariate_heatmap_pearson.png` for the `covariates` preset.
  - `covariate_correlation_analysis_pearson.csv`, `covariate_correlation_analysis_pearson_pvalues.csv`, and `covariate_correlation_heatmap_pearson.png` for the `combined` preset.
- Subject exclusions are applied from `data/config/excluded_subjects.json` by default before correlations are computed.
- The input must contain the configured subject ID column (`subject_id` by default) whenever the exclusion manifest is nonempty; the script fails hard rather than inferring exclusions from row order.
- The script fails if more than 48 subjects remain after exclusions (`--max-subjects 48` by default), so a missed exclusion manifest cannot silently enter the correlation analysis.
- Identifier columns (`subject_id`, the configured subject column, and `homer_subject` when present) are removed from the correlation matrix after exclusions.
- To explore recruitment-order artifacts, pass `--include-subject-id-correlation`; this adds the normalized numeric ID as `recruitment_order_proxy` after exclusions. Treat this as an exploratory diagnostic only, because subject ID is interpretable only if it encodes recruitment/order.

Recommended commands:

```bash
# Covariates-only Pearson correlations + p-values
python covariate_correlation_analysis.py \
  --preset covariates \
  --input path/to/covariates_with_subject_id.csv \
  --out-dir covariate_outputs \
  --excluded-subjects-json data/config/excluded_subjects.json

# Combined dataset correlations + heatmap saved alongside the combined dataset
python covariate_correlation_analysis.py \
  --preset combined \
  --out-dir data/tabular/generated_data \
  --excluded-subjects-json data/config/excluded_subjects.json \
  --max-subjects 48 \
  --include-subject-id-correlation
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

### A5e) Plot channel-vs-ROI beta dynamics for confusing findings

- Python: `plot_beta_discrepancy_dynamics.py`
- Inputs:
  - `data/tabular/generated_data/homer3_betas_plus_combined_sfv_data_inner_join.csv`
  - `data/config/roi_definition.json`
  - `data/config/excluded_subjects.json`

What it does:
- Reads the merged wide beta table and compares one selected channel against one selected ROI.
- Defaults to the current discrepancy discussed in this repo:
  - channel: `S04_D02` `HbR`
  - ROI: `L_DMPFC` `HbO`
- Applies the shared subject-exclusion manifest by default so plotted subjects match the inferential analyses.
- Treats `0` and `NaN` beta values as pruned/missing by default (no beta imputation).
- Computes ROI means from the pre-specified member channels in `data/config/roi_definition.json`.
- Produces a composite figure with:
  - a channel panel across the 4 conditions
  - an ROI-mean panel across the 4 conditions
  - an optional ROI-decomposition panel showing the member channels plus the ROI mean
- Supports selectable styles:
  - `raw_means`
  - `means_only`
  - `raw_only`
- Writes a tidy CSV of the plotted raw values and summary means/intervals alongside the figure for auditability.

Default output directory:
- `data/results/beta_discrepancy_plots/`

Command:

```bash
python plot_beta_discrepancy_dynamics.py
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

### fNIRS (MNE/MNE-NIRS)

- `fnirs_analysis/fnirs_analysis.py`: SNIRF preprocessing + first-level GLM per run; writes per-subject outputs under `glm_results/`
- `fnirs_analysis/homer_betas_qc.py`: subject-level QC from imported Homer beta-wide table (`0`/`NaN` treated as excluded channels); writes subject and cohort QC CSVs
- `r_subject_exclusions.R`: shared R helpers that enforce centralized subject exclusions from JSON across inferential scripts
- `fnirs_analysis/combine_glm_output.py`: combine subject GLM CSVs; writes `combined_glm_long*.csv` + `combined_matrices/`
- `fnirs_analysis/group_analysis_anova.py`: group-level two-stage contrast testing + selective channel inference
- `fnirs_analysis/group_analysis_lme.py`: per-channel mixed-effects + FDR + gated post-hocs
- `fnirs_analysis/FNIRS_TODO.md`: methodological notes and a change list for the fNIRS pipeline

## LaTeX manuscript workflow

The manuscript source lives in `latex/main.tex` with bibliography entries in
`latex/references.bib`. Generated LaTeX artifacts are redirected into
`latex/build/` so the source directory stays stable while you edit.

From the repository root:

```bash
make -C latex pdf
make -C latex watch
make -C latex tidy
```

`make -C latex watch` runs `latexmk -pvc`, which continuously rebuilds
`latex/build/main.pdf` whenever `latex/main.tex` or `latex/references.bib`
changes. Open `latex/build/main.pdf` in a PDF viewer that auto-reloads on file
changes to get live preview while editing. `make -C latex tidy` removes the
legacy flat-layout artifacts that were previously written directly into
`latex/`.

---
