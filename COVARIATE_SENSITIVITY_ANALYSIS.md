# Format x Content LMM Covariate Sensitivity Sweep

Experiment root: `/tmp/sfv_lmm_covariate_sweep.gFL7wX`

## Purpose

Exploratory sensitivity analysis comparing FDR-significant ROI/channel findings across three additive subject-level covariate specifications:

1. Baseline current scripts: `age`
2. Variant: `age + education_years`
3. Variant: `education_years`

The copied temp scripts preserve the existing mixed-model framework (`lmerTest::lmer`, Kenward-Roger inference, BH-FDR). This follows the repository analysis citations for random-effects models (Laird & Ware, 1982), `lme4`/`lmerTest` mixed-model inference (Bates et al., 2015; Kuznetsova et al., 2017), Kenward-Roger inference (Kenward & Roger, 1997; Halekoh & Højsgaard, 2014), and BH-FDR correction (Benjamini & Hochberg, 1995).

## Run Details

- Source repo copied to: `/tmp/sfv_lmm_covariate_sweep.gFL7wX/work`
- Original repository outputs were not used as output targets.
- Scripted subject exclusions removed 3 listed/matched subjects, leaving 48 subjects for all six runs.
- `age` and `education_years` were complete in the merged input.
- Education-years distribution in the merged input: 12 years = 27, 14 years = 3, 16 years = 15, 18 years = 6.
- Significance threshold: `p_fdr < 0.05`.
- Significance sets were compared by `{ROI/channel, chromophore, effect}`.

## Significant ROI Results

All three model specifications produced the same four FDR-significant ROI rows:

| ROI | Chrom | Effect | Baseline q | Age + Education q | Education q |
|---|---|---|---:|---:|---:|
| `R_DLPFC` | HbR | format | 0.0282154562 | 0.0282154562 | 0.0282154562 |
| `L_DLPFC` | HbO | format | 0.0412983650 | 0.0412983663 | 0.0412983666 |
| `M_DMPFC` | HbO | format | 0.0412983650 | 0.0412983663 | 0.0412983666 |
| `L_DMPFC` | HbO | format | 0.0425809828 | 0.0425809828 | 0.0425809840 |

No ROI interaction or content effects were FDR-significant in any run.

## Significant Channel Results

All three model specifications produced the same two FDR-significant channel rows:

| Channel | Chrom | Effect | Baseline q | Age + Education q | Education q |
|---|---|---|---:|---:|---:|
| `S04_D02` | HbR | interaction | 0.0248788432 | 0.0248788436 | 0.0248788436 |
| `S07_D07` | HbR | format | 0.0405758039 | 0.0405758039 | 0.0405758027 |

No channel content effects were FDR-significant in any run.

## Set Differences

Relative to baseline:

- `age + education_years`: 0 added ROI rows, 0 dropped ROI rows; 0 added channel rows, 0 dropped channel rows.
- `education_years`: 0 added ROI rows, 0 dropped ROI rows; 0 added channel rows, 0 dropped channel rows.

The fixed-effect estimates for tested format/content/interaction terms were unchanged across variants in the output tables. Kenward-Roger SE/p/FDR values shifted only at very small numerical levels:

| Unit | Variant vs Baseline | Max absolute uncorrected-p delta | Max absolute FDR-q delta |
|---|---|---:|---:|
| ROI | age + education_years | 2.35e-08 | 2.14e-08 |
| Channel | age + education_years | 2.35e-06 | 2.61e-06 |
| ROI | education_years | 2.26e-08 | 2.05e-08 |
| Channel | education_years | 3.75e-06 | 6.82e-06 |

## Output Files

- `summary/all_tidy_results_with_significance.csv`
- `summary/fdr_significant_results.csv`
- `summary/fdr_significant_counts.csv`
- `summary/set_differences_vs_baseline.csv`
- Per-run CSVs under `results/baseline`, `results/age_education`, and `results/education_only`.
