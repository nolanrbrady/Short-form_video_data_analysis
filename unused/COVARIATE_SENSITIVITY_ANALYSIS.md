# Format x Content LMM Covariate Sensitivity Sweep and Implementation Decision

The analyses were run in isolated temporary workspaces; no production result was overwritten during sensitivity testing.

## Purpose

Exploratory sensitivity analysis comparing retention, engagement, ROI, and channel findings across four additive subject-level covariate specifications:

1. No covariates
2. Baseline at the time of the sweep: `age`
3. `age + education_years`
4. `education_years`

The copied temp scripts preserved the existing mixed-model framework (`lmerTest::lmer`, Kenward-Roger inference and BH-FDR for neural models, Satterthwaite inference and Holm correction for behavioral models). This follows the repository analysis citations for random-effects models (Laird & Ware, 1982), `lme4`/`lmerTest` mixed-model inference (Bates et al., 2015; Kuznetsova et al., 2017), Kenward-Roger inference (Kenward & Roger, 1997; Halekoh & Højsgaard, 2014), Holm correction (Holm, 1979), and BH-FDR correction (Benjamini & Hochberg, 1995).

## Implementation Decision

Following the sensitivity sweep, the production channel, ROI, retention, and engagement omnibus models were updated to include both additive covariates: `age + education_years`. Analyses without an existing covariate-adjustment model, including the correlational and condition-coded post-hoc models, remain unchanged. `education_years` is the documented study-codebook proxy rather than a directly measured continuous education-duration variable.

## Run Details

- Original repository outputs were not used as output targets.
- Scripted subject exclusions removed 3 listed/matched subjects, leaving 48 subjects for all behavioral runs and the same per-unit neural complete-case cohorts across specifications.
- `age` and `education_years` were complete in the merged input.
- Education-years distribution in the merged input: 12 years = 27, 14 years = 3, 16 years = 15, 18 years = 6.
- Significance threshold: `p_fdr < 0.05`.
- Significance sets were compared by `{analysis unit, chromophore, effect}`.

## Behavioral Results

All four specifications produced the same behavioral conclusions. Retention format (Holm-adjusted p = 0.00956) and content (Holm-adjusted p = 0.00545) were significant; the retention interaction was not. Engagement content was significant (Holm-adjusted p = 3.55e-15); engagement format and interaction were not. There were no uncorrected or adjusted significance changes, and the tested within-subject effect estimates were unchanged.

## Significant ROI Results

All four model specifications produced the same four FDR-significant ROI rows:

| ROI | Chrom | Effect | No covariates q | Age baseline q | Age + Education q | Education q |
|---|---|---|---:|---:|---:|---:|
| `R_DLPFC` | HbR | format | 0.0282154562 | 0.0282154562 | 0.0282154562 | 0.0282154562 |
| `L_DLPFC` | HbO | format | 0.0412983632 | 0.0412983650 | 0.0412983663 | 0.0412983666 |
| `M_DMPFC` | HbO | format | 0.0412983632 | 0.0412983650 | 0.0412983663 | 0.0412983666 |
| `L_DMPFC` | HbO | format | 0.0425809840 | 0.0425809828 | 0.0425809828 | 0.0425809840 |

No ROI interaction or content effects were FDR-significant in any run.

## Significant Channel Results

All four model specifications produced the same two FDR-significant channel rows:

| Channel | Chrom | Effect | No covariates q | Age baseline q | Age + Education q | Education q |
|---|---|---|---:|---:|---:|---:|
| `S04_D02` | HbR | interaction | 0.0248788486 | 0.0248788432 | 0.0248788436 | 0.0248788436 |
| `S07_D07` | HbR | format | 0.0405758039 | 0.0405758039 | 0.0405758039 | 0.0405758027 |

No channel content effects were FDR-significant in any run.

## Set Differences

Relative to baseline:

- `no covariates`: 0 added or dropped behavioral, ROI, or channel rows.
- `age + education_years`: 0 added ROI rows, 0 dropped ROI rows; 0 added channel rows, 0 dropped channel rows.
- `education_years`: 0 added ROI rows, 0 dropped ROI rows; 0 added channel rows, 0 dropped channel rows.

The fixed-effect estimates for tested format/content/interaction terms were unchanged across variants in the output tables. Kenward-Roger SE/p/FDR values shifted only at very small numerical levels:

| Unit | Variant vs Baseline | Max absolute uncorrected-p delta | Max absolute FDR-q delta |
|---|---|---:|---:|
| ROI | no covariates | 2.33e-08 | 5.82e-08 |
| Channel | no covariates | 4.01e-06 | 1.07e-05 |
| ROI | age + education_years | 2.35e-08 | 2.14e-08 |
| Channel | age + education_years | 2.35e-06 | 2.61e-06 |
| ROI | education_years | 2.26e-08 | 2.05e-08 |
| Channel | education_years | 3.75e-06 | 6.82e-06 |

The temporary sweep also verified that isolated age-only runs reproduced the then-current saved production outputs exactly before the implementation decision was applied.
