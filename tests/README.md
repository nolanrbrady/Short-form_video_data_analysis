# Tests / Validation

This folder contains lightweight validation harnesses that aim to verify scientific/analysis integrity against `ANALYSIS_SPEC.md`.

## Pipeline C (R): synthetic end-to-end validation

Runs `analyze_format_content_lmm_channelwise.R` on a synthetic dataset that:
- Uses Homer-style subject IDs (`sub_0001`) and combined IDs (`0001`)
- Includes both HbO/HbR and multiple channels
- Includes a pruned beta (`0`) to verify pruned-channel missingness + complete-case logic
- Verifies coefficient recovery (within tolerance) for known ground-truth fixed effects
- Verifies post-hoc gating only triggers for interaction-FDR-significant channels
- Verifies duplicate `subject_id` fails hard

Command:

```bash
Rscript tests/validate_pipeline_c_r.R
```

## Retention Pipeline (R): deterministic + synthetic validation

Runs `analyze_retention_format_content_lmm.R` on generated retention datasets and verifies:
- Deterministic analytic ground truth for Length/Content/Interaction contrasts
- End-to-end coefficient recovery for known generating parameters
- Holm-correction correctness across the 3 omnibus effects
- Interaction-gated post-hoc behavior (on/off cases)
- Zero-as-valid retention handling (not treated as missing)
- Complete-case subject drop only for true `NA`
- Fail-hard integrity checks (duplicate IDs, missing required columns, non-numeric retention values)

Command:

```bash
Rscript tests/validate_retention_pipeline_r.R
```

## Engagement Pipeline (R): deterministic + synthetic validation

Runs `analyze_engagement_format_content_lmm.R` on generated engagement datasets and verifies:
- Deterministic analytic ground truth for Length/Content/Interaction contrasts
- End-to-end coefficient recovery for known generating parameters
- Holm-correction correctness across the 3 omnibus effects
- Interaction-gated post-hoc behavior (on/off cases)
- Zero-as-valid engagement handling (not treated as missing)
- Complete-case subject drop only for true `NA`
- Fail-hard integrity checks (duplicate IDs, missing required columns, non-numeric engagement values)

Command:

```bash
Rscript tests/validate_engagement_pipeline_r.R
```
