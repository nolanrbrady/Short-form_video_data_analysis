## fNIRS Pipeline Review: Items to Address + Supporting Rationale & Literature
- Generated: 12/17/2025

This document summarizes recommended changes to the current fNIRS pipeline:
- Subject-level preprocessing + GLM (`fnirs_analysis.py`)
- Group/second-level inference (`group_analysis_anova.py`, `group_analysis_lme.py`)

It is written to support methodological justification in manuscripts/OSF preregistrations.

---

## Key Study Context (assumed)
- Block design: **120 s task blocks** + **30 s rest**, repeated over ~45 minutes
- **4 conditions**, each presented **4 times**, within-subjects (2×2 structure: Format × Content)
- SNIRF triggers are **onset-only** (no durations stored in file)

---

## Summary of Required Changes (high priority)

### 1) Explicitly model block durations (onset-only triggers)
**What to change**
- Ensure each task annotation is modeled as a **120-second duration** (boxcar) rather than a 0-duration impulse.
- Avoid “double specifying” duration:
  - Either set durations on `raw.annotations` OR pass `stim_dur=120` to design-matrix construction, but not both.
- If you set durations on annotations, apply it **only** to task conditions (not “bad”, “edge”, etc.).

**Why**
- In GLM-based analyses, the event duration defines whether the regressor represents a sustained block vs a momentary impulse. With onset-only triggers, leaving duration at 0 effectively changes the experimental model and can bias amplitude estimates and interpretation.

**Supporting literature**
- Ye, J. C., Tak, S., Jang, K. E., Jung, J., & Jang, J. (2009). *NIRS-SPM: Statistical parametric mapping for near-infrared spectroscopy.* **NeuroImage**.
- Pinti, P., et al. (2020). *A review on the use of the general linear model in functional near-infrared spectroscopy and recommendations for best practice.* (widely cited GLM best-practices review for fNIRS).

---

### 2) Avoid redundant drift removal (preprocessing high-pass + GLM drift terms)
**What to change**
- For slow block designs like 120s/30s, prefer:
  - **Low-pass filtering** in preprocessing (e.g., ~0.2–0.3 Hz), and
  - Handle slow drift within the **GLM drift model** (e.g., cosine basis / high-pass parameter in the design matrix),
  - **Avoid** an extremely low preprocessing FIR high-pass (e.g., ~0.003 Hz) if you are also using GLM drift regressors.

**Why**
- Applying both an ultra-low preprocessing high-pass filter and a GLM drift model can be redundant and may distort low-frequency components relevant to long blocks.
- Very low FIR cutoffs can imply extremely long filters and edge effects; relying on the GLM drift model is usually cleaner for slow paradigms when the design matrix already includes drift regressors.

**Supporting literature**
- Brigadoi, S., & Cooper, R. J. (2015). (review/practical guidance on preprocessing choices and noise sources in fNIRS; commonly cited in fNIRS preprocessing discussions).
- Pinti, P., et al. (2020). (GLM best-practices and careful handling of filtering/drift).
- Ye et al. (2009). (GLM framework and drift/noise handling in fNIRS).

---

### 3) Make bad-channel exclusion explicit at the GLM fit step (SCI)
**What to change**
- After computing SCI and marking channels as bad, explicitly exclude them right before GLM fitting (e.g., `exclude="bads"`).

**Why**
- Channel-quality metrics only improve inference if those channels are actually removed/excluded from the model fit. Relying on implicit downstream behavior can be fragile when objects are copied/transformed (OD → Hb, channel subsetting, etc.).

**Supporting literature**
- Pollonini, L., et al. (2014). *Scalp Coupling Index: a quality metric for near-infrared spectroscopy.* **Neurophotonics**.
- Pinti, P., et al. (2020). (emphasizes QC and robust preprocessing/reporting).

---

### 4) If short-separation channels exist, use them as nuisance regressors (or justify omission)
**What to change**
- If your montage includes short-separation channels:
  - Add short-channel regressors (e.g., mean short HbO/HbR or PCA of short channels) to the design matrix for long-channel GLM.
- If you do not have short channels, document alternative approaches (global signal methods, careful filtering, physiological recording regressors if available).

**Why**
- A central validity challenge in fNIRS is contamination by systemic/superficial physiology. Short-separation regression is widely recommended to improve specificity of cortical effects.

**Supporting literature**
- Gagnon, L., et al. (2014). (short-separation regression / removing superficial/systemic components in fNIRS; canonical reference).
- Yücel, M. A., et al. (2021). *Best practices / guidelines for fNIRS acquisition, analysis, and reporting.* **Neurophotonics**. (best-practice guidance includes attention to systemic physiology and short-separation use where possible).

---

### 5) Second-level LME: enforce complete-case within-subject data per channel (critical correctness fix)
**What to change**
- In `group_analysis_lme.py`, for each channel:
  - Keep only subjects who have **all four conditions** present for that channel before fitting the model and computing contrasts.

**Why**
- Without complete-case enforcement, treatment coding + missing condition levels can produce invalid or misleading contrasts (and the analysis is no longer a clean within-subject comparison for that channel).
- Your docstring claims complete-case filtering is done, but the current code does not enforce it.

**Supporting literature**
- Laird, N. M., & Ware, J. H. (1982). *Random-effects models for longitudinal data.* **Biometrics**. (foundational mixed-effects framework for repeated measures).
- Pinheiro, J. C., & Bates, D. M. (2000). *Mixed-Effects Models in S and S-PLUS.* (standard reference text for practice/assumptions).

---

## Additional Improvements (recommended, not strictly required)

### A) Prevent duplicate aggregation on re-runs
**What to change**
- When appending to an “ALL runs” CSV, deduplicate rows (e.g., by subject/run/channel/condition/chroma).

**Why**
- Re-running pipelines can silently duplicate rows, contaminating group analyses.

**Supporting literature**
- General reproducible data-analysis best practice; no special fNIRS-specific citation needed.

---

### B) Consider propagating first-level uncertainty into group inference (variance weighting)
**What to consider**
- If first-level outputs include standard errors/variances for each beta, consider:
  - inverse-variance weighting, or
  - hierarchical modeling that accounts for first-level uncertainty.

**Why**
- Treating all betas as equally precise ignores heteroskedasticity across channels and subjects and can reduce efficiency or distort inference.

**Supporting literature**
- Friston, K. J., et al. (2002). (summary-statistics / hierarchical inference logic used broadly in neuroimaging).
- Pinti, P., et al. (2020). (discussion of GLM and inference considerations in fNIRS).

---

### C) Multiplicity control is necessary at the channel level (you’re doing this; just document it clearly)
**What to document**
- You correct across channels using FDR in both ANOVA-style and LME pipelines (good).
- Ensure output messages match the actual correction method used (some prints/comments currently refer to “Holm” even when BH-FDR is used).

**Supporting literature**
- Benjamini, Y., & Hochberg, Y. (1995). *Controlling the false discovery rate: a practical and powerful approach to multiple testing.* **JRSSB**.

---

## Notes on Your Two Second-Level Options

### Option 1: Two-stage 2×2 contrast pipeline (`group_analysis_anova.py`)
- Strength: Directly tests the **three planned within-subject effects** (Format, Content, Interaction).
- Caveat: “Global” averaging across channels can be sensitive to which channels survive QC per subject; consider fixing a common channel set or using ROIs.

### Option 2: Per-channel mixed model (`group_analysis_lme.py`)
- Strength: Flexible and standard for repeated measures, with random subject intercepts.
- Critical fix: Must enforce **complete-case** per channel for valid within-subject contrasts.

---

## Minimal “Change List” (implementation checklist)

- [ ] Set **task** annotation durations to 120 s (onset-only triggers); do not double-specify duration.
- [ ] Use **low-pass only** in preprocessing; rely on GLM drift model for slow drift (avoid ultra-low FIR high-pass when using cosine drift).
- [ ] Explicitly exclude **bad (SCI) channels** before GLM fitting.
- [ ] If short-separation channels exist: add them as **nuisance regressors**.
- [ ] In `group_analysis_lme.py`: enforce **complete-case (all 4 conditions) per channel** before LME + contrasts.
- [ ] Deduplicate aggregated CSV outputs on re-run.

---

## References (starter list)
- Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate. **JRSSB**.
- Friston, K. J., et al. (2002). Classical and Bayesian inference in neuroimaging / hierarchical inference (summary-statistics logic). (widely cited neuroimaging inference reference).
- Gagnon, L., et al. (2014). Short-separation regression for systemic physiology removal in fNIRS. **NeuroImage**.
- Laird, N. M., & Ware, J. H. (1982). Random-effects models for longitudinal data. **Biometrics**.
- Pinheiro, J. C., & Bates, D. M. (2000). **Mixed-Effects Models in S and S-PLUS**.
- Pinti, P., et al. (2020). GLM use in fNIRS and best-practice recommendations. (review).
- Pollonini, L., et al. (2014). Scalp Coupling Index (SCI). **Neurophotonics**.
- Ye, J. C., et al. (2009). NIRS-SPM: GLM framework for fNIRS. **NeuroImage**.
- Yücel, M. A., et al. (2021). Best-practice guidelines for fNIRS. **Neurophotonics**.

---