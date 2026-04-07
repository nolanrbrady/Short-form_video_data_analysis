# Preprocessing

- **Scholkmann, F., & Wolf, M. (2013)**
  - **Title:** General equation for the differential pathlength factor of the frontal human head depending on wavelength and age
  - **First Author:** Scholkmann
  - **Year:** 2013
  - **Usage:** General Linear Model / Beer-Lambert Law
  - **Reasoning:** Used as the default for wavelength-dependent partial pathlength factor (PPF) in the Beer-Lambert Law conversion from optical density to haemoglobin concentration when no specific PPF is provided.
  - **Link:** https://doi.org/10.1117/1.JBO.18.10.105004
  - **Source:** `fnirs_analysis/fnirs_analysis.py`

- **Brigadoi, S., & Cooper, R. J. (2015)**
  - **Title:** How short is short? Optimum source–detector distance for short-separation channels in functional near-infrared spectroscopy
  - **First Author:** Brigadoi
  - **Year:** 2015
  - **Usage:** Preprocessing Guidelines
  - **Reasoning:** Referenced for practical guidance on preprocessing choices and identifying noise sources in fNIRS data.
  - **Link:** https://doi.org/10.1117/1.NPh.2.2.025005
  - **Source:** `fnirs_analysis/FNIRS_TODO.md`

- **Gagnon, L., et al. (2014)**
  - **Title:** Further improvement in reducing superficial contamination in NIRS using double short separation measurements
  - **First Author:** Gagnon
  - **Year:** 2014
  - **Usage:** Short-separation Regression
  - **Reasoning:** Cited as the canonical reference for removing superficial and systemic physiological components using short-separation channels.
  - **Link:** https://doi.org/10.1016/j.neuroimage.2013.01.073
  - **Source:** `fnirs_analysis/FNIRS_TODO.md`

- **Yücel, M. A., et al. (2021)**
  - **Title:** Best practices for fNIRS publications
  - **First Author:** Yücel
  - **Year:** 2021
  - **Usage:** Best practices for analysis/reporting, including chromophore transparency and explicit quality handling
  - **Reasoning:** Provides peer-reviewed guidance for transparent fNIRS reporting. Used here to justify explicitly plotting both HbO and HbR and preserving pruned-channel placeholders (`0`/`NaN`) as missing values in subject-level FIR visualizations (no imputation).
  - **Link:** https://doi.org/10.1117/1.NPh.8.1.012101
  - **Source:** `fnirs_analysis/FNIRS_TODO.md`, `homer_fir.py`, `plot_fir_betas_subjects.py`, `plot_beta_discrepancy_dynamics.py`, `README.md`

- **Ye, J. C., Tak, S., Jang, K. E., Jung, J., & Jang, J. (2009)**
  - **Title:** NIRS-SPM: Statistical parametric mapping for near-infrared spectroscopy
  - **First Author:** Ye
  - **Year:** 2009
  - **Usage:** GLM basis-function HRF representation for Homer `idxBasis=1` Gaussian FIR exports
  - **Reasoning:** Peer-reviewed methodological reference for representing the fNIRS HRF as a weighted sum of temporal basis functions in a GLM. Used here to justify reconstructing the latent HRF from Homer-exported Gaussian basis weights before plotting or summarizing, rather than treating the FIR basis coefficients themselves as the hemodynamic response.
  - **Link:** https://doi.org/10.1016/j.neuroimage.2008.08.036
  - **Source:** `homer_fir.py`, `collapse_homer_fir_to_auc.py`, `plot_fir_betas_subjects.py`, `README.md`, `ANALYSIS_SPEC.md`

- **Ning, L., et al. (2024)**
  - **Title:** fNIRS dataset during complex scene analysis with temporal HRF and physiological features
  - **First Author:** Ning
  - **Year:** 2024
  - **Usage:** Time-window AUC summary of reconstructed hemodynamic responses
  - **Reasoning:** Provides peer-reviewed precedent for summarizing fNIRS hemodynamic trajectories with explicit, fixed-window HRF area-under-the-curve features. Used here to justify the choice of a pre-specified task-window AUC summary once the latent HRF has been reconstructed from the Homer FIR basis weights.
  - **Link:** https://doi.org/10.3389/fnhum.2024.1418592
  - **Source:** `homer_fir.py`, `collapse_homer_fir_to_auc.py`, `README.md`, `ANALYSIS_SPEC.md`

# Quality Control

- **Pollonini, L., et al. (2016)**
  - **Title:** PHOEBE: a method for real time mapping of optodes-scalp coupling in functional near-infrared spectroscopy
  - **First Author:** Pollonini
  - **Year:** 2016
  - **Usage:** Scalp Coupling Index (SCI); windowed QC cadence
  - **Reasoning:** Foundational reference for SCI as a coupling-quality metric and supports 10s windowed QC cadence used by PHOEBE for time-local quality assessment.
  - **Link:** https://doi.org/10.1364/BOE.7.005104
  - **Alt Link:** https://opg.optica.org/abstract.cfm?uri=boe-7-12-5104
  - **Alt Link:** https://pubmed.ncbi.nlm.nih.gov/28018728/
  - **Source:** `fnirs_analysis/fnirs_analysis.py`, `fnirs_analysis/qc_check.py`, `fnirs_analysis/fnirs_preprocess_justifications.md`

- **Hocke, Oni, Duszynski, Corrigan, Frederick, Dunn (2018)**
  - **Title:** Automated Processing of fNIRS Data—A Visual Guide to the Pitfalls and Consequences
  - **First Author:** Hocke
  - **Year:** 2018
  - **Usage:** SCI Threshold
  - **Reasoning:** Uses a hard-coded SCI threshold of ≥ 0.7 for channel acceptance, supporting the run-level exclusion criteria.
  - **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC6428450/
  - **Source:** `fnirs_analysis/qc_check.py`, `fnirs_analysis/fnirs_preprocess_justifications.md`

- **Hernandez & Pollonini (2020)**
  - **Title:** NIRSplot: A Tool for Quality Assessment of fNIRS Scans
  - **First Author:** Hernandez
  - **Year:** 2020
  - **Usage:** QC Assessment (SCI + PSP)
  - **Reasoning:** Establishes SCI as a coupling-quality metric and uses peak spectral power (PSP) for time-interval QC in the QT-NIRS/NIRSplot framework.
  - **Link:** https://research.birmingham.ac.uk/en/publications/nirsplot-a-tool-for-quality-assessment-of-fnirs-scans/
  - **Alt Link:** https://scholars.houstonmethodist.org/en/publications/nirsplot-a-tool-for-quality-assessment-of-fnirs-scans/
  - **Alt Link:** https://par.nsf.gov/servlets/purl/10209495
  - **Source:** `fnirs_analysis/qc_check.py`, `fnirs_analysis/fnirs_preprocess_justifications.md`

- **Meier et al. (2025)**
  - **Title:** The effects of protocol factors and participant characteristics on functional near-infrared spectroscopy data quality after stroke
  - **First Author:** Meier
  - **Year:** 2025
  - **Usage:** Signal Quality (SCI + PSP thresholds)
  - **Reasoning:** Applies QT-NIRS defaults (SCI=0.8, PSP=0.1) in practice, supporting threshold choices for windowed trial QC.
  - **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC12489784/
  - **Alt Link:** https://www.sciencedirect.com/science/article/pii/S2666956025000443
  - **Source:** `fnirs_analysis/fnirs_preprocess_justifications.md`

- **Novi et al. (2023)**
  - **Title:** Revealing the spatiotemporal requirements for accurate subject identification with resting-state functional connectivity: a simultaneous fNIRS-fMRI study
  - **First Author:** Novi
  - **Year:** 2023
  - **Usage:** Subject Inclusion Criteria; Beta-table QC thresholding
  - **Reasoning:** Precedent for excluding runs/subjects with < 50% good channels; reused for condition-level channel-sufficiency thresholding in imported Homer beta-table QC.
  - **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC9896013/
  - **Source:** `fnirs_analysis/fnirs_preprocess_justifications.md`, `fnirs_analysis/homer_betas_qc.py`, `README.md`

- **Pinti et al. (2024)**
  - **Title:** Ecological functional near-infrared spectroscopy in mobile children: using short separation channels to correct for systemic contamination during naturalistic neuroimaging
  - **First Author:** Pinti
  - **Year:** 2024
  - **Usage:** Exclusion Criteria; Beta-table QC thresholding
  - **Reasoning:** Supports the exclusion rule of < 50% good quality channels and the minimum requirement of usable blocks (trials), including downstream condition-level channel sufficiency reporting.
  - **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC11460616/
  - **Source:** `fnirs_analysis/fnirs_preprocess_justifications.md`, `fnirs_analysis/homer_betas_qc.py`, `README.md`

- **Dina et al. (2025)**
  - **Title:** Measuring neurodevelopment of inhibitory control in children using naturalistic virtual reality
  - **First Author:** Dina
  - **Year:** 2025
  - **Usage:** Exclusion Criteria; Beta-table QC thresholding
  - **Reasoning:** Replicates the < 50% good-channel exclusion logic and the requirement for a minimum number of blocks per condition; used as precedent for condition-level channel sufficiency summaries.
  - **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC12289916/
  - **Source:** `fnirs_analysis/fnirs_preprocess_justifications.md`, `fnirs_analysis/homer_betas_qc.py`, `README.md`

- **Shiffler, R. E. (1988)**
  - **Title:** Maximum Z Scores and Outliers
  - **First Author:** Shiffler
  - **Year:** 1988
  - **Usage:** Mathematical feasibility bound for sample mean +/- 3 SD outlier screening
  - **Reasoning:** Shows that the largest attainable sample z-score is bounded by `(n - 1) / sqrt(n)`. Used here to justify fail-fast documentation and skipped-column reporting when fewer than 11 observed subjects are available for a `3 SD` screening rule, because such columns cannot produce a detectable outlier under the chosen method.
  - **Link:** https://doi.org/10.1080/00031305.1988.10475530
  - **Source:** `mask_homer_auc_between_subject_outliers.py`, `README.md`, `ANALYSIS_SPEC.md`

- **Leys, C., Ley, C., Klein, O., Bernard, P., & Licata, L. (2013)**
  - **Title:** Detecting outliers: Do not use standard deviation around the mean, use absolute deviation around the median
  - **First Author:** Leys
  - **Year:** 2013
  - **Usage:** Explicit limitation note for mean +/- SD outlier screening
  - **Reasoning:** Documents that mean-and-SD based outlier rules are common but non-robust because the candidate outlier influences both the mean and SD. Used here to transparently describe the limitation of the project-chosen between-subject `mean +/- 3 SD` censoring rule rather than presenting it as a robust default.
  - **Link:** https://doi.org/10.1016/j.jesp.2013.03.013
  - **Source:** `mask_homer_auc_between_subject_outliers.py`, `README.md`, `CITATIONS.md`

- **Fiske et al. (2022)**
  - **Title:** The neural correlates of inhibitory control in 10-month-old infants: A functional near-infrared spectroscopy study
  - **First Author:** Fiske
  - **Year:** 2022
  - **Usage:** Trial Count
  - **Reasoning:** Supports the requirement of at least three blocks per condition for reliable data.
  - **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7616317/
  - **Source:** `fnirs_analysis/fnirs_preprocess_justifications.md`

# Measures (Questionnaires)

- **Kroenke, K., Strine, T. W., Spitzer, R. L., Williams, J. B. W., Berry, J. T., & Mokdad, A. H. (2009)**
  - **Title:** The PHQ-8 as a measure of current depression in the general population
  - **First Author:** Kroenke
  - **Year:** 2009
  - **Usage:** PHQ-8 depressive symptom severity (sum-score construct)
  - **Reasoning:** Peer-reviewed validation reference for the PHQ-8 instrument, cited because this study intentionally administered the 8-item version even though the Qualtrics export still labels those items as “PHQ-9”.
  - **Link:** https://doi.org/10.1007/s00127-009-0113-3
  - **Source:** `process_sociodemographic.py`, `sfv_data_description.md`

- **Spitzer, R. L., Kroenke, K., Williams, J. B. W., & Löwe, B. (2006)**
  - **Title:** A brief measure for assessing generalized anxiety disorder: the GAD-7
  - **First Author:** Spitzer
  - **Year:** 2006
  - **Usage:** GAD-7 generalized anxiety symptom severity (sum-score construct)
  - **Reasoning:** Canonical peer-reviewed reference for the GAD-7 instrument, cited to justify interpretation of the `gad_total` composite score derived from Qualtrics items labeled “GAD” in this study dataset.
  - **Link:** https://doi.org/10.1001/archinte.166.10.1092
  - **Source:** `process_sociodemographic.py`, `sfv_data_description.md`

- **Kessler, R. C., Adler, L., Ames, M., et al. (2005)**
  - **Title:** The World Health Organization Adult ADHD Self-Report Scale (ASRS): A short screening scale for use in the general population
  - **First Author:** Kessler
  - **Year:** 2005
  - **Usage:** WHO ASRS adult ADHD symptom screening (sum-score construct)
  - **Reasoning:** Canonical peer-reviewed reference for the WHO ASRS instrument, cited to justify interpretation of the `asrs_total` composite score derived from Qualtrics items labeled “ASRS” in this study dataset.
  - **Link:** https://doi.org/10.1017/S0033291704002892
  - **Source:** `process_sociodemographic.py`, `sfv_data_description.md`

# Statistical Analysis

- **Laird, N. M., & Ware, J. H. (1982)**
  - **Title:** Random-effects models for longitudinal data
  - **First Author:** Laird
  - **Year:** 1982
  - **Usage:** Random-effects Models
  - **Reasoning:** Foundational framework for the random-intercept mixed models used across the within-subject neural, retention, and engagement analyses, including the additive age-adjusted omnibus specifications.
  - **Link:** https://pubmed.ncbi.nlm.nih.gov/7168798/
  - **Source:** `fnirs_analysis/FNIRS_TODO.md`, `analyze_format_content_lmm_channelwise.R`, `analyze_format_content_lmm_roi.R`, `analyze_retention_format_content_lmm.R`, `analyze_engagement_format_content_lmm.R`, `ANALYSIS_SPEC.md`, `README.md`

- **Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015)**
  - **Title:** Fitting Linear Mixed-Effects Models Using lme4
  - **First Author:** Bates
  - **Year:** 2015
  - **Usage:** Linear mixed-effects modeling
  - **Reasoning:** Primary peer-reviewed reference for the `lme4` framework used by the channelwise, ROI-wise, retention, and engagement LMM analyses in R, including the age-adjusted omnibus fits, the use of simple numerically stabilizing unit changes inside the neural fit path while back-transforming reported estimates into the original beta units, and the practice of explicitly auditing optimizer/convergence warnings alongside returned fits.
  - **Link:** https://doi.org/10.18637/jss.v067.i01
  - **Source:** `analyze_format_content_lmm_channelwise.R`, `r_lmm_convergence_helpers.R`, `tests/validate_pipeline_c_r.R`, `tests/validate_lmm_convergence_helpers_r.R`, `analyze_format_content_lmm_roi.R`, `tests/validate_pipeline_c_roi_r.R`, `analyze_retention_format_content_lmm.R`, `tests/validate_retention_pipeline_r.R`, `analyze_engagement_format_content_lmm.R`, `tests/validate_engagement_pipeline_r.R`, `tests/calibrate_type1_error_r.R`, `tests/calibrate_type2_error_r.R`, `ANALYSIS_SPEC.md`, `README.md`

- **Kuznetsova, A., Brockhoff, P. B., & Christensen, R. H. B. (2017)**
  - **Title:** lmerTest Package: Tests in Linear Mixed Effects Models
  - **First Author:** Kuznetsova
  - **Year:** 2017
  - **Usage:** Fixed-effect tests in linear mixed models
  - **Reasoning:** Justifies the `lmerTest` implementation used for mixed-model fixed-effect inference across analyses in this repo, including the age-adjusted omnibus models and the channelwise/ROI scripts that pair `lmerTest` with Kenward-Roger denominator df calculations.
  - **Link:** https://doi.org/10.18637/jss.v082.i13
  - **Source:** `analyze_format_content_lmm_channelwise.R`, `analyze_format_content_lmm_roi.R`, `analyze_retention_format_content_lmm.R`, `analyze_engagement_format_content_lmm.R`, `ANALYSIS_SPEC.md`, `README.md`

- **Kenward, M. G., & Roger, J. H. (1997)**
  - **Title:** Small Sample Inference for Fixed Effects from Restricted Maximum Likelihood
  - **First Author:** Kenward
  - **Year:** 1997
  - **Usage:** Kenward-Roger denominator degrees of freedom and F-test approximation
  - **Reasoning:** Primary methodological reference for the Kenward-Roger approximation used for fixed-effect tests in the channelwise and ROI format×content LMM scripts, including their age-adjusted omnibus models.
  - **Link:** https://doi.org/10.2307/2533558
  - **Source:** `analyze_format_content_lmm_channelwise.R`, `analyze_format_content_lmm_roi.R`, `ANALYSIS_SPEC.md`

- **Halekoh, U., & Højsgaard, S. (2014)**
  - **Title:** A Kenward-Roger Approximation and Parametric Bootstrap Methods for Tests in Linear Mixed Models - The R Package pbkrtest
  - **First Author:** Halekoh
  - **Year:** 2014
  - **Usage:** R implementation of Kenward-Roger tests via `pbkrtest`
  - **Reasoning:** Justifies the software implementation path used by the R scripts to obtain Kenward-Roger denominator df and p-values for `lmer` models in the age-adjusted channelwise and ROI omnibus analyses.
  - **Link:** https://doi.org/10.18637/jss.v059.i09
  - **Source:** `analyze_format_content_lmm_channelwise.R`, `analyze_format_content_lmm_roi.R`, `ANALYSIS_SPEC.md`

- **Satterthwaite, F. E. (1946)**
  - **Title:** An approximate distribution of estimates of variance components
  - **First Author:** Satterthwaite
  - **Year:** 1946
  - **Usage:** Approximate degrees of freedom (Satterthwaite)
  - **Reasoning:** Foundational reference for the Satterthwaite df approximation used (via `lmerTest`) in the retention and engagement LMM analyses, including the additive age-adjusted omnibus fits.
  - **Link:** https://doi.org/10.2307/3002019
  - **Source:** `analyze_retention_format_content_lmm.R`, `analyze_engagement_format_content_lmm.R`, `ANALYSIS_SPEC.md`, `README.md`

- **Poldrack, R. A. (2007)**
  - **Title:** Region of interest analysis for fMRI
  - **First Author:** Poldrack
  - **Year:** 2007
  - **Usage:** ROI signal extraction strategy
  - **Reasoning:** Provides the methodological framework for extracting a single summary signal from pre-specified ROIs; used here to justify ROI-level summary of channel betas before mixed-model inference and before the post-hoc neural format-effect association analysis.
  - **Link:** https://doi.org/10.1093/scan/nsm006
  - **Source:** `analyze_format_content_lmm_roi.R`, `analyze_correlational_relationships.R`, `analyze_correlational_relationships_roi_means.R`, `plot_beta_discrepancy_dynamics.py`, `ANALYSIS_SPEC.md`, `README.md`

- **Morey, R. D. (2008)**
  - **Title:** Confidence intervals from normalized data: A correction to Cousineau (2005)
  - **First Author:** Morey
  - **Year:** 2008
  - **Usage:** Within-subject confidence intervals for repeated-measures summaries
  - **Reasoning:** Provides the correction factor for Cousineau-normalized confidence intervals in repeated-measures plots. Used here so condition-level beta summaries in the channel-vs-ROI discrepancy visualization reflect the within-subject structure of the design rather than between-subject error bars.
  - **Link:** https://doi.org/10.20982/tqmp.04.2.p061
  - **Source:** `plot_beta_discrepancy_dynamics.py`, `README.md`

- **Kriegeskorte, N., Simmons, W. K., Bellgowan, P. S. F., & Baker, C. I. (2009)**
  - **Title:** Circular analysis in systems neuroscience: the dangers of double dipping
  - **First Author:** Kriegeskorte
  - **Year:** 2009
  - **Usage:** Exploratory interpretation / selective-inference caution for data-informed target selection
  - **Reasoning:** Documents that reusing the same dataset for feature selection and downstream inference can bias apparent evidence strength. Used here to justify labeling the targeted channel/ROI follow-up correlations as exploratory when the target set is selected from this dataset rather than fixed a priori.
  - **Link:** https://doi.org/10.1038/nn.2303
  - **Source:** `analyze_correlational_relationships.R`, `analyze_channel_behavior_relationships.py`, `analyze_behavior_pairwise_correlations.R`, `tests/validate_behavior_pairwise_correlations_r.R`, `README.md`, `ANALYSIS_SPEC.md`, `CITATIONS.md`

- **Spearman, C. (1904)**
  - **Title:** The Proof and Measurement of Association between Two Things
  - **First Author:** Spearman
  - **Year:** 1904
  - **Usage:** Spearman rank correlation for exploratory channel-behavior screening and sensitivity analyses
  - **Reasoning:** Foundational reference for Spearman's rank correlation. Used here because several behavioral outcomes are bounded or only plausibly monotonic with neural effects, making a rank-based sensitivity check more defensible than assuming linear-normal relationships for every association.
  - **Link:** https://doi.org/10.2307/1412159
  - **Source:** `analyze_channel_behavior_relationships.py`, `tests/validate_channel_behavior_relationships_py.py`, `analyze_correlational_relationships.R`, `analyze_correlational_relationships_roi_means.R`, `tests/validate_correlational_relationships_r.R`, `README.md`, `ANALYSIS_SPEC.md`

- **Bonett, D. G. (2020)**
  - **Title:** Point-biserial correlation: Interval estimation, hypothesis testing, meta-analysis, and sample size determination
  - **First Author:** Bonett
  - **Year:** 2020
  - **Usage:** Point-biserial correlation for binary behavioral predictors
  - **Reasoning:** Modern peer-reviewed reference for point-biserial correlation as a standardized effect size and test for dichotomous-vs-continuous associations. Used here for binary behavioral variables such as `pd_status` when screening against continuous channel betas, and to justify interpreting 0/1-coded Pearson behavioral-screen coefficients on the same dichotomous-vs-continuous effect-size scale.
  - **Link:** https://doi.org/10.1111/bmsp.12189
  - **Source:** `analyze_channel_behavior_relationships.py`, `tests/validate_channel_behavior_relationships_py.py`, `analyze_behavior_pairwise_correlations.R`, `tests/validate_behavior_pairwise_correlations_r.R`, `README.md`, `ANALYSIS_SPEC.md`

- **Pearson, K. (1896)**
  - **Title:** Mathematical Contributions to the Theory of Evolution. III. Regression, Heredity, and Panmixia
  - **First Author:** Pearson
  - **Year:** 1896
  - **Usage:** Pearson product-moment correlation
  - **Reasoning:** Foundational source for the product-moment correlation coefficient used in the exploratory post-hoc analysis between continuous pooled long/short neural means and the matching pooled long/short behavioral means, the supplementary raw behavioral task-cell values tested against those same pooled neural means, and the standalone pairwise behavioral screening analysis.
  - **Link:** https://doi.org/10.1098/rsta.1896.0007
  - **Source:** `analyze_correlational_relationships.R`, `analyze_correlational_relationships_roi_means.R`, `analyze_behavior_pairwise_correlations.R`, `tests/validate_correlational_relationships_r.R`, `tests/validate_behavior_pairwise_correlations_r.R`, `README.md`, `ANALYSIS_SPEC.md`

- **Fisher, R. A. (1921)**
  - **Title:** On the "Probable Error" of a Coefficient of Correlation Deduced from a Small Sample
  - **First Author:** Fisher
  - **Year:** 1921
  - **Usage:** Fisher z confidence intervals for Pearson correlation
  - **Reasoning:** Provides the variance-stabilizing transformation underlying the confidence intervals reported for Pearson correlations in the post-hoc pooled long/short analyses, the supplementary raw-behavior versus pooled-neural association analyses, and the standalone pairwise behavioral screening analysis.
  - **Link:** http://hdl.handle.net/2440/15169
  - **Source:** `analyze_correlational_relationships.R`, `analyze_correlational_relationships_roi_means.R`, `analyze_behavior_pairwise_correlations.R`, `tests/validate_correlational_relationships_r.R`, `tests/validate_behavior_pairwise_correlations_r.R`, `README.md`, `ANALYSIS_SPEC.md`

- **Cleveland, W. S. (1979)**
  - **Title:** Robust Locally Weighted Regression and Smoothing Scatterplots
  - **First Author:** Cleveland
  - **Year:** 1979
  - **Usage:** LOESS trend line for Spearman correlation figures
  - **Reasoning:** Provides the classic locally weighted smoothing approach used here to visualize monotonic Spearman relationships without imposing the linear model that underlies Pearson plots. This makes the plot annotation consistent with the rank-based inferential target while still showing the local trend in the original measurement scale.
  - **Link:** https://doi.org/10.1080/01621459.1979.10481038
  - **Source:** `analyze_correlational_relationships_roi_means.R`, `README.md`, `CITATIONS.md`

- **Bender, R., & Lange, S. (2001)**
  - **Title:** Adjusting for multiple testing - when and how?
  - **First Author:** Bender
  - **Year:** 2001
  - **Usage:** Multiple-testing family definition for exploratory/post-hoc association analysis
  - **Reasoning:** Summarizes how multiplicity correction should be matched to the final inferential claim rather than applied mechanically. Used here to justify defining BH families at the analysis-tier x behavior-run x format grouping x association-method level for the pooled long/short primary analysis and the supplementary raw-behavior follow-up rows.
  - **Link:** https://doi.org/10.1016/S0895-4356(00)00314-0
  - **Source:** `analyze_correlational_relationships.R`, `analyze_correlational_relationships_roi_means.R`, `tests/validate_correlational_relationships_r.R`, `README.md`, `ANALYSIS_SPEC.md`

- **Holm, S. (1979)**
  - **Title:** A Simple Sequentially Rejective Multiple Test Procedure
  - **First Author:** Holm
  - **Year:** 1979
  - **Usage:** Multiple testing correction (omnibus effects)
  - **Reasoning:** Justifies the Holm correction used across the three pre-planned omnibus effects in the retention and engagement supporting analyses.
  - **Link:** https://doi.org/10.2307/4615733
  - **Source:** `analyze_retention_format_content_lmm.R`, `tests/validate_retention_pipeline_r.R`, `analyze_engagement_format_content_lmm.R`, `tests/validate_engagement_pipeline_r.R`

- **Burton, A., Altman, D. G., Royston, P., & Holder, R. L. (2006)**
  - **Title:** The design of simulation studies in medical statistics
  - **First Author:** Burton
  - **Year:** 2006
  - **Usage:** Monte Carlo type-I/type-II calibration test design
  - **Reasoning:** Provides guidance for designing simulation studies to evaluate operating characteristics (including false-positive and false-negative behavior) of statistical methods under controlled data-generating mechanisms; used to justify repeated null/non-null calibration harnesses in this repo.
  - **Link:** https://doi.org/10.1002/sim.2673
  - **Source:** `tests/calibrate_type1_error_r.R`, `tests/calibrate_type2_error_r.R`, `tests/README.md`

- **Pinheiro, J. C., & Bates, D. M. (2000)**
  - **Title:** Mixed-Effects Models in S and S-PLUS
  - **First Author:** Pinheiro
  - **Year:** 2000
  - **Usage:** Mixed-Effects Models
  - **Reasoning:** Standard reference text for assumptions and practice in mixed-effects models.
  - **Link:** https://link.springer.com/book/10.1007/b98882
  - **Source:** `fnirs_analysis/FNIRS_TODO.md`

- **Friston, K. J., et al. (2002)**
  - **Title:** Classical and Bayesian Inference in Neuroimaging: Theory
  - **First Author:** Friston
  - **Year:** 2002
  - **Usage:** Inference
  - **Reasoning:** Summary-statistics and hierarchical inference logic used broadly in neuroimaging.
  - **Link:** https://doi.org/10.1006/nimg.2002.1090
  - **Source:** `fnirs_analysis/FNIRS_TODO.md`

- **Benjamini, Y., & Hochberg, Y. (1995)**
  - **Title:** Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing
  - **First Author:** Benjamini
  - **Year:** 1995
  - **Usage:** False Discovery Rate (FDR)
  - **Reasoning:** Practical approach for controlling the false discovery rate in multiple testing.
  - **Link:** https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
  - **Source:** `fnirs_analysis/FNIRS_TODO.md`, `analyze_format_content_lmm_channelwise.R`, `analyze_format_content_lmm_channelwise.py`, `tests/validate_pipeline_c_r.R`, `analyze_format_content_lmm_roi.R`, `tests/validate_pipeline_c_roi_r.R`, `analyze_correlational_relationships.R`, `analyze_channel_behavior_relationships.py`, `tests/validate_correlational_relationships_r.R`, `tests/validate_channel_behavior_relationships_py.py`, `README.md`

- **Lenth, R. V. (2016)**
  - **Title:** Least-Squares Means: The R Package lsmeans
  - **First Author:** Lenth
  - **Year:** 2016
  - **Usage:** Post-hoc contrasts / estimated marginal means
  - **Reasoning:** Peer-reviewed reference for least-squares means / estimated marginal means and their use for post-hoc contrasts in linear (mixed) models; emmeans is the modern successor in R used for post-hoc comparisons in this repo.
  - **Link:** https://doi.org/10.18637/jss.v069.i01
  - **Source:** `analyze_format_content_lmm_channelwise.R`, `tests/validate_pipeline_c_r.R`, `analyze_format_content_lmm_roi.R`, `tests/validate_pipeline_c_roi_r.R`, `analyze_retention_format_content_lmm.R`, `analyze_engagement_format_content_lmm.R`

- **Searle, S. R., Speed, F. M., & Milliken, G. A. (1980)**
  - **Title:** Population marginal means in the linear model: An alternative to least squares means
  - **First Author:** Searle
  - **Year:** 1980
  - **Usage:** Estimated marginal means concept
  - **Reasoning:** Foundational reference for population marginal means / least-squares means terminology and interpretation underlying EMM-based post-hoc summaries and the equal-weight marginal-mean logic used when collapsing the 2x2 design to a long-vs-short main-effect contrast.
  - **Link:** https://doi.org/10.1080/00031305.1980.10483031
  - **Source:** `analyze_format_content_lmm_channelwise.R`, `analyze_format_content_lmm_roi.R`, `analyze_retention_format_content_lmm.R`, `analyze_engagement_format_content_lmm.R`

- **Student (1908)**
  - **Title:** The Probable Error of a Mean
  - **First Author:** Student
  - **Year:** 1908
  - **Usage:** Paired t-tests in post-hoc contrasts (Python)
  - **Reasoning:** Foundational reference for the t-test framework used for within-subject paired contrasts in the Python post-hoc step (uncorrected per study instruction).
  - **Link:** https://doi.org/10.1093/biomet/6.1.1
  - **Source:** `analyze_format_content_lmm_channelwise.py`

- **Ye, J. C., et al. (2009)**
  - **Title:** NIRS-SPM: Statistical parametric mapping for near-infrared spectroscopy
  - **First Author:** Ye
  - **Year:** 2009
  - **Usage:** NIRS-SPM
  - **Reasoning:** Statistical parametric mapping and GLM framework for fNIRS, including drift and noise handling.
  - **Link:** https://doi.org/10.1016/j.neuroimage.2008.08.036
  - **Source:** `fnirs_analysis/FNIRS_TODO.md`

- **Pinti, P., et al. (2019)**
  - **Title:** Current Status and Issues Regarding Pre-processing of fNIRS Neuroimaging Data: An Investigation of Diverse Signal Filtering Methods Within a General Linear Model Framework
  - **First Author:** Pinti
  - **Year:** 2019
  - **Usage:** GLM Best Practices
  - **Reasoning:** Review on the use of the General Linear Model in fNIRS, offering recommendations for filtering, drift handling, and inference.
  - **Link:** https://www.frontiersin.org/articles/10.3389/fnhum.2018.00505/full
  - **Source:** `fnirs_analysis/FNIRS_TODO.md`

# Reproducibility

- **Sandve, G. K., et al. (2013)**
  - **Title:** Ten Simple Rules for Reproducible Computational Research
  - **First Author:** Sandve
  - **Year:** 2013
  - **Usage:** Centralized, auditable workflow controls (exclusions + orchestration + certification + derivation provenance)
  - **Reasoning:** Supports maintaining a single, auditable source of analysis decisions and execution order (participant exclusions, pipeline orchestration, machine-readable certification artifacts, FIR-to-AUC provenance sidecars, outlier-screening audit artifacts, and fail-fast one-row-per-subject ID validation in upstream tabular preprocessing) to prevent script-specific drift across inferential endpoints.
  - **Link:** https://doi.org/10.1371/journal.pcbi.1003285
  - **Source:** `r_subject_exclusions.R`, `pipeline_preprocess_merge.sh`, `collapse_homer_fir_to_auc.py`, `mask_homer_auc_between_subject_outliers.py`, `certify_preprocess_merge_integrity.py`, `validate_homer_fir_auc_conversion.py`, `generate_combined_data.py`, `process_sociodemographic.py`, `README.md`, `analyze_format_content_lmm_roi.R`
