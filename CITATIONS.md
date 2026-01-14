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
  - **Usage:** Best Practices
  - **Reasoning:** Guidelines for fNIRS acquisition, analysis, and reporting, specifically regarding systemic physiology and short-separation channels.
  - **Link:** https://doi.org/10.1117/1.NPh.8.1.012101
  - **Source:** `fnirs_analysis/FNIRS_TODO.md`

# Quality Control

- **Pollonini, L., et al. (2016)**
  - **Title:** PHOEBE: a method for real time mapping of optodes-scalp coupling in functional near-infrared spectroscopy
  - **First Author:** Pollonini
  - **Year:** 2016
  - **Usage:** Scalp Coupling Index (SCI)
  - **Reasoning:** Supports the use of an SCI threshold (e.g., 0.8) to identify and drop channels with poor optode-scalp coupling.
  - **Link:** https://doi.org/10.1364/BOE.7.005104
  - **Source:** `fnirs_analysis/fnirs_analysis.py`, `fnirs_analysis/fnirs_preprocess_justifications.md`

- **Hocke, Oni, Duszynski, Corrigan, Frederick, Dunn (2018)**
  - **Title:** Automated Processing of fNIRS Data—A Visual Guide to the Pitfalls and Consequences
  - **First Author:** Hocke
  - **Year:** 2018
  - **Usage:** SCI Threshold
  - **Reasoning:** Uses a hard-coded SCI threshold of ≥ 0.7 for channel acceptance, supporting the run-level exclusion criteria.
  - **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC6428450/
  - **Source:** `fnirs_analysis/fnirs_preprocess_justifications.md`

- **Hernandez & Pollonini (2020)**
  - **Title:** NIRSplot: A Tool for Quality Assessment of fNIRS Scans
  - **First Author:** Hernandez
  - **Year:** 2020
  - **Usage:** QC Assessment
  - **Reasoning:** Establishes ~0.8 as a reasonable threshold for good optode-scalp coupling in the QT-NIRS/NIRSplot framework.
  - **Link:** https://par.nsf.gov/servlets/purl/10209495
  - **Source:** `fnirs_analysis/fnirs_preprocess_justifications.md`

- **Meier et al. (2025)**
  - **Title:** The effects of protocol factors and participant characteristics on functional near-infrared spectroscopy data quality after stroke
  - **First Author:** Meier
  - **Year:** 2025
  - **Usage:** Signal Quality
  - **Reasoning:** Discusses the effects of protocol factors on fNIRS signal quality and reinforces SCI as a principled data-quality gate.
  - **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC12489784/
  - **Source:** `fnirs_analysis/fnirs_preprocess_justifications.md`

- **Novi et al. (2023)**
  - **Title:** Revealing the spatiotemporal requirements for accurate subject identification with resting-state functional connectivity: a simultaneous fNIRS-fMRI study
  - **First Author:** Novi
  - **Year:** 2023
  - **Usage:** Subject Inclusion Criteria
  - **Reasoning:** Precedent for excluding runs/subjects with < 50% good channels.
  - **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC9896013/
  - **Source:** `fnirs_analysis/fnirs_preprocess_justifications.md`

- **Pinti et al. (2024)**
  - **Title:** Ecological functional near-infrared spectroscopy in mobile children: using short separation channels to correct for systemic contamination during naturalistic neuroimaging
  - **First Author:** Pinti
  - **Year:** 2024
  - **Usage:** Exclusion Criteria
  - **Reasoning:** Supports the exclusion rule of < 50% good quality channels and the minimum requirement of usable blocks (trials).
  - **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC11460616/
  - **Source:** `fnirs_analysis/fnirs_preprocess_justifications.md`

- **Dina et al. (2025)**
  - **Title:** Measuring neurodevelopment of inhibitory control in children using naturalistic virtual reality
  - **First Author:** Dina
  - **Year:** 2025
  - **Usage:** Exclusion Criteria
  - **Reasoning:** Replicates the < 50% good-channel exclusion logic and the requirement for a minimum number of blocks per condition.
  - **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC12289916/
  - **Source:** `fnirs_analysis/fnirs_preprocess_justifications.md`

- **Fiske et al. (2022)**
  - **Title:** The neural correlates of inhibitory control in 10-month-old infants: A functional near-infrared spectroscopy study
  - **First Author:** Fiske
  - **Year:** 2022
  - **Usage:** Trial Count
  - **Reasoning:** Supports the requirement of at least three blocks per condition for reliable data.
  - **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7616317/
  - **Source:** `fnirs_analysis/fnirs_preprocess_justifications.md`

# Statistical Analysis

- **Laird, N. M., & Ware, J. H. (1982)**
  - **Title:** Random-effects models for longitudinal data
  - **First Author:** Laird
  - **Year:** 1982
  - **Usage:** Random-effects Models
  - **Reasoning:** Foundational framework for random-effects models in longitudinal data analysis.
  - **Link:** https://pubmed.ncbi.nlm.nih.gov/7168798/
  - **Source:** `fnirs_analysis/FNIRS_TODO.md`

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
  - **Source:** `fnirs_analysis/FNIRS_TODO.md`

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
