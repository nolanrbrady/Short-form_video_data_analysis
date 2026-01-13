# fNIRS QC / Exclusion Criteria Support (Run → Channel → Subject)

This file lists **peer‑reviewed** precedents that support the three QC choices you committed to for your short‑form video fNIRS study.

---

## 1) Run-level Scalp Coupling Index (SCI) cutoff: **exclude runs with SCI < 0.7**

- **Hocke, Oni, Duszynski, Corrigan, Frederick, Dunn (2018)** — *Automated Processing of fNIRS Data—A Visual Guide to the Pitfalls and Consequences*  
  **Reason for inclusion:** Uses the PHOEBE optode–scalp coupling approach with a **hard-coded SCI threshold of ≥ 0.7** as part of its channel acceptance logic (i.e., channels below this coupling level are treated as poor quality).  
  **Link:** https://www.mdpi.com/1999-4893/11/5/67

- **Hernandez & Pollonini (2020)** — *NIRSplot: A Tool for Quality Assessment of fNIRS Scans*  
  **Reason for inclusion:** Establishes a commonly used “good coupling” prior for SCI, noting **~0.8 as a reasonable threshold** for good optode–scalp coupling in their QT‑NIRS/NIRSplot framework (useful as an upper-end reference; your 0.7 run‑level rule is a slightly more permissive, pragmatic cutoff).  
  **Link:** https://par.nsf.gov/servlets/purl/10209495

- **Meier et al. (2025)** — *The effects of protocol factors and participant characteristics on fNIRS signal quality…*  
  **Reason for inclusion:** Discusses QT‑NIRS and the use of **SCI (with PSP) as a principled data‑quality gate**, reinforcing SCI as a mainstream QC signal rather than an ad‑hoc metric.  
  **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC12489784/

---

## 2) Run/subject inclusion rule: **exclude if < 50% “good” channels**

- **Novi et al. (2023)** — *Revealing the spatiotemporal requirements for accurate subject identification with resting-state functional connectivity: a simultaneous fNIRS–fMRI study*  
  **Reason for inclusion:** Explicitly **excluded individual runs that did not have at least 50% good channels**, and excluded participants who did not have enough good runs after applying that pruning rule. This is a very direct precedent for your “<50% good channels” cutoff.  
  **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC9896013/

- **Pinti et al. (2024)** — *Ecological functional near-infrared spectroscopy in mobile children: using short separation channels to correct for systemic contamination during naturalistic neuroimaging*  
  **Reason for inclusion:** Participants were included only if they had **≥ 50% of channels of good quality** (plus a minimum number of usable blocks). While developmental, it’s a clean, explicit rule that generalizes well to adult fNIRS QC logic.  
  **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC11460616/

- **Dina et al. (2025)** — *Measuring neurodevelopment of inhibitory control in children using naturalistic virtual reality*  
  **Reason for inclusion:** Replicates the same **<50% good‑channel exclusion** logic (and ≥3 blocks) while explicitly stating the cutoff in the methods, providing independent confirmation that “50% good channels” is a defensible inclusion gate.  
  **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC12289916/

---

## 3) Minimum usable trials per condition: **require ≥ 3 of 4 usable trials per condition (subject-level)**

*(Many fNIRS papers are block-based; “≥3 blocks per condition” is the closest structural match to your “≥3 trials per condition” rule.)*

- **Pinti et al. (2024)** — *Ecological functional near-infrared spectroscopy in mobile children…*  
  **Reason for inclusion:** Participants were included only if they had **at least three blocks** meeting performance criteria (in addition to channel‑quality criteria). This supports “3 is the minimum usable unit count” logic for stable estimation.  
  **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC11460616/

- **Dina et al. (2025)** — *Measuring neurodevelopment of inhibitory control…*  
  **Reason for inclusion:** States that participants were excluded if they had **<3 blocks** meeting criteria, reinforcing “≥3” as a practical minimum for condition-level estimation.  
  **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC12289916/

- **Fiske et al. (2022)** — *The neural correlates of inhibitory control in 10-month-old infants: A functional near-infrared spectroscopy study*  
  **Reason for inclusion:** Notes participants were encouraged to complete **at least three blocks of each type** to ensure enough reliable fNIRS data—directly supporting “≥3 per condition” as a defensible minimum.  
  **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7616317/
