# Homer3 ProcStream (Defensible, Literature-Backed) — NIRSport (NIRx) 8×8 PFC, Block Design (120 s task / 30 s rest), **No short-sep channels**

## Study Overview
This is an fNIRS study on 18–32 year-old college students watching four video types (Short-form education, Short-form entertainment, Long-form education, Long-form entertainment). Each trial is **120 s task + 30 s rest**, repeated four times per condition, within-subject. Primary interest is the **2×2 interaction** (Content × Length).

This README captures a **Homer3 ProcStream** that is defensible for your design and explicitly flags what matters most when you **do not have short-sep channels**.

---

## Core constraints you cannot “preprocess away”
### No short-separation regression
Without short channels, **systemic/scalp physiology is harder to remove** and may bias activation estimates. PHOEBE explicitly calls this out:

> “The NIR light… interrogates the cerebral cortex, but to a larger extent also the extracerebral tissue layers… affect the fNIRS signals and result in potential misinterpretation.”  [oai_citation:0‡PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7793571/)

So the pipeline below is “best-possible under constraints,” but you should still interpret results with that limitation in mind (and report it).

---

## Recommended ProcStream (in order)

### 1) `hmrR_PruneChannels`
**Goal:** create a sane `mlActAuto` channel mask **without nuking the montage**.

**Parameters**
- `dRange = [0.1 5000]`  *(keeping your current choice, per your request)*
- `SNRthresh = 2`
- `SDrange = [25 45]` mm

**Why (with paper-grounded rationale)**
- PHOEBE expects you to **explicitly report** channel-quality thresholds used for rejection:

  > “The Methods section may thus include an indication of the SNR threshold… utilized to reject data channels…”  [oai_citation:1‡PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7793571/)

- Adult long channels are typically around **3 cm separation**, so bounding `SDrange` around that is defensible:

  > “Most frequently a source-detector separation of 3 cm is used…”  [oai_citation:2‡PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6428450/)

**Practical QC you should run (non-negotiable)**
- Report `% channels kept` by `mlActAuto` (per subject; also mean±SD). PHOEBE also recommends reporting channel/participant retention to avoid misinterpretation.  [oai_citation:3‡PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7793571/)

---

### 2) `hmrR_Intensity2OD`
**Parameters:** none

**Why**
- Optical density (OD) is the standard space for motion handling and MBLL conversion; Homer pipelines assume this ordering. (This step is foundational and generally uncontroversial.)

---

### 3) `hmrR_MotionArtifactByChannel`  *(detection/masking)*
**Parameters (Brigadoi-style HOMER settings)**
- `tMotion = 1.0`
- `tMask   = 1.0`
- `STDEVthresh = 50`
- `AMPthresh   = 0.4`

**Why**
- Brigadoi et al. explicitly report these HOMER motion-detection parameters:

  > “tMotion = 1; tMask = 1; stdevThreshold = 50; ampThreshold = 0.4.”  [oai_citation:4‡PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3762942/)

**Important sanity check (because these thresholds can be conservative)**
- Verify that motion is actually being flagged for at least some participants (visualize `tIncAuto` / masks). If *nothing* is ever flagged, you’re not “clean,” you’re “blind.”

PHOEBE also explicitly expects motion handling + parameters to be reported:

> “…handling and correction of motion artifacts and related parameters should be reported (e.g., the thresholds for identifying motion…)”  [oai_citation:5‡PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7793571/)

---

### 4) `hmrR_MotionCorrectWavelet`
**Parameters**
- `iqr = 1.5`
- `turnon = 1`

**Why**
- A commonly used wavelet setting is explicitly described as:

  > “wavelet transform (iqr = 1.5…) (based on Molavi and Dumont, 2012)”  [oai_citation:6‡PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6428450/)

Wavelet correction remains a strong default when you expect head motion and want a well-tested HOMER-compatible approach.

---

### 5) `hmrR_BandpassFilt` *(on OD)*
**Parameters**
- `hpf = 0.003`
- `lpf = 0.2`

**Why**
- For **long, sustained blocks**, Pinti et al. recommend a **smaller high-pass cutoff** than the common 0.01 Hz:

  > “For sustained periods of stimulation longer than 100s, smaller cut-off frequency (R\_Tddr\_CompFc\_Low) should be used.”  [oai_citation:7‡Frontiers](https://www.frontiersin.org/journals/human-neuroscience/articles/10.3389/fnhum.2018.00505/full)

  Your block is **120 s**, so `hpf=0.003` is intentionally conservative (to avoid removing very slow task structure).

- A low-pass around **0.2 Hz** has a direct rationale: it targets the higher-frequency physiological bands (respiration/cardiac and instrument noise). One preprocessing guide states:

  > “Band-pass filtering is commonly used to exclude high frequency oscillations (≥0.2 Hz; heartbeat, respiratory rate, instrument noise).”  [oai_citation:8‡PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6428450/)

---

### 6) `hmrR_OD2Conc`
**Parameters**
- `ppf = [DPF(760nm) DPF(850nm)]` *(use your device’s wavelengths; verify in your SNIRF metadata)*

**Why (and what you must report)**
PHOEBE is very explicit that MBLL/DPF choices must be stated:

> “In all cases, the method of choice and relevant parameters (e.g., DPF) should be stated and citations should be provided.”  [oai_citation:9‡PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7793571/)

#### DPF calculation (Scholkmann & Wolf-style equation via published coefficients)
One peer-reviewed preprocessing guide reproduces the age- and wavelength-dependent DPF equation (attributed to Scholkmann & Wolf) with coefficients:

> “DPF = a + b*A^c + d*λ^3 + e*λ^2 + f*λ … a=223.3; b=0.05624; c=0.8493; d=-5.723e−7; e=0.001245; f=−0.9025.”  [oai_citation:10‡PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6428450/)

Using that equation (A in years; λ in nm), the DPF values are:

| Mean age (A) | DPF @ 760 nm | DPF @ 850 nm |
|---:|---:|---:|
| 20 | 6.00 | 4.94 |
| 21 | 6.03 | 4.97 |
| 22 | 6.06 | 5.00 |

**Recommendation for your dataset (mean age ≈ 20):**
- `ppf = [6.00 4.94]` *(assuming wavelengths are 760/850; confirm ordering in SNIRF)*

Best practice (if you have per-subject ages available): compute DPF per subject and apply per subject, then report the method and the age range.

---

### 7) `hmrR_BlockAvg` *(QC / visualization; optional if GLM is primary)*
**Parameters**
- `trange = [-10 140]`

**Why**
Your task is 120 s, and you have ~150 s between onsets. This window captures baseline, full block, and recovery without bleeding into the next trial.

---

### 8) `hmrR_GLM`
This is where two “gotchas” matter: **solver choice** and **making the basis match a 120 s block**.

**Parameters (recommended)**
- `trange = [-10 140]`
- `glmSolveMethod = 2`  *(iWLS / AR-IRLS in Homer3)*
- `driftOrder = 0`  *(required with iWLS in Homer3; see note below)*
- `idxBasis = 2`  *(modified gamma ⊗ square-wave block)*
- `paramsBasis = [0.1 3.0 1.8 3.0]` *(Homer defaults)*
- `rhoSD_ssThresh = 15`
- `flagSSmethod = 0` *(no short-sep regression)*
- `c_vector = 0`

**Why (Homer3 documentation + rationale)**
#### A) `glmSolveMethod`
Homer3 defines the solvers as:

> “glmSolveMethod: 1 - OLS; 2 - iWLS.” 

and includes a critical constraint:

> “NOTE: driftOrder should be set to 0 when using iWLS or may result in spurious results.” 

So if you want AR-IRLS-style robustness (iWLS), you **must** use `driftOrder=0` in Homer3. Drift is then handled by your conservative high-pass (`hpf=0.003`) plus the model intercept.

#### B) `idxBasis` and `paramsBasis`
Homer3 defines the relevant basis exactly:

> “idxBasis: 2 - Modified Gamma function convolved with square-wave of duration T” 

and its parameterization:

> “paramsBasis: [tau1 sigma1 T1 tau2 sigma2 T2]” 

Homer3 also lists default values (for an event-like T=10), which you should **only change in the T slots** to match your block:

> “defaults: paramsBasis = [0.1 3.0 10 1.8 3.0 10]” 

So for your 120 s block, the defensible change is:
- keep tau/sigma at Homer defaults (shape priors)

That yields:
- `paramsBasis = [0.1 3.0 1.8 3.0]`

**AR-IRLS Justification**
Barker et al introduce the iWLS method and demonstrate that it out performs OLS (https://pmc.ncbi.nlm.nih.gov/articles/PMC3756568/https://pmc.ncbi.nlm.nih.gov/articles/PMC3756568/)

---

## How to justify `paramsBasis` values (the logic you put in the Methods)
1) **Pick a basis that matches your paradigm.**  
   For sustained blocks, Homer’s `idxBasis=2` is explicitly “modified gamma ⊗ square-wave of duration T.”   
   That’s exactly what you want: an impulse-like HRF shape “stretched” by the block.

2) **Set the duration parameter(s) T to the true stimulus duration.**  
   Your stimulus is 120 s, so you set `T=120` (for both HbO and HbR entries). This is the *highest-leverage, design-driven* basis choice.

3) **Keep tau/sigma at defaults unless you have a physiologically grounded reason to deviate.**  
   Tau/sigma are shape priors. For a long block, estimates are typically not highly sensitive to small tweaks in tau/sigma compared with getting T right.

4) **Optional (but very defensible) sensitivity check:**  
   Run a tiny subset with `idxBasis=1` (canonical) vs `idxBasis=2` and confirm your key fixed effects don’t flip direction. Then report that basis choice did not qualitatively change conclusions.

---

## Copy/paste ProcStream settings (final)

### `hmrR_PruneChannels`
- `dRange`: `0.1 5000`
- `SNRthresh`: `2`
- `SDrange`: `25 45`

### `hmrR_MotionArtifactByChannel`
- `tMotion`: `1.0`
- `tMask`: `1.0`
- `STDEVthresh`: `50`
- `AMPthresh`: `0.4`

### `hmrR_MotionCorrectWavelet`
- `iqr`: `1.5`
- `turnon`: `1`

### `hmrR_BandpassFilt (OpticalDensity)`
- `hpf`: `0.003`
- `lpf`: `0.2`

### `hmrR_OD2Conc`
- `ppf`: `6.00 4.94`  *(for mean age ≈ 20 and wavelengths 760/850; confirm wavelength order in SNIRF)*

### `hmrR_BlockAvg`
- `trange`: `-10 140`

### `hmrR_GLM`
- `trange`: `-10 140`
- `glmSolveMethod`: `2`
- `driftOrder`: `0`
- `idxBasis`: `2`
- `paramsBasis`: `0.1 3.0 1.8 3.0`
- `rhoSD_ssThresh`: `0`
- `flagSSmethod`: `0`
- `c_vector`: `0`

---

## One last “are we missing anything?” checklist (because no short-seps)
- **Systemic physiology risk** (scalp + systemic): acknowledge and report as a limitation; consider adding a nuisance strategy later (e.g., global/PCA regressors) if reviewers push on it.  [oai_citation:11‡PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7793571/)
- **Verify motion + pruning aren’t silently doing nothing** (or everything).
- **Confirm basis duration matches task** (T=120 is the big one). 