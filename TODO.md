A list of things I need to do for this analysis pipeline.

## Data Quality Evaluation
- [ ] Create a data quality evaluation script that gives subject level quality indications like number of bad channels, SCI, SNR, etc.
- [ ] Look into ways of visualy confirming that the pipeline is successfully preprocessing the data.

## fNIRS Preprocessing
- [ ] Remove the redundant high-pass filtering (Keep low pass filter and GLM Cosine drift model but remove the FIR high-pass filter) (Pinti, P., et al. (2019). The present and future use of functional near-infrared spectroscopy (fNIRS) for cognitive neuroscience. Annals of the New York Academy of Sciences.) 
    * Citations to read for filtering parameters -> (Pinti et al., 2019; Yücel et al., 2021)
- [ ] Consider using a stricter scalp coupling index (from 0.5 to 0.7-0.8) (Pollonini, L., et al. (2014). PHOEBE: a processing system for functional near-infrared spectroscopy data.)
- [ ] Force bad channels to be removed prior to GLM
- [ ] Figure out how to exclude trials if they don't fit the exclusion criteria
- [ ] Check and remove trials that don't meet the quality criteria (SCI < 0.8, Bad Channels > 50%, < 2 usable trials per condition)
- [ ] Consider using boxcar convolved HRF or FIR approach instead of the current `SPM + derivative` approach for HRF_MODEL.
    * SPM + derivative may not be the best fit for the long run tasks.
 

## Analysis Pipeline
- [ ] Change LME from `theta ~ C(condition)` to `theta ~ Length * Content` to properly model the 2x2 interaction. This will require updating the preprocessing logic to create a Length and Content column in the long DF.
    * Consider using `theta ~ Length * Content + (1 + Length * Content | subject_id)` to account for subject-level variation in the effects.
- [ ] Evaluate why we're doing global channel averaging in Step 1 of the `group_analysis_anova.py` script. This is most likely wrong.