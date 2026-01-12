"""
Intention: This file is intended to generate an artifact `fnirs_qc.csv` that contains the QC metrics for each subject within the dataset.
The dataset is comprised of 51 subjects each of whom where exposed to 4 different conditions (short-form education, short-form entertainment, long-form education, long-form entertainment) four sperate times.
The study is a within-subjects design with 4 conditions each condition was presented 4 times.

The QC metrics are as follows:
- Subject ID
- Average SNR
- Average SCI
- List of Bad Channels (using an SCI cut off of 0.8) - These are channels that have an SCI value less than 0.8
- Number of Bad Channels - The number of channels that have an SCI value less than 0.8

The end result should be a CSV file where each participant has a row and each column is a QC metric.

The dataset can be found in `../fNIRs` and the output file should be saved in `../fnirs_analysis/fnirs_qc.csv`.
"""

import mne