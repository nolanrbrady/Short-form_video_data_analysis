#!/usr/bin/env python3
"""
Subject-level fNIRS GLM pipeline (SNIRF → HbO/HbR preprocessing → GLM) using MNE + MNE-NIRS.

Docs referenced:
- read_raw_snirf, OD/MBLL, filtering, SCI, TDDR, GLM examples:
  https://mne.tools/stable/ and https://mne.tools/mne-nirs/stable/

Author: Nolan Brady
"""

from __future__ import annotations
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import mne
from mne.preprocessing import nirs as mne_preprocess
from mne_nirs.experimental_design import make_first_level_design_matrix
from mne_nirs.statistics import run_glm
from mne_nirs.channels import get_short_channels, get_long_channels

# =========================
# User-configurable settings
# =========================

# Paths
DATA_ROOT: str = "../fNIRs"          # root folder to search (recursive)
OUTPUT_ROOT: str = "./glm_results"    # where results will be saved

# Subject ID extraction
# If not None, we will search the file path with this regex and use the first capture group.
# Example patterns: r"(sub-[A-Za-z0-9_]+)"  or  r"/([A-Za-z0-9_-]+)/[^/]*\.snirf$"
SUBJECT_ID_REGEX: Optional[str] = r"(sub-[A-Za-z0-9_]+)"

# Event labeling (use underscores, not slashes, per Nilearn quirk)
# Map your device-specific annotation/trigger names to four condition names.
ANNOTATION_RENAME_MAP: Dict[str, str] = {
    "1": "Short-Form Education",
    "2": "Short-Form Entertainment",
    "3": "Long-Form Entertainment",
    "4": "Long-Form Education",
}
STIMULUS_DURATION_SEC: float = 120.0  # set per your design (uniform duration assumed here)

# Preprocessing toggles
APPLY_TDDR: bool = True                  # Spike/baseline-motion artifact correction (TDDR) on OD
APPLY_FILTER: bool = True                # Band-pass filter on haemoglobin

# A low‑pass in the 0.2–0.3 Hz range is common to suppress respiration (~0.2–0.3 Hz)
LOW_PASS_HZ: float = 0.30                # remove cardiac/respiration

# NOTE: High-pass filter is removed in favor of the GLM drift model.
# The fundamental task frequency is ≈ 1 / (120 + 30) = 0.0067 Hz. A standard high‑pass of 0.01 Hz would attenuate or distort this slow block effect
HIGH_PASS_HZ: float = None               # remove slow drift (was 0.003 Hz)
RESAMPLE_HZ: Optional[float] = None      # e.g., 1.0 for speed; None to keep native

# Beer–Lambert Law
# Wavelength-dependent partial pathlength factor (PPF). If set to a value, this will
# override MNE's default literature-based wavelength-dependent values.
# Set to None to use MNE's recommended defaults (Scholkmann & Wolf, 2013).
BEER_LAMBERT_PPF: Optional[float | Dict[str, float]] = None  # e.g., 6.0; None is recommended

# Optional: drop channels with poor coupling (SCI threshold). None to skip.
# SCI itself is introduced in Pollonini et al. (2016) (PHOEBE), while an ~0.8 “good coupling” reference is used in Hernandez & Pollonini (2020) (NIRSplot).
SCALP_COUPLING_MIN: Optional[float] = 0.8

# GLM design settings
# Adding the derivative makes the amplitude estimate more robust to timing variation.
HRF_MODEL: str = "spm + derivative"                    # 'spm', 'spm + derivative', or 'glover', etc.
DRIFT_MODEL: str = "cosine"               # cosine drift regression
# High-pass cutoff for the GLM drift model. Should be lower than task-related frequencies.
# Rule of thumb: 1 / (2 * (task_duration + inter_stimulus_interval)).
# For 120s blocks, a safe value is e.g., 1/250s = 0.004 Hz.
HIGH_PASS_HZ_DESIGN: float = 0.0033        # choose ~ 1/(2 * max inter-trial interval (task length + rest length)); tune per experiment
NOISE_MODEL: str = "ar1"                  # Nilearn AR(1) noise model
INCLUDE_SHORT_CHANNEL_REGRESSORS: bool = False  # We did not use short channels in this study.

# File discovery
SNIRF_EXTENSIONS: Tuple[str, ...] = (".snirf", ".SNIRF")
IGNORE_SUBJECTS = []

# ==============
# Helper methods
# ==============

def find_snirf_files(root: str) -> List[Path]:
    #TODO: Implement a check for IGNORE_SUBJECTS and make sure they are not included in the SNIRF files listed
    return [p for p in Path(root).rglob("*") if p.suffix in SNIRF_EXTENSIONS]

def extract_subject_id(snirf_path: Path) -> str:
    if SUBJECT_ID_REGEX:
        m = re.search(SUBJECT_ID_REGEX, str(snirf_path))
        if m:
            return m.group(1)
    # fallback: use immediate parent folder name
    return snirf_path.parent.name

def rename_and_set_durations(raw: mne.io.BaseRaw) -> None:
    if ANNOTATION_RENAME_MAP:
        raw.annotations.rename(ANNOTATION_RENAME_MAP)
    if STIMULUS_DURATION_SEC is not None:
        raw.annotations.set_durations(STIMULUS_DURATION_SEC)

def compute_and_drop_poor_sci(raw_od: mne.io.BaseRaw, threshold: float) -> mne.io.BaseRaw:
    """
    Compute Scalp Coupling Index (SCI) on OD and mark channels below threshold as bad.
    """
    sci_vals = mne_preprocess.scalp_coupling_index(raw_od)
    bads = [ch for ch, sci in zip(raw_od.ch_names, sci_vals) if np.isfinite(sci) and sci < threshold]
    if bads:
        # Combine with any existing bad channels and update info
        existing_bads = list(raw_od.info.get("bads", []))
        new_bads = list(sorted(set(existing_bads + bads)))
        raw_od.info["bads"] = new_bads
        print(f"    Marked {len(bads)} channels as bad due to SCI < {threshold}.")
    # Subsequent MNE functions will ignore channels in info["bads"]
    return raw_od

def preprocess_intensity_to_haemoglobin(raw_intensity: mne.io.BaseRaw) -> mne.io.BaseRaw:
    """
    OD -> (optional SCI) -> (optional TDDR) -> BLL -> (optional filter) -> (optional resample)
    Returns haemoglobin Raw.
    """
    # 1. Convert to optical density
    raw_od = mne_preprocess.optical_density(raw_intensity)

    # 2. Optional: Prune channels with poor scalp coupling
    if SCALP_COUPLING_MIN is not None:
        raw_od = compute_and_drop_poor_sci(raw_od, SCALP_COUPLING_MIN)

    # 3. Optional: Correct motion artifacts on OD data
    if APPLY_TDDR:
        raw_od = mne_preprocess.temporal_derivative_distribution_repair(raw_od)

    # 4. Convert to haemoglobin concentration
    # Defaults to MNE's dynamic wavelength-dependent PPF values 
    ppf_arg = {}
    if BEER_LAMBERT_PPF is not None:
        ppf_arg["ppf"] = BEER_LAMBERT_PPF
    raw_haemo = mne_preprocess.beer_lambert_law(raw_od, **ppf_arg)

    # 5. Optional: Band-pass filter the haemoglobin data
    # This is done after HbO conversion, which is the standard practice.
    if APPLY_FILTER:
        raw_haemo = raw_haemo.copy().filter(
            l_freq=HIGH_PASS_HZ, h_freq=LOW_PASS_HZ, method="fir", phase="zero", fir_window="hamming"
        )

    # 6. Optional: Resample data (after filtering to prevent aliasing)
    if RESAMPLE_HZ is not None:
        raw_haemo = raw_haemo.copy().resample(RESAMPLE_HZ, npad="auto")

    return raw_haemo


def fit_glm_and_save(
    subject_id: str,
    run_label: str,
    raw_long: mne.io.BaseRaw,
    design_matrix: pd.DataFrame,
    out_root: Path
) -> Tuple[Path, Path]:
    #TODO: Explicitly remove the channels marked as bad prior to GLM
    # Sometimes MNE doesn't remove the bad channels automatically during the GLM step.
    glm_est = run_glm(raw_long, design_matrix, noise_model=NOISE_MODEL)

    # Save full GLM object (HDF5)
    subj_dir = out_root / subject_id
    subj_dir.mkdir(parents=True, exist_ok=True)
    glm_path = subj_dir / f"{subject_id}_{run_label}_glm.h5"
    glm_est.save(str(glm_path), overwrite=True)

    # Save tidy dataframe for quick inspection / downstream stats
    df = glm_est.to_dataframe()
    df["subject_id"] = subject_id
    df["run"] = run_label
    csv_path = subj_dir / f"{subject_id}_{run_label}_glm_results.csv"
    df.to_csv(csv_path, index=False)

    # Also append/update a subject-level aggregate CSV
    agg_path = subj_dir / f"{subject_id}_ALL_runs_glm_results.csv"
    if agg_path.exists():
        old = pd.read_csv(agg_path)
        new = pd.concat([old, df], axis=0, ignore_index=True)
    else:
        new = df
    new.to_csv(agg_path, index=False)

    return glm_path, csv_path

# ==========
# Main logic
# ==========

def main():
    data_root = Path(DATA_ROOT).expanduser().resolve()
    out_root = Path(OUTPUT_ROOT).expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    snirf_files = find_snirf_files(str(data_root))
    if not snirf_files:
        print(f"No SNIRF files found under {data_root}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(snirf_files)} SNIRF file(s).")

    for snirf_path in sorted(snirf_files):
        subject_id = extract_subject_id(snirf_path)
        run_label = snirf_path.stem  # file name without extension

        print(f"\n>>> Processing {snirf_path} (subject={subject_id}, run={run_label})")

        # 1) Load raw intensity and harmonize annotations
        raw_intensity = mne.io.read_raw_snirf(str(snirf_path), preload=True, verbose=False)
        rename_and_set_durations(raw_intensity)

        # 2) Preprocess to haemoglobin
        raw_haemo = preprocess_intensity_to_haemoglobin(raw_intensity)

        # 3) Split long/short channels for GLM; we analyze long channels
        raw_long = get_long_channels(raw_haemo)

        if raw_long.n_times == 0 or len(raw_long.ch_names) == 0:
            print("No long channels available after preprocessing; skipping.", file=sys.stderr)
            continue

        # 4) Build first-level design matrix
        design_matrix = make_first_level_design_matrix(
            raw_haemo,
            drift_model=DRIFT_MODEL,
            high_pass=HIGH_PASS_HZ_DESIGN,
            hrf_model=HRF_MODEL,
            stim_dur=STIMULUS_DURATION_SEC,
        )

        # 5) Fit GLM and save
        glm_path, csv_path = fit_glm_and_save(
            subject_id=subject_id,
            run_label=run_label,
            raw_long=raw_long,
            design_matrix=design_matrix,
            out_root=out_root
        )

        print(f"Saved GLM:  {glm_path}")
        print(f"Saved CSV:  {csv_path}")

    print("\nDone.")

if __name__ == "__main__":
    main()
