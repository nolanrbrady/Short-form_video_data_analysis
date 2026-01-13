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
- Number of useable trials per condition (Short-Form Education, Short-Form Entertainment, Long-Form Entertainment, Long-Form Education)
  - A "run" is the whole fNIRS scan (one SNIRF file).
  - A "trial" is a 120s segment of the scan starting at a condition trigger annotation ("1".."4").
  - A trial is counted if its trial-window average SCI > 0.8.

The end result should be a CSV file where each participant has a row and each column is a QC metric.

The dataset can be found in `../fNIRs` and the output file should be saved in `../fnirs_analysis/fnirs_qc.csv`.

Author: Nolan Brady
Date: 2026-01-12
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

import mne
from mne.preprocessing import nirs as mne_nirs_preproc


# =======================
# User-configurable values
# =======================
#
# Per your request: all adjustable parameters live here as ALL CAPS constants,
# and the script avoids CLI arguments entirely.

# Root folder containing subject directories like `sub_0001/`, each containing one or more
# `.snirf` runs.
#
# Note: this is interpreted relative to the *current working directory* (same convention as
# `fnirs_analysis.py`). Most users run scripts from the Analysis repo root, where `../fNIRs`
# points to the sibling `fNIRs/` folder.
DATA_ROOT = Path("../fNIRs")

# Output is anchored to this file location (so it doesn't depend on CWD).
OUT_CSV = Path(__file__).resolve().parent / "fnirs_qc.csv"

# The requested SCI cutoff for "bad channels" (per-channel threshold on the per-channel SCI).
SCI_CUTOFF = 0.8

# Trial window duration in seconds (120s per your study design).
TRIAL_DURATION_SEC = 120.0

# Fail-fast safety: do not overwrite output unless explicitly allowed.
OVERWRITE_OUTPUT = True

# Mapping from raw trigger/annotation ID to condition name.
ANNOTATION_RENAME_MAP: Dict[str, str] = {
    "1": "Short-Form Education",
    "2": "Short-Form Entertainment",
    "3": "Long-Form Entertainment",
    "4": "Long-Form Education",
}


@dataclass(frozen=True)
class RunQC:
    """QC metrics computed for a single SNIRF run."""

    subject_id: str
    run_path: Path
    sci: np.ndarray  # per-channel SCI (NaNs allowed)
    snr_db: np.ndarray  # per-channel SNR in dB (NaNs allowed)
    ch_names: List[str]
    # Trial-window average SCI values and bad channel counts for each detected condition trigger.
    # Format: (cond_id, avg_sci, n_bad_channels)
    trial_qc: List[Tuple[str, float, int]]


def _resolve_data_root(user_root: Path) -> Path:
    """
    Resolve the fNIRs data root robustly.

    We try:
    1) The user-provided path resolved from the current working directory
    2) A fallback relative to this script location (common when CWD is `fnirs_analysis/`)
    """
    # 1) User path relative to CWD
    cand1 = user_root.expanduser().resolve()
    if cand1.exists():
        return cand1

    # 2) Fallback: repo layout has `Analysis/` and `fNIRs/` as siblings.
    # qc_check.py is at: <...>/Short Form Video Study Data/Analysis/fnirs_analysis/qc_check.py
    # parents[2] is:      <...>/Short Form Video Study Data
    cand2 = Path(__file__).resolve().parents[2] / "fNIRs"
    if cand2.exists():
        return cand2

    raise FileNotFoundError(
        "Could not locate fNIRs data root. Tried:\n"
        f"- {cand1}\n"
        f"- {cand2}\n"
        "Pass an explicit --data-root pointing to the folder that contains sub_#### directories."
    )


def _find_snirf_files(data_root: Path) -> List[Path]:
    """Find all SNIRF files under data_root."""
    snirf = sorted([p for p in data_root.rglob("*.snirf") if p.is_file()])
    snirf += sorted([p for p in data_root.rglob("*.SNIRF") if p.is_file()])
    # De-dupe in case both globs catch the same file (unlikely, but cheap insurance)
    return sorted(set(snirf))


def _subject_id_from_path(snirf_path: Path) -> str:
    """
    Extract subject ID from the SNIRF path.

    The dataset structure appears as: `.../fNIRs/sub_0001/<run>.snirf`
    so we use the immediate parent directory name.
    """
    return snirf_path.parent.name


def _compute_intensity_snr_db(raw_intensity: mne.io.BaseRaw) -> np.ndarray:
    """
    Compute a simple per-channel intensity SNR in dB.

    There is no single universal SNR definition for fNIRS. For *QC*, a robust and transparent
    summary is:

        SNR_dB = 20 * log10( |mean(intensity)| / std(intensity) )

    computed per channel on the raw intensity timeseries.

    Notes:
    - We avoid custom filtering here; this is intended as a high-level data-quality indicator.
    - Values are computed over the full recording (all samples).
    - Channels with zero/near-zero std are treated as NaN to avoid infinities.
    """
    # Only use fNIRS channels (exclude stim/aux channels if present in the file).
    ch_types = raw_intensity.get_channel_types()
    fnirs_picks = [i for i, t in enumerate(ch_types) if t.startswith("fnirs")]
    if len(fnirs_picks) == 0:
        raise ValueError("No fNIRS channels found in this SNIRF file (cannot compute SNR).")
    picks = fnirs_picks

    data = raw_intensity.get_data(picks=picks)  # shape: (n_channels, n_times)
    mu = np.nanmean(data, axis=1)
    sigma = np.nanstd(data, axis=1, ddof=0)

    # Guard against divide-by-zero / invalid channels
    eps = np.finfo(float).eps
    ratio = np.where(sigma > eps, np.abs(mu) / sigma, np.nan)
    return 20.0 * np.log10(ratio)


def _compute_run_qc(snirf_path: Path) -> RunQC:
    """
    Load a SNIRF run and compute per-channel SCI (on optical density) and SNR (on intensity).
    Also extracts event triggers.
    """
    raw_intensity = mne.io.read_raw_snirf(str(snirf_path), preload=True, verbose=False)

    # Only use fNIRS channels (exclude stim/aux channels if present in the file).
    ch_types = raw_intensity.get_channel_types()
    fnirs_picks = [i for i, t in enumerate(ch_types) if t.startswith("fnirs")]
    if len(fnirs_picks) == 0:
        raise ValueError(f"No fNIRS channels found in this SNIRF file: {snirf_path}")
    picks = fnirs_picks
    ch_names = list(np.array(raw_intensity.ch_names)[picks])

    # SCI is defined on optical density (OD) in MNE.
    raw_od = mne_nirs_preproc.optical_density(raw_intensity)
    sci = mne_nirs_preproc.scalp_coupling_index(raw_od)
    sci = np.asarray(sci, dtype=float)
    sci = sci[picks]

    # SNR on intensity (see function docstring for definition).
    snr_db = _compute_intensity_snr_db(raw_intensity)

    if len(ch_names) != len(sci) or len(ch_names) != len(snr_db):
        raise RuntimeError(
            "Channel count mismatch while computing QC metrics:\n"
            f"- n_ch={len(ch_names)}\n"
            f"- len(sci)={len(sci)}\n"
            f"- len(snr_db)={len(snr_db)}\n"
            f"File: {snirf_path}"
        )

    # Compute trial-window SCI: each trial is a 120s segment starting at a condition trigger.
    # Format: (cond_id, avg_sci, n_bad_channels)
    trial_qc: List[Tuple[str, float, int]] = []
    ann = raw_intensity.annotations
    if ann is not None and len(ann) > 0:
        for desc, onset in zip(ann.description, ann.onset):
            if desc not in ANNOTATION_RENAME_MAP:
                continue

            # Crop OD to the trial window and compute SCI on that segment.
            tmin = float(onset)
            tmax = float(onset) + float(TRIAL_DURATION_SEC)
            # Clamp to available data range; if a run is truncated, we still compute on what exists.
            raw_tmax = float(raw_od.times[-1]) if raw_od.n_times > 0 else 0.0
            tmax = min(tmax, raw_tmax)

            if tmax <= tmin:
                print(
                    "[qc_check] WARNING: trial window is empty after clamping; "
                    f"run={snirf_path} desc={desc} onset={onset} tmin={tmin} tmax={tmax}",
                    file=sys.stderr,
                )
                # Append NaN for SCI and 0 for bad channels (or maybe all channels? let's stick to 0 for now as it's unusable anyway)
                trial_qc.append((desc, float("nan"), 0))
                continue

            raw_trial = raw_od.copy().crop(tmin=tmin, tmax=tmax, include_tmax=False)
            sci_trial = mne_nirs_preproc.scalp_coupling_index(raw_trial)
            sci_trial = np.asarray(sci_trial, dtype=float)[picks]
            
            avg_sci = float(np.nanmean(sci_trial))
            n_bad = int(np.sum(sci_trial < SCI_CUTOFF))
            trial_qc.append((desc, avg_sci, n_bad))

    return RunQC(
        subject_id=_subject_id_from_path(snirf_path),
        run_path=snirf_path,
        sci=sci,
        snr_db=np.asarray(snr_db, dtype=float),
        ch_names=ch_names,
        trial_qc=trial_qc,
    )


def _aggregate_subject_qc(runs: Sequence[RunQC], sci_cutoff: float) -> dict:
    """
    Aggregate run-level QC into subject-level QC.

    We pool per-channel values across runs and compute the mean (ignoring NaNs). We also
    compute the union of "bad channels" (SCI < cutoff) across all runs.
    """
    if len(runs) == 0:
        raise ValueError("Expected at least one run for subject aggregation.")

    subject_id = runs[0].subject_id
    if any(r.subject_id != subject_id for r in runs):
        raise ValueError("Runs passed to _aggregate_subject_qc must all share the same subject_id.")

    # Fail fast: QC is only meaningful if channels are consistent across runs.
    # (Same channel set AND same order, so per-channel aggregation is unambiguous.)
    ref_ch_names = runs[0].ch_names
    for r in runs[1:]:
        if r.ch_names != ref_ch_names:
            raise ValueError(
                f"Channel layout mismatch within subject {subject_id}.\n"
                f"Reference run: {runs[0].run_path}\n"
                f"Mismatched run: {r.run_path}\n"
                "Ref channels (first 10): "
                + ", ".join(ref_ch_names[:10])
                + "\nOther channels (first 10): "
                + ", ".join(r.ch_names[:10])
            )

    # Aggregate per-channel across runs (still per-channel)
    # NOTE: no channel is labeled bad based on any across-channel average).
    sci_stack = np.vstack([r.sci for r in runs])  # (n_runs, n_channels)
    snr_stack = np.vstack([r.snr_db for r in runs])  # (n_runs, n_channels)

    sci_by_channel = np.nanmean(sci_stack, axis=0)
    snr_by_channel = np.nanmean(snr_stack, axis=0)

    # Subject-level summaries (requested outputs)
    avg_sci = float(np.nanmean(sci_by_channel)) if np.any(np.isfinite(sci_by_channel)) else float("nan")
    avg_snr_db = float(np.nanmean(snr_by_channel)) if np.any(np.isfinite(snr_by_channel)) else float("nan")

    # Bad channels are defined *per channel* using the per-channel SCI values.
    bad_channels = [ch for ch, sci in zip(ref_ch_names, sci_by_channel) if np.isfinite(sci) and sci < sci_cutoff]

    # Count useable trials per condition (trial-window average SCI > cutoff).
    valid_trials_counts = {name: 0 for name in ANNOTATION_RENAME_MAP.values()}

    n_channels = len(ref_ch_names)
    for r in runs:
        for cond_id, trial_sci, trial_bad_count in r.trial_qc:
            if cond_id not in ANNOTATION_RENAME_MAP:
                continue
            
            # Trial is useable if:
            # 1. Trial SCI >= cutoff
            # 2. Bad channel count <= 50%
            is_sci_ok = np.isfinite(trial_sci) and float(trial_sci) >= sci_cutoff
            is_ch_count_ok = trial_bad_count <= (n_channels / 2.0)
            
            if is_sci_ok and is_ch_count_ok:
                valid_trials_counts[ANNOTATION_RENAME_MAP[cond_id]] += 1

    # Subject-level exclusion flags
    pass_sci = avg_sci >= sci_cutoff
    pass_ch_count = len(bad_channels) <= (n_channels / 2.0)
    # Check if all conditions have at least 2 useable trials
    pass_trials = all(count >= 2 for count in valid_trials_counts.values())
    
    qc_dict = {
        "subject_id": subject_id,
        "avg_snr_db": avg_snr_db,
        "avg_sci": avg_sci,
        "bad_channels": ", ".join(bad_channels),
        "n_bad_channels": int(len(bad_channels)),
        "n_runs_found": int(len(runs)),
        "exclude_subject": not (pass_sci and pass_ch_count and pass_trials),
        "reason_for_exclusion": (
            ("" if pass_sci else "SCI_low; ") +
            ("" if pass_ch_count else "Too_many_bad_channels; ") +
            ("" if pass_trials else "Insufficient_useable_trials; ")
        ).strip("; "),
    }
    # Add the per-condition usable trial counts
    qc_dict.update(valid_trials_counts)
    return qc_dict


def build_qc_table(data_root: Path, sci_cutoff: float) -> pd.DataFrame:
    """
    Build a per-subject QC table by scanning SNIRF files under data_root.
    """
    # Compute run-level QC
    print(f"Finding SNIRF files in {data_root}...")
    snirf_files = _find_snirf_files(data_root)
    if not snirf_files:
        raise FileNotFoundError(f"No .snirf files found under {data_root}")

    print(f"Found {len(snirf_files)} SNIRF file(s). Computing run-level QC...")
    run_qc: List[RunQC] = []
    for i, p in enumerate(snirf_files, 1):
        print(f"[{i}/{len(snirf_files)}] Computing QC for: {p.name}")
        run_qc.append(_compute_run_qc(p))

    # Group by subject_id
    by_subj: dict[str, List[RunQC]] = {}
    for r in run_qc:
        by_subj.setdefault(r.subject_id, []).append(r)

    print(f"Aggregating QC for {len(by_subj)} subject(s)...")
    rows = []
    for subject_id in sorted(by_subj.keys()):
        print(f"  -> Processing subject: {subject_id}")
        rows.append(_aggregate_subject_qc(by_subj[subject_id], sci_cutoff=sci_cutoff))

    df = pd.DataFrame(rows)

    # Sort so the lowest-quality subjects appear first:
    # - Lower avg_sci is worse (ascending)
    # - More bad channels is worse (descending)
    df = df.sort_values(
        by=["avg_sci", "n_bad_channels", "subject_id"],
        ascending=[True, False, True],
        na_position="last",
    ).reset_index(drop=True)
    return df


def main() -> int:
    try:
        print("Starting fNIRS QC Analysis...")
        data_root = _resolve_data_root(DATA_ROOT)
        df = build_qc_table(data_root=data_root, sci_cutoff=float(SCI_CUTOFF))
        out_path = OUT_CSV.expanduser().resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        if out_path.exists() and not OVERWRITE_OUTPUT:
            raise FileExistsError(
                f"Refusing to overwrite existing output file: {out_path}\n"
                "Set OVERWRITE_OUTPUT = True at the top of qc_check.py if you intend to replace it."
            )
        df.to_csv(out_path, index=False)
        print(f"Wrote QC table with {len(df)} subject(s) to: {out_path}")
        return 0
    except Exception as e:
        print(f"[qc_check] ERROR: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())