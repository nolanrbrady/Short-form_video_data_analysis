"""Plot reconstructed Homer3 FIR HRFs for selected subjects.

This script intentionally uses top-of-file configuration variables (no CLI).
It streams the CSV and only extracts requested subjects plus requested condition
columns to avoid loading the full wide table into memory.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from homer_fir import (
    AUCSummarySettings,
    LatentHRFReconstructor,
    ReconstructionSettings,
    build_condition_feature_map,
    extract_beta_vector,
    load_preprocessing_settings,
    parse_fir_header,
)

# -----------------------------------------------------------------------------
# User-editable configuration (no CLI by request)
# -----------------------------------------------------------------------------

INPUT_CSV = "data/tabular/homer3_glm_betas_wide_fir_pca.csv"
OUTPUT_DIR = "data/results/fir_beta_plots"
SETTINGS_JSON = "data/config/preprocessing_settings.json"
TARGET_SUBJECTS = ["sub_0001", "sub_0002"]
TARGET_CONDITION = "02"
TARGET_CHANNEL: str | None = "S01_D01"  # e.g., "S01_D01"
DPI = 300


def parse_condition_columns(
    header: list[str], condition: str, target_channel: str | None = None
) -> dict[str, dict[str, object]]:
    """Build a channel/chrom map from header columns for one condition."""
    specs = parse_fir_header(header)
    return build_condition_feature_map(specs, condition, target_channel=target_channel)


def extract_target_rows(
    csv_path: Path, target_subjects: Iterable[str]
) -> tuple[list[str], dict[str, list[str]]]:
    """Read header and only rows for requested subjects (streaming)."""
    target_set = set(target_subjects)
    if not target_set:
        raise ValueError("No target subjects configured.")

    rows_by_subject: dict[str, list[str]] = {}

    with csv_path.open("r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if not header:
            raise ValueError("Input CSV is empty.")
        if header[0] != "Subject":
            raise ValueError("Expected first CSV column to be 'Subject'.")

        for row in reader:
            if not row:
                continue
            subject = row[0].strip()
            if subject in target_set:
                if subject in rows_by_subject:
                    raise ValueError(
                        f"Duplicate subject row detected for {subject}; expected one row per subject."
                    )
                rows_by_subject[subject] = row
                if len(rows_by_subject) == len(target_set):
                    break

    missing = sorted(target_set - set(rows_by_subject))
    if missing:
        raise ValueError(f"Missing requested subject(s): {', '.join(missing)}")

    return header, rows_by_subject


def extract_subject_betas(
    row: list[str], cond_col_map: dict[str, dict[str, object]]
) -> tuple[list[int], dict[str, dict[str, np.ndarray]]]:
    """Extract per-channel Gaussian basis weights for HbO/HbR from one row."""
    channels = sorted(cond_col_map.keys())
    first_channel = channels[0]
    basis_idx = list(cond_col_map[first_channel]["HbO"].basis_idx)

    betas: dict[str, dict[str, np.ndarray]] = {}
    for channel in channels:
        betas[channel] = {}
        for chrom in ("HbO", "HbR"):
            betas[channel][chrom] = extract_beta_vector(row, cond_col_map[channel][chrom])
    return basis_idx, betas


def build_gaussian_tbasis(
    trange: tuple[float, float], basis_spacing: float, basis_sigma: float
) -> tuple[np.ndarray, np.ndarray]:
    """Recreate Homer3's idxBasis=1 Gaussian basis matrix."""
    reconstructor = LatentHRFReconstructor(
        ReconstructionSettings(
            trange_start=trange[0],
            trange_stop=trange[1],
            basis_spacing=basis_spacing,
            basis_sigma=basis_sigma,
        ),
        AUCSummarySettings(
            baseline_start=trange[0],
            baseline_stop=trange[0],
            auc_start=trange[0],
            auc_stop=trange[0] + basis_spacing,
        ),
    )
    return reconstructor.time_axis, reconstructor.basis_matrix


def reconstruct_hrf_from_betas(
    basis_idx: list[int],
    beta: np.ndarray,
    trange: tuple[float, float],
    basis_spacing: float,
    basis_sigma: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Reconstruct the Homer idxBasis=1 HRF from Gaussian basis weights."""
    reconstructor = LatentHRFReconstructor(
        ReconstructionSettings(
            trange_start=trange[0],
            trange_stop=trange[1],
            basis_spacing=basis_spacing,
            basis_sigma=basis_sigma,
        ),
        AUCSummarySettings(
            baseline_start=trange[0],
            baseline_stop=trange[0],
            auc_start=trange[0],
            auc_stop=trange[0] + basis_spacing,
        ),
    )
    return reconstructor.reconstruct_from_beta_vector(tuple(basis_idx), beta)


def build_subject_hrf_series(
    row: list[str],
    cond_col_map: dict[str, dict[str, object]],
    reconstructor: LatentHRFReconstructor,
) -> tuple[np.ndarray, dict[str, dict[str, np.ndarray]]]:
    """Extract Gaussian basis weights and reconstruct HRFs for one subject."""
    basis_idx, betas = extract_subject_betas(row, cond_col_map)
    t_hrf: np.ndarray | None = None
    series: dict[str, dict[str, np.ndarray]] = {}

    for channel, chrom_map in betas.items():
        series[channel] = {}
        for chrom, beta in chrom_map.items():
            t_curr, hrf = reconstructor.reconstruct_from_beta_vector(basis_idx, beta)
            if t_hrf is None:
                t_hrf = t_curr
            series[channel][chrom] = hrf

    if t_hrf is None:
        raise ValueError("No HRF series reconstructed.")

    return t_hrf, series


def plot_subject(
    subject: str,
    condition: str,
    t_hrf: np.ndarray,
    series: dict[str, dict[str, np.ndarray]],
    out_dir: Path,
    dpi: int,
) -> Path:
    """Write a single subject figure for one condition (HbO+HbR HRFs)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(14, 7))

    channels = sorted(series.keys())
    cmap = plt.get_cmap("tab20")
    single_channel = len(channels) == 1

    channel_handles: list[Line2D] = []
    x = np.asarray(t_hrf, dtype=float)

    for i, channel in enumerate(channels):
        color = cmap(i % 20)
        y_hbo = series[channel]["HbO"]
        y_hbr = series[channel]["HbR"]

        ax.plot(x, y_hbo, color=color, linestyle="-", linewidth=1.1, alpha=0.85)
        ax.plot(x, y_hbr, color=color, linestyle="--", linewidth=1.1, alpha=0.85)
        channel_handles.append(Line2D([0], [0], color=color, linestyle="-", label=channel))

    style_handles = [
        Line2D([0], [0], color="black", linestyle="-", linewidth=1.6, label="HbO"),
        Line2D([0], [0], color="black", linestyle="--", linewidth=1.6, label="HbR"),
    ]

    if single_channel:
        ax.set_title(
            f"{subject} | Cond{condition} | {channels[0]} Reconstructed HRF (HbO and HbR)"
        )
    else:
        ax.set_title(f"{subject} | Cond{condition} Reconstructed HRF (HbO and HbR)")
    ax.set_xlabel("Time From Onset (s)")
    ax.set_ylabel("HRF")
    ax.grid(True, alpha=0.25)

    style_legend = ax.legend(handles=style_handles, title="Chromophore", loc="upper right")
    ax.add_artist(style_legend)
    if not single_channel:
        ax.legend(
            handles=channel_handles,
            title="Channel",
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            ncol=1,
            fontsize=8,
            title_fontsize=9,
            frameon=False,
        )
        fig.tight_layout(rect=(0, 0, 0.82, 1))
        out_path = out_dir / f"{subject}_cond{condition}_fir_hrf_hbo_hbr.png"
    else:
        fig.tight_layout()
        out_path = out_dir / f"{subject}_cond{condition}_{channels[0]}_fir_hrf_hbo_hbr.png"
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_path


def run_plotting(
    input_csv: str,
    output_dir: str,
    target_subjects: list[str],
    target_condition: str,
    target_channel: str | None = None,
    settings_json: str = SETTINGS_JSON,
    dpi: int = 300,
) -> list[Path]:
    """End-to-end utility used by main() and tests."""
    if target_channel is None:
        raise ValueError(
            "target_channel must be explicitly set (e.g., 'S01_D01'); "
            "refusing to plot all channels by default."
        )

    csv_path = Path(input_csv)
    out_dir = Path(output_dir)
    if not csv_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {csv_path}")

    settings = load_preprocessing_settings(settings_json)
    reconstructor = LatentHRFReconstructor(settings.reconstruction, settings.auc_summary)
    header, rows_by_subject = extract_target_rows(csv_path, target_subjects)
    cond_col_map = parse_condition_columns(
        header, target_condition, target_channel=target_channel
    )

    out_paths: list[Path] = []
    for subject in target_subjects:
        t_hrf, series = build_subject_hrf_series(
            rows_by_subject[subject], cond_col_map, reconstructor
        )
        out_path = plot_subject(subject, target_condition, t_hrf, series, out_dir, dpi)
        out_paths.append(out_path)

    return out_paths


def main() -> None:
    output_paths = run_plotting(
        input_csv=INPUT_CSV,
        output_dir=OUTPUT_DIR,
        target_subjects=TARGET_SUBJECTS,
        target_condition=TARGET_CONDITION,
        target_channel=TARGET_CHANNEL,
        settings_json=SETTINGS_JSON,
        dpi=DPI,
    )
    print(f"[INFO] Wrote {len(output_paths)} plot(s):")
    for p in output_paths:
        print(f"  - {p}")


if __name__ == "__main__":
    main()
