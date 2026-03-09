"""Shared Homer FIR parsing, latent-HRF reconstruction, and AUC summarization.

This module is the single source of truth for:
  - parsing Homer FIR beta-table columns,
  - reconstructing latent HRFs from Gaussian basis weights, and
  - reducing reconstructed HRFs to a pre-registered task-window summary.
"""

from __future__ import annotations

import json
import math
import re
from hashlib import sha256
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path

import numpy as np


FIR_COL_PATTERN = re.compile(r"^(S\d+_D\d+)_Cond(\d{2})_(HbO|HbR)_Basis(\d{3})$")


@dataclass(frozen=True)
class ReconstructionSettings:
    """Settings for reconstructing Homer Gaussian-basis latent HRFs."""

    trange_start: float
    trange_stop: float
    basis_spacing: float
    basis_sigma: float

    @property
    def trange(self) -> tuple[float, float]:
        return (self.trange_start, self.trange_stop)


@dataclass(frozen=True)
class AUCSummarySettings:
    """Settings for baseline correction and task-window AUC summarization."""

    baseline_start: float
    baseline_stop: float
    auc_start: float
    auc_stop: float

    @property
    def baseline_window(self) -> tuple[float, float]:
        return (self.baseline_start, self.baseline_stop)

    @property
    def auc_window(self) -> tuple[float, float]:
        return (self.auc_start, self.auc_stop)


@dataclass(frozen=True)
class PreprocessingSettings:
    """Project-level settings for latent-HRF reconstruction and AUC reduction."""

    reconstruction: ReconstructionSettings
    auc_summary: AUCSummarySettings


@dataclass(frozen=True)
class FIRFeatureSpec:
    """Column specification for one channel/condition/chromophore FIR vector."""

    output_name: str
    channel: str
    condition: str
    chrom: str
    basis_idx: tuple[int, ...]
    column_idx: tuple[int, ...]

    def validate(self) -> None:
        if len(self.basis_idx) != len(self.column_idx):
            raise ValueError(
                f"Schema error for {self.output_name}: basis and column index counts differ."
            )
        expected_basis = tuple(range(1, len(self.basis_idx) + 1))
        if self.basis_idx != expected_basis:
            raise ValueError(
                "Expected contiguous basis numbering starting at 1 for exact Homer "
                f"reconstruction in {self.output_name}; got {list(self.basis_idx[:5])}..."
            )


def load_preprocessing_settings(path: str | Path) -> PreprocessingSettings:
    """Load and validate the shared preprocessing settings JSON."""
    path_obj = Path(path)
    raw = json.loads(path_obj.read_text(encoding="utf-8"))

    try:
        recon = raw["latent_hrf_reconstruction"]
        auc = raw["fir_auc_summary"]
    except KeyError as exc:
        raise ValueError(
            f"Missing required settings section in {path_obj}: {exc.args[0]}"
        ) from exc

    settings = PreprocessingSettings(
        reconstruction=ReconstructionSettings(
            trange_start=float(recon["trange_start"]),
            trange_stop=float(recon["trange_stop"]),
            basis_spacing=float(recon["basis_spacing"]),
            basis_sigma=float(recon["basis_sigma"]),
        ),
        auc_summary=AUCSummarySettings(
            baseline_start=float(auc["baseline_start"]),
            baseline_stop=float(auc["baseline_stop"]),
            auc_start=float(auc["auc_start"]),
            auc_stop=float(auc["auc_stop"]),
        ),
    )
    _validate_settings(settings, path_obj)
    return settings


def settings_to_dict(settings: PreprocessingSettings) -> dict[str, dict[str, float]]:
    """Convert preprocessing settings to a JSON-serializable dictionary."""
    return {
        "latent_hrf_reconstruction": {
            "trange_start": settings.reconstruction.trange_start,
            "trange_stop": settings.reconstruction.trange_stop,
            "basis_spacing": settings.reconstruction.basis_spacing,
            "basis_sigma": settings.reconstruction.basis_sigma,
        },
        "fir_auc_summary": {
            "baseline_start": settings.auc_summary.baseline_start,
            "baseline_stop": settings.auc_summary.baseline_stop,
            "auc_start": settings.auc_summary.auc_start,
            "auc_stop": settings.auc_summary.auc_stop,
        },
    }


def compute_file_sha256(path: str | Path) -> str:
    """Compute a stable SHA-256 hash for a file."""
    digest = sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def default_auc_provenance_path(auc_csv: str | Path) -> Path:
    """Return the default sidecar provenance path for a derived AUC CSV."""
    return Path(auc_csv).with_suffix(".provenance.json")


def _validate_settings(settings: PreprocessingSettings, path_obj: Path) -> None:
    recon = settings.reconstruction
    auc = settings.auc_summary
    if recon.trange_stop <= recon.trange_start:
        raise ValueError(f"Invalid reconstruction trange in {path_obj}: stop must be > start.")
    if recon.basis_spacing <= 0 or recon.basis_sigma <= 0:
        raise ValueError(
            f"Invalid reconstruction spacing/sigma in {path_obj}: values must be > 0."
        )
    if auc.baseline_stop < auc.baseline_start:
        raise ValueError(
            f"Invalid baseline window in {path_obj}: stop must be >= start."
        )
    if auc.auc_stop <= auc.auc_start:
        raise ValueError(f"Invalid AUC window in {path_obj}: stop must be > start.")
    if auc.baseline_start < recon.trange_start or auc.baseline_stop > recon.trange_stop:
        raise ValueError(
            f"Baseline window {auc.baseline_window} falls outside reconstruction trange "
            f"{recon.trange} in {path_obj}."
        )
    if auc.auc_start < recon.trange_start or auc.auc_stop > recon.trange_stop:
        raise ValueError(
            f"AUC window {auc.auc_window} falls outside reconstruction trange "
            f"{recon.trange} in {path_obj}."
        )


def _parse_float_or_nan(raw: str) -> float:
    """Parse a numeric string to float without coercing legitimate zeros."""
    try:
        val = float(raw)
    except (TypeError, ValueError):
        return float("nan")
    if math.isnan(val):
        return float("nan")
    return val


def normalize_excluded_channel_betas(beta: np.ndarray) -> np.ndarray:
    """Convert excluded channels to all-NaN and reject partial missingness.

    Project data-integrity policy for imported Homer betas:
    - excluded channels may appear as an all-zero basis-weight vector or an
      all-NaN basis-weight vector,
    - individual zero-valued basis weights remain valid coefficients, and
    - partial NaNs inside one basis vector are treated as a hard data error.
    """
    if np.isnan(beta).all():
        return np.full(beta.shape, np.nan, dtype=float)

    finite_mask = np.isfinite(beta)
    if finite_mask.all() and np.allclose(beta, 0.0):
        return np.full(beta.shape, np.nan, dtype=float)

    if np.isnan(beta).any():
        raise ValueError(
            "Unexpected partial NaNs within one channel's basis weights. "
            "Excluded channels must be entirely all-zero or all-NaN."
        )

    return beta


def parse_fir_header(header: list[str]) -> OrderedDict[str, FIRFeatureSpec]:
    """Parse FIR feature columns from a Homer wide table header."""
    if not header:
        raise ValueError("Input CSV is empty.")
    if header[0] != "Subject":
        raise ValueError("Expected first CSV column to be 'Subject'.")

    grouped: OrderedDict[str, dict[str, object]] = OrderedDict()
    for col_idx, col_name in enumerate(header[1:], start=1):
        match = FIR_COL_PATTERN.match(col_name)
        if not match:
            continue
        channel, cond, chrom, basis = match.groups()
        output_name = f"{channel}_Cond{cond}_{chrom}"
        if output_name not in grouped:
            grouped[output_name] = {
                "channel": channel,
                "condition": cond,
                "chrom": chrom,
                "basis_idx": [],
                "column_idx": [],
            }
        grouped[output_name]["basis_idx"].append(int(basis))
        grouped[output_name]["column_idx"].append(col_idx)

    if not grouped:
        raise ValueError(
            "No FIR columns matched pattern like 'S01_D01_Cond01_HbO_Basis001'."
        )

    specs: OrderedDict[str, FIRFeatureSpec] = OrderedDict()
    for output_name, meta in grouped.items():
        pairs = sorted(
            zip(meta["basis_idx"], meta["column_idx"], strict=True), key=lambda x: x[0]
        )
        spec = FIRFeatureSpec(
            output_name=output_name,
            channel=str(meta["channel"]),
            condition=str(meta["condition"]),
            chrom=str(meta["chrom"]),
            basis_idx=tuple(pair[0] for pair in pairs),
            column_idx=tuple(pair[1] for pair in pairs),
        )
        spec.validate()
        specs[output_name] = spec
    return specs


def build_condition_feature_map(
    specs: OrderedDict[str, FIRFeatureSpec],
    condition: str,
    target_channel: str | None = None,
) -> dict[str, dict[str, FIRFeatureSpec]]:
    """Build a channel/chrom map for one condition from parsed FIR specs."""
    condition_specs: dict[str, dict[str, FIRFeatureSpec]] = {}
    channels_in_condition: set[str] = set()

    for spec in specs.values():
        if spec.condition != condition:
            continue
        channels_in_condition.add(spec.channel)
        if target_channel is not None and spec.channel != target_channel:
            continue
        if spec.channel not in condition_specs:
            condition_specs[spec.channel] = {}
        condition_specs[spec.channel][spec.chrom] = spec

    if target_channel is not None and target_channel not in channels_in_condition:
        available = ", ".join(sorted(channels_in_condition)) or "none"
        raise ValueError(
            f"Requested target channel {target_channel} not found in Cond{condition}. "
            f"Available channels: {available}"
        )

    if not condition_specs:
        msg = f"No FIR columns found for condition Cond{condition}."
        if target_channel is not None:
            msg = (
                f"No FIR columns found for target channel {target_channel} in "
                f"condition Cond{condition}."
            )
        raise ValueError(msg)

    for channel, chrom_map in condition_specs.items():
        for chrom in ("HbO", "HbR"):
            if chrom not in chrom_map:
                raise ValueError(
                    f"Schema error: channel {channel} missing {chrom} columns in Cond{condition}."
                )

        hbo_basis = chrom_map["HbO"].basis_idx
        hbr_basis = chrom_map["HbR"].basis_idx
        if hbo_basis != hbr_basis:
            raise ValueError(
                f"Schema error: basis mismatch between HbO/HbR for channel {channel} in Cond{condition}."
            )

    return condition_specs


def extract_beta_vector(row: list[str], spec: FIRFeatureSpec) -> np.ndarray:
    """Extract one channel/condition/chrom FIR basis-weight vector from a CSV row."""
    values: list[float] = []
    for col_idx in spec.column_idx:
        if col_idx >= len(row):
            raise ValueError(
                f"Row is shorter than expected while reading {spec.output_name}; "
                f"expected column index {col_idx}, row length {len(row)}."
            )
        values.append(_parse_float_or_nan(row[col_idx]))
    return normalize_excluded_channel_betas(np.array(values, dtype=float))


class LatentHRFReconstructor:
    """Reconstruct latent HRFs and compute baseline-corrected AUC summaries.

    Citation: Ye et al. (2009), NeuroImage 44(2), 428-447.
    The fNIRS HRF can be represented as a weighted sum of temporal basis
    functions in a GLM. This reconstructor uses the exported Homer Gaussian
    basis weights to recover that latent HRF estimate.

    Citation: Yucel et al. (2021), Neurophotonics 8(1), 012101.
    Best-practice reporting for fNIRS requires explicit, auditable disclosure
    of preprocessing/modeling assumptions, which is why reconstruction and
    summary windows are centralized in a shared config and validated here.
    """

    def __init__(self, settings: ReconstructionSettings, auc_settings: AUCSummarySettings):
        self.settings = settings
        self.auc_settings = auc_settings
        self.time_axis, self.basis_matrix = self._build_gaussian_basis()
        self._validate_window_alignment(self.auc_settings.baseline_window, label="baseline")
        self._validate_window_alignment(self.auc_settings.auc_window, label="AUC")

    def _build_gaussian_basis(self) -> tuple[np.ndarray, np.ndarray]:
        start, stop = self.settings.trange
        spacing = self.settings.basis_spacing
        sigma = self.settings.basis_sigma

        t_hrf = np.arange(start, stop, spacing, dtype=float)
        n_basis = int(np.floor((stop - start) / spacing) - 1)
        if n_basis <= 0:
            raise ValueError(
                f"Invalid HRF basis configuration for trange={self.settings.trange}, "
                f"basis_spacing={spacing}."
            )

        tbasis = np.zeros((len(t_hrf), n_basis), dtype=float)
        for basis_number in range(1, n_basis + 1):
            center = start + basis_number * spacing
            gaussian = np.exp(-((t_hrf - center) ** 2) / (2 * sigma**2))
            tbasis[:, basis_number - 1] = gaussian / gaussian.max()
        return t_hrf, tbasis

    def _validate_window_alignment(self, window: tuple[float, float], *, label: str) -> None:
        start, stop = window
        if not np.any(np.isclose(self.time_axis, start)):
            raise ValueError(
                f"{label} window start {start} is not represented on the HRF time axis. "
                "Refusing to interpolate new samples."
            )
        if not np.any(np.isclose(self.time_axis, stop)):
            raise ValueError(
                f"{label} window stop {stop} is not represented on the HRF time axis. "
                "Refusing to interpolate new samples."
            )

    def reconstruct_from_beta_vector(
        self,
        basis_idx: tuple[int, ...] | list[int],
        beta: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Reconstruct the Homer latent HRF from Gaussian basis weights."""
        expected_basis = list(range(1, len(basis_idx) + 1))
        if list(basis_idx) != expected_basis:
            raise ValueError(
                "Expected contiguous basis numbering starting at 1 for exact "
                f"Homer reconstruction; got {list(basis_idx)[:5]}..."
            )
        if self.basis_matrix.shape[1] != len(beta):
            raise ValueError(
                f"Basis count mismatch: reconstructed basis matrix has {self.basis_matrix.shape[1]} "
                f"columns but beta vector has {len(beta)} values."
            )
        if np.isnan(beta).all():
            return self.time_axis, np.full(len(self.time_axis), np.nan, dtype=float)
        if np.isnan(beta).any():
            raise ValueError(
                "Cannot reconstruct HRF exactly because some Gaussian basis weights "
                "are missing/pruned. Refusing to substitute missing weights."
            )
        return self.time_axis, self.basis_matrix @ beta

    def summarize_auc(self, hrf: np.ndarray) -> float:
        """Compute baseline-corrected task-window AUC by trapezoidal integration.

        Citation: Pinti et al. (2024), Frontiers in Human Neuroscience 18:1418592.
        This implementation follows the common practice of summarizing an fNIRS
        hemodynamic trajectory over a prespecified time window with an area-
        under-the-curve feature, with all window choices kept explicit and
        auditable rather than data-driven.
        """
        if np.isnan(hrf).all():
            return float("nan")
        if np.isnan(hrf).any():
            raise ValueError("AUC summary received a partially missing HRF series.")

        baseline_mask = self._window_mask(self.auc_settings.baseline_window)
        auc_mask = self._window_mask(self.auc_settings.auc_window)
        baseline_value = float(np.mean(hrf[baseline_mask]))
        corrected = hrf - baseline_value
        return float(np.trapezoid(corrected[auc_mask], self.time_axis[auc_mask]))

    def _window_mask(self, window: tuple[float, float]) -> np.ndarray:
        start, stop = window
        mask = np.isclose(self.time_axis, start) | np.isclose(self.time_axis, stop)
        mask = mask | ((self.time_axis > start) & (self.time_axis < stop))
        if mask.sum() < 2:
            raise ValueError(
                f"Window {window} contains fewer than 2 HRF samples on the configured time axis."
            )
        return mask
