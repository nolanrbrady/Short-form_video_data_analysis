"""Validate the FIR-to-AUC conversion against excluded-channel policy.

Validation defaults remain explicit in this file, while CLI path overrides let
the pipeline configure raw/AUC/settings locations from one place.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from homer_fir import (
    compute_file_sha256,
    default_auc_provenance_path,
    extract_beta_vector,
    load_preprocessing_settings,
    parse_fir_header,
    settings_to_dict,
)


RAW_FIR_CSV = "data/tabular/homer3_glm_betas_wide_fir_pca.csv"
AUC_CSV = "data/tabular/generated_data/homer3_glm_betas_wide_auc.csv"
SETTINGS_JSON = "data/config/preprocessing_settings.json"
BASIS_FAMILY = "Homer3 idxBasis=1 Gaussian"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate that the FIR-to-AUC conversion preserves excluded-channel policy."
    )
    parser.add_argument(
        "--raw-fir-csv",
        default=RAW_FIR_CSV,
        help="Path to the raw Homer FIR basis-weight CSV.",
    )
    parser.add_argument(
        "--auc-csv",
        default=AUC_CSV,
        help="Path to the derived single-beta AUC CSV.",
    )
    parser.add_argument(
        "--settings-json",
        default=SETTINGS_JSON,
        help="Path to the shared preprocessing settings JSON.",
    )
    parser.add_argument(
        "--provenance-json",
        default=None,
        help="Optional override for the FIR-to-AUC provenance JSON path.",
    )
    return parser.parse_args()


def _parse_auc_value(raw: str) -> float:
    try:
        value = float(raw)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Encountered non-numeric AUC value token: {raw!r}") from exc
    if np.isnan(value):
        return float("nan")
    return value


def validate_homer_fir_auc_conversion(
    raw_fir_csv: str = RAW_FIR_CSV,
    auc_csv: str = AUC_CSV,
    settings_json: str = SETTINGS_JSON,
    provenance_json: str | None = None,
) -> None:
    """Fail hard if excluded FIR vectors are not represented correctly in AUC output.

    Citation: Sandve et al. (2013), PLoS Computational Biology 9(10), e1003285.
    Reproducible computational research benefits from explicit, automated checks
    on derived artifacts so that critical data-integrity invariants fail loudly
    rather than relying on manual inspection.
    """
    raw_path = Path(raw_fir_csv)
    auc_path = Path(auc_csv)
    settings_path = Path(settings_json)
    provenance_path = (
        Path(provenance_json) if provenance_json is not None else default_auc_provenance_path(auc_path)
    )
    if not raw_path.exists():
        raise FileNotFoundError(f"Raw FIR CSV not found: {raw_path}")
    if not auc_path.exists():
        raise FileNotFoundError(f"AUC CSV not found: {auc_path}")
    if not settings_path.exists():
        raise FileNotFoundError(f"Settings JSON not found: {settings_path}")
    if not provenance_path.exists():
        raise FileNotFoundError(f"FIR-to-AUC provenance JSON not found: {provenance_path}")

    settings = load_preprocessing_settings(settings_path)
    expected_provenance = {
        "analysis_step": "collapse_homer_fir_to_auc",
        "basis_family": BASIS_FAMILY,
        "input_csv": str(raw_path),
        "input_csv_sha256": compute_file_sha256(raw_path),
        "settings_json": str(settings_path),
        "settings_json_sha256": compute_file_sha256(settings_path),
        "output_csv": str(auc_path),
        "output_csv_sha256": compute_file_sha256(auc_path),
        "settings": settings_to_dict(settings),
    }
    payload = json.loads(provenance_path.read_text(encoding="utf-8"))
    for key, expected_value in expected_provenance.items():
        observed_value = payload.get(key)
        if observed_value != expected_value:
            raise ValueError(
                f"Provenance mismatch for {key}: expected {expected_value!r}, "
                f"found {observed_value!r} in {provenance_path}."
            )

    with raw_path.open("r", newline="") as raw_f, auc_path.open("r", newline="") as auc_f:
        raw_reader = csv.reader(raw_f)
        auc_reader = csv.reader(auc_f)
        raw_header = next(raw_reader, None)
        auc_header = next(auc_reader, None)
        if raw_header is None:
            raise ValueError("Raw FIR CSV is empty.")
        if auc_header is None:
            raise ValueError("AUC CSV is empty.")
        if raw_header[0] != "Subject" or auc_header[0] != "Subject":
            raise ValueError("Expected first column to be 'Subject' in both CSVs.")

        specs = parse_fir_header(raw_header)
        if payload.get("feature_count") != len(specs):
            raise ValueError(
                f"Provenance feature_count mismatch: expected {len(specs)}, "
                f"found {payload.get('feature_count')!r}."
            )
        expected_auc_header = ["Subject", *specs.keys()]
        if auc_header != expected_auc_header:
            raise ValueError(
                "AUC CSV header does not match the FIR-derived feature order. "
                "Expected columns derived from the raw FIR schema."
            )
        basis_count = len(next(iter(specs.values())).basis_idx)
        if payload.get("basis_count_per_feature") != basis_count:
            raise ValueError(
                f"Provenance basis_count_per_feature mismatch: expected {basis_count}, "
                f"found {payload.get('basis_count_per_feature')!r}."
            )

        row_count = 0
        for raw_row_number, (raw_row, auc_row) in enumerate(
            zip(raw_reader, auc_reader, strict=True), start=2
        ):
            if not raw_row or not auc_row:
                raise ValueError(
                    "Encountered blank row while validating FIR-to-AUC conversion; "
                    "raw and AUC files must be row-aligned and non-empty."
                )

            raw_subject = raw_row[0].strip()
            auc_subject = auc_row[0].strip()
            if raw_subject != auc_subject:
                raise ValueError(
                    "Subject mismatch between raw FIR and AUC CSVs at row "
                    f"{raw_row_number}: raw={raw_subject!r}, auc={auc_subject!r}."
                )

            for auc_col_index, spec in enumerate(specs.values(), start=1):
                beta = extract_beta_vector(raw_row, spec)
                source_is_excluded = np.isnan(beta).all()
                auc_value = _parse_auc_value(auc_row[auc_col_index])
                auc_is_nan = np.isnan(auc_value)

                if source_is_excluded and not auc_is_nan:
                    raise ValueError(
                        f"Excluded FIR basis vector for subject {raw_subject}, "
                        f"{spec.output_name} should map to NaN in the AUC CSV "
                        f"but found {auc_row[auc_col_index]!r}."
                    )
                if auc_is_nan and not source_is_excluded:
                    raise ValueError(
                        f"AUC CSV contains NaN for subject {raw_subject}, {spec.output_name} "
                        "even though the source FIR basis vector was not an all-zero/all-NaN "
                        "excluded channel."
                    )
            row_count += 1

        if row_count == 0:
            raise ValueError("Validation found no data rows in the FIR/AUC CSVs.")
        if payload.get("row_count") != row_count:
            raise ValueError(
                f"Provenance row_count mismatch: expected {row_count}, "
                f"found {payload.get('row_count')!r}."
            )


def main() -> None:
    args = parse_args()
    validate_homer_fir_auc_conversion(
        raw_fir_csv=args.raw_fir_csv,
        auc_csv=args.auc_csv,
        settings_json=args.settings_json,
        provenance_json=args.provenance_json,
    )
    print("[PASS] FIR-to-AUC conversion lint passed")


if __name__ == "__main__":
    main()
