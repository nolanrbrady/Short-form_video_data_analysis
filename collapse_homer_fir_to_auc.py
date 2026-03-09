"""Collapse Homer FIR basis-weight exports to single beta values via AUC.

This script intentionally uses top-of-file configuration variables instead of a
CLI so the methodological assumptions remain easy to inspect and edit.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from homer_fir import (
    LatentHRFReconstructor,
    compute_file_sha256,
    default_auc_provenance_path,
    extract_beta_vector,
    load_preprocessing_settings,
    parse_fir_header,
    settings_to_dict,
)


INPUT_CSV = "data/tabular/homer3_glm_betas_wide_fir.csv"
OUTPUT_CSV = "data/tabular/generated_data/homer3_glm_betas_wide_auc.csv"
SETTINGS_JSON = "data/config/preprocessing_settings.json"
BASIS_FAMILY = "Homer3 idxBasis=1 Gaussian"


def _format_output_value(value: float) -> str:
    if np.isnan(value):
        return "NaN"
    return f"{value:.17g}"


def collapse_homer_fir_to_auc(
    input_csv: str = INPUT_CSV,
    output_csv: str = OUTPUT_CSV,
    settings_json: str = SETTINGS_JSON,
) -> Path:
    """Convert FIR basis-weight columns into single beta columns via AUC."""
    settings = load_preprocessing_settings(settings_json)
    reconstructor = LatentHRFReconstructor(settings.reconstruction, settings.auc_summary)

    input_path = Path(input_csv)
    output_path = Path(output_csv)
    settings_path = Path(settings_json)
    provenance_path = default_auc_provenance_path(output_path)
    if not input_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_path}")
    if not settings_path.exists():
        raise FileNotFoundError(f"Settings JSON not found: {settings_path}")

    with input_path.open("r", newline="") as in_f:
        reader = csv.reader(in_f)
        header = next(reader, None)
        if header is None:
            raise ValueError("Input CSV is empty.")
        specs = parse_fir_header(header)
        output_header = ["Subject", *specs.keys()]
        row_count = 0

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", newline="") as out_f:
            writer = csv.writer(out_f)
            writer.writerow(output_header)

            seen_subjects: set[str] = set()
            for row_number, row in enumerate(reader, start=2):
                if not row:
                    continue
                subject = row[0].strip()
                if not subject:
                    raise ValueError(f"Missing Subject value at input row {row_number}.")
                if subject in seen_subjects:
                    raise ValueError(
                        f"Duplicate subject row detected for {subject}; expected one row per subject."
                    )
                seen_subjects.add(subject)

                out_row = [subject]
                for spec in specs.values():
                    beta = extract_beta_vector(row, spec)
                    _, hrf = reconstructor.reconstruct_from_beta_vector(spec.basis_idx, beta)
                    auc = reconstructor.summarize_auc(hrf)
                    out_row.append(_format_output_value(auc))
                writer.writerow(out_row)
                row_count += 1

    provenance_payload = {
        "analysis_step": "collapse_homer_fir_to_auc",
        "basis_family": BASIS_FAMILY,
        "input_csv": str(input_path),
        "input_csv_sha256": compute_file_sha256(input_path),
        "settings_json": str(settings_path),
        "settings_json_sha256": compute_file_sha256(settings_path),
        "output_csv": str(output_path),
        "output_csv_sha256": compute_file_sha256(output_path),
        "row_count": row_count,
        "feature_count": len(specs),
        "basis_count_per_feature": reconstructor.basis_matrix.shape[1],
        "settings": settings_to_dict(settings),
    }
    provenance_path.write_text(
        json.dumps(provenance_payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    return output_path


def main() -> None:
    output_path = collapse_homer_fir_to_auc()
    print(f"[INFO] Wrote FIR-to-AUC beta table: {output_path}")
    print(f"[INFO] Wrote FIR-to-AUC provenance: {default_auc_provenance_path(output_path)}")


if __name__ == "__main__":
    main()
