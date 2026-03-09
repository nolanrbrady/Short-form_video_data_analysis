#!/usr/bin/env python
"""
Certification checks for preprocess+merge outputs.

This script verifies that merged data integrity is preserved and writes
machine-readable audit artifacts for reproducibility.

Reproducibility note:
Centralized, explicit assertions and persistent audit outputs reduce silent
analysis drift (Sandve et al., 2013; see CITATIONS.md).
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import pandas as pd


BETA_COL_PATTERN = re.compile(r"^S\d+_D\d+_Cond\d{2}_(HbO|HbR)$")


@dataclass(frozen=True)
class AssertionResult:
    name: str
    passed: bool
    detail: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Certify preprocess+merge integrity for SFV tabular + Homer3 join outputs."
    )
    parser.add_argument(
        "--combined_csv",
        default="data/tabular/combined_sfv_data.csv",
        help="Path to combined tabular CSV.",
    )
    parser.add_argument(
        "--homer_csv",
        default="data/tabular/homer3_glm_betas_wide_auc.csv",
        help="Path to Homer3 wide beta CSV.",
    )
    parser.add_argument(
        "--merged_csv",
        default="data/tabular/homer3_betas_plus_combined_sfv_data_inner_join.csv",
        help="Path to merged output CSV.",
    )
    parser.add_argument(
        "--out_json",
        default="data/results/preprocess_merge_certification.json",
        help="Path to output certification JSON.",
    )
    parser.add_argument(
        "--out_audit_csv",
        default="data/results/preprocess_merge_id_audit.csv",
        help="Path to output ID audit CSV.",
    )
    parser.add_argument(
        "--out_dropped_csv",
        default="data/results/preprocess_merge_dropped_ids.csv",
        help="Path to output dropped-ID CSV.",
    )
    return parser.parse_args()


def require_path_exists(path: Path) -> None:
    if not path.exists():
        raise SystemExit(f"[FAIL] Required input does not exist: {path}")


def normalize_subject_ids(raw: pd.Series, *, col_name: str, table_name: str) -> pd.Series:
    # Strict ID normalization: extract numeric component only; fail on unparseable values.
    raw_str = raw.astype("string")
    extracted = raw_str.str.extract(r"(\d+)", expand=False)
    bad_mask = extracted.isna()
    if bad_mask.any():
        bad_examples = raw_str[bad_mask].dropna().unique().tolist()[:10]
        raise SystemExit(
            "[FAIL] Unparseable subject IDs in "
            f"{table_name}.{col_name}. Examples: {bad_examples}"
        )
    return extracted.astype("int64")


def assert_no_duplicates(ids: pd.Series, *, table_name: str) -> AssertionResult:
    dupes = ids[ids.duplicated()].unique().tolist()
    return AssertionResult(
        name=f"{table_name}: no duplicate normalized subject IDs",
        passed=(len(dupes) == 0),
        detail="OK" if len(dupes) == 0 else f"Duplicates found: {dupes[:10]}",
    )


def assert_non_empty(ids: Iterable[int], *, label: str) -> AssertionResult:
    ids_list = list(ids)
    return AssertionResult(
        name=f"{label}: non-empty subject set",
        passed=(len(ids_list) > 0),
        detail=f"n={len(ids_list)}",
    )


def write_csv(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def main() -> None:
    args = parse_args()

    combined_path = Path(args.combined_csv)
    homer_path = Path(args.homer_csv)
    merged_path = Path(args.merged_csv)
    out_json_path = Path(args.out_json)
    out_audit_path = Path(args.out_audit_csv)
    out_dropped_path = Path(args.out_dropped_csv)

    require_path_exists(combined_path)
    require_path_exists(homer_path)
    require_path_exists(merged_path)

    combined = pd.read_csv(combined_path)
    homer = pd.read_csv(homer_path)
    merged = pd.read_csv(merged_path)

    required_columns = [
        (combined, "subject_id", str(combined_path)),
        (homer, "Subject", str(homer_path)),
        (merged, "subject_id", str(merged_path)),
    ]
    for df, col_name, table_name in required_columns:
        if col_name not in df.columns:
            raise SystemExit(f"[FAIL] Missing required column '{col_name}' in {table_name}")

    beta_cols = [c for c in merged.columns if BETA_COL_PATTERN.match(c)]
    if not beta_cols:
        raise SystemExit(
            "[FAIL] Merged CSV has no beta columns matching pattern "
            "'S##_D##_Cond##_HbO/HbR'."
        )

    combined_ids = normalize_subject_ids(combined["subject_id"], col_name="subject_id", table_name="combined")
    homer_ids = normalize_subject_ids(homer["Subject"], col_name="Subject", table_name="homer")
    merged_ids = normalize_subject_ids(merged["subject_id"], col_name="subject_id", table_name="merged")

    combined_set = set(combined_ids.tolist())
    homer_set = set(homer_ids.tolist())
    merged_set = set(merged_ids.tolist())
    intersection_set = combined_set.intersection(homer_set)
    combined_only = sorted(combined_set - merged_set)
    homer_only = sorted(homer_set - merged_set)

    assertions: list[AssertionResult] = []
    assertions.append(assert_non_empty(combined_set, label="combined"))
    assertions.append(assert_non_empty(homer_set, label="homer"))
    assertions.append(assert_non_empty(merged_set, label="merged"))
    assertions.append(assert_no_duplicates(combined_ids, table_name="combined"))
    assertions.append(assert_no_duplicates(homer_ids, table_name="homer"))
    assertions.append(assert_no_duplicates(merged_ids, table_name="merged"))
    assertions.append(
        AssertionResult(
            name="merged set equals intersection(combined, homer)",
            passed=(merged_set == intersection_set),
            detail=f"merged_n={len(merged_set)}, intersection_n={len(intersection_set)}",
        )
    )
    assertions.append(
        AssertionResult(
            name="merged row count equals intersection size",
            passed=(len(merged) == len(intersection_set)),
            detail=f"merged_rows={len(merged)}, intersection_n={len(intersection_set)}",
        )
    )
    assertions.append(
        AssertionResult(
            name="no subject loss from combined after merge (strict)",
            passed=(len(combined_only) == 0),
            detail="OK" if len(combined_only) == 0 else f"combined_only={combined_only[:10]}",
        )
    )
    assertions.append(
        AssertionResult(
            name="no subject loss from homer after merge (strict)",
            passed=(len(homer_only) == 0),
            detail="OK" if len(homer_only) == 0 else f"homer_only={homer_only[:10]}",
        )
    )

    universe = sorted(combined_set.union(homer_set).union(merged_set))
    id_audit = pd.DataFrame(
        {
            "subject_id_norm": universe,
            "in_combined": [sid in combined_set for sid in universe],
            "in_homer": [sid in homer_set for sid in universe],
            "in_merged": [sid in merged_set for sid in universe],
        }
    )
    id_audit["status"] = id_audit.apply(
        lambda row: (
            "intersection"
            if row["in_combined"] and row["in_homer"] and row["in_merged"]
            else "combined_only"
            if row["in_combined"] and (not row["in_homer"]) and (not row["in_merged"])
            else "homer_only"
            if row["in_homer"] and (not row["in_combined"]) and (not row["in_merged"])
            else "unexpected_membership"
        ),
        axis=1,
    )

    dropped_rows = []
    for sid in combined_only:
        dropped_rows.append({"subject_id_norm": sid, "source": "combined_only"})
    for sid in homer_only:
        dropped_rows.append({"subject_id_norm": sid, "source": "homer_only"})
    dropped_df = pd.DataFrame(dropped_rows, columns=["subject_id_norm", "source"])

    write_csv(out_audit_path, id_audit)
    write_csv(out_dropped_path, dropped_df)

    payload = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "combined_csv": str(combined_path),
            "homer_csv": str(homer_path),
            "merged_csv": str(merged_path),
        },
        "outputs": {
            "id_audit_csv": str(out_audit_path),
            "dropped_ids_csv": str(out_dropped_path),
        },
        "counts": {
            "combined_rows": int(len(combined)),
            "homer_rows": int(len(homer)),
            "merged_rows": int(len(merged)),
            "combined_unique_ids": int(len(combined_set)),
            "homer_unique_ids": int(len(homer_set)),
            "merged_unique_ids": int(len(merged_set)),
            "intersection_size": int(len(intersection_set)),
            "combined_only_count": int(len(combined_only)),
            "homer_only_count": int(len(homer_only)),
            "beta_columns_in_merged": int(len(beta_cols)),
        },
        "assertions": [
            {"name": a.name, "passed": a.passed, "detail": a.detail} for a in assertions
        ],
    }
    payload["passed"] = all(a.passed for a in assertions)
    write_json(out_json_path, payload)

    for a in assertions:
        status = "PASS" if a.passed else "FAIL"
        print(f"[{status}] {a.name} :: {a.detail}")
    print(f"[INFO] wrote: {out_json_path}")
    print(f"[INFO] wrote: {out_audit_path}")
    print(f"[INFO] wrote: {out_dropped_path}")

    if not payload["passed"]:
        raise SystemExit("[FAIL] Certification failed. See JSON/CSV artifacts for details.")

    print("[OK] Certification passed.")


if __name__ == "__main__":
    main()
