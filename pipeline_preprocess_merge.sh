#!/usr/bin/env bash
set -euo pipefail

# Single entry-point orchestrator for preprocessing + merge.
#
# Reproducibility note:
# Centralizing execution order and hard-fail checks reduces drift across runs
# and preserves an auditable workflow (Sandve et al., 2013; see CITATIONS.md).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

STEP_INDEX=0
STEP_TOTAL=5

log() {
  echo "[pipeline] $*"
}

die() {
  echo "[pipeline][error] $*" >&2
  exit 1
}

need_file() {
  local p="$1"
  [[ -f "$p" ]] || die "Required file not found: $p"
}

need_cmd() {
  local c="$1"
  command -v "$c" >/dev/null 2>&1 || die "Required command not found in PATH: $c"
}

run_step() {
  local title="$1"
  shift
  STEP_INDEX=$((STEP_INDEX + 1))
  log "[STEP ${STEP_INDEX}/${STEP_TOTAL}] ${title}"
  "$@"
  log "[OK] ${title}"
}

log "Validating runtime dependencies"
need_cmd python
need_cmd Rscript

log "Validating required input files"
need_file "demographic/combined_engagement_data.csv"
need_file "data/tabular/recall_assessment_score_diffs.csv"
need_file "qualtrics/final_SF_demographic_data.csv"
need_file "data/tabular/homer3_glm_betas_wide.csv"

log "Cleaning data/results (unconditional)"
mkdir -p "data/results"
find "data/results" -mindepth 1 -maxdepth 1 -exec rm -rf {} +

run_step "Preprocess engagement ratings" \
  python process_engagement.py

run_step "Preprocess socio-demographic covariates" \
  python process_sociodemographic.py

run_step "Build combined tabular dataset" \
  python generate_combined_data.py

run_step "Merge Homer3 betas with combined tabular data" \
  Rscript merge_homer3_betas_with_combined_data.R \
    --homer_csv data/tabular/homer3_glm_betas_wide.csv \
    --combined_csv data/tabular/combined_sfv_data.csv \
    --out_csv data/tabular/homer3_betas_plus_combined_sfv_data_inner_join.csv

run_step "Certify preprocess+merge integrity" \
  python certify_preprocess_merge_integrity.py \
    --combined_csv data/tabular/combined_sfv_data.csv \
    --homer_csv data/tabular/homer3_glm_betas_wide.csv \
    --merged_csv data/tabular/homer3_betas_plus_combined_sfv_data_inner_join.csv \
    --out_json data/results/preprocess_merge_certification.json \
    --out_audit_csv data/results/preprocess_merge_id_audit.csv \
    --out_dropped_csv data/results/preprocess_merge_dropped_ids.csv

log "Pipeline completed successfully."
