#!/usr/bin/env bash
set -euo pipefail

# Single entry-point orchestrator for preprocessing + merge.
#
# Reproducibility note:
# Centralizing execution order and hard-fail checks reduces drift across runs
# and preserves an auditable workflow (Sandve et al., 2013; see CITATIONS.md).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

GENERATED_TABULAR_DIR="data/tabular/generated_data"
ENGAGEMENT_INPUT_CSV="demographic/combined_engagement_data.csv"
QUALTRICS_INPUT_CSV="qualtrics/final_SF_demographic_data.csv"
HOMER_RAW_FIR_CSV="data/tabular/homer3_glm_betas_wide_fir_pca.csv"
PREPROCESS_SETTINGS_JSON="data/config/preprocessing_settings.json"
PREPROCESS_CERT_JSON="data/results/preprocess_merge_certification.json"
PREPROCESS_ID_AUDIT_CSV="data/results/preprocess_merge_id_audit.csv"
PREPROCESS_DROPPED_IDS_CSV="data/results/preprocess_merge_dropped_ids.csv"

# Generated Files
RECALL_DIFFS_CSV="${GENERATED_TABULAR_DIR}/recall_assessment_score_diffs.csv"
HOMER_AUC_CSV="${GENERATED_TABULAR_DIR}/homer3_glm_betas_wide_auc.csv"
COMBINED_TABULAR_CSV="${GENERATED_TABULAR_DIR}/combined_sfv_data.csv"
MERGED_TABULAR_CSV="${GENERATED_TABULAR_DIR}/homer3_betas_plus_combined_sfv_data_inner_join.csv"

STEP_INDEX=0
STEP_TOTAL=7

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
need_file "$ENGAGEMENT_INPUT_CSV"
need_file "$RECALL_DIFFS_CSV"
need_file "$QUALTRICS_INPUT_CSV"
need_file "$HOMER_RAW_FIR_CSV"
need_file "$PREPROCESS_SETTINGS_JSON"

log "Cleaning data/results (unconditional)"
mkdir -p "data/results"
find "data/results" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
mkdir -p "$GENERATED_TABULAR_DIR"

run_step "Preprocess engagement ratings" \
  python process_engagement.py

run_step "Preprocess socio-demographic covariates" \
  python process_sociodemographic.py

run_step "Build combined tabular dataset" \
  python generate_combined_data.py

run_step "Collapse Homer FIR betas to baseline-corrected task-window AUC betas" \
  python collapse_homer_fir_to_auc.py \
    --input-csv "$HOMER_RAW_FIR_CSV" \
    --output-csv "$HOMER_AUC_CSV" \
    --settings-json "$PREPROCESS_SETTINGS_JSON"

run_step "Lint FIR-to-AUC conversion against excluded-channel policy" \
  python validate_homer_fir_auc_conversion.py \
    --raw-fir-csv "$HOMER_RAW_FIR_CSV" \
    --auc-csv "$HOMER_AUC_CSV" \
    --settings-json "$PREPROCESS_SETTINGS_JSON"

run_step "Merge Homer3 betas with combined tabular data" \
  Rscript merge_homer3_betas_with_combined_data.R \
    --homer_csv "$HOMER_AUC_CSV" \
    --combined_csv "$COMBINED_TABULAR_CSV" \
    --out_csv "$MERGED_TABULAR_CSV"

run_step "Certify preprocess+merge integrity" \
  python certify_preprocess_merge_integrity.py \
    --combined_csv "$COMBINED_TABULAR_CSV" \
    --homer_csv "$HOMER_AUC_CSV" \
    --merged_csv "$MERGED_TABULAR_CSV" \
    --out_json "$PREPROCESS_CERT_JSON" \
    --out_audit_csv "$PREPROCESS_ID_AUDIT_CSV" \
    --out_dropped_csv "$PREPROCESS_DROPPED_IDS_CSV"

log "Pipeline completed successfully."
