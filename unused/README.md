# Archived material

These files are not used to generate the current manuscript results. They were
moved here on 2026-09-30 to distinguish the active workflow from historical or
alternative analyses. Nothing was deleted and no numerical analysis rule was
changed. Preserving this audit trail follows Sandve et al. (2013),
https://doi.org/10.1371/journal.pcbi.1003285; see the root `CITATIONS.md`.

## Inventory and recovery

`move_manifest.csv` records each original path, archived path, SHA-256 checksum,
and reason for moving it. All moved files preserve their original bytes.

The archive contains:

- The older demographic analysis, broad channel-behavior screen, broader
  pooled/raw correlation follow-up, and covariate correlation diagnostics.
- Their dedicated configuration/tests and saved covariate diagnostic outputs.
- The older study presentation, MNE planning list, and historical
  covariate-sensitivity decision record. The latter explains the model-choice
  history; its old numerical findings are not current results.
- Historical README/specification/test instructions in `LEGACY_WORKFLOWS.md`,
  `LEGACY_ANALYSIS_SPEC.md`, and `LEGACY_TESTS.md`. Some referenced optional
  scripts were already absent before this cleanup.

Archived scripts/tests are preserved for recovery, not maintained as runnable
alternatives in their new locations. Some resolve inputs or imports relative to
their original locations. Restore the manifest's original layout in an isolated
workspace before revisiting them; do not run old commands against production
outputs. If restoring in this repository, first ensure the original destination
does not exist, then move the corresponding archived file back. Do not overwrite
newer files.

## Active workflow

The root README describes the retained preprocessing pipeline, four primary
mixed-model analyses, behavioral pairwise correlations, selected pooled-mean
correlations, demographics table, and publication plots. Their shared helpers,
raw inputs, current configuration, result files, and validation tests remain in
place. Recall auditing, FIR plots, retention sensitivity, and citation checking
remain at the root as optional QC/robustness tools, as requested.

Ordinary pytest discovery is restricted to `tests/`; archived tests are not part
of current validation. This prevents retired analyses from being accidentally
reintroduced as dependencies or requirements for reproducing current results.
