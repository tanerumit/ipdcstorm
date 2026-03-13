## 1. Goal

Identify and reduce residual positive hindcast bias under the restored legacy validation path by adding internal diagnostics that separate frequency, intensity, interaction, `R34_source`, and event-retention effects for Saba, Statia, and St_Martin.

## 2. Scope

Included:
- internal hindcast diagnostic helpers and attribution-grid outputs
- `R34_source` stratification, retention summaries, and wind-mode retention comparison
- targeted tests for decomposition, retention, stratification, and metadata capture
- brief note documenting the near-threshold assumption

Not included:
- any public API changes
- any default change away from raw lambda, legacy wind fallback, or legacy sampler
- empirical retuning or model recalibration

## 3. Problem solved

After restoring legacy defaults, residual hindcast bias at Saba and Statia still needed to be localized to specific mechanisms instead of being treated as a single aggregate error. The package did not yet expose internal diagnostics to separate effective frequency, intensity, `R34_source`, and threshold-retention behavior under the restored baseline.

## 4. Summary

Added internal hindcast diagnostics in `R/hazard_validation.R` that now emit:
- frequency / intensity / interaction decomposition for each site
- per-site `R34_source` summaries across `observed`, `partial`, `climo`, and `none`
- annual event-retention summaries over train/test splits, including zero-event years, TS years, HUR years, and near-threshold years
- deterministic run metadata on all new diagnostic tables
- same-seed wind-mode retention comparisons for legacy vs diagnostic wind handling

Workspace diagnostics under the restored baseline (`legacy` wind, `raw` rates, `legacy` sampler) showed:
- Saba residual RL bias stayed positive at RP 5/10/25/50: `+6.3%`, `+8.9%`, `+10.2%`, `+10.8%`
- Statia residual RL bias stayed more strongly positive: `+10.9%`, `+17.3%`, `+23.7%`, `+28.2%`
- St_Martin was near-neutral to slightly low at shorter RP: `-6.5%`, `-2.1%`, `-0.0%`, `+1.0%`
- decomposition under the restored baseline was intensity-dominated at Saba and Statia, with frequency contribution negative at both sites
- baseline retained events at Saba and Statia were concentrated in `observed`, `partial`, and `climo`, with `none` contributing zero retained TS+/HUR years
- switching to the diagnostic wind mode mainly redistributed retained TS+ events from `partial` to `observed`, with no HUR-count change in the comparison tables
- near-threshold years were present but modest: Saba `2` train / `2` test, Statia `3` train / `2` test

## 5. Files changed

- `R/hazard_validation.R`
  - approx. `+502 / -10` diff lines
  - added internal metadata, decomposition, `R34_source`, retention, and wind-comparison helpers
  - threaded diagnostics through `.validate_hindcast_all()` and `.run_hindcast_attribution_grid()`
- `tests/testthat/test-hazard-validation-exported.R`
  - approx. `+128 / -0` diff lines
  - added focused tests for decomposition, retention, `R34_source` summaries, wind comparison, and attribution-grid outputs
- `dev/notes/hindcast-near-threshold-band.md`
  - new `8` line note
  - documents the diagnostic near-threshold definition

## 6. Commands run

- `Rscript -e "parse(file='R/hazard_validation.R')"`
- `Rscript -e "parse(file='tests/testthat/test-hazard-validation-exported.R')"`
- `Rscript -e "devtools::test(filter = 'hazard[-_ ](validation|outer_cutoff|rmw)|hindcast|wind|bias|retention')"`
- workspace baseline/grid diagnostic run via `pkgload::load_all('.')` and `ipdcstorm:::.run_hindcast_attribution_grid(...)`
- focused workspace diagnostic run via `pkgload::load_all('.')` and `ipdcstorm:::.run_hindcast_attribution_case(...)` for:
  - `wind_field_mode = 'legacy', annual_rate_mode = 'raw', sampler_mode = 'legacy'`
  - `wind_field_mode = 'diagnostic_new', annual_rate_mode = 'raw', sampler_mode = 'legacy'`

## 7. Test results

- Parse checks: passed
- Targeted test filter: passed
  - `0` failures
  - `0` warnings
  - `104` passes
- Workspace diagnostics: completed from workspace code path with fixed seeds

## 8. Behavior changes

- No exported behavior changed
- No package defaults changed
- Internal hindcast attribution results now include:
  - `baseline_case_id`
  - `baseline_diagnostics`
  - `case_diagnostics`
  - `wind_retention_comparison`
- Internal per-site hindcast diagnostics now carry deterministic metadata fields:
  - `model_seed`
  - `validation_seed`
  - `data_id`
  - `parameter_id`
  - `lambda_scaler_id`

## 9. Follow-ups/risks

- The dominant residual bias at Saba and Statia under restored defaults appears to be intensity / tail-shape, not effective frequency.
- `partial` and `observed` pathways dominate retained TS+ events at Saba and Statia; `none` is not overrepresented in retained biased years under the legacy baseline.
- The diagnostic wind comparison suggests fallback geometry mainly reclassifies TS+ retention between `partial` and `observed` pathways rather than creating new HUR retention.
- Near-threshold retention exists but does not by itself look like the primary driver of the remaining positive RL bias at Saba and Statia.
- The next patch, if any, is better justified as a narrow fallback-geometry / tail-intensity investigation than another frequency or exposure-rate adjustment.
