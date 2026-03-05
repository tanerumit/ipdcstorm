Goal
- Run the real validation workflow with Tier 1A + Tier 1B on the packaged IBTrACS input, verify deterministic Tier 1B outputs across repeated runs, confirm fold usability for Miami/Saba/Statia, document the blocked-CV definition, and add a brief methods note.

Scope
- Validation-layer verification plus minimal follow-up edits for Tier 1B messaging and documentation.
- No hazard physics, simulation logic, exported signatures, dependencies, or scoring defaults were changed.

Summary
- Ran the production validation path from `inst/extcode/01-pipeline-baseline.R` using:
  - `cfg = make_hazard_cfg(data_path = "inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv", search_radius_km = 800, start_year = 1970L, n_sim_years = 2000L)`
  - the standard 5-site `targets`
  - `validate_hazard_model(..., validation_cfg = make_validation_cfg(), severities = c("TS", "HUR"), sst_cfg = NULL)`
- Executed the real workflow twice back-to-back and confirmed deterministic Tier 1B per-site metrics and deterministic CSV outputs.
- Verified Tier 1A and Tier 1B both print in the real workflow console. A direct log capture via base `sink()` missed `message()` output, but the shell console output from both runs showed both sections.
- Confirmed Miami, Saba, and Statia each retain `folds_used = 4` out of `folds_total = 6`; no default change to `block_size` or `test_n_pos_years >= 2` was needed.
- Confirmed fold construction is index-based: contiguous groups of observed annual-max years, not fixed calendar decades.
- Added a Tier 1B header note that states the block definition and defaults.
- Added manuscript-table metadata attributes and a short methods paragraph in the validation tutorial.
- Tightened the Tier 1B required-columns error message to explicitly mention accepted aliases.

Files Changed
- `R/hazard_validation.R`
- `inst/tutorials/tutorial_2_multisite_validation_climate.qmd`

Commands Run
- `Get-Content -Path AGENTS.md`
- `Get-Content -Path inst\extcode\01-pipeline-baseline.R`
- `rg -n run_validation_suite .`
- `rg -n validate_hazard_model .`
- Real workflow run twice via `@' ... '@ | Rscript -` using the production `validate_hazard_model(...)` call and `make_validation_cfg()`
- `Rscript -e "parse(file='R/hazard_validation.R')"`
- Real workflow run once more after the follow-up edits via `@' ... '@ | Rscript -`
- `Get-Date -Format "yyyy-MM-dd__HHmm"`

Test Results
- `Rscript -e "parse(file='R/hazard_validation.R')"`: passed
- Real workflow (run 1): completed
  - Tier 1A printed
  - Tier 1B printed
  - Tier 1B site summary:
    - Miami: `6 / 4`, `mean_test_p0 = 0.65`, `mean_test_q90 = 48.6`, `max_test_q95 = 95.9`
    - Puerto_Rico: `6 / 4`, `0.60`, `55.2`, `76.4`
    - Saba: `6 / 4`, `0.375`, `58.6`, `107`
    - St_Martin: `6 / 6`, `0.483`, `60.5`, `113`
    - Statia: `6 / 4`, `0.55`, `52.3`, `89.5`
- Real workflow (run 2): completed
  - Tier 1A printed
  - Tier 1B printed
  - Tier 1B site summary matched run 1 exactly
- Determinism checks:
  - Tier 1B metrics identical across run 1 and run 2: `TRUE`
  - CSV hashes identical across run 1 and run 2: `TRUE`
  - Matching CSV hashes were confirmed for:
    - `output/validation/modern_blocked_cv_summary.csv`
    - `output/validation/modern_blocked_cv_compact.csv`
    - `output/validation/validation_summary.csv`
    - `output/validation/hindcast_return_levels.csv`
- Real workflow (post-edit verification): completed
  - Tier 1B header now prints:
    - `Block definition: contiguous index-based groups of 10 observed annual-max years (not fixed calendar decades).`
    - `Defaults: era >= 1970, TS+ storms (\`Vmax >= 34 kt\`); scored folds require test_n_pos_years >= 2.`
  - Paper-table metadata attributes:
    - `era_label = "1970+"`
    - `block_size_years = 10`
    - `block_definition = "contiguous index-based groups of observed annual-max years"`
    - `storm_vmax_min = 34`
    - `threshold_kt = 34`
    - `scoring_filter = "test_n_pos_years >= 2"`

Behavior Changes
- Tier 1B now states its block construction and scoring defaults in the console output.
- The compact manuscript table now carries in-memory metadata attributes describing era, block definition, thresholds, and scoring filter.
- The Tier 1B fail-fast error message now explicitly states the stable required inputs and accepted aliases:
  - `iso_time`
  - site wind from `V_site_kt` or `site_wind_kt`
  - storm wind from `Vmax_kt` or `storm_wind_kt`
- The validation tutorial now includes a short Tier 1B methods note consistent with the implemented defaults.

Follow-ups / Risks
- Tier 1A still reports holdout insufficiency for Miami, Puerto Rico, and Saba under the current single-split default because those sites have fewer than 20 annual years after event filtering; this is unchanged and remains a separate limitation from Tier 1B.
- The repeated warning from `plot_cdf_comparison()` (`scale_y_log10` infinite values) still appears during the full validation workflow; this was observed but not changed because it is outside the requested Tier 1B scope.
