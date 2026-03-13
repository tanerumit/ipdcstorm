# Goal
Reduce positive hindcast bias by tightening wind-field inflation around R34 handling, constraining hindcast intensity generation, and making Tier 1A use adjusted lambda by default.

# Scope
- `R/hazard_core.R`
- `R/hazard_validation.R`
- `man/dot-estimate_site_wind_holland.Rd`
- `man/make_validation_cfg.Rd`
- `tests/testthat/test-hazard-validation-exported.R`
- `tests/testthat/test-hazard_outer_cutoff.R`
- `tests/testthat/test-rmw_estimation.R`

# Summary
- Added explicit R34 provenance handling in `compute_site_winds_full()` and `.estimate_site_wind_holland()` with `observed`, `partial`, `climo`, and `none` pathways.
- Replaced the broad mean-radius fallback for incomplete directional R34 with a tighter partial fallback based on the minimum observed quadrant, and propagated the chosen source into diagnostics and solver behavior.
- Tightened fallback outer-cutoff handling and disabled R34 calibration inflation outside the fully observed R34 path.
- Changed `make_validation_cfg()` so Tier 1A uses adjusted rates by default (`hindcast_use_raw_rates = FALSE`).
- Constrained hindcast intensity sampling so draws are bounded by the observed pool support instead of freely extending to the class upper bound.

# Files Changed
- `R/hazard_core.R`
- `R/hazard_validation.R`
- `man/dot-estimate_site_wind_holland.Rd`
- `man/make_validation_cfg.Rd`
- `tests/testthat/test-hazard-validation-exported.R`
- `tests/testthat/test-hazard_outer_cutoff.R`
- `tests/testthat/test-rmw_estimation.R`

# Commands Run
- `Rscript -e "parse(file='R/hazard_core.R'); parse(file='R/hazard_validation.R')"`
- `Rscript -e "parse(file='tests/testthat/test-hazard_outer_cutoff.R'); parse(file='tests/testthat/test-rmw_estimation.R'); parse(file='tests/testthat/test-hazard-validation-exported.R')"`
- `Rscript -e "pkgload::load_all('.'); testthat::test_dir('tests/testthat', reporter='summary', filter='hazard[-_ ](core|validation)|hindcast|wind|rmw|outer_cutoff')"`
- `Rscript -e "devtools::document()"`
- `Rscript inst/extcode/pipelines/01-baseline-validation.R`
- `Rscript -e "pkgload::load_all('.'); ... run_hazard_model(...); run_validation_suite(...)"` for five-site benchmark reruns using workspace code

# Test Results
- Parse checks passed for touched R files and touched test files.
- Scoped tests passed under `pkgload::load_all('.')`.
- `devtools::document()` completed and regenerated the two touched Rd files.
- The packaged baseline pipeline run completed after creating `output/baseline/comparison`, but that script loads the installed package, not the workspace patch.
- The five-site benchmark rerun under `pkgload::load_all('.')` completed, but the hindcast acceptance criteria were not met.

# Behavior Changes
- Site-wind output now records the R34 source pathway and uses stricter fallback geometry when directional radii are incomplete or climatology-filled.
- Tier 1A hindcast now defaults to adjusted lambda.
- Hindcast intensity draws are bounded to the observed training-pool range.

# Results Delta
- Relative to the archived snapshot in `output/validation/old/hindcast_return_levels.csv`, the patched workspace benchmark moved in the wrong direction for the key Caribbean sites.
- Saba RL bias worsened from about `+11/+11/+12/+12%` to `+44/+39/+39/+41%` at RP `5/10/25/50`.
- Statia RL bias worsened from about `+13/+16/+21/+25%` to `+40/+40/+47/+54%`.
- St. Martin improved slightly at long return periods but still stayed positive, roughly `+18/+9/+5/+4%`.
- Simulated GEV shape also remained too heavy-tailed under the workspace run: `sim xi` was `-0.16` at Saba and `+0.02` at Statia versus observed `-0.30` and `-0.26`.
- The adjusted-lambda default appears to be a secondary upward driver, while the current wind/intensity changes did not offset it enough.

# Follow-ups / Risks
- Current workspace benchmark suggests this patch should not be treated as accepted bias reduction work yet.
- The main unresolved risk is systematic upper-tail overestimation at Saba and Statia under the patched workspace code.
- Before landing, the next pass should isolate the hindcast sampler and adjusted-lambda default separately from the wind-field changes to identify which component is dominating the tail inflation.
