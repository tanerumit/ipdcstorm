## 1. Goal
Determine whether residual positive hindcast bias at Saba and Statia is concentrated in specific `R34_source` pathways under the restored legacy validation path, without changing defaults or operational behavior.

## 2. Scope
Included:
- Internal-only hindcast diagnostics in `R/hazard_validation.R`
- Targeted tests for `R34_source` stratification, deterministic completion, metadata attachment, and top-tail reproducibility
- Workspace diagnostic run for `Saba`, `Statia`, and `St_Martin` with fixed seeds on the legacy path

Excluded:
- Any package-default changes
- New fallback rules
- Wind-field physics retuning, lambda recalibration, or sampler changes
- Public API / export changes

## 3. Problem solved
Search-radius sensitivity and prior attribution left a residual upper-tail inflation signal at Saba and Statia. Existing diagnostics showed event retention by `R34_source`, but did not isolate which pathways were actually driving annual-max tail years and threshold exceedances.

## 4. Summary
Extended the existing hindcast retention diagnostics to stratify retained events, annual-max contributors, threshold exceedances, and top annual maxima by `R34_source`, while attaching run metadata to every table. Added an annual-max comparison summary that keeps observed vs simulated distribution diagnostics alongside dominant observed/top-tail `R34_source` labels.

Workspace diagnostics on the restored legacy path show:
- `Saba`: top-tail and annual-max exceedances above 64 kt and 96 kt are entirely `climo` in training years; top-5 annual maxima are 4 `climo`, 1 `observed`.
- `Statia`: top-tail and annual-max exceedances above 64 kt and 96 kt are entirely `climo` in training years; top-5 annual maxima are 4 `climo`, 1 `observed`.
- `St_Martin`: upper tail is mixed; `climo` dominates the train-period top tail, but `observed` still contributes materially, including the strongest 2017 event.

## 5. Files changed
- `R/hazard_validation.R` — `+232 / -10`
  Added internal `R34_source` tail helpers, expanded hindcast diagnostics, and propagated new diagnostic tables through the hindcast attribution collection path.
- `tests/testthat/test-hazard-validation-exported.R` — `+64 / -0`
  Added targeted tests for four-pathway stratification, missing-path completion, metadata propagation, top-tail tables, and reproducibility.

## 6. Commands run
```powershell
Rscript -e "parse(file='R/hazard_validation.R')"
Rscript -e "parse(file='tests/testthat/test-hazard-validation-exported.R')"
Rscript -e "devtools::test(filter = 'hazard[-_ ](validation|rmw)|hindcast|wind|R34|tail')"
```

```powershell
@'
pkgload::load_all('.')
library(dplyr)
library(tibble)

targets <- tibble::tribble(
  ~location, ~lat, ~lon,
  'St_Martin', 18.0708, -63.0501,
  'Saba', 17.6350, -63.2300,
  'Statia', 17.4890, -62.9740
)

cfg <- make_hazard_cfg(
  data_path = 'inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv',
  search_radius_km = 600,
  historical_start_year = 1970L,
  simulation_years = 2000L
)
validation_cfg <- make_validation_cfg(
  holdout_years = 10L,
  n_sim = 2000L,
  return_periods = c(5, 10, 25, 50),
  conf_level = 0.90,
  seed = 42L,
  save_plots = FALSE,
  save_tables = FALSE
)
old_opts <- options(ipdcstorm.wind_field_mode = 'legacy', ipdcstorm.hindcast_sampler_mode = 'legacy')
on.exit(options(old_opts), add = TRUE)
out <- run_hazard_model(cfg = cfg, targets = targets, seed = 2026L, climate = make_climate_cfg(scenario = 'stationary'), verbose = FALSE)
hc <- ipdcstorm:::.validate_hindcast_all(
  out = out,
  holdout_years = validation_cfg$holdout_years,
  n_sim = validation_cfg$n_sim,
  return_periods = validation_cfg$return_periods,
  conf_level = validation_cfg$conf_level,
  seed = validation_cfg$seed,
  beta_sst = 0,
  gamma_intensity = 0,
  use_raw_rates = TRUE,
  xi_bounds = validation_cfg$advanced$xi_bounds,
  n_boot = 0L
)
'@ | Rscript -
```

## 7. Test results
- `Rscript -e "devtools::test(filter = 'hazard[-_ ](validation|rmw)|hindcast|wind|R34|tail')"`: PASS
- Parse checks for touched files: PASS

## 8. Behavior changes
- No public behavior changed.
- Hindcast diagnostics now expose additional internal tables:
  - `threshold_exceedance`
  - `top_annual_max`
  - `annual_max_comparison`
- Existing `r34_source_summary` now includes retained-event counts, affected-year counts, annual-max sample shares, top-tail shares, and upper quantiles.

## 9. Follow-ups/risks
- The new attribution is observational/diagnostic: it localizes which historical `R34_source` pathways dominate the tail, but does not by itself change the hindcast generator.
- For Saba and Statia, evidence supports a next patch focused on incomplete / climatological `R34` handling rather than fully observed cases.
- A future follow-up could add a compact printed/reporting helper for these diagnostics if repeated workspace use becomes common, but that was intentionally not promoted into exported API here.
