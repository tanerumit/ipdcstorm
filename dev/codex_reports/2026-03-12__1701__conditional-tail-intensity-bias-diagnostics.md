## 1. Goal
Identify which retained high-end hindcast events drive residual upper-tail overestimation at Saba and Statia under the restored legacy path, split the diagnosis by `R34_source`, and determine whether the next patch should target incomplete-`R34` fallback geometry or fully observed conditional wind decay/taper behavior.

## 2. Scope
Included:
- Internal-only hindcast tail diagnostics in `R/hazard_validation.R`
- Deterministic event-level top-tail extraction for observed test years and simulated annual maxima
- Pathway-split tail summaries and metadata capture
- Targeted test coverage in `tests/testthat/test-hazard-validation-exported.R`

Excluded:
- Any package default changes
- Any new fallback rule
- Any sampler, lambda, or physics retuning
- Any public API changes

## 3. Problem solved
The prior attribution isolated residual bias as intensity-dominated but did not identify which retained events and which `R34_source` pathways actually populate the conditional tail at Saba and Statia. This patch adds reproducible internal diagnostics to localize that tail behavior before any surgical model change.

## 4. Summary
Added internal helpers to:
- carry peak-event geometry and `R34_source` provenance into hindcast retention diagnostics
- extract top observed test-period annual-max contributor events with closest approach, RMW, quadrant, and directional `R34` completeness
- summarize pathway shares in the top tail versus the overall retained-event pool
- attach deterministic run metadata to all new diagnostic tables

Fixed one sparse-data edge case so threshold-share columns remain present even when no test annual maxima exceed a requested threshold.

Workspace legacy-path diagnostics with fixed seeds indicate:
- Saba test-period top tail is dominated by `observed` events, not `climo`
- Statia test-period top tail is also dominated by `observed` events, not `climo`
- `partial` contributes materially but is secondary in the top-5 annual maxima at both sites
- This supports a next patch focused on conditional wind decay/taper for retained intense events rather than incomplete-`R34` fallback geometry alone

## 5. Files changed
- `R/hazard_validation.R` — internal diagnostic helpers and hindcast diagnostic plumbing; large internal-only edit
- `tests/testthat/test-hazard-validation-exported.R` — targeted deterministic tests for tail extraction, pathway summaries, missing-pathway handling, and metadata capture

Approximate diff:
- `R/hazard_validation.R`: `+607 / -?` lines in current diff stat
- `tests/testthat/test-hazard-validation-exported.R`: `+105` lines in current diff stat

## 6. Commands run
```powershell
Rscript -e "parse(file='R/hazard_validation.R')"
Rscript -e "parse(file='tests/testthat/test-hazard-validation-exported.R')"
Rscript -e "devtools::test(filter = 'hazard[-_ ](validation|rmw)|hindcast|wind|R34|tail')"
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
out <- run_hazard_model(
  cfg = cfg,
  targets = targets,
  seed = 2026L,
  climate = make_climate_cfg(scenario = 'stationary'),
  verbose = FALSE
)
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
for (loc in c('Saba', 'Statia', 'St_Martin')) {
  hi <- hc$per_island[[loc]]
  comp <- hc$comparison %>% filter(location == loc) %>%
    transmute(return_period, obs_test_rl, sim_rl, rl_bias_pct = 100 * (sim_rl - obs_test_rl) / pmax(obs_test_rl, 1))
  xi_tbl <- tibble(
    location = loc,
    sim_xi = hi$gev_fit$gev_fit$xi %||% NA_real_,
    obs_xi = hi$obs_gev$gev_fit$xi %||% NA_real_
  )
  overall_share <- hi$diagnostics$tail_pathway_summary %>%
    select(peak_r34_source, n_overall_retained_events, overall_retained_event_share,
           n_top_tail_events, top_n_annual_max_share,
           tidyselect::starts_with('n_threshold_exceedances_'),
           tidyselect::starts_with('threshold_exceedance_share_'))
  top_events <- hi$diagnostics$tail_event_detail %>%
    filter(sample == 'observed_test') %>%
    select(location, year, storm_id, storm_class, annual_max_rank,
           simulated_site_wind_kt, observed_site_year_annual_max_kt,
           closest_approach_km, peak_rmw_used_km, peak_r34_source,
           peak_r34_completeness, peak_quadrant)
  sim_top <- hi$diagnostics$tail_event_detail %>%
    filter(sample == 'simulated_annual_max') %>%
    select(location, sample_year, annual_max_rank, simulated_site_wind_kt)
  cmp <- hi$diagnostics$tail_pathway_comparison

  cat('\n=== ', loc, ' ===\n', sep = '')
  cat('RL bias (%):\n')
  print(as.data.frame(comp), row.names = FALSE)
  cat('GEV xi:\n')
  print(as.data.frame(xi_tbl), row.names = FALSE)
  cat('Tail comparison:\n')
  print(as.data.frame(cmp), row.names = FALSE)
  cat('Pathway summary:\n')
  print(as.data.frame(overall_share), row.names = FALSE)
  cat('Observed test top-tail events:\n')
  print(as.data.frame(top_events), row.names = FALSE)
  cat('Simulated annual-max top tail:\n')
  print(as.data.frame(sim_top), row.names = FALSE)
}
'@ | Rscript -
```

## 7. Test results
- `Rscript -e "devtools::test(filter = 'hazard[-_ ](validation|rmw)|hindcast|wind|R34|tail')"`: passed
- Result: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 125 ]`
- Workspace diagnostic run: passed after fixing sparse threshold-column handling in `tail_pathway_summary`

## 8. Behavior changes
No public behavior change. New outputs remain internal-only and are attached to hindcast diagnostics as additional internal tables:
- `tail_event_detail`
- `tail_pathway_summary`
- `tail_pathway_comparison`

## 9. Follow-ups/risks
- The new decision signal is strongest for observed test-period annual-max contributors; simulated annual-max top rows still do not carry a source-specific synthetic event ID because the legacy hindcast sampler does not preserve that provenance.
- For the current fixed-seed validation run, Saba and Statia top-5 observed test annual maxima are dominated by `observed` cases. A focused next patch should therefore prioritize conditional wind decay/taper behavior for retained intense events, especially where closest approach is several multiples of RMW.
- Incomplete-`R34` handling is still present in the retained pool and may remain a secondary contributor, but this run does not support treating it as the primary driver at Saba or Statia.
