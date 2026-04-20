# Run hazard model and validate in one step

End-to-end convenience wrapper: runs
[`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md)
then
[`run_validation_suite()`](https://tanerumit.github.io/ipdcstorm/reference/run_validation_suite.md).
This is the recommended entry point for standard validation workflows.

The validation uses a three-tier framework:

- **Tier 1 (Hindcast)**: Hold out the last N years of the historical
  record, simulate synthetic annual maxima from the training period, and
  compare return levels. The simulated return level should fall within
  the observed return level's confidence interval (computed via
  delta-method or parametric bootstrap of the hurdle-GEV model).

- **Tier 2 (Rate check)**: Compare model storm rates (lambda) against
  published HURDAT2/IBTrACS climatologies for the target region.

- **Tier 3 (Wind field)**: Spot-check site-level winds from the Holland
  wind profile model against historical station observations.

## Usage

``` r
validate_hazard_model(
  cfg,
  targets,
  validation_cfg = make_validation_cfg(),
  storm_classes = c("TS", "HUR")
)
```

## Arguments

- cfg:

  A `hazard_cfg` object from
  [`make_hazard_cfg()`](https://tanerumit.github.io/ipdcstorm/reference/make_hazard_cfg.md).
  Controls basin, year range, search radius, and core model parameters.

- targets:

  Data frame of target locations with columns:

  - `location` (character): Site name (e.g., "Saba")

  - `lat` / `lon` (numeric): Coordinates in decimal degrees

  - `search_radius_km` (numeric): Radius for storm proximity filtering

- validation_cfg:

  A `validation_cfg` object from
  [`make_validation_cfg()`](https://tanerumit.github.io/ipdcstorm/reference/make_validation_cfg.md).
  Controls holdout period, simulation size, return periods, confidence
  level, and output paths. Default uses 90% CI with 10-year holdout and
  inherited synthetic years from the hazard model output.

- storm_classes:

  Character vector of storm classes to include (default:
  `c("TS", "HUR")`). Must match IBTrACS classification codes. Use
  `"HUR"` alone for hurricane-only analysis.

## Value

A list with three elements:

- `out`:

  Full hazard model output from
  [`run_hazard_model()`](https://tanerumit.github.io/ipdcstorm/reference/run_hazard_model.md),
  containing `events`, `fit`, `annual_max`, and `config` components.

- `val`:

  Validation results from
  [`run_validation_suite()`](https://tanerumit.github.io/ipdcstorm/reference/run_validation_suite.md): -
  `hindcast`: Per-island return level comparison table with CI
  membership flags. Access via `val$hindcast$comparison`. -
  `rate_check`: Rate sanity check results with OK/FLAG status. -
  `wind_field`: Wind field spot-check results with bias estimates. -
  `summary`: Tibble with compact hindcast tail-diagnostic columns.

- `artifacts`:

  Named list of saved file paths: - `plots`: Named character vector of
  PNG paths. - `tables`: Named character vector of CSV/MD paths.

## Examples

``` r
# --- Basic validation with defaults (90% CI) ---
targets_df <- data.frame(
  location         = c("Saba", "St. Eustatius", "St. Martin"),
  lat              = c(17.63, 17.49, 18.04),
  lon              = c(-63.23, -62.97, -63.07),
  search_radius_km = c(200, 200, 200)
)

result <- validate_hazard_model(
  cfg     = make_hazard_cfg(simulation_years = 2000),
  targets = targets_df
)
#> 
#> ── Loading data ──────────────────────────────────────────────
#>    Input file       : ibtracs_1970.csv
#>    Raw NA input     : 49,905 track points (1970–2025)
#>    Model start year : 1970
#> 
#> ── Target/location filtering ─────────────────────────────────
#>    Location          TS    HUR
#>    Saba              26      3
#>    St. Eustatius     27      4
#>    St. Martin        29      7
#>    Model sample     : 704 target-event records (1970–2025) across 3 locations
#> 
#> ── Rate calibration ──────────────────────────────────────────
#>    Adjustments      : 6 location/class adjustments (mode=target, clamped=0, missing_ref=4)
#>    Some location/class rates were increased to match reference rates; set advanced$lambda_scaling_mode='down_only' to prevent upscaling.
#> 
#> ── Climate ───────────────────────────────────────────────────
#>    Climate mode     : baseline (delta_sst = 0)
#>    Climate input    : scenario helper
#>    Climate scenario : stationary
#>    SST baseline     : 1991–2020
#>    Sensitivity mode : fixed
#>    Rate effect      : beta_0=0.60 -> beta=0.60 (+82% per +1°C) [guardrail]
#>    Count regime     : baseline | raw 1.000x -> basin 1.000x
#>    Intensity effect : gamma_0=0.244 -> gamma=0.244 (+24% per +1°C)
#>    Local response   : redistribution 0.000 | range 1.000x–1.000x
#>    Perturbation     : disabled
#> 
#> ── Simulation ────────────────────────────────────────────────
#>    Synthetic years  : 2,000
#>    Random seed      : 1540460832
#>    Rate scaling     : 1.000x
#>    ✔ Done
#> 
#> ── Run metadata ──────────────────────────────────────────────
#>    seed=1540460832
#>    lambda_mode=target
#> ════════════════════════════════════════════════════════════════════════
#>   HAZARD MODEL VALIDATION SUITE
#>   Holdout: 10 yr  |  Sim: 2,000 yr  |  CI: 95%  |  Seed: 42
#> ════════════════════════════════════════════════════════════════════════
#> 
#> ── TIER 1  HINDCAST (return levels) ────────────────────────────────────
#>   Simulated return levels are compared against observed GEV estimates;
#>   ✓ means the simulated value falls within the observed 95% confidence interval.
#> 
#>   Training periods and GEV fits:
#>     Saba           1970–2025   obs: ξ=-0.30, n=19      sim: ξ=-0.25, n=740
#>     St. Eustatius  1970–2024   obs: ξ=-0.25, n=20      sim: ξ=-0.07, n=771
#>     St. Martin     1970–2025   obs: ξ=-0.30, n=23      sim: ξ=-0.26, n=713
#> 
#>   Return-level comparison (kt):
#>   Location       RP   Obs     Sim     95% CI (obs)   Bias        
#>   ──────────────────────────────────────────────────────────────
#>   Saba           5    44      49      [ 35,  52]     +11%    ✓   
#>   Saba           10   53      58      [ 46,  65]     +10%    ✓   
#>   Saba           25   60      66      [ 55,  91]     +11%    ✓   
#>   Saba           50   64      71      [ 61, 116]     +11%    ✓   
#> 
#>   St. Eustatius  5    45      50      [ 35,  51]     +12%    ✓   
#>   St. Eustatius  10   52      61      [ 46,  64]     +16%    ✓   
#>   St. Eustatius  25   59      72      [ 53,  87]     +22%    ✓   
#>   St. Eustatius  50   63      79      [ 59, 110]     +26%    ✓   
#> 
#>   St. Martin     5    49      47      [ 36,  59]     -3%     ✓   
#>   St. Martin     10   59      57      [ 50,  79]     -3%     ✓   
#>   St. Martin     25   68      66      [ 61, 116]     -3%     ✓   
#>   St. Martin     50   73      71      [ 70, 154]     -3%     ✓   
#> 
#> ── TIER 2  STORM RATE CHECK ────────────────────────────────────────────
#>   Calibrated annual passage rates are compared against published
#>   HURDAT2/IBTrACS climatology bounds; rates outside range are flagged.
#> 
#>   Location       Class  Fitted   Ref      Ratio    Expected Status  
#>   ─────────────────────────────────────────────────────────────────
#>   Saba           TS     0.489    1.450    0.34     0.64     ✓       
#> 
#>   St. Eustatius  TS     0.500    NA       NA       NA       ✗       
#> 
#>   St. Martin     TS     0.393    NA       NA       NA       ✗       
#> 
#>   Saba           HUR    0.064    0.550    0.12     0.30     ✓       
#> 
#>   St. Eustatius  HUR    0.087    NA       NA       NA       ✗       
#> 
#>   St. Martin     HUR    0.125    NA       NA       NA       ✗       
#> 
#>   ▸ 4 rate(s) outside reference bounds — consider adjusting search_radius_km or historical_start_year.
#> 
#> ── TIER 3  WIND FIELD SPOT-CHECKS ──────────────────────────────────────
#>   Holland wind-profile estimates are compared against historical station
#>   observations; ✓ means the model wind is within the observational tolerance.
#> 
#>   Mean bias: -16.9 kt   MAE: 16.9 kt   RMSE: 16.9 kt   (1 storms)
#> 
#>   Storm      Location       Obs     Model   Bias     Tolerance     
#>   ───────────────────────────────────────────────────────────────
#>   IRMA       Saba           80 kt   63 kt   -21%     C ±35%    ✓   
#> 
#>   ▸ Quality tiers — A: direct station (±15%)  B: 10-min converted (±25%)  C: estimated (±35%)
#> 
#> ════════════════════════════════════════════════════════════════════════
#>   SUMMARY
#>   A model with failing tiers should be recalibrated before drawing
#>   conclusions from climate scenario analyses.
#> 
#>   Tier 1  Hindcast      12 / 12  within 95% CI      ✓
#>   Tier 2  Rate check     2 /  6  within bounds       ✗
#>   Tier 3  Wind field     1 /  1  within tolerance    ✓
#> 
#>   Climate  counts      1.006x basin (raw 1.000x)  baseline
#>            HUR share   -0.10 pp
#>            redistrib.  RMSE 0.49 pp  max 0.69 pp
#> 
#> ════════════════════════════════════════════════════════════════════════

# Inspect results
result$val$summary                # Compact hindcast summary diagnostics
#> # A tibble: 3 × 4
#>   location      storm_class delta_top1_p50 delta_overall_p99
#>   <chr>         <chr>                <dbl>             <dbl>
#> 1 Saba          TS                    5.78              16.3
#> 2 St. Eustatius TS                    2.52              13.1
#> 3 St. Martin    TS                  -29.9              -15.4
result$val$hindcast$comparison    # Detailed return-level table
#> NULL
result$artifacts$plots            # Paths to saved figures
#> NULL

# --- Stricter validation with 95% CI ---
result_95 <- validate_hazard_model(
  cfg            = make_hazard_cfg(simulation_years = 4000),
  targets        = targets_df,
  validation_cfg = make_validation_cfg(
    conf_level     = 0.95,
    n_sim          = 10000,
    return_periods = c(10, 25, 50, 100)
  )
)
#> 
#> ── Loading data ──────────────────────────────────────────────
#>    Input file       : ibtracs_1970.csv
#>    Raw NA input     : 49,905 track points (1970–2025)
#>    Model start year : 1970
#> 
#> ── Target/location filtering ─────────────────────────────────
#>    Location          TS    HUR
#>    Saba              26      3
#>    St. Eustatius     27      4
#>    St. Martin        29      7
#>    Model sample     : 704 target-event records (1970–2025) across 3 locations
#> 
#> ── Rate calibration ──────────────────────────────────────────
#>    Adjustments      : 6 location/class adjustments (mode=target, clamped=0, missing_ref=4)
#>    Some location/class rates were increased to match reference rates; set advanced$lambda_scaling_mode='down_only' to prevent upscaling.
#> 
#> ── Climate ───────────────────────────────────────────────────
#>    Climate mode     : baseline (delta_sst = 0)
#>    Climate input    : scenario helper
#>    Climate scenario : stationary
#>    SST baseline     : 1991–2020
#>    Sensitivity mode : fixed
#>    Rate effect      : beta_0=0.60 -> beta=0.60 (+82% per +1°C) [guardrail]
#>    Count regime     : baseline | raw 1.000x -> basin 1.000x
#>    Intensity effect : gamma_0=0.244 -> gamma=0.244 (+24% per +1°C)
#>    Local response   : redistribution 0.000 | range 1.000x–1.000x
#>    Perturbation     : disabled
#> 
#> ── Simulation ────────────────────────────────────────────────
#>    Synthetic years  : 4,000
#>    Random seed      : 7117163
#>    Rate scaling     : 1.000x
#>    ✔ Done
#> 
#> ── Run metadata ──────────────────────────────────────────────
#>    seed=7117163
#>    lambda_mode=target
#> ════════════════════════════════════════════════════════════════════════
#>   HAZARD MODEL VALIDATION SUITE
#>   Holdout: 10 yr  |  Sim: 10,000 yr  |  CI: 95%  |  Seed: 42
#> ════════════════════════════════════════════════════════════════════════
#> 
#> ── TIER 1  HINDCAST (return levels) ────────────────────────────────────
#>   Simulated return levels are compared against observed GEV estimates;
#>   ✓ means the simulated value falls within the observed 95% confidence interval.
#> 
#>   Training periods and GEV fits:
#>     Saba           1970–2025   obs: ξ=-0.30, n=19      sim: ξ=-0.23, n=3723
#>     St. Eustatius  1970–2024   obs: ξ=-0.25, n=20      sim: ξ=-0.08, n=3935
#>     St. Martin     1970–2025   obs: ξ=-0.30, n=23      sim: ξ=-0.23, n=3688
#> 
#>   Return-level comparison (kt):
#>   Location       RP   Obs     Sim     95% CI (obs)   Bias        
#>   ──────────────────────────────────────────────────────────────
#>   Saba           10   53      60      [ 46,  65]     +14%    ✓   
#>   Saba           25   60      69      [ 55,  91]     +16%    ✓   
#>   Saba           50   64      74      [ 61, 116]     +17%    ✓   
#>   Saba           100  67      79      [ 66, 150]     +18%    ✓   
#> 
#>   St. Eustatius  10   52      63      [ 46,  64]     +20%    ✓   
#>   St. Eustatius  25   59      74      [ 53,  87]     +26%    ✓   
#>   St. Eustatius  50   63      82      [ 59, 110]     +30%    ✓   
#>   St. Eustatius  100  66      89      [ 63, 141]     +35%    ✓   
#> 
#>   St. Martin     10   59      60      [ 50,  79]     +3%     ✓   
#>   St. Martin     25   68      70      [ 61, 116]     +4%     ✓   
#>   St. Martin     50   73      76      [ 70, 154]     +5%     ✓   
#>   St. Martin     100  76      81      [ 78, 185]     +6%     ✓   
#> 
#> ── TIER 2  STORM RATE CHECK ────────────────────────────────────────────
#>   Calibrated annual passage rates are compared against published
#>   HURDAT2/IBTrACS climatology bounds; rates outside range are flagged.
#> 
#>   Location       Class  Fitted   Ref      Ratio    Expected Status  
#>   ─────────────────────────────────────────────────────────────────
#>   Saba           TS     0.489    1.450    0.34     0.64     ✓       
#> 
#>   St. Eustatius  TS     0.500    NA       NA       NA       ✗       
#> 
#>   St. Martin     TS     0.393    NA       NA       NA       ✗       
#> 
#>   Saba           HUR    0.064    0.550    0.12     0.30     ✓       
#> 
#>   St. Eustatius  HUR    0.087    NA       NA       NA       ✗       
#> 
#>   St. Martin     HUR    0.125    NA       NA       NA       ✗       
#> 
#>   ▸ 4 rate(s) outside reference bounds — consider adjusting search_radius_km or historical_start_year.
#> 
#> ── TIER 3  WIND FIELD SPOT-CHECKS ──────────────────────────────────────
#>   Holland wind-profile estimates are compared against historical station
#>   observations; ✓ means the model wind is within the observational tolerance.
#> 
#>   Mean bias: -16.9 kt   MAE: 16.9 kt   RMSE: 16.9 kt   (1 storms)
#> 
#>   Storm      Location       Obs     Model   Bias     Tolerance     
#>   ───────────────────────────────────────────────────────────────
#>   IRMA       Saba           80 kt   63 kt   -21%     C ±35%    ✓   
#> 
#>   ▸ Quality tiers — A: direct station (±15%)  B: 10-min converted (±25%)  C: estimated (±35%)
#> 
#> ════════════════════════════════════════════════════════════════════════
#>   SUMMARY
#>   A model with failing tiers should be recalibrated before drawing
#>   conclusions from climate scenario analyses.
#> 
#>   Tier 1  Hindcast      12 / 12  within 95% CI      ✓
#>   Tier 2  Rate check     2 /  6  within bounds       ✗
#>   Tier 3  Wind field     1 /  1  within tolerance    ✓
#> 
#>   Climate  counts      0.989x basin (raw 1.000x)  baseline
#>            HUR share   -0.18 pp
#>            redistrib.  RMSE 0.79 pp  max 1.01 pp
#> 
#> ════════════════════════════════════════════════════════════════════════

# --- With climate conditioning ---
result_climate <- validate_hazard_model(
  cfg = make_hazard_cfg(
    climate = make_climate_cfg(scenario = "ssp245", target_year = 2050)
  ),
  targets = targets_df
)
#> 
#> ── Loading data ──────────────────────────────────────────────
#>    Input file       : ibtracs_1970.csv
#>    Raw NA input     : 49,905 track points (1970–2025)
#>    Model start year : 1970
#> 
#> ── Target/location filtering ─────────────────────────────────
#>    Location          TS    HUR
#>    Saba              26      3
#>    St. Eustatius     27      4
#>    St. Martin        29      7
#>    Model sample     : 704 target-event records (1970–2025) across 3 locations
#> 
#> ── Rate calibration ──────────────────────────────────────────
#>    Adjustments      : 6 location/class adjustments (mode=target, clamped=0, missing_ref=4)
#>    Some location/class rates were increased to match reference rates; set advanced$lambda_scaling_mode='down_only' to prevent upscaling.
#> 
#> ── Climate ───────────────────────────────────────────────────
#>    Climate mode     : future (delta_sst = +0.50°C)
#>    Climate input    : scenario helper
#>    Climate scenario : SSP245
#>    Target year      : 2050.0
#>    SST baseline     : 1991–2020
#>    Sensitivity mode : fixed
#>    Rate effect      : beta_0=0.60 -> beta=0.60 (+82% per +1°C) [guardrail]
#>    Count regime     : delta_only_fixed_guardrail | raw 1.350x -> basin 1.028x
#>    Intensity effect : gamma_0=0.244 -> gamma=0.244 (+24% per +1°C)
#>    Local response   : redistribution 0.060 | range 0.977x–1.060x
#>    Perturbation     : disabled
#> 
#> ── Simulation ────────────────────────────────────────────────
#>    Synthetic years  : 1,000
#>    Random seed      : 1566150152
#>    Rate scaling     : 1.028x
#>    ✔ Done
#> 
#> ── Run metadata ──────────────────────────────────────────────
#>    seed=1566150152
#>    lambda_mode=target
#> ════════════════════════════════════════════════════════════════════════
#>   HAZARD MODEL VALIDATION SUITE
#>   Holdout: 10 yr  |  Sim: 1,000 yr  |  CI: 95%  |  Seed: 42
#> ════════════════════════════════════════════════════════════════════════
#> 
#> ── TIER 1  HINDCAST (return levels) ────────────────────────────────────
#>   Simulated return levels are compared against observed GEV estimates;
#>   ✓ means the simulated value falls within the observed 95% confidence interval.
#> 
#>   Training periods and GEV fits:
#>     Saba           1970–2025   obs: ξ=-0.30, n=19      sim: ξ=-0.24, n=381
#>     St. Eustatius  1970–2024   obs: ξ=-0.25, n=20      sim: ξ=-0.05, n=402
#>     St. Martin     1970–2025   obs: ξ=-0.30, n=23      sim: ξ=-0.21, n=379
#> 
#>   Return-level comparison (kt):
#>   Location       RP   Obs     Sim     95% CI (obs)   Bias        
#>   ──────────────────────────────────────────────────────────────
#>   Saba           5    44      52      [ 35,  52]     +16%    ✓   
#>   Saba           10   53      61      [ 46,  65]     +16%    ✓   
#>   Saba           25   60      70      [ 55,  91]     +17%    ✓   
#>   Saba           50   64      75      [ 61, 116]     +18%    ✓   
#> 
#>   St. Eustatius  5    45      53      [ 35,  51]     +18%    ✗   
#>   St. Eustatius  10   52      64      [ 46,  64]     +22%    ✓   
#>   St. Eustatius  25   59      76      [ 53,  87]     +29%    ✓   
#>   St. Eustatius  50   63      84      [ 59, 110]     +34%    ✓   
#> 
#>   St. Martin     5    49      51      [ 36,  59]     +4%     ✓   
#>   St. Martin     10   59      63      [ 50,  79]     +7%     ✓   
#>   St. Martin     25   68      74      [ 61, 116]     +9%     ✓   
#>   St. Martin     50   73      80      [ 70, 154]     +11%    ✓   
#> 
#> ── TIER 2  STORM RATE CHECK ────────────────────────────────────────────
#>   Calibrated annual passage rates are compared against published
#>   HURDAT2/IBTrACS climatology bounds; rates outside range are flagged.
#> 
#>   Location       Class  Fitted   Ref      Ratio    Expected Status  
#>   ─────────────────────────────────────────────────────────────────
#>   Saba           TS     0.489    1.450    0.34     0.64     ✓       
#> 
#>   St. Eustatius  TS     0.500    NA       NA       NA       ✗       
#> 
#>   St. Martin     TS     0.393    NA       NA       NA       ✗       
#> 
#>   Saba           HUR    0.064    0.550    0.12     0.30     ✓       
#> 
#>   St. Eustatius  HUR    0.087    NA       NA       NA       ✗       
#> 
#>   St. Martin     HUR    0.125    NA       NA       NA       ✗       
#> 
#>   ▸ 4 rate(s) outside reference bounds — consider adjusting search_radius_km or historical_start_year.
#> 
#> ── TIER 3  WIND FIELD SPOT-CHECKS ──────────────────────────────────────
#>   Holland wind-profile estimates are compared against historical station
#>   observations; ✓ means the model wind is within the observational tolerance.
#> 
#>   Mean bias: -16.9 kt   MAE: 16.9 kt   RMSE: 16.9 kt   (1 storms)
#> 
#>   Storm      Location       Obs     Model   Bias     Tolerance     
#>   ───────────────────────────────────────────────────────────────
#>   IRMA       Saba           80 kt   63 kt   -21%     C ±35%    ✓   
#> 
#>   ▸ Quality tiers — A: direct station (±15%)  B: 10-min converted (±25%)  C: estimated (±35%)
#> 
#> ════════════════════════════════════════════════════════════════════════
#>   SUMMARY
#>   A model with failing tiers should be recalibrated before drawing
#>   conclusions from climate scenario analyses.
#> 
#>   Tier 1  Hindcast      11 / 12  within 95% CI      ✓
#>   Tier 2  Rate check     2 /  6  within bounds       ✗
#>   Tier 3  Wind field     1 /  1  within tolerance    ✓
#> 
#>   Climate  counts      1.017x basin (raw 1.350x)  delta_only_fixed_guardrail
#>            HUR share   +2.30 pp
#>            redistrib.  RMSE 1.97 pp  max 2.64 pp
#> 
#> ════════════════════════════════════════════════════════════════════════
```
