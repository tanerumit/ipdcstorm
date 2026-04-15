# Model Validation Report

**Package**: ipdcstorm
**Generated**: 2026-04-15
**Output directory**: output/validation/

## 1. Executive Summary

Tier 1 passes overall (14 of 15 return-level checks within the observed 95% CI; 93.3%), Tier 2 passes with no flagged rate checks, and Tier 3 flags a wind-field issue because St_Martin shows a +22.4 kt bias for IRMA. On that basis, the baseline model is close to usable, but it is not yet cleanly validated for climate scenario analysis because all three tiers do not pass simultaneously. The main concern is not the storm-rate fit; it is near-eye wind overprediction at St_Martin for a major event.

From a technical standpoint, this pattern suggests that the stochastic event generation and annual-rate calibration are broadly coherent, while the remaining weakness sits in the physical wind translation from storm properties to site exposure. Scientifically, that matters because climate-scenario analysis is most credible only when both the frequency structure and the local intensity mapping are defensible.

## 2. Tier 1 — Hindcast Return Levels

Tier 1 meets the package pass threshold of 80%. Saba passes all five return periods with moderate positive bias, suggesting a mild tendency to overpredict site wind return levels by about 9.9% to 11.8%. Statia also passes all five return periods, but the bias is consistently positive and increases with return period, from 12.5% at 5 years to 30.8% at 100 years, which points to an increasingly heavy simulated upper tail rather than a rate-fit problem. St_Martin is the only marginal location: four of five return periods pass, and the only miss is the 100-year return level, where the simulated value falls just below the observed CI. That pattern is more consistent with tail-intensity fitting than with a broad calibration failure.

Technically, Tier 1 is the most relevant check on whether the fitted stationary model can recreate the observed severity distribution when asked to predict outside its training window. Scientifically, the main caution is that wide bootstrap confidence intervals at long return periods can hide moderate tail bias, so a formal pass here should be interpreted as evidence of broad consistency rather than proof that the extreme-value tail is perfectly specified.

| Location | Checks in CI | Pass rate | Mean bias (%) | Bias range (%) | Tier status |
| --- | ---: | ---: | ---: | ---: | --- |
| Saba | 5 / 5 | 100.0% | 10.8 | 9.9 to 11.8 | ✓ |
| St_Martin | 4 / 5 | 80.0% | -2.7 | -3.4 to -2.2 | ✓ marginal |
| Statia | 5 / 5 | 100.0% | 21.5 | 12.5 to 30.8 | ✓ |

## 3. Tier 2 — Storm Rate Check

Tier 2 passes cleanly. All six fitted annual passage rates are marked `OK`, so there is no evidence that `search_radius_km = 800` or `historical_start_year = 1970` is producing obvious over- or under-counting relative to the package’s reference climatology. Because the rate fit is already inside the accepted bounds, there is no immediate reason to change the search radius or truncate the record further for this baseline configuration.

Technically, this tier isolates the count process from the intensity model, which is useful because it tests whether the event-selection geometry and historical sampling window are reasonable before interpreting any return-level result. Scientifically, a clean Tier 2 pass strengthens the case that the package is sampling the correct regional storm climate, so the main residual uncertainty shifts away from event frequency and toward intensity and wind-field structure.

| Location | Storm class | Fitted λ | Reference λ | Target λ | Flag |
| --- | --- | ---: | ---: | ---: | --- |
| St_Martin | TS | 0.393 | 0.800 | 0.720 | OK |
| Saba | TS | 0.489 | 1.450 | 0.935 | OK |
| Statia | TS | 0.500 | 1.450 | 0.935 | OK |
| St_Martin | HUR | 0.125 | 0.400 | 0.180 | OK |
| Saba | HUR | 0.064 | 0.550 | 0.165 | OK |
| Statia | HUR | 0.087 | 0.550 | 0.165 | OK |

## 4. Tier 3 — Wind Field Spot-Check

Tier 3 does not pass. Only two rows contain valid model-versus-observed comparisons, so this tier is data-sparse, but one of those two is a material miss: IRMA at St_Martin is overpredicted by 22.4 kt, which exceeds the skill threshold for an event-level flag. Saba shows a -16.4 kt bias for IRMA, which is tolerable but still notable. No valid wind-field comparison is available for Statia in the saved CSV, so wind attenuation there remains unverified.

Technically, Tier 3 is a structural check on the Holland-profile wind model rather than on the stochastic calibration, and it is therefore the best diagnostic for whether site winds are physically plausible once a storm track is given. Scientifically, the sparse sample limits confidence, but the St_Martin IRMA miss is still important because high-impact hazard studies are often dominated by a few extreme close-passage events, exactly where wind-profile assumptions matter most.

| Location | Valid events | Mean bias (kt) | RMSE (kt) | Tier status |
| --- | ---: | ---: | ---: | --- |
| Saba | 1 | -16.4 | 16.4 | Tolerable |
| St_Martin | 1 | 22.4 | 22.4 | ✗ |
| Statia | 0 | NA | NA | — |

The flagged event is IRMA at St_Martin (+22.4 kt). Bias direction is mixed across the two available sites, so this is not yet evidence of a basin-wide bias; it is evidence that the near-eye wind-field formulation is too strong for at least one critical local case.

## 5. Calibration Recommendations

**Issue**: Tier 3 — St_Martin wind-field spot-check is flagged high for IRMA.  
**Likely cause**: the Holland-profile configuration used in validation is overpredicting near-eye site wind at very short distance, even though the statistical rate calibration is acceptable.  
**Recommended action**: keep `make_hazard_cfg(search_radius_km = 800, historical_start_year = 1970)` unchanged for now, and re-run validation with wind-field diagnostics enabled to test whether the St_Martin miss is sensitive to the package’s wind-field mode. In this codebase that means re-running validation after setting `options(ipdcstorm.wind_field_mode = "diagnostic_new")` before `run_validation_suite()`.

**Issue**: Tier 1 — Statia passes but shows a strong positive bias gradient at long return periods.  
**Likely cause**: simulated tail intensity is somewhat heavy relative to the observed fit, even though the values remain inside the wide observed CI.  
**Recommended action**: keep the current rate configuration, but review the highest-intensity event library and repeat the validation with the same `make_validation_cfg(holdout_years = 10, n_sim = 2000, return_periods = c(5, 10, 25, 50, 100), conf_level = 0.95)` to confirm the long-return-period bias is stable rather than seed-specific.

From a technical perspective, the calibration priority should stay on the component that actually failed: the wind-field translation, not the passage-rate fit. Scientifically, that is the more defensible sequence because adjusting storm counts to fix what is currently a local wind-bias problem would risk degrading a Tier 2 component that already agrees with reference climatology.

## 6. Configuration Used

The script `inst/extcode/02-model-validation.R` used `SEED = 42`, `N_SIM = 2000`, `HIST_START_YEAR = 1970`, `SEARCH_RADIUS_KM = 800`, `HOLDOUT_YEARS = 10`, `RETURN_PERIODS = c(5, 10, 25, 50, 100)`, and `CONF_LEVEL = 0.95`, with outputs written to `output/validation/`. The detailed tier CSVs were available, but `validation_summary.csv` did not contain the compact pass/fail schema expected by the skill, so the tier conclusions above were derived from `hindcast_return_levels.csv`, `rate_check.csv`, and `wind_field.csv`.

This configuration is scientifically reasonable for a baseline stationary validation because it separates calibration from climate forcing and uses a holdout long enough to test out-of-sample behavior. The main technical limitation is not the chosen settings themselves, but the small number of valid Tier 3 observational comparisons available under this configuration.
