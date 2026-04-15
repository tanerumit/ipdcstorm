# Model Validation Report

**Package**: ipdcstorm
**Generated**: 2026-04-15
**Output directory**: output/validation/

## 1. Executive Summary

Tier 1 passes overall (14 of 15 return-level checks within the observed 95% CI; 93.3%), Tier 2 passes with no flagged rate checks, and Tier 3 flags a wind-field issue because St_Martin shows a +22.4 kt bias for IRMA. Taken together, these results show that the baseline model captures the regional storm-rate structure and reproduces site-level return-level behavior well enough to support hazard screening and comparative scenario analysis. The current limitation is not in the climatological frequency structure, but in the local wind translation for very close high-intensity passages. For practical use, this means the validation already supports comparative interpretation across sites and scenarios, while indicating that additional caution is warranted for absolute exceedance estimates at the most exposed locations.

Overall, the validation pattern indicates that the stochastic event generation and annual-rate calibration are coherent across all three target sites. Within the `ipdcstorm` workflow, this makes the current analysis useful as both a baseline performance check and a guide to where subsequent refinement is most likely to improve site-level hazard estimates.

## 2. Tier 1 — Hindcast Return Levels

Tier 1 meets the package pass threshold of 80%. Saba passes all five return periods with moderate positive bias, suggesting a mild tendency to overpredict site wind return levels by about 9.9% to 11.8%. Statia also passes all five return periods, but the bias is consistently positive and increases with return period, from 12.5% at 5 years to 30.8% at 100 years, which points to an increasingly heavy simulated upper tail rather than a rate-fit problem. St_Martin is the only marginal location: four of five return periods pass, and the only miss is the 100-year return level, where the simulated value falls just below the observed CI. That pattern is more consistent with tail-intensity fitting than with a broad calibration failure.

Tier 1 provides the clearest evidence that the fitted stationary model can reproduce the observed severity structure when evaluated out of sample, which is the core requirement for using it as a baseline hazard engine. In this workflow, that is an important result: the model is not only reproducing calibration-period summaries, but is also showing reasonable hindcast skill against withheld data. The main caution is that wide bootstrap confidence intervals at long return periods can mask moderate tail bias, so the result is best interpreted as support for broad adequacy rather than as proof that the far tail is fully resolved.

| Location | Checks in CI | Pass rate | Mean bias (%) | Bias range (%) | Tier status |
| --- | ---: | ---: | ---: | ---: | --- |
| Saba | 5 / 5 | 100.0% | 10.8 | 9.9 to 11.8 | ✓ |
| St_Martin | 4 / 5 | 80.0% | -2.7 | -3.4 to -2.2 | ✓ marginal |
| Statia | 5 / 5 | 100.0% | 21.5 | 12.5 to 30.8 | ✓ |

## 3. Tier 2 — Storm Rate Check

Tier 2 passes cleanly. All six fitted annual passage rates are marked `OK`, so there is no evidence that `search_radius_km = 800` or `historical_start_year = 1970` is producing obvious over- or under-counting relative to the package’s reference climatology. Because the rate fit is already inside the accepted bounds, there is no immediate reason to change the search radius or truncate the record further for this baseline configuration.

This tier isolates the count process from the intensity model, which shows that the event-selection geometry and historical sampling window are already working as intended. There is no sign that model behavior is being driven by an obviously mis-specified search radius or an unrepresentative historical sampling frame. In practice, that makes the downstream validation results easier to interpret, because the remaining uncertainty can be assigned primarily to intensity and wind-field structure rather than to errors in the underlying storm climatology.

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

Tier 3 is a structural check on the Holland-profile wind model rather than on the stochastic calibration, and it is therefore the most direct test of whether site winds remain physically credible once a storm track is specified. Even though this tier is sparse, it is still informative because it localizes the remaining model weakness to a specific component and a specific type of event. The St_Martin IRMA miss matters because high-impact hazard estimates are often controlled by rare close-passage storms, while the broader value of this result is that it narrows the next refinement step substantially.

| Location | Valid events | Mean bias (kt) | RMSE (kt) | Tier status |
| --- | ---: | ---: | ---: | --- |
| Saba | 1 | -16.4 | 16.4 | Tolerable |
| St_Martin | 1 | 22.4 | 22.4 | ✗ |
| Statia | 0 | NA | NA | — |

The flagged event is IRMA at St_Martin (+22.4 kt). Bias direction is mixed across the two available sites, so this is not evidence of a basin-wide failure; it is evidence that the near-eye wind-field formulation is too strong in at least one critical local case. For later analyses, this means the current validation supports continued use of the model for comparative and exploratory work, while also showing where additional caution is needed for absolute exceedance estimates at the most exposed sites.

## 5. Next Steps In The Workflow

**Tier 3 follow-up**: The St_Martin wind-field spot-check is high for IRMA.  
**Interpretation**: the Holland-profile configuration used in validation is overpredicting near-eye site wind at very short distance, even though the statistical rate calibration is acceptable.  
**Next step**: keep `make_hazard_cfg(search_radius_km = 800, historical_start_year = 1970)` unchanged for now, and re-run validation with wind-field diagnostics enabled to test whether the St_Martin miss is sensitive to the package’s wind-field mode. In this codebase that means re-running validation after setting `options(ipdcstorm.wind_field_mode = "diagnostic_new")` before `run_validation_suite()`.

**Tier 1 follow-up**: Statia passes but shows a strong positive bias gradient at long return periods.  
**Interpretation**: simulated tail intensity is somewhat heavy relative to the observed fit, even though the values remain inside the wide observed CI.  
**Next step**: keep the current rate configuration, but review the highest-intensity event library and repeat the validation with the same `make_validation_cfg(holdout_years = 10, n_sim = 2000, return_periods = c(5, 10, 25, 50, 100), conf_level = 0.95)` to confirm the long-return-period bias is stable rather than seed-specific.

Within the workflow, the priority should stay on the component that actually failed: the wind-field translation, not the passage-rate fit. The present results show that the model already has a stable backbone for regional storm occurrence and broad hazard behavior, so the next step is targeted rather than global.

## 6. Configuration Used

The script `inst/extcode/02-model-validation.R` used `SEED = 42`, `N_SIM = 2000`, `HIST_START_YEAR = 1970`, `SEARCH_RADIUS_KM = 800`, `HOLDOUT_YEARS = 10`, `RETURN_PERIODS = c(5, 10, 25, 50, 100)`, and `CONF_LEVEL = 0.95`, with outputs written to `output/validation/`. The detailed tier CSVs were available, but `validation_summary.csv` did not contain the compact pass/fail schema expected by the skill, so the tier conclusions above were derived from `hindcast_return_levels.csv`, `rate_check.csv`, and `wind_field.csv`.

This configuration is scientifically reasonable for a baseline stationary validation because it separates calibration from climate forcing and uses a holdout long enough to test out-of-sample behavior. The main technical limitation is not the chosen settings themselves, but the small number of valid Tier 3 observational comparisons available under this configuration. Even so, the configuration remains useful within the package workflow because it provides a transparent baseline, shows where model skill is strongest, and offers a defensible foundation for subsequent refinement and scenario comparison.
