# Climate frequency guardrail note

## Cause

Future count inflation was traced to climate-rate calibration using target-stacked annual counts in `run_hazard_model()`. When the same basin storm affected multiple targets, those target-conditioned annual counts were summed during `beta_0` calibration, making the fitted SST sensitivity depend on the target set instead of basin activity.

## Fix

- `run_hazard_model()` now derives climate-calibration counts from basin-consistent, de-duplicated storm-year records via `compute_annual_counts(events_all, ...)`, which removes target duplication from the annual activity series.
- `estimate_beta_sst()` now rejects annual-count inputs that still contain a `location` column, so target-conditioned series cannot silently enter the climate calibrator.
- The exact annual total count series used for `beta_0` is returned as `annual_count_series` with provenance `annual_count_source = "basin_unique_storm_year_counts"`.

## Guardrails

- `beta_0` plausibility guardrail: if the shrunk estimate exceeds `1.2 1/degC`, it falls back to the literature prior `0.6 1/degC`.
- Shrinkage is strengthened by both overlap length and coefficient uncertainty, so large-standard-error MLEs are pulled harder toward the prior.
- Future SST rate multiplier guardrail: if `sst_scale = exp(beta_sst * delta_sst)` exceeds `4x`, `beta_sst` is reduced to the corresponding bound `log(4) / delta_sst`.

## Rationale

- Units remain explicit: `delta_sst` in degC, `beta_sst` in `1/degC`, annual counts in storms/year.
- Basin-scale climate sensitivity should be invariant to the number or geometry of targets; only local hazard rates should vary with target filtering.
- The guardrails are intended to block implausible order-of-magnitude frequency inflation while preserving stationary behavior and leaving the wind/intensity physics unchanged.
