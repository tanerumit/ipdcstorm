# Unified climate workflow

- Date: 2026-03-12
- Scope: hazard config construction, climate resolution, `run_hazard_model()`, and direct validation callers.

## API change

Climate configuration is now embedded inside `hazard_cfg` and resolved only from `cfg$climate`.
`run_hazard_model()` no longer accepts a separate `climate` argument.
Baseline runs must use `make_climate_cfg(scenario = "stationary")`.

## Why baseline and future now share one pathway

The previous interface had two baseline representations:

- `climate = NULL`, which bypassed climate resolution completely.
- `make_climate_cfg(scenario = "stationary")`, which ran the resolver with `delta_sst = 0`.

That created two metadata shapes and two execution paths for conceptually similar runs.
The refactor removes the bypass so both stationary and future runs calibrate baseline sensitivities, resolve scenario metadata, and return climate metadata through the same code path.

## Expected result changes

Legacy `climate = NULL` runs can now differ from the unified stationary baseline because the old off-path hard-set `beta_sst = 0`, `gamma = 0`, and omitted resolved climate metadata.
The new stationary baseline resolves `beta_0` and `gamma_0` through the same calibration workflow used for future scenarios, while still keeping `delta_sst = 0`.

These result changes are expected and acceptable because they are the direct consequence of removing the legacy climate-off mode and making the full hazard configuration explicit in `cfg`.
