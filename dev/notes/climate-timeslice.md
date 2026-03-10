## Climate workflow: stationary time-slice delta shift

The hazard model climate workflow now uses a stationary time-slice SST shift instead of a transient annual forcing path.

- `delta_sst` is a single scalar in degC representing the mean SST over a future period minus the mean SST over a baseline period, with both periods referenced to the same climatology.
- Estimation is now a pre-step. `prepare_sst_data()`, `estimate_beta_sst()`, and `estimate_gamma_intensity()` remain standalone helpers and are not called inside `run_hazard_model()`.
- Scenario lookup is now a separate pre-step through `get_scenario_delta()`, which interpolates a scalar shift at the midpoint of a requested future period.
- Model response is now a separate pre-step through `make_climate_response()`, which packages `beta_sst`, `gamma`, and optional Level 3 perturbation factors.
- Runner input is now explicit through `make_climate_input()`, which combines a scalar shift and response package for a stationary hazard run.

Methodologically, all simulation years are repeated draws from one shifted climate state. Year-to-year variability comes from the Gamma activity factor and Poisson/Binomial sampling, not from time-varying SST forcing.

Reproducibility also changed:

- `run_hazard_model()` now sets the RNG seed once at entry.
- If `seed = NULL`, the runner generates a seed, uses it immediately, and records it in `run_metadata$seed`.
- No transient SST trajectory is generated or consumed in the runner or simulator path.
