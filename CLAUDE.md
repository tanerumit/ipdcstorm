# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`ipdcstorm` is an R package for hurricane hazard assessment targeting the Dutch Caribbean islands (Saba, St. Eustatius, St. Martin). It uses IBTrACS North Atlantic best-track data to estimate site-specific tropical cyclone risks via a Holland-type wind profile model with climate change extensions.

## Commands

```bash
# Run all tests
R -q -e "devtools::test()"

# Run a scoped subset (preferred — substitute <REGEX> for the filter)
R -q -e "devtools::test(filter = '<REGEX>')"

# Run a single test file directly
Rscript -e "testthat::test_file('tests/testthat/test-smoke.R')"

# Parse-check a touched source file
Rscript -e "parse(file='R/hazard_core.R')"

# Regenerate documentation (Rd + NAMESPACE)
R -q -e "devtools::document()"

# Full package check
R -q -e "devtools::check()"

# Load package for interactive development
R -q -e "devtools::load_all()"
```

## Architecture

The package is structured as a sequential pipeline:

1. **`R/hazard_ibtracs.R`** — IBTrACS ingestion: `read_ibtracs_clean()` reads/cleans IBTrACS CSV, standardizes USA best-track fields (wind, pressure, R34/R50/R64 directional radii, RMW), and applies basin/season filters with strict fail-fast validation.

2. **`R/hazard_core.R`** — Core wind field computations: geometry helpers (great-circle distance, bearing, quadrant), Holland-type radial wind profile (`.estimate_site_wind_holland()`), RMW estimation (Knaff & Zehr 2007), R34 climatological infill (Knaff et al. 2015), forward-motion asymmetry (`.add_forward_motion_asymmetry()`). Entry point: `compute_site_winds_full(df, target_lat, target_lon)`.

3. **`R/hazard_run.R`** — Model orchestration: `make_hazard_cfg()` (user-facing config constructor), `run_hazard_model(cfg, targets, seed)` (end-to-end pipeline across multiple target locations). Returns a list with `sim`, `events`, `rates`, `fit`, `cfg`, and `run_metadata`.

4. **`R/hazard_climate.R`** — Climate change workflow. Three-tier hierarchy:
   - **L1** (`make_sst_cfg` + SST activity scaling): Poisson/NegBin regression of MDR SST → annual activity factor
   - **L2** (gamma intensity): TS/HUR split adjustment via SST-conditioned `estimate_gamma_intensity()`
   - **L3** (`perturb_event()`): direct storm property perturbation (optional expert extension)
   Built-in MDR SST (ERSST v5, 1970–2024) via `get_mdr_sst_builtin()`.

5. **`R/hazard_downscale.R`** — Temporal downscaling: `build_event_library()` → `generate_daily_hazard_impact()`. Produces daily wind + damage-forcing time series from stochastic event sampling.

6. **`R/hazard_validation.R`** — Three-tier validation: hindcast hold-out RL comparison, rate sanity vs HURDAT2 climatologies, wind field spot-checks. Entry points: `make_validation_cfg()`, `run_validation_suite()`, `validate_hazard_model()`.

7. **`R/hazard_viz.R`** — Plotting functions for tracks, wind profiles, return levels, and diagnostics.

8. **`R/hazard_utils.R`** — Internal utilities shared across modules.

## Key Scientific Conventions

- **Internal units**: wind in kt (converted to m/s at output boundaries only), distances in km, pressure in hPa. Unit conversions at interfaces: 1 nm = 1.852 km.
- **RMW precedence** (deterministic): observed USA_RMW → calibrated mapping from mean wind radii → Knaff fallback.
- **R34 quality gate**: directional R34 used when ≥3 quadrants present; mean radius when ≥2; climatological infill otherwise.
- **Holland B parameterization**: Vickery & Wadhera (2008) style; bounded [1.0, 2.5].
- **Storm classes**: TD (<34 kt), TS (34–63 kt), HUR (≥64 kt).

## Branch and Commit Conventions

- **Never commit to `master` directly.** All work on the `agent` branch.
- Commit message format: `<type>: <what>` (max 72 chars). Types: `fix`, `feat`, `refactor`, `docs`, `test`, `chore`.
- One logical change per commit.

## Completion Reports

Every task must produce exactly one report in `tools/codex_reports/` named:
`YYYY-MM-DD__HHMM__<task-slug>.md`

Required sections: goal, scope, summary, files changed, commands run, test results, behavior changes, follow-ups/risks.

## Allowed Dependencies

Core (Imports in DESCRIPTION): `dplyr`, `tidyr`, `ggplot2`, `patchwork`, `readr`, `lubridate`, `geosphere`, `tibble`, `purrr`, `magrittr`, `rlang`, `MASS`, `ncdf4`, `scales`, `stringr`.
Do not add new dependencies without explicit approval.

## Pipeline Scripts

Entry-point pipeline scripts live in `inst/extcode/pipelines/`. Sample data and demo scripts are in `inst/extcode/`.

## Additional Rules

The `.agents/skills/r-coding/SKILL.md` file defines R coding constraints (API stability, minimal diff, error handling, NAMESPACE hygiene) that always apply alongside this file.
