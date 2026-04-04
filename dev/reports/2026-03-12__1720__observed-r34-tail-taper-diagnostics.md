## 1. Goal
Diagnose the normalized-radius regime associated with observed-`R34` tail contributors at Saba and Statia, implement one narrow observed-only decay/taper candidate, and validate whether that change reduces upper-tail hindcast bias without degrading other sites.

## 2. Scope
Included:
- Internal normalized-radius diagnostics for hindcast tail contributors
- Observed-`R34` tail detail and radius-bin summaries with deterministic metadata
- One internal observed-only wind-field candidate mode in `.estimate_site_wind_holland()`
- Targeted tests and a fixed-seed workspace validation rerun

Excluded:
- Any package default change
- Any rate-mode, search-radius, sampler, or fallback-rule change
- Any partial/climo/none pathway behavior change
- Any exported API change

## 3. Problem solved
The prior `R34_source` split showed the residual tail issue at Saba and Statia is mainly in the `observed` branch. This patch adds the missing normalized-radius diagnostics needed to localize those contributors and tests an observed-only taper without changing operational defaults.

## 4. Summary
Added internal diagnostics in `R/hazard_validation.R` to:
- carry `peak_normalized_radius = closest_approach_km / peak_rmw_used_km`
- extract observed-`R34` tail contributors
- bin retained observed-`R34` events by normalized radius
- attach deterministic metadata to all new tables

Added one internal comparison mode, `ipdcstorm.wind_field_mode = "observed_r34_adjusted"`, that applies a localized observed-only post-calibration taper for intense storms in roughly `1.5-2.75 x RMW`.

Validation outcome:
- Diagnosed strongest observed-path tail regime at Saba and Statia: the main high-end contributors sit near `1.88-1.95 x RMW`, within the `[1.75, 2.5)` bin.
- The candidate taper reduces those specific event winds, but it does **not** improve hindcast RL bias at Saba or Statia.
- `xi` improves slightly at Saba and Statia, but the RL bias gate fails because the historical hindcast tail is reduced at least as much as the simulated tail.
- Because the mode remains internal-only and defaults stay on legacy, no operational behavior changed.

## 5. Files changed
- `R/hazard_core.R` — added internal observed-only comparison mode for moderate-radius tapering
- `R/hazard_validation.R` — added normalized-radius provenance, observed-`R34` tail detail, radius-bin summaries, and metadata plumbing
- `tests/testthat/test-rmw_estimation.R` — added isolation tests for the observed-only taper
- `tests/testthat/test-hazard-validation-exported.R` — added deterministic normalized-radius diagnostic tests
- `man/dot-estimate_site_wind_holland.Rd` — updated internal method note for the new comparison mode

## 6. Commands run
```powershell
Rscript -e "parse(file='R/hazard_core.R'); parse(file='R/hazard_validation.R'); parse(file='tests/testthat/test-rmw_estimation.R'); parse(file='tests/testthat/test-hazard-validation-exported.R')"
Rscript -e "devtools::test(filter = 'hazard[-_ ](validation|rmw)|hindcast|wind|R34|tail')"
```

Workspace observed-`R34` diagnostics and candidate validation were run with fixed seeds via `pkgload::load_all('.')`, using:
- baseline: `ipdcstorm.wind_field_mode = "legacy"`
- candidate: `ipdcstorm.wind_field_mode = "observed_r34_adjusted"`
- shared settings: raw lambda, legacy sampler, `model_seed = 2026`, `validation_seed = 42`

## 7. Test results
- Parse checks: passed
- Scoped tests: passed
- Result: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 138 ]`

## 8. Behavior changes
No default behavior changed.

New internal diagnostic tables:
- `observed_r34_tail_detail`
- `observed_r34_radius_summary`
- `observed_r34_radius_comparison`

New internal comparison mode:
- `ipdcstorm.wind_field_mode = "observed_r34_adjusted"`

## 9. Follow-ups/risks
- The normalized-radius diagnosis is clear: Saba and Statia observed-path high-end contributors cluster near `~1.9 x RMW`.
- The candidate taper fails the main acceptance question. RL bias changed as follows:
  - Saba RP `5/10/25/50`: `+10.83/+9.94/+10.47/+11.12%` to `+11.14/+10.39/+10.93/+11.52%`
  - Statia RP `5/10/25/50`: `+12.47/+16.10/+21.74/+26.20%` to `+13.45/+18.29/+25.10/+30.31%`
- `xi` moved slightly toward observed:
  - Saba: `-0.245` to `-0.252` vs observed `-0.300`
  - Statia: `-0.068` to `-0.079` vs observed `-0.248`
- Non-target sites did not degrade materially:
  - Miami improved modestly
  - Puerto Rico was effectively unchanged
  - St_Martin was unchanged
- Next step should **not** continue on this exact observed-`R34` decay/taper mechanism as-is. The internal diagnostics suggest the residual problem is localized in observed-path moderate-radius contributors, but this direct taper reduces the historical hindcast tail too symmetrically to improve RL bias. A different mechanism is needed.
