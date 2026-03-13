## 1. Goal
Diagnose the pre-calibration tail mechanism for observed-`R34` storms at Saba and Statia, add internal diagnostics focused on observed-`RMW` behavior in the `~1.75-2.5 x RMW` regime, test one narrow pre-calibration candidate, and keep only changes that pass the rollback gate.

## 2. Scope
Included:
- Internal forensic diagnostics for observed-`R34` tail contributors
- Explicit `RMW` provenance carried through trackpoint wind diagnostics
- Intermediate pre-calibration radial-response terms in hindcast diagnostics
- Focused pre-calibration radius-band summaries for observed-`R34` plus observed-`RMW` events
- Fixed-seed workspace reruns comparing legacy, calibration-adjusted, and pre-calibration candidate paths
- Targeted deterministic tests for the new forensic outputs

Excluded:
- Any package default change
- Any public API change
- Any rate-mode, search-radius, sampler, or fallback-rule change
- Any retained operational candidate after it failed acceptance

## 3. Problem solved
The prior forensic pass established that the dominant Saba/Statia tail bias was entering before `R34` calibration, but the package did not yet expose enough internal structure to isolate the pre-calibration response by radius band and `RMW` provenance. This patch adds that forensic visibility and records the result of a narrow candidate trial that was subsequently rolled back because it failed the hindcast acceptance criteria.

## 4. Summary
The new diagnostics confirm the dominant excess is in observed-`R34` plus observed-`RMW` events before calibration. For the key 2017 IRMA-like contributors:
- Saba: `2017242N16333`, `63.621 kt`, `r/RMW = 1.882`, `peak_rmw_provenance = observed`, pre-cal `62.480 kt`, post-cal `62.401 kt`
- Statia: `2017242N16333`, `62.277 kt`, `r/RMW = 1.951`, `peak_rmw_provenance = observed`, pre-cal `61.934 kt`, post-cal `61.851 kt`

The added forensic fields show:
- `peak_gradient_factor`
- `peak_surface_factor`
- `peak_steepening_factor`
- `peak_precal_response_factor`
- pre/post-calibration winds
- forward-motion increment

The focused observed-`RMW` summaries now split the tail into:
- `<1.75`
- `[1.75,2.5)`
- `>=2.5`

and provide deterministic summaries and comparisons for the observed-`R34` plus observed-`RMW` subset.

I tested one narrow internal pre-calibration candidate during workspace validation, restricted to observed-`R34` plus observed-`RMW` points in `[1.75,2.5)`. It reduced the key top-event pre-calibration winds:
- Saba 2017 top event: `63.621 -> 59.018 kt`
- Statia 2017 top event: `62.277 -> 56.769 kt`

But it failed the hindcast gate overall and was rolled back from the final diff:
- Saba RL bias RP `5/10/25/50`: `+10.83/+9.94/+10.47/+11.12%` to `+10.93/+10.08/+10.55/+11.11%`
- Statia RL bias RP `5/10/25/50`: `+12.47/+16.10/+21.74/+26.20%` to `+13.05/+17.47/+23.84/+28.76%`
- Compared with the rejected calibration-adjusted candidate, the pre-calibration candidate was worse at Statia and not materially better at Saba.

Interpretation:
- The pre-calibration candidate did reduce excess in the target `[1.75,2.5)` regime.
- That reduction did not translate into better Saba/Statia hindcast RL bias.
- The failure mode matches the rollback rule: local cosmetic suppression without a net hindcast improvement.
- The next step should move deeper into the underlying radial wind formulation rather than add another localized modifier.

## 5. Files changed
- `R/hazard_core.R` — 65 added lines, `RMW` provenance plumbing for diagnostics
- `R/hazard_validation.R` — 1350 added / 11 removed lines, extended forensic diagnostics, focused pre-calibration radius-band summaries, metadata plumbing
- `tests/testthat/test-rmw_estimation.R` — 60 added lines, deterministic provenance coverage
- `tests/testthat/test-hazard-validation-exported.R` — 164 added / 1 removed lines, deterministic forensic summary coverage
- `man/dot-estimate_site_wind_holland.Rd` — 6 added lines, internal method note update

## 6. Commands run
```powershell
Rscript -e "parse(file='R/hazard_core.R'); parse(file='R/hazard_validation.R'); parse(file='tests/testthat/test-rmw_estimation.R'); parse(file='tests/testthat/test-hazard-validation-exported.R')"
Rscript -e "devtools::test(filter = 'hazard[-_ ](validation|rmw)|hindcast|wind|R34|tail')"
Rscript -e "tools::checkRd('man/dot-estimate_site_wind_holland.Rd'); cat('Rd OK\n')"
```

Workspace forensic and comparison reruns were executed from `pkgload::load_all('.')` with:
- `model_seed = 2026`
- `validation_seed = 42`
- raw hindcast lambda
- legacy sampler
- baseline `ipdcstorm.wind_field_mode = 'legacy'`
- comparison `ipdcstorm.wind_field_mode = 'observed_r34_calibration_adjusted'`
- trial-only `ipdcstorm.wind_field_mode = 'observed_rmw_precal_adjusted'` before rollback

## 7. Test results
- Parse checks: passed
- Scoped tests: passed
- Result: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 155 ]`
- Touched internal Rd check: passed for `man/dot-estimate_site_wind_holland.Rd`

## 8. Behavior changes
No default behavior changed.

The final diff adds diagnostics only:
- `RMW_provenance` is now carried explicitly through trackpoint wind outputs
- observed-`R34` forensic tables now include intermediate radial-response terms
- new internal summaries isolate observed-`R34` plus observed-`RMW` tail behavior by pre-calibration radius band

The failed pre-calibration candidate was not kept in the final code.

## 9. Follow-ups/risks
- The dominant biased Saba/Statia events are still observed-`R34` plus observed-`RMW` cases near `~1.9 x RMW`.
- The excess is clearly entering before calibration; the calibration step is near-neutral or slightly damping for the key events.
- A narrow modifier can reduce the target-band event winds, but the hindcast validation objective does not improve, so more modifiers of that kind are unlikely to be the right next move.
- Non-target sites did not degrade under the trial candidate in a catastrophic way, but the candidate still failed because Statia worsened materially and Saba did not improve enough.
- Next step should move deeper into the radial wind formulation itself, not another localized post-hoc multiplier.
