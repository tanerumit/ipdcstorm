## 1. Goal
Replace `damage_method` + `damage_params` in the daily hazard/impact downscaling API with a single `damage` argument that carries both method selection and method-specific parameters, without changing the underlying damage formulations.

## 2. Scope
Included:
- `generate_daily_hazard_impact()` and its direct internal helper in [R/hazard_downscale.R](C:/Users/taner/WS/ipdcstorm/R/hazard_downscale.R)
- Targeted daily hazard/downscale tests in [tests/testthat/test-hazard-downscale-exported.R](C:/Users/taner/WS/ipdcstorm/tests/testthat/test-hazard-downscale-exported.R)
- Regenerated Rd files for the touched roxygen blocks

Excluded:
- `pulse_shape`
- Broader hazard/downscaling API redesign
- Unrelated exported functions and unrelated dirty worktree changes

## 3. Problem solved
The old split between `damage_method` and `damage_params` made the public API more ambiguous and allowed conceptually incompatible fields to be mixed. The new `damage` list keeps the method and its parameters together and rejects unsupported fields up front.

## 4. Summary
Added a strict internal `.validate_damage_spec()` helper that requires a named list with `damage$method`, fills method defaults, rejects unknown fields, and validates provided numeric scalars. `generate_daily_hazard_impact()` now accepts `damage = list(method = ...)` and dispatches both the intensity and power-law paths from that unified object. `damage_intensity` is still produced for both methods; for `powerlaw` it remains a generic bounded wind-intensity index using `V0 = 34`, `V1 = 120`, and the shared exponent `p`.

## 5. Files changed
- [R/hazard_downscale.R](C:/Users/taner/WS/ipdcstorm/R/hazard_downscale.R): +504/-142 lines in diff; API refactor, strict damage validation, roxygen updates, internal dispatch updates
- [tests/testthat/test-hazard-downscale-exported.R](C:/Users/taner/WS/ipdcstorm/tests/testthat/test-hazard-downscale-exported.R): +130/-0 lines in diff; new default/override/validation coverage for `damage`
- [man/generate_daily_hazard_impact.Rd](C:/Users/taner/WS/ipdcstorm/man/generate_daily_hazard_impact.Rd): +229/-27 lines in diff; regenerated public docs for `damage`
- [man/dot-generate_daily_hazard_impact_single.Rd](C:/Users/taner/WS/ipdcstorm/man/dot-generate_daily_hazard_impact_single.Rd): +9/-11 lines in diff; regenerated internal worker docs
- [man/dot-validate_damage_spec.Rd](C:/Users/taner/WS/ipdcstorm/man/dot-validate_damage_spec.Rd): new 18-line generated Rd file for the internal validator

## 6. Commands run
```powershell
git branch --show-current
Get-Content -Path C:\Users\taner\WS\shared-tools\01_skills\r-coding\SKILL.md -TotalCount 220
rg -n "generate_daily_hazard_impact|damage_method|damage_params|damage_intensity|powerlaw|pulse_shape" R tests/testthat README* vignettes man
$i=1; Get-Content -Path R\hazard_downscale.R | ForEach-Object { '{0,4}: {1}' -f $i, $_; $i++ } | Select-Object -Skip 930 -First 430
$i=1; Get-Content -Path tests\testthat\test-hazard-downscale-exported.R | ForEach-Object { '{0,4}: {1}' -f $i, $_; $i++ }
git status --short
@'
source('R/hazard_downscale.R')
... baseline fixture and old API calls ...
saveRDS(list(intensity = res_int, powerlaw = res_pow), 'C:/Users/taner/AppData/Local/Temp/damage_api_baseline.rds')
'@ | Rscript -
Rscript -e "parse(file='R/hazard_downscale.R')"
Rscript -e "parse(file='tests/testthat/test-hazard-downscale-exported.R')"
Rscript -e "devtools::test(filter = 'daily|damage|downscale')"
Rscript -e "devtools::document()"
Rscript -e "parse(file='R/hazard_downscale.R'); parse(file='tests/testthat/test-hazard-downscale-exported.R')"
Rscript -e "devtools::test(filter = 'daily|damage|downscale')"
@'
source('R/hazard_downscale.R')
... baseline fixture and new API calls ...
print(report)
'@ | Rscript -
git diff --numstat -- R/hazard_downscale.R tests/testthat/test-hazard-downscale-exported.R man/generate_daily_hazard_impact.Rd man/dot-generate_daily_hazard_impact_single.Rd man/dot-validate_damage_spec.Rd
git status --short -- R/hazard_downscale.R tests/testthat/test-hazard-downscale-exported.R man/generate_daily_hazard_impact.Rd man/dot-generate_daily_hazard_impact_single.Rd man/dot-validate_damage_spec.Rd
```

## 7. Test results
- `Rscript -e "parse(file='R/hazard_downscale.R')"`: pass
- `Rscript -e "parse(file='tests/testthat/test-hazard-downscale-exported.R')"`: pass
- `Rscript -e "devtools::test(filter = 'daily|damage|downscale')"`: pass, 38 tests
- `Rscript -e "devtools::document()"`: pass
- Fixed-seed old-vs-new comparison on `wind_kt`, `wind_gust_kt`, `damage_intensity`, `damage_rate`, `cum_damage`: exact equality for both `intensity` and `powerlaw`

## 8. Behavior changes
- `generate_daily_hazard_impact()` now accepts `damage` and no longer accepts `damage_method` or `damage_params`
- Invalid `damage` inputs now fail fast when the list is unnamed, the method is unsupported, fields are incompatible with the selected method, or provided parameter values are malformed
- `damage_intensity` remains in the output for both methods; for `powerlaw` it is explicitly documented as a generic bounded wind-intensity index rather than the power-law damage formula itself

## 9. Follow-ups/risks
- No compatibility shim was added for the removed `damage_method` / `damage_params` arguments; external callers must update immediately
- The repo had unrelated pre-existing dirty changes, especially under `man/` and `tools/codex_reports/`; those were left untouched except for the files listed above
- Results delta: none. The fixed-seed baseline comparison was exactly equal on the requested output fields for representative intensity and power-law calls
