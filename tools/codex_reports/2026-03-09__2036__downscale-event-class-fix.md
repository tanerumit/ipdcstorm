Goal
- Fix the downscaled daily/event `event_class` mismatch so class labels match the final realized event peak used in generated daily outputs.

Scope
- Patched `R/hazard_downscale.R`.
- Added one targeted test in `tests/testthat/`.
- Did not modify core canonical event-table logic or visualization code.

Summary
- Added an internal helper to classify realized downscaled event peaks using canonical TS/HUR thresholds in knots, reusing `classify_severity()` when available and falling back to the local identical helper otherwise.
- Updated `generate_daily_year_extended()` to classify each event from the final realized in-year pulse peak `max(contrib)` before writing `event_class` into daily rows.
- Added an internal post-generation consistency check to stop if any daily event label disagrees with its realized daily peak.
- Added a focused test proving mismatched carried labels are corrected to `HUR`, `TS`, and `TD` based on realized pulse peaks.

Files Changed
- `R/hazard_downscale.R`
- `tests/testthat/test-downscale-classification.R`

Commands Run
- `Rscript -e "parse(file='R/hazard_downscale.R')"`
- `Rscript -e "library(dplyr); library(tidyr); library(ggplot2); library(magrittr); source('R/hazard_core.R'); source('R/hazard_downscale.R')"`
- `Rscript -e "devtools::test(filter = 'downscale|hazard|classification')"`
- `git diff -- R/hazard_downscale.R tests/testthat/test-downscale-classification.R`
- `git diff --no-index -- NUL tests/testthat/test-downscale-classification.R`

Test Results
- Parse: passed.
- Source/load: passed.
- `devtools::test(filter = 'downscale|hazard|classification')`: passed (11 tests total in matched files, including new downscale classification test).
- Unrelated warning during test startup: `plot_summary_panel` listed as exported but not present in namespace.

Behavior Changes
- Downscaled daily `event_class` is now derived from the final realized event peak in the generated daily pulse rather than trusting the sampled incoming label.
- This can now yield `TD` for sub-34 kt realized peaks in downscaled daily outputs, which is the threshold-consistent result.

Follow-ups/Risks
- `sample_events_for_year_extended()` still constructs a provisional `event_class`, but `generate_daily_year_extended()` now overwrites daily output labels from realized pulse peaks, which is the earliest safe point once pulse shape and year truncation are known.
- If any other code consumes sampled-event tables directly before daily generation, it will still see provisional labels; that was not changed because the final realized peak is not known until pulse generation.
