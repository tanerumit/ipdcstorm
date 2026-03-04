# Diagnostics A-D for Saba and Statia

Run date: 2026-03-04

This summary is based on one validation pass executed from the repo root with the packaged IBTrACS file at `inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv` (the default `data/ibtracs/...` path is not present in this workspace). The printed A-D diagnostics now come directly from `validate_hazard_model()`.

## Repro command

```powershell
@'
pkgload::load_all(".")
targets <- tibble::tribble(
  ~name, ~lat, ~lon,
  "St_Martin", 18.0708, -63.0501,
  "Saba", 17.6350, -63.2300,
  "Statia", 17.4890, -62.9740,
  "Puerto_Rico", 18.2208, -66.5901,
  "Miami", 25.7617, -80.1918
)
cfg <- make_hazard_cfg(data_path = "inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv")
val_cfg <- make_validation_cfg(save_plots = FALSE, save_tables = FALSE)
validate_hazard_model(
  cfg = cfg,
  targets = targets,
  validation_cfg = val_cfg
)
'@ | Rscript -
```

## Findings

### Saba

- A) R34 calibration inflation: `delta_top1_p50 = 0.00 kt`, `delta_top1_p90 = 0.00 kt`, `delta_overall_p99 = +0.0819 kt`.
  Interpretation: disabling calibration does not reduce the upper tail. Hypothesis A is not supported in this run.
- B) Outer cutoff hit rate: `frac_beyond_all = 0.8258`, `frac_beyond_top1 = 0.0000`.
  Interpretation: top-1% points are not beyond the outer cutoff, so the cutoff cannot explain tail bias.
- C) Radii / geometry: `R34_p50_all = 166.68 km`, `R34_p90_all = 277.80 km`, `R34_p50_top1 = 166.68 km`, `R34_p90_top1 = 166.68 km`, `dist_p50_top1 = 88.33 km`, `dist/R34_p50_top1 = 0.321`.
  Interpretation: top-1% points sit well inside R34, which is consistent with a permissive geometry / envelope explanation.
- D) Forward-motion asymmetry: `boost_p50 = 0.036 kt`, `boost_p90 = 2.415 kt`, `boost_p99 = 4.788 kt`, `frac_boost_gt10 = 0.000106`.
  Interpretation: asymmetry boosts are modest and do not support asymmetry as the main tail inflator.

### Statia

- A) R34 calibration inflation: `delta_top1_p50 = 0.00 kt`, `delta_top1_p90 = 0.00 kt`, `delta_overall_p99 = +0.0096 kt`.
  Interpretation: disabling calibration does not reduce the upper tail. Hypothesis A is not supported in this run.
- B) Outer cutoff hit rate: `frac_beyond_all = 0.8232`, `frac_beyond_top1 = 0.0000`.
  Interpretation: top-1% points are not beyond the outer cutoff, so the cutoff cannot explain tail bias.
- C) Radii / geometry: `R34_p50_all = 166.68 km`, `R34_p90_all = 277.80 km`, `R34_p50_top1 = NA`, `R34_p90_top1 = NA`, `dist_p50_top1 = 88.40 km`, `dist/R34_p50_top1 = NA`.
  Interpretation: the top-1% subset has no finite `dist/R34` values, so the direct ratio test is inconclusive for Statia in this run.
- D) Forward-motion asymmetry: `boost_p50 = 0.044 kt`, `boost_p90 = 2.435 kt`, `boost_p99 = 4.802 kt`, `frac_boost_gt10 = 0.000000`.
  Interpretation: asymmetry boosts are modest and do not support asymmetry as the main tail inflator.

## Overall read

- Strongest support: C for Saba, because the top tail sits deep inside the implied R34 envelope.
- Not supported: A, B, and D for both sites in this run.
- Statia needs a follow-up geometry check on the rows with missing `R34_km` in the top-1% subset, because the direct `dist/R34` statistic is unavailable there.
