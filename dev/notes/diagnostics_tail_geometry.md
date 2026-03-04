# Tail Geometry Diagnostics (TS+ q99)

Run date: 2026-03-04

This run used the packaged IBTrACS file at `inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv`, because the default `data/ibtracs/...` path is not present in this workspace.

## Run command used

```powershell
@'
pkgload::load_all(".")
cfg <- make_hazard_cfg(data_path = "inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv")
targets <- tibble::tribble(
  ~name, ~lat, ~lon,
  "St_Martin", 18.0708, -63.0501,
  "Saba", 17.6350, -63.2300,
  "Statia", 17.4890, -62.9740,
  "Puerto_Rico", 18.2208, -66.5901,
  "Miami", 25.7617, -80.1918
)
out <- run_hazard_model(cfg = cfg, targets = targets, verbose = FALSE)
.run_validation_diagnostics_tail_geometry(out$trackpoints, targets)
'@ | Rscript -
```

## Per-site table

| site | n_pos | top_n | R34_p50 (km) | R34_p90 (km) | R34_unique | dist_p50_top (km) | dist_over_R34_p50_top | dist_over_R34mean_p50_top |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Saba | 27 | 1 | 203.72 | 277.80 | 12 | 53.56 | 0.321 | 0.231 |
| Statia | 35 | 1 | 203.72 | 305.58 | 16 | 54.20 | 0.325 | 0.234 |
| St_Martin | 31 | 1 | 240.76 | 277.80 | 13 | 26.64 | 0.131 | 0.113 |

## Rounded R34 frequency check (TS+ subset)

- Saba top rounded `R34_km` values: `277.80` (6), `111.12` (4), `166.68` (3), `203.72` (3), `129.64` (2)
- Statia top rounded `R34_km` values: `277.80` (7), `111.12` (5), `203.72` (5), `166.68` (4), `74.08` (2)
- St_Martin top rounded `R34_km` values: `277.80` (8), `111.12` (5), `185.20` (3), `240.76` (3), `129.64` (2)

## Interpretation

- There is clear radius quantization in all three TS+ subsets: `R34_unique` is only `12-16` distinct rounded values despite `27-35` rows. That points to repeated discrete radius values rather than a smooth continuum.
- The specific `90 nm` value (`166.68 km`) does recur, but it is not the median in this TS+ restricted diagnostic. Both Saba and Statia have `R34_p50 = 203.72 km` (`110 nm`), so the earlier `166.68 km` median from the all-winds diagnostic is not reproduced in the TS+ subset.
- Saba and Statia show a material directional-vs-mean geometry gap in the tail-driving point:
  - Saba: `dist/R34 = 0.321` vs `dist/R34_mean = 0.231`
  - Statia: `dist/R34 = 0.325` vs `dist/R34_mean = 0.234`
  This is an absolute shift of about `0.09` and a relative drop of about `28%`, which means the storm-wide mean radius is materially larger than the directional radius for the tail point. That is large enough to matter if any upstream logic is mixing directional and mean radii.
- St_Martin also shows a gap, but smaller: `0.131` vs `0.113` (about `14%` relative), so the directional-vs-mean mismatch is stronger for Saba and Statia than for the control site in this run.

## Most likely upstream causes to inspect

- `.get_directional_radius()` in `R/hazard_core.R`: confirms how quadrant-specific radii are selected and whether repeated exact quadrant values are expected from the source fields.
- The raw quadrant radius ingestion / conversion path that populates `r34_*_nm`: repeated multiples of `1.852` strongly suggest discrete nautical-mile source radii are being preserved without smoothing, which can create quantized `R34_km`.
- `mean_radius_nm()` inside `compute_site_winds_full()` in `R/hazard_core.R`: this is where the storm-wide mean radius is formed; any directional-vs-mean mismatch for tail points will show up here.
- Bearing/quadrant assignment (`calculate_bearing()` and `.get_quadrant()` in `R/hazard_core.R`): if the chosen quadrant flips or is not the physically relevant one for a near-miss point, the directional radius used at the site can diverge from the mean radius in a way that changes `dist/R34`.
