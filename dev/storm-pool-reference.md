# Storm pool reference — stress-test anchor hurricanes

Reference for the six pinned-focal SIDs used by the stress-test experiment
(`inst/extcode/05-stress-test-experiment.R`, `PINNED_SID_POOL`). All six are
Category-4 or Category-5 Atlantic hurricanes that passed within ~300 km of the
Saba / Statia mid-point and are present in `inst/extdata/ibtracs_1970.csv`.

Last updated: 2026-04-19.

---

## 1. Meteorological comparison (IBTrACS v04r01, NA basin)

Values below are lifetime extrema from IBTrACS (per-storm maxima/minima across
the full track); "—" indicates the field is not populated in the packaged
subset for that storm. RMW (radius of maximum wind) and R34 (34-kt wind radii)
are under-populated in pre-2000 records.

| SID | Name | Season | Peak wind (kt) | Peak wind (mph) | Saffir–Simpson | Min pres (hPa) | Peak RMW (km) | Max R34 (nm) | Mean fwd speed (kt) | Track obs |
|---|---|---|---|---|---|---|---|---|---|---|
| 2017242N16333 | **Irma**    | 2017 | 155 | 178 | Cat 5 | 914 | 9.3  | 360 | 12.0 | 123 |
| 2017260N12310 | **Maria**   | 2017 | 150 | 173 | Cat 5 | 908 | 9.3  | 220 | 12.8 | 131 |
| 1989254N13340 | **Hugo**    | 1989 | 140 | 161 | Cat 5 | 918 | 18.5 | —   | 17.0 | 124 |
| 1998259N10335 | **Georges** | 1998 | 135 | 155 | Cat 4 | 937 | 37.0 | —   | 11.6 | 134 |
| 1999318N17278 | **Lenny**   | 1999 | 135 | 155 | Cat 4 | 933 | —    | —   |  8.0 |  77 |
| 1995241N11333 | **Luis**    | 1995 | 130 | 150 | Cat 4 | 935 | —    | —   | 14.3 | 121 |

Saffir–Simpson category is inferred from peak sustained wind. Lifetime peaks
occur at sea and are generally NOT equal to winds felt on the Dutch Caribbean
islands — see Section 2.

---

## 2. Closest approach to Dutch Caribbean targets

Minimum great-circle distance (Haversine) between each storm's track and each
target, with the IBTrACS wind value at the closest-approach fix. "Wind @
closest" is the storm's over-water intensity at the closest point; it is not a
site-level (overland, downscaled) wind.

### Sint Maarten  (18.071 °N, 63.050 °W)

| Storm | Min dist (km) | Wind @ closest (kt) | Date of closest approach |
|---|---|---|---|
| Irma    |   6 | 155 | 2017-09-06 11:15 Z |
| Lenny   |   6 |  85 | 1999-11-18 (hovering) |
| Luis    |  37 | 115 | 1995-09-05 |
| Georges |  95 |  95 | 1998-09-21 12:00 Z |
| Hugo    | 139 | 125 | 1989-09-17 18:00 Z |
| Maria   | 168 | 148 | 2017-09-19 21:00 Z |

### Saba  (17.635 °N, 63.230 °W)

| Storm | Min dist (km) | Wind @ closest (kt) | Date of closest approach |
|---|---|---|---|
| Lenny   |  41 | 110 | 1999-11-18 18:00 Z |
| Georges |  47 |  95 | 1998-09-21 12:00 Z |
| Irma    |  52 | 155 | 2017-09-06 12:00 Z |
| Luis    |  84 | 115 | 1995-09-05 21:00 Z |
| Hugo    |  87 | 125 | 1989-09-17 18:00 Z |
| Maria   | 117 | 148 | 2017-09-19 21:00 Z |

### Sint Eustatius / Statia  (17.489 °N, 62.974 °W)

| Storm | Min dist (km) | Wind @ closest (kt) | Date of closest approach |
|---|---|---|---|
| Georges |  33 |  99 | 1998-09-21 09:00 Z |
| Lenny   |  45 |  65 | 1999-11-19 15:00 Z |
| Irma    |  54 | 155 | 2017-09-06 09:00 Z |
| Luis    |  83 | 115 | 1995-09-05 18:00 Z |
| Hugo    |  86 | 125 | 1989-09-17 18:00 Z |
| Maria   | 114 | 145 | 2017-09-19 18:00 Z |

**Pattern.** Irma is the only member of the pool whose eyewall went directly
over a Dutch target (Sint Maarten, within 6 km at Cat-5). Lenny also passes
within 6 km of Sint Maarten but as a Cat-2 by that point (its Cat-4 peak was
reached further west). Georges, Statia's closest approach, was at Cat-2/3 (95
kt). Maria, despite its Cat-5 peak, stayed >110 km south of every island.

---

## 3. Recorded impacts on the Dutch Caribbean

Short summaries focused on wind, surge, damage, and fatalities specifically
for Sint Maarten, Saba, and Sint Eustatius. Sources given per storm.

### Irma (2017-09-06)
- **Track geometry** — Category-5 eyewall passed directly over Sint Maarten
  early 6 September 2017; Saba and Statia were 50–55 km south, inside the
  southern eyewall.
- **Winds / surge** — Sustained winds estimated ~155 kt (~287 km/h) at
  landfall; gusts exceeded 200 km/h before instruments at Princess Juliana
  Airport failed. Storm surge 2–3 m on north-facing Sint Maarten coasts.
- **Damage** — ~70% of Sint Maarten structures damaged or destroyed; Princess
  Juliana Airport heavily damaged; 4 fatalities on Sint Maarten (Dutch side).
  Saba and Statia reported widespread roof damage and utility outages but
  no fatalities.
- **USD damage** — Sint Maarten (Dutch side): **~US$2.5 billion (2017)**
  per the Sint Maarten National Recovery and Resilience Plan; Netherlands
  Antilles whole-kingdom exposure via the SXM Trust Fund.

Sources:
- [NHC Tropical Cyclone Report AL112017 — Hurricane Irma](https://www.nhc.noaa.gov/data/tcr/AL112017_Irma.pdf)
- [Sint Maarten National Recovery and Resilience Plan (government.sx)](https://www.government.sx/Portals/0/Recovery%20Plan/NRRP_final.pdf)
- [Sint Maarten — Hurricane Irma Situation Report (ReliefWeb/UN OCHA)](https://reliefweb.int/report/sint-maarten-kingdom-netherlands/sint-maarten-hurricane-irma-situation-report-no-1)

### Maria (2017-09-19/20)
- **Track geometry** — Passed ~115 km south of Saba/Statia on 19 Sep 2017 as
  a Category-5 and moved northwest toward Puerto Rico; Sint Maarten was ~170
  km north of the track.
- **Winds / surge** — Dutch Caribbean islands experienced tropical-storm to
  low-end hurricane-force gusts (~50–70 kt / 90–130 km/h); no significant
  surge relative to Irma. Major Maria surge and wind damage occurred on
  Dominica and Puerto Rico, not the Dutch islands.
- **Damage** — Compounding effect on Irma-damaged structures (tarps stripped,
  additional roof damage, power-restoration setbacks). No Maria-attributed
  fatalities on the Dutch Caribbean.
- **USD damage** — No standalone Maria figure published for the Dutch
  Caribbean; incremental losses bundled into the Irma PDNA. Regional Maria
  losses were dominated by Puerto Rico (~US$90 billion, 2017) and Dominica
  (~US$1.3 billion, 2017) per NHC TCR.

Sources:
- [NHC Tropical Cyclone Report AL152017 — Hurricane Maria](https://www.nhc.noaa.gov/data/tcr/AL152017_Maria.pdf)
- [NWS San Juan — Major Hurricane Maria (September 20, 2017)](https://www.weather.gov/sju/maria2017)

### Hugo (1989-09-17/18)
- **Track geometry** — Eye tracked over Guadeloupe on 17 Sep 1989 then
  northwest across the northeastern Caribbean; closest approach to Sint
  Maarten was ~137 km to the southwest; Saba and Statia were ~86–87 km north
  of the eye.
- **Winds / surge** — Sint Maarten station recorded peak sustained winds
  ~46 mph (~40 kt) and peak gusts ~78 mph (~68 kt). Saba and Statia
  experienced hurricane-force gusts; Saba and Statia "lost much of their
  vegetation" per the NHC preliminary report.
- **Damage** — Sint Maarten: roof damage, downed power lines, ~25 sailboats
  damaged (one lost with 4 people). Saba and Statia: "many homes, piers, and
  public buildings" sustained severe damage.
- **USD damage / fatalities** — **US$50 million (1989)** and **11 fatalities**
  across the Netherlands Antilles per the NHC preliminary report.

Sources:
- [NHC Preliminary Report — Hurricane Hugo (via NWS ILM archive)](https://www.weather.gov/media/ilm/climate/Hugo/NHC_report_Hugo.pdf)
- [NHC — TPC NHC Hurricane Hugo (1989) archive](https://www.nhc.noaa.gov/1989hugo.html)

### Luis (1995-09-05/06)
- **Track geometry** — Eye passed ~37 km north of Sint Maarten on 5 Sep 1995
  as a strong Category-4; Sint Maarten was in the southern eyewall; Saba and
  Statia were 83–84 km south-southwest and experienced outer eyewall /
  strong rainbands.
- **Winds / surge** — Luis produced Sint Maarten's most catastrophic wind and
  surge event prior to Irma. Sustained winds approaching 120 kt, with reports
  of gusts >150 kt at Princess Juliana. Offshore significant wave heights
  exceeded 10 m; surge 2–4 m on north-facing coasts.
- **Damage** — Philipsburg ~70% destroyed; ~15% of residences rendered
  uninhabitable; 1,300 of ~1,500 boats sheltered in Simpson Bay Lagoon sunk
  or grounded. Airport terminal, hotels, churches, and schools damaged.
- **USD damage / fatalities** — **~US$1.8 billion (1995)** on the Dutch side
  of Sint Maarten; **8 fatalities** on Sint Maarten per the NHC preliminary
  report.

Sources:
- [NHC Preliminary Report — Hurricane Luis (AL131995)](https://www.nhc.noaa.gov/data/tcr/AL131995_Luis.pdf)
- [Caribbean — Hurricane Luis UN DHA Situation Reports 1–10 (ReliefWeb)](https://reliefweb.int/report/antigua-and-barbuda/caribbean-hurricane-luis-sep-1995-un-dha-situation-reports-1-10)
- [ECLAC — Effects of Hurricane Luis on St. Kitts & Nevis / Netherlands Antilles](https://www.cepal.org/en/publications/25592-st-kitts-nevis-anguilla-sint-maarten-effects-hurricane-luis)

### Lenny (1999-11-18/19)
- **Track geometry** — Unusual west-to-east ("wrong-way") track. The eye
  stalled ~6 km south of Sint Maarten on 18 Nov 1999; Saba was ~41 km south
  of the eye, Statia ~45 km. The hurricane lingered near the islands for
  roughly 24 hours before moving east-northeast.
- **Winds / surge** — Sint Maarten reported winds from the west-southwest,
  the unusual direction producing surge and swell on normally-leeward
  west-facing coasts. Winds in the 80–90 kt range with gusts ~110 kt.
  Philipsburg rainfall totaled 27.56 inches (~700 mm).
- **Damage** — Severe west-coast damage on Sint Maarten (Simpson Bay, marinas
  and hotels typically sheltered from easterly trade-wind weather). Saba's
  Fort Bay harbor heavily damaged; Statia's oil terminal and coastal roads
  damaged. Containers floated ashore on Sint Maarten.
- **USD damage / fatalities** — Region-wide ~US$330 million (1999) per NHC;
  Netherlands Antilles share ~US$60–70 million (1999) per ECLAC. Fatalities
  in the Dutch Caribbean were limited (no deaths on Sint Maarten are
  attributed to Lenny in the NHC report; regional total ~13 across the
  northeastern Caribbean).

Sources:
- [NHC Preliminary Report — Hurricane Lenny (AL161999)](https://www.nhc.noaa.gov/data/tcr/AL161999_Lenny.pdf)
- [ECLAC — Effects of Hurricane Lenny on Netherlands Antilles and neighbours](https://www.cepal.org/en/publications/25548-saint-kitts-nevis-antigua-barbuda-netherlands-antilles-effects-hurricane-lenny)
- [Caribbean — Hurricane Lenny OCHA Situation Report (ReliefWeb)](https://reliefweb.int/report/antigua-and-barbuda/caribbean-hurricane-lenny-ocha-situation-report-no-5)

### Georges (1998-09-20/21)
- **Track geometry** — Passed over Antigua and St. Kitts on 20–21 Sep 1998
  as Category-3/4, then ~33 km south of Statia and ~47 km south of Saba;
  Sint Maarten was ~95 km north of the track.
- **Winds / surge** — On Saba and Statia: tropical-storm-force sustained
  winds with hurricane-force gusts (~60–75 kt gusts); surge <1.5 m. Sint
  Maarten experienced tropical-storm conditions and minor coastal flooding.
- **Damage** — Limited on the Dutch Caribbean relative to other pool storms:
  some roof damage, crop losses, and power disruption on Saba/Statia; only
  tropical-storm conditions on Sint Maarten. Main regional impacts were on
  St. Kitts, Puerto Rico, the Dominican Republic, and Haiti.
- **USD damage / fatalities** — No separately itemised Dutch Caribbean figure
  in the NHC TCR or ECLAC assessments; losses were minor relative to the
  region-wide ~US$9.7 billion (1998). No Dutch Caribbean fatalities reported.

Sources:
- [NHC Preliminary Report — Hurricane Georges (AL071998)](https://www.nhc.noaa.gov/data/tcr/AL071998_Georges.pdf)
- [Atlantic Hurricane Season of 1998 — NHC Monthly Weather Review](https://www.nhc.noaa.gov/data/mwreview/1998.pdf)

---

## 4. Consolidated impact summary (Dutch Caribbean only)

| Storm | Primary island affected | Peak on-island wind | Surge | Fatalities (Dutch Carib.) | Damage (Dutch Carib., year USD) |
|---|---|---|---|---|---|
| **Irma** 2017 | Sint Maarten (direct Cat-5 hit) | ~155 kt sust. | 2–3 m | 4 | ~US$2.5 B (2017) |
| **Luis** 1995 | Sint Maarten (southern eyewall) | ~120 kt sust. | 2–4 m | 8 | ~US$1.8 B (1995) |
| **Hugo** 1989 | Saba / Statia / Sint Maarten | ~40 kt sust, 68 kt gust (SXM) | <2 m | 11 (N. Antilles) | ~US$50 M (1989, N. Antilles) |
| **Lenny** 1999 | Sint Maarten (west-coast surge) | 80–90 kt sust. | 2–4 m W-coast | 0 (SXM) | ~US$60–70 M (1999, N. Antilles) |
| **Georges** 1998 | Statia / Saba (passing south) | TS sust., 60–75 kt gusts | <1.5 m | 0 | Not itemised (minor) |
| **Maria** 2017 | Compounding Irma (south of islands) | 50–70 kt gusts | minor | 0 | Bundled into Irma PDNA |

**Stress-test relevance.** The pool captures a deliberate gradient:
- **Extreme:** Irma, Luis — direct or near-direct Cat-4/5 hits on Sint Maarten
  with site-level sustained winds ≥120 kt and >US$1 B damage each.
- **Severe:** Lenny — Cat-4 peak, Cat-2 at islands, but unusual westerly
  approach producing leeward-coast surge not captured by easterly-track
  hurricanes; Hugo — Cat-4/5 peak at sea, Cat-1-equivalent gusts at SXM with
  major vegetation/marine damage on Saba/Statia.
- **Moderate:** Georges — passing-south Cat-4 that clipped Saba/Statia with
  TS-force sustained winds.
- **Compounding:** Maria — Cat-5 south of the islands that added damage to
  Irma-compromised infrastructure without itself generating major
  Dutch-Caribbean on-island winds.

The stress-test pool therefore spans the full range of geometries the Dutch
Caribbean actually sees: direct eyewall hits, outer-eyewall clips,
westerly-surge events, and compounding follow-ons.
