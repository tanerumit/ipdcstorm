# Spatial Coherence in Multi-Location Sampling — Design Options

**Context:** The current model samples storm events independently for each location.
As a result, the same historical storm (e.g. Irma, `2017242N16333`) can appear in
Saba's simulated year but not St. Martin's — even though both were hit simultaneously
in reality. This breaks spatial coherence and undermines portfolio-level stress-year
selection and co-occurrence analysis.

---

## Option 1 — Shared event pool (storm-first sampling) ✓ *selected*

Sample storms once at the basin level per simulated year, then assign each
sampled storm to all locations whose event library contains it. Co-occurrence
is enforced by construction: if Irma is drawn it appears simultaneously at
every location with Irma in its library, each with its own pre-computed
site-level wind profile.

**Implementation change:** restructure the per-year sampling loop from
*location → draw* to *draw → assign to locations*. Annual storm counts move
to basin level; per-location rates are used as weights when resolving how many
events to draw from each severity class.

**Pros:** exact spatial coherence, no new data needed, wind profiles unchanged.  
**Cons:** Poisson rates must be recalibrated at basin level; a single draw can
produce very different per-location counts in sparse years.

---

## Option 2 — Conditional resampling from a reference location

Keep per-location libraries but add a conditioning step: sample the anchor
location first, then sample each remaining location conditional on the anchor's
draw, using historical pairwise co-occurrence probabilities estimated from
`out$events`.

**Pros:** minimal architectural change, graftable onto current model.  
**Cons:** asymmetric (anchor choice matters), pairwise conditioning does not
guarantee trivariate or higher-order coherence, co-occurrence priors are noisy
for rare events.

---

## Option 3 — Copula coupling of annual metrics

Leave the sampling engine unchanged. Post-process the independent per-location
annual metric series by replacing each location's annual ranks with a joint
draw from a Gaussian or t-copula calibrated to historical inter-location
correlations (available from `04-spatial-correlation.R`).

**Pros:** completely decoupled from sampling machinery, straightforward to add
as a post-processing step.  
**Cons:** reshuffles years without guaranteeing shared `event_id` values across
locations — track-based and aftermath queries are unaffected. More a statistical
fix than a physical one; best suited to annual-metric outputs only.

---

## Option 4 — Full synthetic track model

Simulate complete storm tracks through the basin (beta-advection or parametric
model), apply a wind field model (Holland profile) to compute simultaneous wind
at all locations for each track.

**Pros:** gold standard, physically grounded, naturally handles all pairwise and
higher-order spatial relationships.  
**Cons:** substantial additional complexity (track model, wind field model,
landfall interaction); appropriate for a future v2 architecture rather than an
incremental fix.

---

*Written: 2026-04-12*
