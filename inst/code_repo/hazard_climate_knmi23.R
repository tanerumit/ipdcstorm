################################################################################
# hazard_climate_knmi23.R
# KNMI'23 Climate Scenario Support for Dutch Caribbean Hurricane Hazard Model
#
# Adds four KNMI'23 scenario modes (Ld, Ln, Hd, Hn) to generate_sst_scenario()
# and provides helper/metadata functions.
#
# The d/n (dry/wet) distinction does NOT change SST trajectories — it modifies
# Level 3 precipitation scaling only. Both Hd and Hn share the same SST path;
# likewise for Ld and Ln.
#
# References:
#   - Van der Wiel et al. (2024), Earth's Future, KNMI'23 methodology
#   - KNMI'23 BES climate scenarios (Van Dorland et al., 2023)
#   - Witteveen+Bos (2024), Climate change and adaptation efforts BES islands
#   - IPCC AR6 WG1 Ch. 4 & 9 (constrained GSAT and tropical SST scaling)
################################################################################


# =============================================================================
# 1) KNMI'23 SCENARIO METADATA
# =============================================================================

#' KNMI'23 scenario reference table
#'
#' @description
#' Returns a tibble describing the four KNMI'23 climate scenarios for the
#' Dutch Caribbean, including their SSP equivalents, MDR SST anomaly targets,
#' and recommended Level 3 precipitation scaling modifiers.
#'
#' The temperature axis (H/L) determines the emissions pathway and SST warming.
#' The moisture axis (d/n) affects precipitation response only (Level 3).
#'
#' MDR SST anomalies are relative to the 1991-2020 baseline, derived from:
#'   - KNMI'23 global warming levels (IPCC AR6 constrained estimates)
#'   - Tropical Atlantic SST scaling factor beta ~ 0.71 K/K
#'   - Cross-validated against Hibbert et al. (2025) Caribbean CMIP6 SST study
#'
#' @return Tibble with columns: scenario, ssp, variant, delta_sst_2050,
#'   delta_sst_2100, precip_scale, air_temp_2050, air_temp_2100, description.
#' @export
knmi_scenario_info <- function() {
  tibble::tibble(
    scenario       = c("knmi_Ld", "knmi_Ln", "knmi_Hd", "knmi_Hn"),
    ssp            = c("SSP1-2.6", "SSP1-2.6", "SSP5-8.5", "SSP5-8.5"),
    variant        = c("dry", "wet", "dry", "wet"),
    # MDR SST anomalies (degC relative to 1991-2020 baseline)
    delta_sst_2050 = c(0.35, 0.35, 0.80, 0.80),
    delta_sst_2100 = c(0.35, 0.35, 2.60, 2.60),
    # Level 3 precipitation scaling (per degC)
    precip_scale   = c(0.05, 0.10, 0.05, 0.10),
    # BES air temperature projections for reference (degC above 1991-2020)
    air_temp_2050  = c(0.8, 0.8, 1.3, 1.3),
    air_temp_2100  = c(0.8, 0.8, 3.4, 3.4),
    description    = c(
      "Low warming + dry: Paris-aligned, reduced precipitation",
      "Low warming + wet: Paris-aligned, marginal precipitation increase",
      "High warming + dry: fossil-fueled, substantially drier",
      "High warming + wet: fossil-fueled, wetter precipitation response"
    )
  )
}


#' Get KNMI'23 Level 3 cc_params adjusted for scenario variant
#'
#' @description
#' Returns the default cc_params list with the precipitation scaling factor
#' adjusted for the dry (d) or wet (n) variant of the KNMI'23 scenarios.
#'
#' @param scenario Character; one of "knmi_Ld", "knmi_Ln", "knmi_Hd", "knmi_Hn".
#' @param base_params Optional named list; base cc_params to modify.
#'   If NULL, uses `default_cc_params()`.
#'
#' @return Named list of cc_params with adjusted precip_scale.
#' @export
knmi_cc_params <- function(scenario, base_params = NULL) {
  info <- knmi_scenario_info()
  row <- info[info$scenario == scenario, ]
  if (nrow(row) == 0) {
    stop("Unknown KNMI scenario: ", scenario,
         ". Must be one of: ", paste(info$scenario, collapse = ", "))
  }

  if (is.null(base_params)) base_params <- default_cc_params()
  base_params$precip_scale <- row$precip_scale
  base_params
}


# =============================================================================
# 2) UPDATED generate_sst_scenario() — drop-in replacement
# =============================================================================

#' Generate SST anomaly scenarios for simulation (with KNMI'23 support)
#'
#' @description
#' Creates a time series of SST anomalies (delta_SST) for use in non-stationary
#' simulation. Extends the original function with four KNMI'23 scenario modes.
#'
#' Supported modes:
#'
#' 1. **stationary**: Constant anomaly (default: 0).
#' 2. **trend**: Linear ramp from delta_start to delta_end.
#' 3. **ssp245**: CMIP6 SSP2-4.5 trajectory (~+0.5C by 2050, +1.0C by 2100).
#' 4. **ssp585**: CMIP6 SSP5-8.5 trajectory (~+1.0C by 2050, +2.5C by 2100).
#' 5. **historical_resample**: Random resampling of historical anomalies.
#' 6. **knmi_Ld**: KNMI'23 Low warming, dry (SSP1-2.6).
#' 7. **knmi_Ln**: KNMI'23 Low warming, wet (SSP1-2.6).
#' 8. **knmi_Hd**: KNMI'23 High warming, dry (SSP5-8.5).
#' 9. **knmi_Hn**: KNMI'23 High warming, wet (SSP5-8.5).
#'
#' KNMI'23 notes:
#'   - Ld/Ln share the same SST trajectory (SSP1-2.6 stabilizes ~2050).
#'   - Hd/Hn share the same SST trajectory (SSP5-8.5 continuous warming).
#'   - The d/n distinction is metadata only — use `knmi_cc_params()` for L3.
#'   - SST targets derived from KNMI'23 global warming levels × tropical
#'     Atlantic scaling factor (beta ~ 0.71 K/K).
#'
#' @param n_years Integer; number of simulation years.
#' @param mode Character; scenario mode (see above).
#' @param start_year Integer; first calendar year of the simulation.
#' @param delta_start Numeric; starting delta_SST for "trend" mode.
#' @param delta_end Numeric; ending delta_SST for "trend" mode.
#' @param sst_hist_df Optional tibble for "historical_resample" mode.
#' @param baseline_years Baseline period (default 1991:2020).
#'
#' @return Tibble with columns: sim_year, calendar_year, sst_anomaly, scenario.
#' @export
generate_sst_scenario <- function(n_years,
                                  mode = c("stationary", "trend",
                                           "ssp245", "ssp585",
                                           "historical_resample",
                                           "knmi_Ld", "knmi_Ln",
                                           "knmi_Hd", "knmi_Hn"),
                                  start_year = 2025L,
                                  delta_start = 0,
                                  delta_end = 1.0,
                                  sst_hist_df = NULL,
                                  baseline_years = 1991L:2020L) {
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Package `tibble` is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` is required.")

  mode <- match.arg(mode)

  n_years <- as.integer(n_years)
  if (!is.finite(n_years) || n_years <= 0) stop("n_years must be a positive integer.")
  start_year <- as.integer(start_year)

  sim_year <- seq_len(n_years)
  cal_year <- start_year + sim_year - 1L

  # helper: robust linear interpolation (with clamping)
  .lerp <- function(x, x0, x1, y0, y1) {
    if (x1 == x0) return(rep(y1, length(x)))
    t <- (x - x0) / (x1 - x0)
    t <- pmin(1, pmax(0, t))
    y0 + (y1 - y0) * t
  }

  # determine a realistic starting anomaly from historical (builtin) if possible
  .get_start_anom <- function(year0) {
    sst_hist <- get_mdr_sst_builtin() |>
      compute_sst_anomaly(baseline_years = baseline_years)
    a <- sst_hist$sst_anomaly[match(year0, sst_hist$year)]
    if (length(a) == 1 && is.finite(a)) return(as.numeric(a))
    # fallback: use the last available anomaly (2024) if start_year beyond range
    a2 <- sst_hist$sst_anomaly[match(max(sst_hist$year), sst_hist$year)]
    if (length(a2) == 1 && is.finite(a2)) return(as.numeric(a2))
    0
  }

  # ---- KNMI'23 scenario SST targets ----
  # MDR SST anomalies (degC relative to 1991-2020)
  # L (SSP1-2.6): stabilizes at +0.35C around 2050 (Paris-aligned)
  # H (SSP5-8.5): +0.8C by 2050, +2.6C by 2100 (continuous warming)
  .knmi_targets <- list(
    knmi_Ld = list(a_2050 = 0.35, a_2100 = 0.35),
    knmi_Ln = list(a_2050 = 0.35, a_2100 = 0.35),
    knmi_Hd = list(a_2050 = 0.80, a_2100 = 2.60),
    knmi_Hn = list(a_2050 = 0.80, a_2100 = 2.60)
  )

  sst_anom <- switch(
    mode,

    stationary = rep(as.numeric(delta_start), n_years),

    trend = {
      seq(as.numeric(delta_start), as.numeric(delta_end), length.out = n_years)
    },

    ssp245 = {
      a0 <- .get_start_anom(start_year)
      y0 <- start_year; y1 <- 2050L; y2 <- 2100L
      a1 <- 0.5; a2 <- 1.0
      a_t <- ifelse(
        cal_year <= y1,
        .lerp(cal_year, y0, y1, a0, a1),
        ifelse(cal_year <= y2,
               .lerp(cal_year, y1, y2, a1, a2),
               a2)
      )
      as.numeric(a_t)
    },

    ssp585 = {
      a0 <- .get_start_anom(start_year)
      y0 <- start_year; y1 <- 2050L; y2 <- 2100L
      a1 <- 1.0; a2 <- 2.5
      a_t <- ifelse(
        cal_year <= y1,
        .lerp(cal_year, y0, y1, a0, a1),
        ifelse(cal_year <= y2,
               .lerp(cal_year, y1, y2, a1, a2),
               a2)
      )
      as.numeric(a_t)
    },

    historical_resample = {
      if (is.null(sst_hist_df)) {
        sst_hist_df <- get_mdr_sst_builtin() |>
          compute_sst_anomaly(baseline_years = baseline_years)
      }
      if (!("sst_anomaly" %in% names(sst_hist_df))) {
        stop("sst_hist_df must contain 'sst_anomaly'. Run compute_sst_anomaly() first.")
      }
      pool <- sst_hist_df$sst_anomaly[is.finite(sst_hist_df$sst_anomaly)]
      if (length(pool) == 0) stop("No valid SST anomalies to resample.")
      sample(pool, n_years, replace = TRUE)
    },

    # ---- KNMI'23 scenarios ----
    # All four KNMI modes use the same piecewise-linear logic with
    # scenario-specific targets. d/n share SST; distinction is L3 only.
    knmi_Ld = , knmi_Ln = , knmi_Hd = , knmi_Hn = {
      targets <- .knmi_targets[[mode]]
      a0 <- .get_start_anom(start_year)
      y0 <- start_year; y1 <- 2050L; y2 <- 2100L
      a1 <- targets$a_2050
      a2 <- targets$a_2100

      a_t <- ifelse(
        cal_year <= y1,
        .lerp(cal_year, y0, y1, a0, a1),
        ifelse(cal_year <= y2,
               .lerp(cal_year, y1, y2, a1, a2),
               a2)
      )
      as.numeric(a_t)
    }
  )

  tibble::tibble(
    sim_year = sim_year,
    calendar_year = as.integer(cal_year),
    sst_anomaly = as.numeric(sst_anom),
    scenario = mode
  )
}


# =============================================================================
# 3) CONVENIENCE: make_sst_cfg() with KNMI scenario support
# =============================================================================

#' Build a climate configuration object (extended for KNMI'23)
#'
#' @description
#' Extended version of make_sst_cfg() that accepts KNMI'23 scenario names
#' and auto-adjusts Level 3 cc_params for the d/n precipitation variant.
#'
#' When scenario is a KNMI'23 mode and cc_params is an empty list(),
#' the precipitation scaling is automatically set based on the d/n variant.
#'
#' @inheritParams make_sst_cfg
#' @export
make_sst_cfg <- function(enabled = TRUE,
                         sst_source = c("builtin", "csv", "ersst_nc"),
                         sst_path = NULL,
                         baseline_years = 1991L:2020L,
                         beta_sst = NULL,
                         beta_prior = 0.6,
                         gamma_intensity = NULL,
                         gamma_prior = 0.065,
                         scenario = "stationary",
                         scenario_start_year = 2025L,
                         cc_params = NULL) {
  sst_source <- match.arg(sst_source)

  # Resolve Level 3 cc_params (NULL = disabled; list() = use defaults)
  if (!is.null(cc_params) && !is.list(cc_params)) {
    stop("cc_params must be NULL (disabled) or a named list of scaling factors.")
  }

  # Auto-adjust cc_params for KNMI scenarios when cc_params = list()
  is_knmi <- grepl("^knmi_", scenario)
  if (is_knmi && !is.null(cc_params) && length(cc_params) == 0) {
    cc_params <- knmi_cc_params(scenario)
  }

  cfg <- list(
    enabled = enabled,
    sst_source = sst_source,
    sst_path = sst_path,
    baseline_years = baseline_years,
    beta_sst = beta_sst,
    beta_prior = beta_prior,
    gamma_intensity = gamma_intensity,
    gamma_prior = gamma_prior,
    scenario = scenario,
    scenario_start_year = scenario_start_year,
    cc_params = cc_params
  )
  class(cfg) <- c("sst_cfg", "list")
  cfg
}
