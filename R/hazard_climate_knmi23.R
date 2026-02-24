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
# NOTE
# =============================================================================
# Core computation functions (generate_sst_scenario(), make_sst_cfg(), etc.) live in
# hazard_climate.R. This file provides KNMI'23 climate-information sources only.
