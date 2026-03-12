################################################################################
# hazard_climate.R
# Climate workflow for hazard simulations
#
# Default climate-adjusted workflow:
#   - Rate effect via SST-conditioned activity scaling (beta_SST)
#   - Intensity effect via SST-conditioned TS/HUR split adjustment (gamma)
#
# Optional expert extension:
#   - Storm perturbation of sampled event properties
#
# Modulates the annual activity factor using observed/projected SST anomalies
# in the Main Development Region (MDR: 10a?"20A?N, 80a?"20A?W):
#
#   A_t ~ Gamma(k, k) A- exp(I2_SST A. I"SST_t)
#
# Components:
#   1) MDR SST data ingestion (ERSST v5 CSV/NetCDF + built-in fallback)
#   2) Anomaly computation relative to 1991a?"2020 climatological mean
#   3) I2_SST estimation via Poisson/NegBin regression on historical record
#   4) SST scenario generation (historical + CMIP6 SSP projections)
#   5) Non-stationary count simulation
#   6) Intensity effect (gamma estimation)
#   7) Optional storm perturbation (perturb_event())
#
# Data sources:
#   - Historical: NOAA ERSST v5 (https://www.ncei.noaa.gov/products/extended-reconstructed-sst)
#   - Projections: CMIP6 `tos` variable (SSP2-4.5 and SSP5-8.5)
################################################################################


# =============================================================================
# 1) BUILT-IN MDR SST REFERENCE DATA
# =============================================================================

#' Built-in MDR SST annual means from NOAA ERSST v5
#'
#' @description
#' Returns a tibble of annual-mean SST (degC) averaged over the Main Development
#' Region (MDR: 10-20N, 80-20W) for the Atlantic hurricane season (Aug-Oct,
#' the peak months most predictive of activity).
#'
#' Values are derived from NOAA ERSST v5 monthly data, spatially averaged over
#' the MDR box, then temporally averaged over ASO (Aug-Sep-Oct) each year.
#'
#' These are provided as a built-in fallback so the model can run without
#' requiring users to download and process NetCDF files.
#'
#' Source: NOAA/NCEI ERSST v5, accessed 2024.
#' Reference: Huang et al. (2017), J. Climate, 30, 8179-8205.
#'
#' @return Tibble with columns:
#' \itemize{
#'   \item `year`: integer calendar year.
#'   \item `sst_mdr_aso`: MDR ASO seasonal mean SST (degC).
#' }
#' @export
get_mdr_sst_builtin <- function() {
  # MDR (10-20N, 80-20W) ASO mean SST from ERSST v5
  # These values represent the Aug-Sep-Oct average over the MDR box.
  tibble::tibble(
    year = 1970L:2024L,
    sst_mdr_aso = c(
      # 1970-1979
      26.86, 26.64, 26.72, 26.80, 26.62, 26.68, 26.58, 26.82, 26.80, 26.78,
      # 1980-1989
      26.92, 26.74, 26.66, 26.84, 26.76, 26.80, 26.72, 26.88, 26.72, 26.82,
      # 1990-1999
      26.88, 26.74, 26.60, 26.78, 26.90, 27.00, 27.04, 26.82, 27.18, 26.96,
      # 2000-2009
      27.04, 27.10, 27.08, 27.16, 27.22, 27.34, 27.16, 27.10, 27.14, 27.04,
      # 2010-2019
      27.38, 27.06, 27.18, 27.14, 27.08, 27.22, 27.26, 27.28, 27.20, 27.22,
      # 2020-2024
      27.30, 27.18, 27.28, 27.48, 27.62
    )
  )
}


# =============================================================================
# 2) SST DATA INGESTION
# =============================================================================

#' Read MDR SST from a CSV file
#'
#' @description
#' Reads a user-supplied CSV containing annual MDR SST values. The CSV must
#' have at minimum a `year` column and an SST column (name specified by
#' `sst_col`). This allows users to supply their own ERSST v5 extractions
#' or alternative SST products.
#'
#' @param csv_path Character; path to the CSV file.
#' @param sst_col Character; name of the SST column (default: "sst_mdr_aso").
#' @param year_col Character; name of the year column (default: "year").
#'
#' @return Tibble with columns: year, sst_mdr_aso.
#' @export
read_mdr_sst_csv <- function(csv_path,
                             sst_col = "sst_mdr_aso",
                             year_col = "year") {
  if (!file.exists(csv_path)) stop("SST CSV not found: ", csv_path)
  if (!requireNamespace("readr", quietly = TRUE)) stop("Package `readr` is required.")

  df <- readr::read_csv(csv_path, show_col_types = FALSE)

  if (!(year_col %in% names(df))) stop("Year column '", year_col, "' not found in CSV.")
  if (!(sst_col %in% names(df))) stop("SST column '", sst_col, "' not found in CSV.")

  tibble::tibble(
    year = as.integer(df[[year_col]]),
    sst_mdr_aso = as.numeric(df[[sst_col]])
  ) |>
    dplyr::filter(is.finite(.data$year), is.finite(.data$sst_mdr_aso)) |>
    dplyr::arrange(.data$year)
}


#' Read MDR SST from ERSST v5 NetCDF (optional)
#'
#' @description
#' Reads monthly ERSST v5 NetCDF data, subsets to the MDR box (10-20N, 80-20W),
#' averages spatially, then computes ASO (Aug-Sep-Oct) seasonal means per year.
#'
#' Requires the `ncdf4` package. If not available, falls back to the built-in
#' reference data with a warning.
#'
#' @param nc_path Character; path to ERSST v5 NetCDF file (e.g., "sst.mnmean.nc").
#' @param mdr_lat Range of latitudes for MDR (default: c(10, 20)).
#' @param mdr_lon Range of longitudes for MDR (default: c(-80, -20)).
#'   Note: ERSST uses 0-360 longitude convention; this function converts internally.
#' @param aso_months Integer vector of months for seasonal average (default: 8:10).
#'
#' @return Tibble with columns: year, sst_mdr_aso.
#' @export
read_mdr_sst_ersst <- function(nc_path,
                               mdr_lat = c(10, 20),
                               mdr_lon = c(-80, -20),
                               aso_months = 8L:10L) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    warning("Package `ncdf4` not available. Falling back to built-in MDR SST data.")
    return(get_mdr_sst_builtin())
  }
  if (!file.exists(nc_path)) stop("ERSST NetCDF not found: ", nc_path)

  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  lon_nc <- ncdf4::ncvar_get(nc, "lon")   # 0-360 convention in ERSST
  lat_nc <- ncdf4::ncvar_get(nc, "lat")
  time_nc <- ncdf4::ncvar_get(nc, "time")

  # Convert MDR lon to 0-360 if needed
  mdr_lon_360 <- ifelse(mdr_lon < 0, mdr_lon + 360, mdr_lon)

  lon_idx <- which(lon_nc >= min(mdr_lon_360) & lon_nc <= max(mdr_lon_360))
  lat_idx <- which(lat_nc >= min(mdr_lat) & lat_nc <= max(mdr_lat))

  if (length(lon_idx) == 0 || length(lat_idx) == 0) {
    warning("No grid cells found in MDR box. Falling back to built-in data.")
    return(get_mdr_sst_builtin())
  }

  # Read SST for MDR subset
  sst <- ncdf4::ncvar_get(nc, "sst",
                           start = c(min(lon_idx), min(lat_idx), 1),
                           count = c(length(lon_idx), length(lat_idx), -1))

  # Time a?' dates (ERSST time is days since 1800-01-01)
  time_origin <- as.Date("1800-01-01")
  dates <- time_origin + as.integer(time_nc)
  years <- as.integer(format(dates, "%Y"))
  months <- as.integer(format(dates, "%m"))

  # Spatial average (area-weighted by cos(lat))
  lat_sub <- lat_nc[lat_idx]
  weights <- cos(lat_sub * pi / 180)
  weights <- weights / sum(weights)

  n_time <- dim(sst)[3]
  sst_mdr <- rep(NA_real_, n_time)

  for (t in seq_len(n_time)) {
    slice <- sst[, , t]
    if (all(is.na(slice))) next
    # Weighted average: mean over lon, then weighted mean over lat
    lon_means <- colMeans(slice, na.rm = TRUE)
    sst_mdr[t] <- sum(lon_means * weights, na.rm = TRUE)
  }

  # ASO seasonal means
  df <- tibble::tibble(year = years, month = months, sst = sst_mdr) |>
    dplyr::filter(.data$month %in% aso_months, is.finite(.data$sst)) |>
    dplyr::group_by(.data$year) |>
    dplyr::summarise(sst_mdr_aso = mean(.data$sst, na.rm = TRUE), .groups = "drop") |>
    dplyr::filter(is.finite(.data$sst_mdr_aso)) |>
    dplyr::arrange(.data$year)

  df
}


# =============================================================================
# 3) SST ANOMALY COMPUTATION
# =============================================================================

#' Compute SST anomalies relative to a climatological baseline
#'
#' @description
#' Computes `delta_SST_t = SST_t - SST_clim`, where `SST_clim` is the mean SST over
#' the specified baseline period (default: 1991-2020, the current WMO standard
#' climatological normal).
#'
#' @param sst_df Tibble with columns: year, sst_mdr_aso.
#' @param baseline_years Integer vector of years defining the climatological
#'   baseline (default: 1991:2020).
#'
#' @return The input tibble with added columns:
#'   \item{sst_clim}{Climatological mean SST (degC) over the baseline period.}
#'   \item{sst_anomaly}{SST anomaly (degC) relative to baseline.}
#'
#' @examples
#' sst <- get_mdr_sst_builtin()
#' sst_anom <- compute_sst_anomaly(sst, baseline_years = 1991:2020)
#' head(sst_anom)
#'
#' @export
compute_sst_anomaly <- function(sst_df, baseline_years = 1991L:2020L) {
  if (!all(c("year", "sst_mdr_aso") %in% names(sst_df))) {
    stop("sst_df must contain columns: year, sst_mdr_aso")
  }

  baseline_sst <- sst_df |>
    dplyr::filter(.data$year %in% baseline_years) |>
    dplyr::pull(.data$sst_mdr_aso)

  if (length(baseline_sst) < 10) {
    warning("Fewer than 10 years overlap with baseline period (",
            min(baseline_years), "-", max(baseline_years),
            "). Anomalies may be unreliable.")
  }

  sst_clim <- mean(baseline_sst, na.rm = TRUE)

  sst_df |>
    dplyr::mutate(
      sst_clim = sst_clim,
      sst_anomaly = .data$sst_mdr_aso - sst_clim
    )
}


# =============================================================================
# 4) I2_SST ESTIMATION
# =============================================================================

#' Estimate the SST-activity scaling coefficient I2_SST
#'
#' @description
#' Fits a Poisson (or negative binomial) GLM of annual TS counts on MDR SST
#' anomaly to estimate I2_SST in:
#'
#'   `E[N_t] = exp(alpha + beta_SST * sst_anomaly_t)`
#'
#' The coefficient I2_SST represents the log-linear sensitivity of annual TS
#' activity to SST anomalies. Typical values from the literature are around
#' 0.4-0.8 per degC for the North Atlantic (Villarini et al. 2011;
#' Vecchi et al. 2021).
#'
#' If the `MASS` package is available, a negative binomial GLM is preferred
#' (accounts for overdispersion). Otherwise, a quasi-Poisson GLM is used.
#'
#' @param annual_counts Tibble from `compute_annual_counts()` with columns:
#'   year, storm_class, n_events.
#' @param sst_df Tibble with columns: year, sst_anomaly.
#' @param min_year Integer; earliest year to include (default: cfg$start_year).
#' @param beta_prior Optional numeric; if provided, shrinks the estimate toward
#'   this prior value using a simple Bayesian-inspired weighted average:
#'   beta_final = w * beta_mle + (1 - w) * beta_prior, where w = min(1, n_years/50).
#'   This stabilizes estimates when the historical record is short. Climate-rate
#'   calibration must use basin-consistent annual counts in storms/year; passing
#'   target-conditioned counts with a `location` column is rejected because that
#'   duplicates annual activity across targets.
#' @param verbose Logical; print diagnostic output.
#'
#' @return A list with:
#'   \item{beta_sst}{Estimated (or shrunk) I2_SST coefficient.}
#'   \item{beta_se}{Standard error of the MLE I2_SST.}
#'   \item{beta_mle}{Raw MLE estimate before shrinkage.}
#'   \item{alpha}{Intercept (log baseline rate).}
#'   \item{method}{Character; "negbin", "quasipoisson", or "literature_fallback".}
#'   \item{n_years}{Number of years used in estimation.}
#'   \item{r_squared_dev}{Deviance-based pseudo-RA2 (proportion of deviance explained).}
#'   \item{aic}{AIC of the fitted model (NA for quasi-Poisson).}
#'   \item{annual_count_series}{Tibble of the basin-consistent annual total count
#'     series used for calibration, in storms/year.}
#'   \item{annual_count_source}{Character label describing the provenance of
#'     `annual_count_series`.}
#'   \item{guardrail}{List describing any regularization or fallback applied to
#'     keep `beta_sst` within documented plausibility bounds.}
#'   \item{fit_data}{Tibble of the joined data used for fitting.}
#'   \item{glm_fit}{The fitted GLM object.}
#'
#' @examples
#' # Small example workflow:
#' sst_df <- get_mdr_sst_builtin() |>
#'   compute_sst_anomaly()
#' # annual_counts should come from compute_annual_counts()
#' # beta_info <- estimate_beta_sst(annual_counts, sst_df)
#'
#' @export
estimate_beta_sst <- function(annual_counts,
                              sst_df,
                              min_year = 1970L,
                              beta_prior = NULL,
                              verbose = TRUE) {
  beta_guardrail_max <- 1.2
  beta_guardrail_se_ref <- 0.35
  beta_guardrail_fallback <- if (is.null(beta_prior) || !is.finite(beta_prior)) 0.6 else as.numeric(beta_prior)

  if ("location" %in% names(annual_counts)) {
    stop(
      "annual_counts used for climate calibration must be basin-consistent and must not contain a location column. ",
      "Pass counts from de-duplicated basin events.",
      call. = FALSE
    )
  }

  # Aggregate to total annual counts (across severities)
  annual_count_series <- annual_counts |>
    dplyr::filter(.data$year >= min_year) |>
    dplyr::group_by(.data$year) |>
    dplyr::summarise(N = sum(.data$n_events), .groups = "drop")

  # Require sst_anomaly column
  if (!("sst_anomaly" %in% names(sst_df))) {
    stop("sst_df must contain 'sst_anomaly'. Run compute_sst_anomaly() first.")
  }

  # Join
  fit_data <- annual_count_series |>
    dplyr::inner_join(
      sst_df |> dplyr::select("year", "sst_anomaly"),
      by = "year"
    ) |>
    dplyr::filter(is.finite(.data$N), is.finite(.data$sst_anomaly))

  n_years <- nrow(fit_data)

  if (n_years < 15) {
    warning("Only ", n_years, " years of overlapping data. ",
            "I2_SST estimate will be unreliable. Consider using beta_prior.")
  }

  if (n_years < 5) {
    message("[I2_SST] Insufficient data (", n_years, " years). Using literature fallback: I2 = 0.6")
    return(list(
      beta_sst = 0.6,
      beta_se = NA_real_,
      beta_mle = NA_real_,
      alpha = NA_real_,
      method = "literature_fallback",
      n_years = n_years,
      r_squared_dev = NA_real_,
      aic = NA_real_,
      annual_count_series = annual_count_series,
      annual_count_source = "basin_unique_storm_year_counts",
      guardrail = list(
        triggered = TRUE,
        reason = "insufficient_overlap",
        beta_max = beta_guardrail_max,
        beta_fallback = 0.6
      ),
      fit_data = fit_data,
      glm_fit = NULL
    ))
  }

  # Fit GLM: N ~ exp(I? + I2A.I"SST)
  # Prefer negative binomial (handles overdispersion); fall back to quasi-Poisson
  glm_fit <- NULL
  method <- "quasipoisson"

  if (requireNamespace("MASS", quietly = TRUE)) {
    glm_fit <- tryCatch({
      MASS::glm.nb(N ~ sst_anomaly, data = fit_data)
    }, error = function(e) NULL)
    if (!is.null(glm_fit)) method <- "negbin"
  }

  if (is.null(glm_fit)) {
    glm_fit <- stats::glm(N ~ sst_anomaly, data = fit_data,
                           family = stats::quasipoisson(link = "log"))
    method <- "quasipoisson"
  }

  coefs <- stats::coef(glm_fit)
  alpha_hat <- coefs["(Intercept)"]
  beta_mle <- coefs["sst_anomaly"]

  # Standard error
  se <- tryCatch(
    stats::summary.glm(glm_fit)$coefficients["sst_anomaly", "Std. Error"],
    error = function(e) NA_real_
  )

  # Deviance pseudo-RA2
  dev_null <- glm_fit$null.deviance
  dev_resid <- glm_fit$deviance
  r2_dev <- if (is.finite(dev_null) && dev_null > 0) {
    1 - dev_resid / dev_null
  } else NA_real_

  # AIC (not available for quasi-Poisson)
  aic_val <- tryCatch(stats::AIC(glm_fit), error = function(e) NA_real_)

  # Optional shrinkage toward prior
  beta_final <- beta_mle
  if (!is.null(beta_prior) && is.finite(beta_prior)) {
    w_n <- min(1.0, n_years / 50)
    w_se <- if (is.finite(se) && se > 0) min(1.0, beta_guardrail_se_ref / se) else 0
    w <- w_n * w_se
    beta_final <- w * beta_mle + (1 - w) * beta_prior
    if (verbose) {
      message(sprintf("[I2_SST] Shrinkage: I2_MLE=%.3f, I2_prior=%.3f, w=%.2f a?' I2_final=%.3f",
                      beta_mle, beta_prior, w, beta_final))
    }
  }

  guardrail <- list(
    triggered = FALSE,
    reason = "none",
    beta_max = beta_guardrail_max,
    beta_fallback = beta_guardrail_fallback
  )
  if (!is.finite(beta_final)) {
    beta_final <- beta_guardrail_fallback
    guardrail$triggered <- TRUE
    guardrail$reason <- "non_finite_beta"
  }
  if (is.finite(beta_final) && beta_final > beta_guardrail_max) {
    if (verbose) {
      message(sprintf(
        "[I2_SST] Guardrail: beta_final=%.3f exceeds max %.3f 1/degC; using fallback %.3f.",
        beta_final, beta_guardrail_max, beta_guardrail_fallback
      ))
    }
    beta_final <- beta_guardrail_fallback
    guardrail$triggered <- TRUE
    guardrail$reason <- "beta_above_plausibility_limit"
  }

  if (verbose) {
    message(sprintf("[I2_SST] Method: %s | n_years: %d | I2_SST: %.3f (SE: %.3f) | pseudo-RA2: %.3f",
                    method, n_years, beta_final, se, r2_dev))
    message(sprintf("[I2_SST] Interpretation: +1A?C SST a?' %.0f%% change in annual rate",
                    100 * (exp(beta_final) - 1)))
  }

  list(
    beta_sst = as.numeric(beta_final),
    beta_se = as.numeric(se),
    beta_mle = as.numeric(beta_mle),
    alpha = as.numeric(alpha_hat),
    method = method,
    n_years = n_years,
    r_squared_dev = r2_dev,
    aic = aic_val,
    annual_count_series = annual_count_series,
    annual_count_source = "basin_unique_storm_year_counts",
    guardrail = guardrail,
    fit_data = fit_data,
    glm_fit = glm_fit
  )
}


# =============================================================================
# 4b) INTENSITY DISTRIBUTION SHIFT
# =============================================================================

#' Extract TS/HUR lambdas from a lambda table
#'
#' @param lambda_table Tibble from \code{compute_lambda_table()} with
#'   severities "TS" and "HUR".
#'
#' @return List with elements \code{ts}, \code{hur}, \code{total}.
#' @keywords internal
.extract_lambdas <- function(lambda_table) {
  lam_ts  <- lambda_table$lambda[lambda_table$storm_class == "TS"]
  lam_hur <- lambda_table$lambda[lambda_table$storm_class == "HUR"]
  if (length(lam_ts) == 0)  lam_ts  <- 0
  if (length(lam_hur) == 0) lam_hur <- 0
  list(
    ts = as.numeric(lam_ts),
    hur = as.numeric(lam_hur),
    total = as.numeric(lam_ts) + as.numeric(lam_hur)
  )
}

#' Estimate the historical hurricane fraction from annual counts
#'
#' @description
#' Computes p_HUR_base = I>_HUR / (I>_TS + I>_HUR) from the historical record.
#' This is the baseline probability that an event reaching TS intensity or
#' above becomes a hurricane (>=64 kt).
#'
#' @param lambda_table Tibble from `compute_lambda_table()` with severities
#'   "TS" and "HUR".
#'
#' @return Numeric scalar: baseline hurricane fraction p_HUR_base.
#' @export
compute_p_hur_base <- function(lambda_table) {
  lam <- .extract_lambdas(lambda_table)
  if (lam$total <= 0) return(0.5)  # safeguard
  as.numeric(lam$hur / lam$total)
}


#' Estimate I3 (intensity shift coefficient) from historical data
#'
#' @description
#' Fits a logistic regression of annual hurricane fraction on SST anomaly
#' to estimate the intensification trend coefficient I3 in:
#'
#'   `p_HUR(t) = p_HUR_base * (1 + gamma * sst_anomaly_t)`
#'
#' The coefficient I3 captures how the fraction of storms reaching hurricane
#' intensity shifts with SST warming. Literature estimates suggest roughly
#' 5-8% increase in Cat 4-5 fraction per degC of tropical SST warming
#' (Knutson et al. 2020;
#' Kossin et al. 2020).
#'
#' Uses a binomial GLM: cbind(n_HUR, n_TS) ~ sst_anomaly, then converts the
#' logistic coefficient to the linear I3 parameterization.
#'
#' @param annual_counts Tibble with columns: year, storm_class, n_events.
#' @param sst_df Tibble with columns: year, sst_anomaly.
#' @param min_year Integer; earliest year to include.
#' @param gamma_prior Optional numeric; prior value for shrinkage (default: 0.065,
#'   i.e. 6.5% increase in HUR fraction per degC).
#' @param verbose Logical.
#'
#' @return A list with:
#'   \item{gamma}{Estimated (or shrunk) I3 coefficient.}
#'   \item{gamma_mle}{Raw MLE from logistic regression (converted to linear scale).}
#'   \item{gamma_se}{Approximate standard error on I3.}
#'   \item{p_hur_base}{Baseline hurricane fraction.}
#'   \item{method}{Character; estimation method used.}
#'   \item{n_years}{Number of years in fit.}
#'   \item{fit_data}{Tibble used for fitting.}
#'
#' @export
estimate_gamma_intensity <- function(annual_counts,
                                     sst_df,
                                     min_year = 1970L,
                                     gamma_prior = 0.065,
                                     verbose = TRUE) {
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package `tidyr` is required for gamma estimation.")
  if (!("sst_anomaly" %in% names(sst_df))) {
    stop("sst_df must contain 'sst_anomaly'. Run compute_sst_anomaly() first.")
  }

  # Pivot annual counts to wide: year, n_TS, n_HUR
  # If multi-island, aggregate to basin-wide totals first
  ac <- annual_counts |>
    dplyr::filter(.data$year >= min_year)

  if ("location" %in% names(ac)) {
    ac <- ac |>
      dplyr::group_by(.data$year, .data$storm_class) |>
      dplyr::summarise(n_events = sum(.data$n_events), .groups = "drop")
  }

  wide <- ac |>
    tidyr::pivot_wider(
      id_cols = "year",
      names_from = "storm_class",
      values_from = "n_events",
      values_fill = 0L
    )

  # Ensure both columns exist
  if (!("TS" %in% names(wide))) wide$TS <- 0L
  if (!("HUR" %in% names(wide))) wide$HUR <- 0L

  # Join with SST
  fit_data <- wide |>
    dplyr::inner_join(
      sst_df |> dplyr::select("year", "sst_anomaly"),
      by = "year"
    ) |>
    dplyr::filter(is.finite(.data$sst_anomaly)) |>
    dplyr::mutate(
      n_total = .data$TS + .data$HUR,
      p_hur = ifelse(.data$n_total > 0, .data$HUR / .data$n_total, NA_real_)
    ) |>
    dplyr::filter(.data$n_total > 0)  # need at least 1 event to estimate fraction

  n_years <- nrow(fit_data)

  if (n_years < 10) {
    if (verbose) message(sprintf("[gamma] Only %d years with events. Using prior: gamma=%.3f", n_years, gamma_prior))
    return(list(
      gamma = gamma_prior,
      gamma_mle = NA_real_,
      gamma_se = NA_real_,
      p_hur_base = mean(fit_data$p_hur, na.rm = TRUE),
      method = "prior_fallback",
      n_years = n_years,
      fit_data = fit_data
    ))
  }

  # Baseline hurricane fraction (overall)
  p_hur_base <- sum(fit_data$HUR) / sum(fit_data$n_total)

  # Binomial GLM: logit(p_HUR) ~ sst_anomaly
  glm_fit <- tryCatch({
    stats::glm(
      cbind(HUR, TS) ~ sst_anomaly,
      data = fit_data,
      family = stats::binomial(link = "logit")
    )
  }, error = function(e) NULL)

  if (is.null(glm_fit)) {
    if (verbose) message("[gamma] Binomial GLM failed. Using prior.")
    return(list(
      gamma = gamma_prior,
      gamma_mle = NA_real_,
      gamma_se = NA_real_,
      p_hur_base = p_hur_base,
      method = "prior_fallback",
      n_years = n_years,
      fit_data = fit_data
    ))
  }

  # Extract logistic coefficient and convert to linear I3
  # In logistic model: logit(p) = I? + I2_logistic A. I"SST
  # At the baseline: p_base = logistic(I?)
  # Linear approximation: dp/dSST a?^ p_base(1-p_base) A. I2_logistic
  # So: I3 = dp/dSST / p_base = (1-p_base) A. I2_logistic
  beta_logistic <- stats::coef(glm_fit)["sst_anomaly"]
  gamma_mle <- as.numeric((1 - p_hur_base) * beta_logistic)

  # Standard error (delta method)
  se_logistic <- tryCatch(
    summary(glm_fit)$coefficients["sst_anomaly", "Std. Error"],
    error = function(e) NA_real_
  )
  gamma_se <- as.numeric((1 - p_hur_base) * se_logistic)

  # Shrinkage toward prior
  w <- min(1.0, n_years / 50)
  gamma_final <- w * gamma_mle + (1 - w) * gamma_prior
  gamma_nonneg <- max(0, gamma_final)

  if (verbose) {
    message(sprintf("[gamma] Method: binomial_glm | n_years: %d | gamma_MLE: %.4f (SE: %.4f)",
                    n_years, gamma_mle, gamma_se))
    message(sprintf("[gamma] Shrinkage: w=%.2f -> gamma_final: %.4f", w, gamma_final))
    if (gamma_final < 0) {
      message(sprintf("[gamma] Constrained to non-negative: %.4f -> 0.0000", gamma_final))
    }
    message(sprintf("[gamma] p_HUR_base: %.3f | +1C SST -> p_HUR changes by %+.1f%%",
                    p_hur_base, 100 * gamma_nonneg))
  }

  list(
    gamma = as.numeric(gamma_nonneg),
    gamma_mle = as.numeric(gamma_mle),
    gamma_se = as.numeric(gamma_se),
    p_hur_base = as.numeric(p_hur_base),
    method = "binomial_glm",
    n_years = n_years,
    fit_data = fit_data
  )
}


#' Compute time-varying storm-class split
#'
#' @description
#' Computes per-year hurricane probability and corresponding TS/HUR rates:
#'
#'   `p_HUR(t) = clamp(p_HUR_base * (1 + gamma * sst_anomaly_t), 0.01, 0.99)`
#'   `lambda_HUR(t) = lambda_total * p_HUR(t)`
#'   `lambda_TS(t)  = lambda_total * (1 - p_HUR(t))`
#'
#' @param lambda_table Tibble from `compute_lambda_table()`.
#' @param sst_anomaly Numeric vector of SST anomalies per simulation year (degC).
#' @param gamma Numeric; intensity shift coefficient.
#' @param p_hur_base Optional numeric; if NULL, computed from lambda_table.
#'
#' @return Tibble with columns: sim_year, sst_anomaly, p_hur, lam_TS, lam_HUR.
#' @export
compute_severity_split <- function(lambda_table,
                                   sst_anomaly,
                                   gamma = 0,
                                   p_hur_base = NULL) {
  lam <- .extract_lambdas(lambda_table)
  lam_total <- lam$total

  if (is.null(p_hur_base)) {
    p_hur_base <- compute_p_hur_base(lambda_table)
  }

  n <- length(sst_anomaly)

  # Time-varying p_HUR, clamped to (0.01, 0.99)
  p_hur_t <- pmin(0.99, pmax(0.01, p_hur_base * (1 + gamma * sst_anomaly)))

  tibble::tibble(
    sim_year = seq_len(n),
    sst_anomaly = sst_anomaly,
    p_hur_base = p_hur_base,
    p_hur = p_hur_t,
    lam_total = lam_total,
    lam_TS = lam_total * (1 - p_hur_t),
    lam_HUR = lam_total * p_hur_t
  )
}


# =============================================================================
# 5A) SST SCENARIO METADATA (CLIMATE INFORMATION SOURCES)
# =============================================================================

#' IPCC AR6 / CMIP6 SSP MDR SST anomaly targets
#'
#' @description
#' Returns a tibble describing default MDR SST anomaly targets (degC relative to
#' the 1991-2020 baseline) used by the built-in SSP scenario generator.
#'
#' These are deliberately simple, stakeholder-facing targets used for the
#' piecewise-linear SST anomaly generator (to 2050 and 2100, then held constant).
#' They are not intended to reproduce a particular CMIP6 model member.
#'
#' @return Tibble with columns: scenario, source, delta_sst_2050, delta_sst_2100, description.
#' @export
ipcc_ar6_sst_scenario_info <- function() {
  tibble::tibble(
    scenario       = c("ssp126", "ssp245", "ssp585"),
    source         = "ipcc_ar6",
    delta_sst_2050 = c(0.3, 0.5, 1.0),
    delta_sst_2100 = c(0.4, 1.0, 2.5),
    description    = c(
      "Low emissions pathway (order-of-magnitude AR6)",
      "Intermediate emissions pathway (order-of-magnitude AR6)",
      "High emissions pathway (order-of-magnitude AR6)"
    )
  )
}

#' Retrieve available SST scenario definitions
#'
#' @description
#' Returns a combined scenario table for all available climate-information sources.
#' If KNMI'23 support is available (knmi_scenario_info()), those scenarios are
#' included as well.
#'
#' @param source Character; one of "all", "ipcc_ar6", "knmi23".
#'
#' @return Tibble with at least: scenario, source, delta_sst_2050, delta_sst_2100.
#' @export
sst_scenario_info <- function(source = c("all", "ipcc_ar6", "knmi23")) {

  source <- match.arg(source)
  ipcc <- ipcc_ar6_sst_scenario_info()

  knmi <- NULL
  if (exists("knmi_scenario_info", mode = "function")) {
    kn <- knmi_scenario_info()
    # normalize to shared column set
    knmi <- kn |>
      dplyr::transmute(
        scenario = .data$scenario,
        source = "knmi23",
        delta_sst_2050 = .data$delta_sst_2050,
        delta_sst_2100 = .data$delta_sst_2100,
        description = if ("description" %in% names(kn)) as.character(.data$description) else ""
      )
  } else {
    knmi <- tibble::tibble(
      scenario = character(0),
      source = character(0),
      delta_sst_2050 = numeric(0),
      delta_sst_2100 = numeric(0),
      description = character(0)
    )
  }

  out <- dplyr::bind_rows(ipcc, knmi)

  out <- switch(
    source,
    all = out,
    ipcc_ar6 = out[out$source == "ipcc_ar6", , drop = FALSE],
    knmi23   = out[out$source == "knmi23", , drop = FALSE]
  )

  out
}

# =============================================================================
# 5) TIME-SLICE CLIMATE INPUTS
# =============================================================================

.validate_baseline_years <- function(baseline_years, arg = "baseline_years") {
  baseline_years <- as.integer(baseline_years)
  if (length(baseline_years) < 10L || any(!is.finite(baseline_years))) {
    stop(arg, " must be an integer vector of length >= 10.", call. = FALSE)
  }
  baseline_years
}

#' Look up a time-slice scenario SST shift
#'
#' @description
#' Returns a scalar `delta_sst` for a future time slice by interpolating the
#' scenario targets in `sst_scenario_info("all")` at the midpoint of
#' `future_period`.
#'
#' @param scenario Character scalar naming a scenario in `sst_scenario_info("all")`.
#' @param future_period Integer vector describing the future period of interest.
#' @param baseline_years Integer vector of climatological reference years.
#'
#' @return Numeric scalar `delta_sst` in degC.
#' @examples
#' get_scenario_delta("ssp585", future_period = 2035:2065)
#' @export
get_scenario_delta <- function(scenario,
                               future_period,
                               baseline_years = 1991L:2020L) {

  .validate_baseline_years(baseline_years)

  if (!is.character(scenario) || length(scenario) != 1L || !nzchar(scenario)) {
    stop("scenario must be a single non-empty character value.", call. = FALSE)
  }

  future_period <- as.integer(future_period)
  if (length(future_period) == 0L || any(!is.finite(future_period))) {
    stop("future_period must contain finite integer years.", call. = FALSE)
  }

  info <- sst_scenario_info("all")

  row <- info[info$scenario == scenario, , drop = FALSE]
  if (nrow(row) == 0L) {
    stop(
      "Unknown SST scenario: ", scenario,
      ". Available: ", paste(info$scenario, collapse = ", "),
      call. = FALSE
    )
  }

  midpoint <- mean(future_period)
  y <- if (midpoint <= 2025) {
    0
  } else if (midpoint <= 2050) {
    (midpoint - 2025) / (2050 - 2025) * as.numeric(row$delta_sst_2050[[1]])
  } else if (midpoint <= 2100) {
    frac <- (midpoint - 2050) / (2100 - 2050)
    as.numeric(row$delta_sst_2050[[1]]) +
      frac * (as.numeric(row$delta_sst_2100[[1]]) - as.numeric(row$delta_sst_2050[[1]]))
  } else {
    as.numeric(row$delta_sst_2100[[1]])
  }

  if (!is.finite(y)) {
    stop("Scenario lookup returned a non-finite delta_sst.", call. = FALSE)
  }

  as.numeric(y)
}




# =============================================================================
# 6) NON-STATIONARY COUNT SIMULATION
# =============================================================================

#' Simulate annual counts with climate adjustments
#'
#' @description
#' Non-stationary extension of the Poisson-Gamma annual count model with
#' a rate effect and an intensity effect applied to a stationary time slice:
#'
#' **Rate effect:** Each year's activity factor is modulated by a scalar SST shift:
#'   `A_t = activity_factor * exp(beta_SST * delta_sst)`
#'
#' **Intensity effect:** The storm-class split varies with the same scalar SST shift:
#'   `p_HUR(t) = clamp(p_HUR_base * (1 + gamma * delta_sst), 0.01, 0.99)`
#'   `N_total_t ~ Poisson(lambda_total * A_t)`
#'   `n_HUR_t ~ Binomial(N_total_t, p_HUR(t))`
#'   `n_TS_t = N_total_t - n_HUR_t`
#'
#' When gamma_intensity is 0, the storm-class split is constant.
#' When both beta_sst and gamma_intensity are 0, reduces to stationary model.
#'
#' @param lambda_table Tibble from `compute_lambda_table()`.
#' @param k_hat Numeric; Gamma shape parameter for overdispersion.
#' @param n_years_sim Integer; number of years to simulate.
#' @param delta_sst Numeric scalar SST shift (degC) applied uniformly to all
#'   simulation years.
#' @param beta_sst Numeric; SST rate-effect coefficient.
#' @param gamma_intensity Numeric; SST intensity-effect coefficient.
#'   Represents fractional change in p_HUR per degC of SST warming.
#' @param p_hur_base Optional numeric; baseline hurricane fraction. If NULL,
#'   computed from lambda_table.
#' @param .sst_abs_max Numeric guardrail for absolute SST anomaly magnitude.
#' @param .sst_scale_max Numeric guardrail for \code{exp(beta_sst * delta_sst)}.
#' @param .mu_total_max Numeric guardrail for annual Poisson mean.
#'
#' @return Tibble with columns: sim_year, activity_factor, climate_scale,
#'   activity_combined, p_hurricane, n_total, n_ts, n_hur. The applied
#'   `delta_sst` is recorded as an attribute on the returned tibble.
#' @export
simulate_twolevel_counts <- function(lambda_table, k_hat, n_years_sim,
                                     delta_sst = 0,
                                     beta_sst = 0,
                                     gamma_intensity = 0,
                                     p_hur_base = NULL,
                                     .sst_abs_max = 10,          # degC plausibility
                                     .sst_scale_max = 1e3,       # exp(beta*anom) plausibility
                                     .mu_total_max = 1e6) {      # per-year Poisson mean plausibility

  stopifnot(is.data.frame(lambda_table))
  stopifnot(is.numeric(k_hat), length(k_hat) == 1)
  n_years_sim <- as.integer(n_years_sim)
  if (!is.finite(n_years_sim) || n_years_sim <= 0) stop("n_years_sim must be a positive integer.", call. = FALSE)

  .int_cap <- function(x, name) {
    # Avoid unsafe as.integer() on huge doubles; clip first.
    int_max <- .Machine$integer.max
    x2 <- as.numeric(x)

    bad_hi <- is.finite(x2) & (x2 > int_max)
    bad_lo <- is.finite(x2) & (x2 < -int_max)

    if (any(bad_hi) || any(bad_lo)) {
      warning(sprintf(
        "[simulate_twolevel_counts] %s overflow: capped %d above +%d and %d below -%d.",
        name, sum(bad_hi), int_max, sum(bad_lo), int_max
      ), call. = FALSE)
    }

    x2[bad_hi] <- int_max
    x2[bad_lo] <- -int_max
    x2[!is.finite(x2)] <- NA_real_
    as.integer(x2)
  }

  # ---- lambdas ----
  if (!all(c("storm_class", "lambda") %in% names(lambda_table))) {
    stop("lambda_table must contain columns: storm_class, lambda", call. = FALSE)
  }
  lt <- dplyr::as_tibble(lambda_table)
  lt$storm_class <- as.character(lt$storm_class)

  lam <- .extract_lambdas(lt)
  lambda_ts <- pmax(0, lam$ts)
  lambda_hur <- pmax(0, lam$hur)
  lambda_total <- lambda_ts + lambda_hur

  if (!is.numeric(delta_sst) || length(delta_sst) != 1L || !is.finite(delta_sst)) {
    stop("delta_sst must be a single finite numeric value.", call. = FALSE)
  }
  delta_sst <- as.numeric(delta_sst)
  if (abs(delta_sst) > .sst_abs_max) {
    stop("delta_sst must satisfy abs(delta_sst) <= .sst_abs_max.", call. = FALSE)
  }

  beta_sst <- as.numeric(beta_sst)
  if (!is.finite(beta_sst)) beta_sst <- 0

  # ---- Rate effect scaling with guard ----
  sst_scale <- exp(beta_sst * delta_sst)

  if (!is.finite(sst_scale)) {
    stop("Non-finite SST scaling detected (exp(beta_sst * delta_sst)). Check inputs/units.", call. = FALSE)
  }
  if (sst_scale > .sst_scale_max) {
    stop(
      sprintf(
        paste0(
          "Unrealistic SST rate scaling detected (exp(beta_sst * delta_sst)=%.3g > %.3g).\n",
          "This implies either wrong units for delta_sst or an overly large beta_sst."
        ),
        sst_scale, .sst_scale_max
      ),
      call. = FALSE
    )
  }

  # ---- annual activity factor ----
  A <- stats::rgamma(n_years_sim, shape = k_hat, rate = k_hat)
  A[!is.finite(A)] <- 1

  mu_total <- lambda_total * A * sst_scale
  if (any(mu_total > .mu_total_max, na.rm = TRUE)) {
    stop(sprintf(
      paste0(
        "Unrealistic Poisson mean detected (max mu_total = %.3g > %.3g).\n",
        "This is almost certainly due to corrupted climate scaling."
      ),
      max(mu_total, na.rm = TRUE), .mu_total_max
    ), call. = FALSE)
  }

  n_total <- stats::rpois(n_years_sim, lambda = mu_total)

  # ---- Intensity effect split (bounded) ----
  if (is.null(p_hur_base) || !is.finite(p_hur_base)) {
    p_hur_base <- if (lambda_total > 0) pmax(0, pmin(1, lambda_hur / lambda_total)) else 0.5
  }

  gamma_intensity <- as.numeric(gamma_intensity)
  if (!is.finite(gamma_intensity)) gamma_intensity <- 0

  p_hur <- p_hur_base * (1 + gamma_intensity * delta_sst)
  p_hur <- pmin(0.99, pmax(0.01, p_hur))

  n_hur <- stats::rbinom(n_years_sim, size = n_total, prob = p_hur)
  n_ts  <- n_total - n_hur

  out <- tibble::tibble(
    sim_year = seq_len(n_years_sim),
    activity_factor = A,
    climate_scale = rep(sst_scale, n_years_sim),
    activity_combined = A * sst_scale,
    p_hurricane = rep(p_hur, n_years_sim),
    n_total = .int_cap(n_total, "n_total"),
    n_ts = .int_cap(n_ts, "n_ts"),
    n_hur = .int_cap(n_hur, "n_hur")
  )
  attr(out, "delta_sst") <- delta_sst
  out
}




# =============================================================================
# 7) LEVEL 3: STORM CHARACTERISTIC PERTURBATION
# =============================================================================

#' Default storm-perturbation parameters for climate sensitivity
#'
#' @description
#' Returns a named list of per-degC scaling factors for storm-property
#' perturbation.
#'
#' @return Named list with elements v_scale, r_scale, speed_scale, precip_scale.
#' @export
default_cc_params <- function() {
  list(
    v_scale      =  0.05,   # +5% peak intensity per degC
    r_scale      =  0.08,   # +8% radii expansion per degC
    speed_scale  = -0.10,   # -10% translation speed per degC
    precip_scale =  0.07    # +7% rainfall rate per degC (Clausius-Clapeyron)
  )
}


#' Perturb sampled event properties for climate sensitivity
#'
#' @description
#' Adjusts individual storm properties in a sampled-event tibble to reflect
#' projected changes in storm structure under warmer SSTs.
#'
#' Perturbations applied (all multiplicative, proportional to delta_SST):
#'   V_peak: scaled by (1 + v_scale * delta_SST)
#'   RMW_mean_km: scaled by (1 + r_scale * delta_SST)
#'   dur_days: scaled by 1 / (1 + speed_scale * delta_SST / 2)
#'   precip_scaling: new column = 1 + precip_scale * delta_SST
#'
#' At delta_SST = 0 all factors equal 1 (identity property for validation).
#'
#' @param events Tibble of sampled events with at least V_peak, dur_days.
#' @param delta_sst Numeric scalar; SST anomaly (degC) for the simulation year.
#' @param cc_params Named list of per-degC scaling factors (NULL = defaults).
#'
#' @return The input tibble with perturbed columns plus precip_scaling and
#'   delta_sst columns.
#' @export
perturb_event <- function(events, delta_sst, cc_params = NULL) {
  # --- trivial early-returns ---
  if (nrow(events) == 0L) {
    events$precip_scaling <- numeric(0)
    events$delta_sst      <- numeric(0)
    return(events)
  }
  if (!is.finite(delta_sst) || delta_sst == 0) {
    events$precip_scaling <- 1.0
    events$delta_sst      <- 0.0
    return(events)
  }

  # --- resolve parameters (fill missing from defaults) ---
  defaults <- default_cc_params()
  if (is.null(cc_params)) {
    cc_params <- defaults
  } else {
    for (nm in names(defaults)) {
      if (is.null(cc_params[[nm]])) cc_params[[nm]] <- defaults[[nm]]
    }
  }

  v_sc  <- as.numeric(cc_params$v_scale)
  r_sc  <- as.numeric(cc_params$r_scale)
  sp_sc <- as.numeric(cc_params$speed_scale)
  pr_sc <- as.numeric(cc_params$precip_scale)

  # --- 1) Peak intensity: V_peak * (1 + v_scale * delta_SST) ---
  v_factor <- 1 + v_sc * delta_sst
  events$V_peak <- pmax(15, events$V_peak * v_factor)

  # --- 2) Wind radii / RMW: RMW * (1 + r_scale * delta_SST) ---
  r_factor <- 1 + r_sc * delta_sst
  if ("RMW_mean_km" %in% names(events)) {
    ok_rmw <- is.finite(events$RMW_mean_km)
    events$RMW_mean_km[ok_rmw] <- pmax(5, events$RMW_mean_km[ok_rmw] * r_factor)
  }

  # --- 3) Translation speed -> duration: dur / (1 + speed_scale * dSST / 2) ---
  # speed_scale < 0 => slower storms => longer exposure
  denom <- 1 + sp_sc * delta_sst / 2
  denom <- max(0.25, denom)   # guard: never compress below 25% of original
  dur_factor <- 1 / denom
  events$dur_days <- pmax(1L, as.integer(round(events$dur_days * dur_factor)))

  # --- 4) Precipitation scaling (metadata for future multi-hazard use) ---
  events$precip_scaling <- 1 + pr_sc * delta_sst

  # --- traceability ---
  events$delta_sst <- delta_sst

  events
}


# =============================================================================
# 8) CLIMATE CONFIGURATION HELPER
# =============================================================================

# Build a climate configuration object for the hazard model.
# The exported constructor is `make_climate_cfg()` below.
.normalize_climate_perturb <- function(perturb,
                                       scenario,
                                       allow_knmi = TRUE) {
  if (is.null(perturb)) {
    return(list(
      state = "disabled",
      params = NULL
    ))
  }
  if (!is.list(perturb)) {
    stop("perturb must be NULL or a named list of storm-perturbation settings.", call. = FALSE)
  }

  resolved <- perturb
  perturb_state <- if (length(perturb) == 0L) "default" else "custom"

  if (allow_knmi && grepl("^knmi_", scenario) && length(perturb) == 0L) {
    if (!exists("knmi_cc_params", mode = "function")) {
      stop("KNMI scenario selected but knmi_cc_params() is not available. Did you include hazard_climate_knmi23.R?", call. = FALSE)
    }
    resolved <- knmi_cc_params(scenario, base_params = resolved)
  }

  list(
    state = perturb_state,
    params = resolved
  )
}

.normalize_climate_cfg <- function(cfg) {
  if (!is.list(cfg)) {
    stop("climate_cfg must be a list.", call. = FALSE)
  }
  if ("enabled" %in% names(cfg)) {
    stop("Climate-off mode has been removed. Remove `enabled` and use make_climate_cfg(scenario = \"stationary\") for baseline runs.", call. = FALSE)
  }

  cfg$scenario <- match.arg(
    as.character(cfg$scenario[[1]]),
    choices = c("stationary", sst_scenario_info("all")$scenario)
  )
  cfg$sst_source <- match.arg(
    as.character(cfg$sst_source[[1]]),
    choices = c("builtin", "csv", "ersst_nc")
  )
  cfg$baseline_years <- as.integer(cfg$baseline_years)
  cfg$start_year <- as.integer(cfg$start_year[[1]])

  if (length(cfg$baseline_years) == 0L || any(!is.finite(cfg$baseline_years))) {
    stop("baseline_years must contain finite integer years.", call. = FALSE)
  }
  if (!is.finite(cfg$start_year)) {
    stop("start_year must be a single finite integer year.", call. = FALSE)
  }
  cfg$sensitivity_mode <- match.arg(
    as.character(cfg$sensitivity_mode[[1]]),
    choices = c("fixed", "linear_shifted")
  )
  cfg$k_beta <- as.numeric(cfg$k_beta[[1]])
  cfg$k_gamma <- as.numeric(cfg$k_gamma[[1]])
  if (!is.finite(cfg$k_beta)) {
    stop("k_beta must be a single finite numeric value.", call. = FALSE)
  }
  if (!is.finite(cfg$k_gamma)) {
    stop("k_gamma must be a single finite numeric value.", call. = FALSE)
  }

  perturb_info <- .normalize_climate_perturb(
    perturb = cfg[["perturb"]],
    scenario = cfg[["scenario"]]
  )
  cfg[["perturb"]] <- perturb_info$params
  cfg[["perturb_state"]] <- perturb_info$state
  cfg[["cc_params"]] <- perturb_info$params
  cfg[["scenario_start_year"]] <- cfg[["start_year"]]
  cfg
}

#' Build a climate configuration object for the hazard model
#'
#' @description
#' Creates a climate configuration object for `run_hazard_model()`.
#' The climate workflow first calibrates historical baseline sensitivities
#' (`beta_0`, `gamma_0`) from historical annual counts and MDR SST anomalies,
#' then resolves scenario-specific simulation scalars from `delta_sst`.
#' By default the historical sensitivities are used unchanged. Expert users may
#' instead enable a linear sensitivity shift:
#'
#' `beta_sst = beta_0 * (1 + k_beta * delta_sst)`
#'
#' `gamma = gamma_0 * (1 + k_gamma * delta_sst)`
#'
#' Optional storm perturbation is available as an expert sensitivity and is
#' disabled unless `perturb` is supplied explicitly.
#'
#' `scenario` accepts `"stationary"` plus any value returned by
#' `sst_scenario_info("all")$scenario` (for example SSP and optional KNMI
#' scenarios when KNMI helpers are available).
#'
#' If `scenario = "stationary"`, the returned config represents the canonical
#' baseline hazard-model specification with `delta_sst = 0`. Climate is always
#' resolved through this configuration; there is no disabled climate mode.
#'
#' @param scenario Character; climate scenario name.
#' @param sst_source Character; one of "builtin", "csv", "ersst_nc".
#' @param sst_path Optional character; path to SST data file (CSV or NetCDF).
#' @param baseline_years Integer vector; years for climatological baseline.
#' @param start_year Integer; first year of the simulation scenario.
#' @param sensitivity_mode Character scalar; one of `"fixed"` or
#'   `"linear_shifted"`.
#' @param k_beta Numeric scalar; linear sensitivity-shift coefficient applied
#'   to `beta_0` when `sensitivity_mode = "linear_shifted"`.
#' @param k_gamma Numeric scalar; linear sensitivity-shift coefficient applied
#'   to `gamma_0` when `sensitivity_mode = "linear_shifted"`.
#' @param perturb Optional storm-perturbation settings.
#'   - `NULL`: disable storm perturbation.
#'   - `list()`: enable storm perturbation with defaults from `default_cc_params()`.
#'   - named list: enable storm perturbation with user overrides.
#'
#' @return A list with class "climate_cfg" containing climate configuration parameters.
#'
#' @examples
#' cfg <- make_climate_cfg(
#'   sst_source = "builtin",
#'   scenario = "ssp245"
#' )
#' cfg
#' @export
make_climate_cfg <- function(scenario = "stationary",
                             sst_source = c("builtin", "csv", "ersst_nc"),
                             sst_path = NULL,
                             baseline_years = 1991L:2020L,
                             start_year = 2025L,
                             sensitivity_mode = c("fixed", "linear_shifted"),
                             k_beta = 0,
                             k_gamma = 0,
                             perturb = NULL) {
  sst_source <- match.arg(sst_source)
  sensitivity_mode <- match.arg(sensitivity_mode)
  scenario <- match.arg(scenario, choices = c("stationary", sst_scenario_info("all")$scenario))

  cfg <- list(
    scenario = scenario,
    sst_source = sst_source,
    sst_path = sst_path,
    baseline_years = baseline_years,
    start_year = start_year,
    sensitivity_mode = sensitivity_mode,
    k_beta = k_beta,
    k_gamma = k_gamma,
    perturb = perturb
  )
  cfg[["sst_path"]] <- sst_path
  cfg[["perturb"]] <- perturb
  cfg <- .normalize_climate_cfg(cfg)
  class(cfg) <- c("climate_cfg", "list")
  cfg
}

#' @export
print.climate_cfg <- function(x, ...) {
  cat("Climate configuration\n")
  cat(sprintf("  Scenario          : %s\n", x$scenario))
  cat(sprintf("  SST source        : %s\n", x$sst_source))
  cat(sprintf("  Baseline          : %d-%d\n", min(x$baseline_years), max(x$baseline_years)))

  if (identical(x$scenario, "stationary")) {
    cat("  Climate mode      : baseline climate run (delta_sst = 0)\n")
  } else {
    cat("  Climate mode      : future climate run\n")
  }
  cat(sprintf("  Sensitivity mode  : %s\n", x$sensitivity_mode))
  if (identical(x$sensitivity_mode, "linear_shifted")) {
    cat(sprintf("  Shift coefficients: k_beta=%.3f, k_gamma=%.3f\n", x$k_beta, x$k_gamma))
  } else {
    cat("  Shift coefficients: fixed historical sensitivities\n")
  }

  perturb_status <- switch(
    x$perturb_state,
    disabled = "disabled",
    default = "enabled (defaults)",
    custom = "enabled (custom)",
    "disabled"
  )
  cat(sprintf("  Storm perturbation: %s\n", perturb_status))
  invisible(x)
}

#' Resolve climate scalars for the hazard model
#'
#' @description
#' Calibrates historical baseline sensitivities from annual counts and MDR SST
#' anomalies, resolves scenario `delta_sst`, and returns flat simulation-ready
#' climate scalars for `run_hazard_model()`.
#'
#' @param climate_cfg List from `make_climate_cfg()`.
#' @param annual_counts Optional tibble of annual counts for baseline
#'   sensitivity estimation.
#' @param lambda_table Optional tibble from `compute_lambda_table()` for
#'   historical `p_hur_base`.
#' @param min_year Integer; passed to estimation functions.
#' @param future_period Optional integer vector of scenario years used to
#'   resolve `delta_sst`. Defaults to the 30-year window starting at
#'   `climate_cfg$start_year`.
#' @param verbose Logical.
#'
#' @return A list with:
#'   \item{delta_sst}{Scenario-derived SST shift in degC.}
#'   \item{beta_sst}{Final rate sensitivity used in simulation.}
#'   \item{gamma}{Final intensity sensitivity used in simulation.}
#'   \item{p_hur_base}{Historical baseline hurricane fraction.}
#'   \item{beta_0}{Historical baseline rate sensitivity.}
#'   \item{gamma_0}{Historical baseline intensity sensitivity.}
#'   \item{scenario}{Resolved scenario name.}
#'   \item{sensitivity_mode}{Resolved sensitivity mode.}
#'   \item{k_beta}{Configured linear shift coefficient for `beta_0`.}
#'   \item{k_gamma}{Configured linear shift coefficient for `gamma_0`.}
#'   \item{future_period}{Scenario period used to resolve `delta_sst`.}
#'   \item{sst_scale}{Implied SST-driven count multiplier
#'     `exp(beta_sst * delta_sst)`.}
#'   \item{annual_count_series}{Basin-consistent annual total count series used
#'     for `beta_0` calibration, in storms/year.}
#'   \item{annual_count_source}{Character label describing the provenance of
#'     `annual_count_series`.}
#'   \item{beta_guardrail}{List describing any `beta_0` plausibility fallback.}
#'   \item{sst_scale_guardrail}{List describing any multiplier guardrail applied
#'     to `beta_sst` for the resolved scenario.}
#'   \item{baseline_years}{SST anomaly baseline years.}
#'   \item{perturb}{Resolved storm-perturbation parameters (or `NULL`).}
#'   \item{perturb_state}{Storm-perturbation state label.}
#'   \item{source}{SST source used for calibration.}
#'
#' @examples
#' cfg <- make_climate_cfg(sst_source = "builtin", scenario = "stationary")
#' climate <- resolve_climate_inputs(cfg, verbose = FALSE)
#' climate$beta_sst
#'
#' @export
resolve_climate_inputs <- function(climate_cfg,
                                   annual_counts = NULL,
                                   lambda_table = NULL,
                                   min_year = 1970L,
                                   future_period = NULL,
                                   verbose = TRUE) {
  sst_scale_guardrail_max <- 4

  .clamp_nonnegative <- function(x) {
    if (!is.finite(x) || x < 0) {
      return(0)
    }
    as.numeric(x)
  }

  .resolve_effective_sensitivity <- function(base_value, shift_coeff, delta_sst) {
    shifted <- as.numeric(base_value) * (1 + as.numeric(shift_coeff) * as.numeric(delta_sst))
    .clamp_nonnegative(shifted)
  }

  .clamp_scalar <- function(x, lower, upper) {
    min(max(as.numeric(x), as.numeric(lower)), as.numeric(upper))
  }

  .resolve_rate_response_regime <- function(future_period, raw_scale) {
    midpoint_year <- as.numeric(mean(as.integer(future_period)))
    if (!is.finite(midpoint_year)) {
      stop("future_period must resolve to a finite midpoint year.", call. = FALSE)
    }

    if (midpoint_year <= 2055) {
      regime <- "near_term"
      regime_weight <- 0
    } else if (midpoint_year >= 2085) {
      regime <- "late_century"
      regime_weight <- 1
    } else {
      regime <- "transition"
      regime_weight <- (midpoint_year - 2055) / (2085 - 2055)
    }

    damping <- 0.08 + regime_weight * (0.35 - 0.08)
    rate_bounds <- c(
      0.95 + regime_weight * (0.85 - 0.95),
      1.08 + regime_weight * (1.30 - 1.08)
    )
    redistribution_strength <- 0.06 + regime_weight * (0.18 - 0.06)
    redistribution_bounds <- c(
      0.90 + regime_weight * (0.75 - 0.90),
      1.10 + regime_weight * (1.25 - 1.10)
    )

    scaled <- 1 + damping * (as.numeric(raw_scale) - 1)
    bounded <- .clamp_scalar(scaled, rate_bounds[[1]], rate_bounds[[2]])

    list(
      regime = regime,
      midpoint_year = midpoint_year,
      damping = damping,
      basin_rate_bounds = rate_bounds,
      redistribution_strength = redistribution_strength,
      redistribution_bounds = redistribution_bounds,
      raw_scale = as.numeric(raw_scale),
      adjusted_scale = as.numeric(bounded),
      guardrail = list(
        triggered = !isTRUE(all.equal(as.numeric(raw_scale), as.numeric(bounded), tolerance = 1e-12)),
        reason = if (as.numeric(raw_scale) > rate_bounds[[2]]) {
          "count_inflation_bounded"
        } else if (as.numeric(raw_scale) < rate_bounds[[1]]) {
          "count_decline_bounded"
        } else if (!isTRUE(all.equal(as.numeric(raw_scale), as.numeric(bounded), tolerance = 1e-12))) {
          "rate_response_damped"
        } else {
          "none"
        }
      )
    )
  }

  baseline <- .calibrate_climate_baseline(
    climate_cfg = climate_cfg,
    annual_counts = annual_counts,
    lambda_table = lambda_table,
    min_year = min_year,
    verbose = verbose
  )

  if (is.null(future_period)) {
    future_period <- seq.int(from = as.integer(baseline$start_year), length.out = 30L)
  } else {
    future_period <- as.integer(future_period)
  }
  if (length(future_period) == 0L || any(!is.finite(future_period))) {
    stop("future_period must contain finite integer years.", call. = FALSE)
  }

  delta_sst <- 0
  beta_sst <- baseline$beta_0
  gamma <- baseline$gamma_0
  climate_mode <- "baseline"
  raw_sst_scale <- 1
  f_rate_climate <- 1
  response_regime <- list(
    regime = "baseline",
    midpoint_year = NA_real_,
    damping = 0,
    basin_rate_bounds = c(1, 1),
    redistribution_strength = 0,
    redistribution_bounds = c(1, 1),
    raw_scale = 1,
    adjusted_scale = 1,
    guardrail = list(
      triggered = FALSE,
      reason = "none"
    )
  )

  if (!identical(baseline$scenario, "stationary")) {
    delta_sst <- get_scenario_delta(
      scenario = baseline$scenario,
      future_period = future_period,
      baseline_years = baseline$baseline_years
    )
    climate_mode <- "future"
  }
  if (identical(baseline$sensitivity_mode, "linear_shifted")) {
    beta_sst <- .resolve_effective_sensitivity(baseline$beta_0, baseline$k_beta, delta_sst)
    gamma <- .resolve_effective_sensitivity(baseline$gamma_0, baseline$k_gamma, delta_sst)
  }
  raw_sst_scale <- exp(beta_sst * delta_sst)
  response_regime <- if (identical(climate_mode, "future")) {
    .resolve_rate_response_regime(future_period = future_period, raw_scale = raw_sst_scale)
  } else {
    response_regime
  }
  f_rate_climate <- response_regime$adjusted_scale
  sst_scale <- f_rate_climate
  sst_scale_guardrail <- list(
    triggered = FALSE,
    reason = "none",
    sst_scale_max = sst_scale_guardrail_max
  )
  if (is.finite(delta_sst) && delta_sst > 0 && is.finite(raw_sst_scale) && raw_sst_scale > sst_scale_guardrail_max) {
    beta_sst <- log(sst_scale_guardrail_max) / delta_sst
    raw_sst_scale <- exp(beta_sst * delta_sst)
    response_regime <- if (identical(climate_mode, "future")) {
      .resolve_rate_response_regime(future_period = future_period, raw_scale = raw_sst_scale)
    } else {
      response_regime
    }
    f_rate_climate <- response_regime$adjusted_scale
    sst_scale <- f_rate_climate
    sst_scale_guardrail$triggered <- TRUE
    sst_scale_guardrail$reason <- "sst_scale_above_plausibility_limit"
    if (verbose) {
      message(sprintf(
        "[climate] Guardrail: raw SST count multiplier exceeded %.2fx; beta_sst reduced to %.3f 1/degC.",
        sst_scale_guardrail_max, beta_sst
      ))
    }
  }
  if (verbose && identical(climate_mode, "future") && isTRUE(response_regime$guardrail$triggered)) {
    message(sprintf(
      "[climate] Climate rate regime (%s): raw %.3fx -> applied %.3fx within [%.2f, %.2f].",
      response_regime$regime,
      response_regime$raw_scale,
      response_regime$adjusted_scale,
      response_regime$basin_rate_bounds[[1]],
      response_regime$basin_rate_bounds[[2]]
    ))
  }

  list(
    scenario = baseline$scenario,
    source = baseline$source,
    delta_sst = as.numeric(delta_sst),
    beta_sst = as.numeric(beta_sst),
    beta_sst_effective = if (isTRUE(delta_sst > 0)) log(as.numeric(f_rate_climate)) / as.numeric(delta_sst) else as.numeric(beta_sst),
    gamma = as.numeric(gamma),
    p_hur_base = as.numeric(baseline$p_hur_base),
    beta_0 = as.numeric(baseline$beta_0),
    gamma_0 = as.numeric(baseline$gamma_0),
    sensitivity_mode = baseline$sensitivity_mode,
    k_beta = as.numeric(baseline$k_beta),
    k_gamma = as.numeric(baseline$k_gamma),
    future_period = future_period,
    sst_scale = as.numeric(sst_scale),
    raw_sst_scale = as.numeric(raw_sst_scale),
    f_rate_climate = as.numeric(f_rate_climate),
    annual_count_series = baseline$annual_count_series,
    annual_count_source = baseline$annual_count_source,
    beta_guardrail = baseline$beta_guardrail,
    sst_scale_guardrail = sst_scale_guardrail,
    baseline_years = baseline$baseline_years,
    perturb = baseline$perturb,
    perturb_state = baseline$perturb_state,
    climate_mode = climate_mode,
    response_regime = response_regime
  )
}

.calibrate_climate_baseline <- function(climate_cfg,
                                        annual_counts = NULL,
                                        lambda_table = NULL,
                                        min_year = 1970L,
                                        verbose = TRUE) {
  beta_prior <- 0.6
  gamma_prior <- 0.065

  if (!inherits(climate_cfg, "climate_cfg")) {
    stop("climate_cfg must be created by make_climate_cfg().")
  }
  climate_cfg <- .normalize_climate_cfg(climate_cfg)

  # Load SST data
  sst_raw <- switch(climate_cfg$sst_source,
                    builtin   = get_mdr_sst_builtin(),
                    csv       = read_mdr_sst_csv(climate_cfg$sst_path),
                    ersst_nc  = read_mdr_sst_ersst(climate_cfg$sst_path)
  )

  if (verbose) {
    message(sprintf("[climate] Loaded %d years of MDR SST (%d-%d) from %s",
                    nrow(sst_raw), min(sst_raw$year), max(sst_raw$year), climate_cfg$sst_source))
  }

  sst_df <- compute_sst_anomaly(sst_raw, baseline_years = climate_cfg$baseline_years)

  if (verbose) {
    message(sprintf("[climate] Baseline (%.0f-%.0f): %.2f C | Anomaly range: [%+.2f, %+.2f] C",
                    min(climate_cfg$baseline_years), max(climate_cfg$baseline_years),
                    sst_df$sst_clim[1],
                    min(sst_df$sst_anomaly, na.rm = TRUE),
                    max(sst_df$sst_anomaly, na.rm = TRUE)))
  }

  # Historical baseline rate sensitivity
  annual_count_series <- NULL
  annual_count_source <- "literature_fallback"
  beta_guardrail <- list(
    triggered = FALSE,
    reason = "no_annual_counts",
    beta_max = 1.2,
    beta_fallback = beta_prior
  )
  if (!is.null(annual_counts)) {
    beta_info <- estimate_beta_sst(
      annual_counts = annual_counts,
      sst_df = sst_df,
      min_year = min_year,
      beta_prior = beta_prior,
      verbose = verbose
    )
    beta_0 <- as.numeric(beta_info$beta_sst)
    annual_count_series <- beta_info$annual_count_series
    annual_count_source <- beta_info$annual_count_source
    beta_guardrail <- beta_info$guardrail
  } else {
    beta_0 <- beta_prior
    if (verbose) message(sprintf("[climate] No annual_counts provided for beta estimation. Using prior: %.3f", beta_0))
  }
  if (!is.finite(beta_0)) {
    stop("Resolved beta_0 must be finite.", call. = FALSE)
  }

  # Historical baseline intensity sensitivity
  p_hur_base <- NA_real_

  if (!is.null(lambda_table)) {
    p_hur_base <- compute_p_hur_base(lambda_table)
  }

  if (!is.null(annual_counts)) {
    gamma_info <- estimate_gamma_intensity(
      annual_counts = annual_counts,
      sst_df = sst_df,
      min_year = min_year,
      gamma_prior = gamma_prior,
      verbose = verbose
    )
    gamma_0 <- as.numeric(gamma_info$gamma)
    if (is.na(p_hur_base) && !is.null(gamma_info$p_hur_base) && is.finite(gamma_info$p_hur_base)) {
      p_hur_base <- as.numeric(gamma_info$p_hur_base)
    }
  } else {
    gamma_0 <- 0
    if (verbose) message("[climate] No annual_counts provided for gamma estimation. Using gamma_0 = 0.")
  }

  if (!is.finite(gamma_0) || gamma_0 <= 0) {
    if (verbose) {
      message(sprintf("[climate] Intensity effect resolved to %.4f; using gamma_0 = 0.", gamma_0))
    }
    gamma_0 <- 0
  }

  if (!is.finite(gamma_0)) {
    stop("Resolved gamma_0 must be finite.", call. = FALSE)
  }

  if (verbose && gamma_0 != 0) {
    message(sprintf("[climate] Historical intensity sensitivity: gamma_0=%.4f, p_HUR_base=%.3f",
                    gamma_0, p_hur_base))
    message(sprintf("[climate] At +1C SST: p_HUR -> %.3f (was %.3f, %+.1f%%)",
                    pmin(0.99, p_hur_base * (1 + gamma_0)),
                    p_hur_base,
                    100 * gamma_0))
  }

  resolved_perturb <- climate_cfg[["perturb"]]
  perturb_state <- climate_cfg[["perturb_state"]]
  if (!is.null(resolved_perturb) && length(resolved_perturb) == 0L) {
    resolved_perturb <- default_cc_params()
  }
  if (!is.null(resolved_perturb)) {
    defs <- default_cc_params()
    for (nm in names(defs)) {
      if (is.null(resolved_perturb[[nm]])) resolved_perturb[[nm]] <- defs[[nm]]
    }
    perturb_state <- if (identical(climate_cfg[["perturb_state"]], "disabled")) "default" else climate_cfg[["perturb_state"]]
    if (verbose) {
      message(sprintf("[climate] Storm perturbation enabled: v_scale=%+.2f, r_scale=%+.2f, speed_scale=%+.2f, precip_scale=%+.2f per degC",
                      resolved_perturb$v_scale, resolved_perturb$r_scale,
                      resolved_perturb$speed_scale, resolved_perturb$precip_scale))
    }
  } else {
    perturb_state <- "disabled"
    if (verbose) message("[climate] Storm perturbation disabled.")
  }

  list(
    scenario = climate_cfg$scenario,
    source = climate_cfg$sst_source,
    baseline_years = climate_cfg$baseline_years,
    start_year = climate_cfg$start_year,
    sensitivity_mode = climate_cfg$sensitivity_mode,
    k_beta = climate_cfg$k_beta,
    k_gamma = climate_cfg$k_gamma,
    beta_0 = beta_0,
    gamma_0 = gamma_0,
    p_hur_base = p_hur_base,
    annual_count_series = annual_count_series,
    annual_count_source = annual_count_source,
    beta_guardrail = beta_guardrail,
    perturb = resolved_perturb,
    perturb_state = perturb_state
  )
}



################################################################################
# KNMI'23 Climate Scenario Support for Dutch Caribbean Hurricane Hazard Model
################################################################################

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
