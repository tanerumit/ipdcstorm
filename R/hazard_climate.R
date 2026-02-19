################################################################################
# hazard_climate.R
# Climate Change Modifications (Levels 1-3)
#
# Level 1: SST-Conditioned Rate Scaling (beta_SST)
# Level 2: Intensity Distribution Shift (gamma)
# Level 3: Storm Characteristic Perturbation (V_peak, RMW, duration, precip)
#
# Modulates the annual activity factor using observed/projected SST anomalies
# in the Main Development Region (MDR: 10â€“20Â°N, 80â€“20Â°W):
#
#   A_t ~ Gamma(k, k) Ã— exp(Î²_SST Â· Î”SST_t)
#
# Components:
#   1) MDR SST data ingestion (ERSST v5 CSV/NetCDF + built-in fallback)
#   2) Anomaly computation relative to 1991â€“2020 climatological mean
#   3) Î²_SST estimation via Poisson/NegBin regression on historical record
#   4) SST scenario generation (historical + CMIP6 SSP projections)
#   5) Non-stationary count simulation
#   6) Intensity distribution shift (Level 2: gamma estimation)
#   7) Storm characteristic perturbation (Level 3: perturb_event())
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
#' Returns a tibble of annual-mean SST (Â°C) averaged over the Main Development
#' Region (MDR: 10â€“20Â°N, 80â€“20Â°W) for the Atlantic hurricane season (Augâ€“Oct,
#' the peak months most predictive of activity).
#'
#' Values are derived from NOAA ERSST v5 monthly data, spatially averaged over
#' the MDR box, then temporally averaged over ASO (Aug-Sep-Oct) each year.
#'
#' These are provided as a built-in fallback so the model can run without
#' requiring users to download and process NetCDF files.
#'
#' Source: NOAA/NCEI ERSST v5, accessed 2024.
#' Reference: Huang et al. (2017), J. Climate, 30, 8179â€“8205.
#'
#' @return Tibble with columns: year (integer), sst_mdr_aso (numeric, Â°C).
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
#' Reads monthly ERSST v5 NetCDF data, subsets to the MDR box (10â€“20Â°N, 80â€“20Â°W),
#' averages spatially, then computes ASO (Augâ€“Sepâ€“Oct) seasonal means per year.
#'
#' Requires the `ncdf4` package. If not available, falls back to the built-in
#' reference data with a warning.
#'
#' @param nc_path Character; path to ERSST v5 NetCDF file (e.g., "sst.mnmean.nc").
#' @param mdr_lat Range of latitudes for MDR (default: c(10, 20)).
#' @param mdr_lon Range of longitudes for MDR (default: c(-80, -20)).
#'   Note: ERSST uses 0â€“360 longitude convention; this function converts internally.
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

  # Time â†’ dates (ERSST time is days since 1800-01-01)
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
#' Computes Î”SST_t = SST_t âˆ’ SST_clim, where SST_clim is the mean SST over
#' the specified baseline period (default: 1991â€“2020, the current WMO standard
#' climatological normal).
#'
#' @param sst_df Tibble with columns: year, sst_mdr_aso.
#' @param baseline_years Integer vector of years defining the climatological
#'   baseline (default: 1991:2020).
#'
#' @return The input tibble with added columns:
#'   \item{sst_clim}{Climatological mean SST (Â°C) over the baseline period.}
#'   \item{sst_anomaly}{Î”SST (Â°C) relative to baseline.}
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
# 4) Î²_SST ESTIMATION
# =============================================================================

#' Estimate the SSTâ€“activity scaling coefficient Î²_SST
#'
#' @description
#' Fits a Poisson (or negative binomial) GLM of annual TC counts on MDR SST
#' anomaly to estimate Î²_SST in:
#'
#'   `E[N_t] = exp(Î± + Î²_SST Â· Î”SST_t)`
#'
#' The coefficient Î²_SST represents the log-linear sensitivity of annual TC
#' activity to SST anomalies. Typical values from the literature are 0.4â€“0.8
#' per Â°C for the North Atlantic (Villarini et al. 2011; Vecchi et al. 2021).
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
#'   Î²_final = wÂ·Î²_MLE + (1-w)Â·Î²_prior, where w = min(1, n_years/50).
#'   This stabilizes estimates when the historical record is short.
#' @param verbose Logical; print diagnostic output.
#'
#' @return A list with:
#'   \item{beta_sst}{Estimated (or shrunk) Î²_SST coefficient.}
#'   \item{beta_se}{Standard error of the MLE Î²_SST.}
#'   \item{beta_mle}{Raw MLE estimate before shrinkage.}
#'   \item{alpha}{Intercept (log baseline rate).}
#'   \item{method}{Character; "negbin", "quasipoisson", or "literature_fallback".}
#'   \item{n_years}{Number of years used in estimation.}
#'   \item{r_squared_dev}{Deviance-based pseudo-RÂ² (proportion of deviance explained).}
#'   \item{aic}{AIC of the fitted model (NA for quasi-Poisson).}
#'   \item{fit_data}{Tibble of the joined data used for fitting.}
#'   \item{glm_fit}{The fitted GLM object.}
#'
#' @export
estimate_beta_sst <- function(annual_counts,
                              sst_df,
                              min_year = 1970L,
                              beta_prior = NULL,
                              verbose = TRUE) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` is required.")
  if (!requireNamespace("stats", quietly = TRUE)) stop("Package `stats` is required.")

  # Aggregate to total annual counts (across severities)
  totals <- annual_counts |>
    dplyr::filter(.data$year >= min_year) |>
    dplyr::group_by(.data$year) |>
    dplyr::summarise(N = sum(.data$n_events), .groups = "drop")

  # Require sst_anomaly column
  if (!("sst_anomaly" %in% names(sst_df))) {
    stop("sst_df must contain 'sst_anomaly'. Run compute_sst_anomaly() first.")
  }

  # Join
  fit_data <- totals |>
    dplyr::inner_join(
      sst_df |> dplyr::select("year", "sst_anomaly"),
      by = "year"
    ) |>
    dplyr::filter(is.finite(.data$N), is.finite(.data$sst_anomaly))

  n_years <- nrow(fit_data)

  if (n_years < 15) {
    warning("Only ", n_years, " years of overlapping data. ",
            "Î²_SST estimate will be unreliable. Consider using beta_prior.")
  }

  if (n_years < 5) {
    message("[Î²_SST] Insufficient data (", n_years, " years). Using literature fallback: Î² = 0.6")
    return(list(
      beta_sst = 0.6,
      beta_se = NA_real_,
      beta_mle = NA_real_,
      alpha = NA_real_,
      method = "literature_fallback",
      n_years = n_years,
      r_squared_dev = NA_real_,
      aic = NA_real_,
      fit_data = fit_data,
      glm_fit = NULL
    ))
  }

  # Fit GLM: N ~ exp(Î± + Î²Â·Î”SST)
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

  # Deviance pseudo-RÂ²
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
    w <- min(1.0, n_years / 50)
    beta_final <- w * beta_mle + (1 - w) * beta_prior
    if (verbose) {
      message(sprintf("[Î²_SST] Shrinkage: Î²_MLE=%.3f, Î²_prior=%.3f, w=%.2f â†’ Î²_final=%.3f",
                      beta_mle, beta_prior, w, beta_final))
    }
  }

  if (verbose) {
    message(sprintf("[Î²_SST] Method: %s | n_years: %d | Î²_SST: %.3f (SE: %.3f) | pseudo-RÂ²: %.3f",
                    method, n_years, beta_final, se, r2_dev))
    message(sprintf("[Î²_SST] Interpretation: +1Â°C SST â†’ %.0f%% change in annual rate",
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
    fit_data = fit_data,
    glm_fit = glm_fit
  )
}


# =============================================================================
# 4b) INTENSITY DISTRIBUTION SHIFT (Level 2 Climate Modification)
# =============================================================================

#' Estimate the historical hurricane fraction from annual counts
#'
#' @description
#' Computes p_HUR_base = Î»_HUR / (Î»_TS + Î»_HUR) from the historical record.
#' This is the baseline probability that an event reaching TS intensity or
#' above becomes a hurricane (â‰¥64 kt).
#'
#' @param lambda_table Tibble from `compute_lambda_table()` with severities
#'   "TS" and "HUR64plus".
#'
#' @return Numeric scalar: baseline hurricane fraction p_HUR_base.
#' @export
compute_p_hur_base <- function(lambda_table) {
  lam_TS  <- lambda_table$lambda[lambda_table$storm_class == "TS"]
  lam_HUR <- lambda_table$lambda[lambda_table$storm_class == "HUR64plus"]
  if (length(lam_TS) == 0)  lam_TS  <- 0
  if (length(lam_HUR) == 0) lam_HUR <- 0
  lam_total <- lam_TS + lam_HUR
  if (lam_total <= 0) return(0.5)  # safeguard
  as.numeric(lam_HUR / lam_total)
}


#' Estimate Î³ (intensity shift coefficient) from historical data
#'
#' @description
#' Fits a logistic regression of annual hurricane fraction on SST anomaly
#' to estimate the intensification trend coefficient Î³ in:
#'
#'   p_HUR(t) = p_HUR_base Ã— (1 + Î³ Â· Î”SST_t)
#'
#' The coefficient Î³ captures how the fraction of storms reaching hurricane
#' intensity shifts with SST warming. Literature estimates suggest 5â€“8% increase
#' in Cat 4â€“5 fraction per Â°C of tropical SST warming (Knutson et al. 2020;
#' Kossin et al. 2020).
#'
#' Uses a binomial GLM: cbind(n_HUR, n_TS) ~ sst_anomaly, then converts the
#' logistic coefficient to the linear Î³ parameterization.
#'
#' @param annual_counts Tibble with columns: year, storm_class, n_events.
#' @param sst_df Tibble with columns: year, sst_anomaly.
#' @param min_year Integer; earliest year to include.
#' @param gamma_prior Optional numeric; prior value for shrinkage (default: 0.065,
#'   i.e. 6.5% increase in HUR fraction per Â°C).
#' @param verbose Logical.
#'
#' @return A list with:
#'   \item{gamma}{Estimated (or shrunk) Î³ coefficient.}
#'   \item{gamma_mle}{Raw MLE from logistic regression (converted to linear scale).}
#'   \item{gamma_se}{Approximate standard error on Î³.}
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
  if (!("HUR64plus" %in% names(wide))) wide$HUR64plus <- 0L

  # Join with SST
  fit_data <- wide |>
    dplyr::inner_join(
      sst_df |> dplyr::select("year", "sst_anomaly"),
      by = "year"
    ) |>
    dplyr::filter(is.finite(.data$sst_anomaly)) |>
    dplyr::mutate(
      n_total = .data$TS + .data$HUR64plus,
      p_hur = ifelse(.data$n_total > 0, .data$HUR64plus / .data$n_total, NA_real_)
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
  p_hur_base <- sum(fit_data$HUR64plus) / sum(fit_data$n_total)

  # Binomial GLM: logit(p_HUR) ~ sst_anomaly
  glm_fit <- tryCatch({
    stats::glm(
      cbind(HUR64plus, TS) ~ sst_anomaly,
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

  # Extract logistic coefficient and convert to linear Î³
  # In logistic model: logit(p) = Î± + Î²_logistic Â· Î”SST
  # At the baseline: p_base = logistic(Î±)
  # Linear approximation: dp/dSST â‰ˆ p_base(1-p_base) Â· Î²_logistic
  # So: Î³ = dp/dSST / p_base = (1-p_base) Â· Î²_logistic
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

  if (verbose) {
    message(sprintf("[gamma] Method: binomial_glm | n_years: %d | gamma_MLE: %.4f (SE: %.4f)",
                    n_years, gamma_mle, gamma_se))
    message(sprintf("[gamma] Shrinkage: w=%.2f -> gamma_final: %.4f", w, gamma_final))
    message(sprintf("[gamma] p_HUR_base: %.3f | +1C SST -> p_HUR changes by %+.1f%%",
                    p_hur_base, 100 * gamma_final))
  }

  list(
    gamma = as.numeric(gamma_final),
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
#'   p_HUR(t) = clamp(p_HUR_base Ã— (1 + Î³ Â· Î”SST_t), 0.01, 0.99)
#'   Î»_HUR(t) = Î»_total Ã— p_HUR(t)
#'   Î»_TS(t)  = Î»_total Ã— (1 âˆ’ p_HUR(t))
#'
#' @param lambda_table Tibble from `compute_lambda_table()`.
#' @param sst_anomaly Numeric vector of Î”SST values per simulation year.
#' @param gamma Numeric; intensity shift coefficient.
#' @param p_hur_base Optional numeric; if NULL, computed from lambda_table.
#'
#' @return Tibble with columns: sim_year, sst_anomaly, p_hur, lam_TS, lam_HUR.
#' @export
compute_severity_split <- function(lambda_table,
                                   sst_anomaly,
                                   gamma = 0,
                                   p_hur_base = NULL) {
  lam_TS  <- lambda_table$lambda[lambda_table$storm_class == "TS"]
  lam_HUR <- lambda_table$lambda[lambda_table$storm_class == "HUR64plus"]
  if (length(lam_TS) == 0)  lam_TS  <- 0
  if (length(lam_HUR) == 0) lam_HUR <- 0
  lam_total <- lam_TS + lam_HUR

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
# 5) SST SCENARIO GENERATION
# =============================================================================

#' Generate SST anomaly scenarios for simulation
#'
#' @description
#' Generates a per-year SST anomaly (degC relative to baseline_years).
#' For SSP scenarios, uses piecewise-linear warming to 2050 and 2100 targets
#' and holds constant thereafter to avoid unrealistic multi-century extrapolation.
#'
#' @param n_years Integer; number of simulation years.
#' @param mode Character; one of "stationary", "trend", "ssp245", "ssp585",
#'   or "historical_resample".
#' @param start_year Integer; first calendar year of the simulation (for labelling).
#' @param delta_start Numeric; starting Î”SST (Â°C) for "trend" mode.
#' @param delta_end Numeric; ending Î”SST (Â°C) for "trend" mode.
#' @param sst_hist_df Optional tibble with historical SST anomalies (year, sst_anomaly)
#'   for "historical_resample" mode. If NULL, uses built-in data.
#' @param baseline_years Baseline period for anomaly computation (default 1991:2020).
#'
#' @return Tibble with columns: sim_year, calendar_year, sst_anomaly, scenario.
#' @export
generate_sst_scenario <- function(n_years,
                                  mode = c("stationary", "trend",
                                           "ssp126", "ssp245", "ssp585",
                                           "historical_resample"),
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
  
  sst_anom <- switch(
    mode,
    
    stationary = rep(as.numeric(delta_start), n_years),
    
    trend = {
      seq(as.numeric(delta_start), as.numeric(delta_end), length.out = n_years)
    },
    
    ssp126 = {
      # Targets relative to 1991-2020:
      # ~+0.3C by 2050, ~+0.4C by 2100 (low-emissions pathway)
      a0 <- .get_start_anom(start_year)
      y0 <- start_year
      y1 <- 2050L
      y2 <- 2100L
      a1 <- 0.3
      a2 <- 0.4
      
      a_t <- ifelse(
        cal_year <= y1,
        .lerp(cal_year, y0, y1, a0, a1),
        ifelse(cal_year <= y2,
               .lerp(cal_year, y1, y2, a1, a2),
               a2)
      )
      as.numeric(a_t)
    },
    
    ssp245 = {
      # Targets relative to 1991â€“2020 (typical AR6 order of magnitude):
      # ~+0.5C by 2050, ~+1.0C by 2100 (median-ish)
      a0 <- .get_start_anom(start_year)
      y0 <- start_year
      y1 <- 2050L
      y2 <- 2100L
      a1 <- 0.5
      a2 <- 1.0
      
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
      # Targets relative to 1991â€“2020:
      # ~+1.0C by 2050, ~+2.5C by 2100 (median-ish)
      a0 <- .get_start_anom(start_year)
      y0 <- start_year
      y1 <- 2050L
      y2 <- 2100L
      a1 <- 1.0
      a2 <- 2.5
      
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
# 6) NON-STATIONARY COUNT SIMULATION
# =============================================================================

#' Simulate two-level annual counts with climate modifications
#'
#' @description
#' Non-stationary extension of the Poisson-Gamma annual count model with
#' two levels of climate modification:
#'
#' **Level 1 (Rate scaling):** Each year's activity factor is modulated by SST:
#'   A_t = activity_factor Ã— exp(Î²_SST Â· Î”SST_t)
#'
#' **Level 2 (Intensity shift):** The storm-class split varies with SST:
#'   p_HUR(t) = clamp(p_HUR_base Ã— (1 + Î³ Â· Î”SST_t), 0.01, 0.99)
#'   N_total_t ~ Poisson(Î»_total Ã— A_t)
#'   n_HUR_t ~ Binomial(N_total_t, p_HUR(t))
#'   n_ts_t = N_total_t âˆ’ n_hur_t
#'
#' When gamma_intensity is 0, class split is constant (Level 1 only).
#' When both beta_sst and gamma_intensity are 0, reduces to stationary model.
#'
#' @param lambda_table Tibble from `compute_lambda_table()`.
#' @param k_hat Numeric; Gamma shape parameter for overdispersion.
#' @param n_years_sim Integer; number of years to simulate.
#' @param sst_anomaly Optional numeric vector of Î”SST values (Â°C) per year.
#' @param beta_sst Numeric; SST rate scaling coefficient (Level 1).
#' @param gamma_intensity Numeric; intensity shift coefficient (Level 2).
#'   Represents fractional change in p_HUR per Â°C of SST warming.
#' @param p_hur_base Optional numeric; baseline hurricane fraction. If NULL,
#'   computed from lambda_table.
#'
#' @return Tibble with columns: sim_year, activity_factor, sst_anomaly,
#'   climate_scale, activity_combined, p_hurricane, n_total, n_ts, n_hur.
#' @export
simulate_twolevel_counts <- function(lambda_table, k_hat, n_years_sim,
                                     sst_anomaly = NULL,
                                     beta_sst = 0,
                                     gamma_intensity = 0,
                                     p_hur_base = NULL,
                                     .sst_abs_max = 10,          # degC plausibility
                                     .sst_scale_max = 1e3,       # exp(beta*anom) plausibility
                                     .mu_total_max = 1e6) {      # per-year Poisson mean plausibility
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` is required.")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Package `tibble` is required.")
  
  stopifnot(is.data.frame(lambda_table))
  stopifnot(is.numeric(k_hat), length(k_hat) == 1)
  n_years_sim <- as.integer(n_years_sim)
  if (!is.finite(n_years_sim) || n_years_sim <= 0) stop("n_years_sim must be a positive integer.", call. = FALSE)
  
  # ---- helpers ----
  .is_dateish <- function(x) inherits(x, c("Date", "POSIXct", "POSIXt"))
  
  .coerce_sst_vec <- function(x, n) {
    # Allow: NULL, numeric vector, Date/POSIX vector (but treated as corruption), or data.frame with sst_anomaly column.
    if (is.null(x)) return(rep(0, n))
    
    if (is.data.frame(x)) {
      nm <- names(x)
      if (!("sst_anomaly" %in% nm)) {
        stop(
          "sst_anomaly was passed as a data.frame/tibble but has no `sst_anomaly` column.\n",
          "Fix caller: pass sst_scenario$sst_anomaly (numeric) OR a tibble with column `sst_anomaly`.",
          call. = FALSE
        )
      }
      x <- x[["sst_anomaly"]]
    }
    
    if (.is_dateish(x)) {
      # Do not silently convert dates to numeric; that always means the wrong thing.
      stop(
        "sst_anomaly is Date/POSIX. You passed a time column, not anomalies.\n",
        "Fix caller: pass numeric Î”SST in Â°C (typically within about Â±0â€“3).",
        call. = FALSE
      )
    }
    
    x <- suppressWarnings(as.numeric(x))
    if (length(x) != n) {
      stop(sprintf("sst_anomaly length (%d) must equal n_years_sim (%d).", length(x), n), call. = FALSE)
    }
    x
  }
  
  .diagnose_corruption <- function(x) {
    fin <- is.finite(x)
    if (!any(fin)) return("sst_anomaly has no finite values.")
    rng <- range(x[fin])
    mx  <- max(abs(x[fin]))
    
    msg <- sprintf("max |sst_anomaly| = %.2f; range = [%.2f, %.2f]", mx, rng[1], rng[2])
    
    # Heuristics for the common wrong-columns
    if (mx > 1500) {
      return(paste0(
        msg, "\n",
        "This looks like a Date numeric (days since 1970) or other timestamp-derived number."
      ))
    }
    if (mx > 500 && mx < 5000) {
      return(paste0(
        msg, "\n",
        "This looks like a calendar year vector (e.g., 2025..3024) or similar."
      ))
    }
    if (mx > 20 && mx < 300) {
      return(paste0(
        msg, "\n",
        "This looks like a 'years since baseline' index (e.g., year-2005) or another trend index, not Â°C anomalies."
      ))
    }
    paste0(msg, "\nLikely wrong units/column (absolute SST, Kelvin, Fahrenheit, etc.).")
  }
  
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
  
  lambda_ts  <- lt$lambda[match("TS", lt$storm_class)]
  lambda_hur <- lt$lambda[match("HUR64plus", lt$storm_class)]
  if (!is.finite(lambda_ts))  lambda_ts <- 0
  if (!is.finite(lambda_hur)) lambda_hur <- 0
  lambda_total <- pmax(0, lambda_ts) + pmax(0, lambda_hur)
  
  # ---- SST anomalies (validated) ----
  sst_vec <- .coerce_sst_vec(sst_anomaly, n_years_sim)
  
  fin <- is.finite(sst_vec)
  if (!any(fin)) stop("sst_anomaly has no finite values.", call. = FALSE)
  
  abs_max <- max(abs(sst_vec[fin]))
  if (is.finite(abs_max) && abs_max > .sst_abs_max) {
    stop(
      "SST anomaly looks corrupt. Expected Î”SST in Â°C (typically within Â±0â€“3; hard guard Â±10).\n",
      .diagnose_corruption(sst_vec), "\n",
      "Fix caller: pass `sst_scenario$sst_anomaly` (NOT calendar_year/sim_year/years_from_2005).",
      call. = FALSE
    )
  }
  
  beta_sst <- as.numeric(beta_sst)
  if (!is.finite(beta_sst)) beta_sst <- 0
  
  # ---- L1 scaling with guard ----
  lin <- beta_sst * sst_vec
  sst_scale <- exp(lin)
  
  if (any(!is.finite(sst_scale))) {
    stop("Non-finite SST scaling detected (exp(beta_sst * sst_anomaly)). Check inputs/units.", call. = FALSE)
  }
  if (max(sst_scale, na.rm = TRUE) > .sst_scale_max) {
    stop(sprintf(
      paste0(
        "Unrealistic SST rate scaling detected (max exp(beta*anom)=%.3g > %.3g).\n",
        "This implies either wrong units for sst_anomaly or an overly large beta_sst."
      ),
      max(sst_scale, na.rm = TRUE), .sst_scale_max
    ), call. = FALSE)
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
  
  # ---- L2 severity split (bounded) ----
  if (is.null(p_hur_base) || !is.finite(p_hur_base)) {
    p_hur_base <- if (lambda_total > 0) pmax(0, pmin(1, lambda_hur / lambda_total)) else 0.5
  }
  
  gamma_intensity <- as.numeric(gamma_intensity)
  if (!is.finite(gamma_intensity)) gamma_intensity <- 0
  
  p_hur <- p_hur_base * (1 + gamma_intensity * sst_vec)
  p_hur <- pmin(0.99, pmax(0.01, p_hur))
  
  n_hur <- stats::rbinom(n_years_sim, size = n_total, prob = p_hur)
  n_ts  <- n_total - n_hur
  
  tibble::tibble(
    sim_year = seq_len(n_years_sim),
    activity_factor = A,
    sst_anomaly = sst_vec,
    climate_scale = sst_scale,
    activity_combined = A * sst_scale,
    p_hurricane = p_hur,
    n_total = .int_cap(n_total, "n_total"),
    n_ts = .int_cap(n_ts, "n_ts"),
    n_hur = .int_cap(n_hur, "n_hur")
  )
}




# =============================================================================
# 7) LEVEL 3: STORM CHARACTERISTIC PERTURBATION
# =============================================================================

#' Default climate-change perturbation parameters
#'
#' @description
#' Returns a named list of per-degC scaling factors for Level 3 storm property
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


#' Perturb sampled event properties for climate change (Level 3)
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
# 8) SST CONFIGURATION HELPER
# =============================================================================

#' Build a climate configuration object for the hazard model
#'
#' @description
#' Creates a climate configuration object for `run_hazard_model()` with
#' user-facing inputs first and expert knobs grouped in `advanced`.
#'
#' Scenarios follow the model's accepted values:
#' `"stationary"`, `"ssp126"`, `"ssp245"`, `"ssp585"`.
#'
#' When `scenario = "stationary"`, climate conditioning is disabled by default
#' and behaves like a stationary run.
#'
#' @param enabled Logical; whether to enable SST-conditioned modifications.
#' @param sst_source Character; one of "builtin", "csv", "ersst_nc".
#' @param sst_path Optional character; path to SST data file (CSV or NetCDF).
#' @param baseline_years Integer vector; years for climatological baseline.
#' @param scenario Character; SST scenario for projection years.
#' @param scenario_start_year Integer; first year of the simulation scenario.
#' @param advanced Optional named list of expert parameters. Most users should
#'   leave this as `NULL`. Supported names:
#'   `beta_sst`, `beta_prior`, `gamma_intensity`, `gamma_prior`, `cc_params`.
#'
#' @return A list with class "sst_cfg" containing all climate configuration parameters.
#' @export
make_sst_cfg <- function(enabled = TRUE,
                         sst_source = c("builtin", "csv", "ersst_nc"),
                         sst_path = NULL,
                         baseline_years = 1991L:2020L,
                         scenario = "stationary",
                         scenario_start_year = 2025L,
                         advanced = NULL) {
  sst_source <- match.arg(sst_source)
  scenario <- match.arg(scenario, choices = c("stationary", "ssp126", "ssp245", "ssp585"))

  defaults <- list(
    beta_sst = NULL,
    beta_prior = 0.6,
    gamma_intensity = NULL,
    gamma_prior = 0.065,
    cc_params = NULL
  )

  if (is.null(advanced)) {
    advanced <- defaults
  } else {
    if (!is.list(advanced)) {
      stop("advanced must be NULL or a named list.", call. = FALSE)
    }
    unknown <- setdiff(names(advanced), names(defaults))
    if (length(unknown) > 0) {
      stop("Unknown names in advanced: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    advanced <- utils::modifyList(defaults, advanced)
  }

  if (!is.null(advanced$cc_params) && !is.list(advanced$cc_params)) {
    stop("cc_params must be NULL (disabled) or a named list of scaling factors.")
  }

  resolved_enabled <- isTRUE(enabled) && !identical(scenario, "stationary")

  cfg <- list(
    enabled = resolved_enabled,
    sst_source = sst_source,
    sst_path = sst_path,
    baseline_years = baseline_years,
    beta_sst = advanced$beta_sst,
    beta_prior = advanced$beta_prior,
    gamma_intensity = advanced$gamma_intensity,
    gamma_prior = advanced$gamma_prior,
    scenario = scenario,
    scenario_start_year = scenario_start_year,
    cc_params = advanced$cc_params,
    advanced = advanced
  )
  class(cfg) <- c("sst_cfg", "list")
  cfg
}

#' @export
print.sst_cfg <- function(x, ...) {
  cat("Climate configuration\n")
  cat(sprintf("  Scenario      : %s\n", x$scenario))
  cat(sprintf("  SST source    : %s\n", x$sst_source))
  cat(sprintf("  Baseline      : %d-%d\n", min(x$baseline_years), max(x$baseline_years)))

  if (isTRUE(x$enabled)) {
    beta_mode <- if (!is.null(x$beta_sst) && is.finite(x$beta_sst)) {
      sprintf("fixed at %.3f", x$beta_sst)
    } else {
      sprintf("estimated from data (prior: %.3f)", x$beta_prior)
    }
    gamma_mode <- if (!is.null(x$gamma_intensity) && is.finite(x$gamma_intensity)) {
      sprintf("fixed at %.4f", x$gamma_intensity)
    } else {
      sprintf("estimated from data (prior: %.4f)", x$gamma_prior)
    }
    cat(sprintf("  beta_SST      : %s\n", beta_mode))
    cat(sprintf("  gamma intensity: %s\n", gamma_mode))
  } else {
    cat("  Climate mode  : disabled (stationary)\n")
  }

  l3_status <- if (is.null(x$cc_params)) "disabled" else "enabled"
  cat(sprintf("  Level 3 perturb: %s\n", l3_status))
  invisible(x)
}

#' Load and prepare SST data based on configuration
#'
#' @description
#' Reads SST data from the configured source, computes anomalies, and
#' optionally estimates Î²_SST (Level 1) and Î³ (Level 2). Returns a
#' processed climate object ready for use in the hazard model.
#'
#' @param sst_cfg List from `make_sst_cfg()`.
#' @param annual_counts Optional tibble of annual counts for Î² and Î³ estimation.
#' @param lambda_table Optional tibble from `compute_lambda_table()` for p_HUR_base.
#' @param min_year Integer; passed to estimation functions.
#' @param verbose Logical.
#'
#' @return A list with:
#'   \item{sst_df}{Tibble of historical SST with anomalies.}
#'   \item{beta_sst}{Estimated or user-supplied Î²_SST.}
#'   \item{beta_info}{Full output from `estimate_beta_sst()` (or NULL).}
#'   \item{gamma_intensity}{Estimated or user-supplied Î³.}
#'   \item{gamma_info}{Full output from `estimate_gamma_intensity()` (or NULL).}
#'   \item{p_hur_base}{Baseline hurricane fraction.}
#'   \item{sst_cfg}{The input configuration.}
#'
#' @export
prepare_sst_data <- function(sst_cfg,
                             annual_counts = NULL,
                             lambda_table = NULL,
                             min_year = 1970L,
                             verbose = TRUE) {
  if (!inherits(sst_cfg, "sst_cfg")) {
    stop("sst_cfg must be created by make_sst_cfg().")
  }
  
  if (!sst_cfg$enabled) {
    if (verbose) message("[SST] Climate conditioning disabled.")
    return(list(
      sst_df = NULL, beta_sst = 0, beta_info = NULL,
      gamma_intensity = 0, gamma_info = NULL, p_hur_base = NA_real_,
      cc_params = NULL,
      sst_cfg = sst_cfg
    ))
  }
  
  # Load SST data
  sst_raw <- switch(sst_cfg$sst_source,
                    builtin   = get_mdr_sst_builtin(),
                    csv       = read_mdr_sst_csv(sst_cfg$sst_path),
                    ersst_nc  = read_mdr_sst_ersst(sst_cfg$sst_path)
  )
  
  if (verbose) {
    message(sprintf("[SST] Loaded %d years of MDR SST (%d-%d) from %s",
                    nrow(sst_raw), min(sst_raw$year), max(sst_raw$year), sst_cfg$sst_source))
  }
  
  # Compute anomalies (compute_sst_anomaly should be strict; see next patch)
  sst_df <- compute_sst_anomaly(sst_raw, baseline_years = sst_cfg$baseline_years)
  
  if (verbose) {
    message(sprintf("[SST] Baseline (%.0f-%.0f): %.2f C | Anomaly range: [%+.2f, %+.2f] C",
                    min(sst_cfg$baseline_years), max(sst_cfg$baseline_years),
                    sst_df$sst_clim[1],
                    min(sst_df$sst_anomaly, na.rm = TRUE),
                    max(sst_df$sst_anomaly, na.rm = TRUE)))
  }
  
  # L1: beta
  beta_info <- NULL
  if (!is.null(sst_cfg$beta_sst) && is.finite(sst_cfg$beta_sst)) {
    beta_sst <- as.numeric(sst_cfg$beta_sst)
    if (verbose) message(sprintf("[L1] Using user-supplied beta_SST = %.3f", beta_sst))
  } else if (!is.null(annual_counts)) {
    beta_info <- estimate_beta_sst(
      annual_counts = annual_counts,
      sst_df = sst_df,
      min_year = min_year,
      beta_prior = sst_cfg$beta_prior,
      verbose = verbose
    )
    beta_sst <- as.numeric(beta_info$beta_sst)
  } else {
    beta_sst <- if (!is.null(sst_cfg$beta_prior)) as.numeric(sst_cfg$beta_prior) else 0
    if (verbose) message(sprintf("[L1] No annual_counts provided for beta estimation. Using prior: %.3f", beta_sst))
  }
  
  # L2: gamma (IMPORTANT FIX)
  gamma_info <- NULL
  p_hur_base <- NA_real_
  
  if (!is.null(lambda_table)) {
    p_hur_base <- compute_p_hur_base(lambda_table)
  }
  
  if (!is.null(sst_cfg$gamma_intensity) && is.finite(sst_cfg$gamma_intensity)) {
    gamma_intensity <- as.numeric(sst_cfg$gamma_intensity)
    if (verbose) message(sprintf("[L2] Using user-supplied gamma = %.4f", gamma_intensity))
  } else if (!is.null(annual_counts)) {
    gamma_info <- estimate_gamma_intensity(
      annual_counts = annual_counts,
      sst_df = sst_df,
      min_year = min_year,
      gamma_prior = sst_cfg$gamma_prior,
      verbose = verbose
    )
    gamma_intensity <- as.numeric(gamma_info$gamma)
    if (is.na(p_hur_base) && !is.null(gamma_info$p_hur_base) && is.finite(gamma_info$p_hur_base)) {
      p_hur_base <- as.numeric(gamma_info$p_hur_base)
    }
  } else {
    # FIX: without data AND without explicit user value, do NOT activate L2.
    gamma_intensity <- 0
    if (verbose) message("[L2] No annual_counts for gamma estimation and no user gamma provided. Using gamma = 0.")
  }
  
  if (verbose && is.finite(gamma_intensity) && gamma_intensity != 0) {
    message(sprintf("[L2] Intensity shift: gamma=%.4f, p_HUR_base=%.3f",
                    gamma_intensity, p_hur_base))
    message(sprintf("[L2] At +1C SST: p_HUR -> %.3f (was %.3f, %+.1f%%)",
                    pmin(0.99, p_hur_base * (1 + gamma_intensity)),
                    p_hur_base,
                    100 * gamma_intensity))
  }
  
  # L3 cc_params passthrough
  resolved_cc_params <- sst_cfg$cc_params
  if (!is.null(resolved_cc_params) && length(resolved_cc_params) == 0) {
    resolved_cc_params <- default_cc_params()
  }
  if (!is.null(resolved_cc_params)) {
    defs <- default_cc_params()
    for (nm in names(defs)) {
      if (is.null(resolved_cc_params[[nm]])) resolved_cc_params[[nm]] <- defs[[nm]]
    }
    if (verbose) {
      message(sprintf("[L3] Storm perturbation enabled: v_scale=%+.2f, r_scale=%+.2f, speed_scale=%+.2f, precip_scale=%+.2f per degC",
                      resolved_cc_params$v_scale, resolved_cc_params$r_scale,
                      resolved_cc_params$speed_scale, resolved_cc_params$precip_scale))
    }
  } else {
    if (verbose) message("[L3] Storm perturbation disabled (cc_params = NULL).")
  }

  list(
    sst_df = sst_df,
    beta_sst = beta_sst,
    beta_info = beta_info,
    gamma_intensity = gamma_intensity,
    gamma_info = gamma_info,
    p_hur_base = p_hur_base,
    cc_params = resolved_cc_params,
    sst_cfg = sst_cfg
  )
}




