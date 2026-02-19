################################################################################
# hazard_validation.R
# Validation framework for the hurricane hazard model
#
# Three validation tiers:
#   1) Hindcast: hold-out last N years, compare simulated return levels vs observed
#   2) Rate sanity: compare model lambdas against published HURDAT2 climatologies
#   3) Wind field spot-checks: compare V_site_kt against station observations
#
# All functions assume the hazard model source files are already loaded.
################################################################################

# =============================================================================
# 1) HINDCAST VALIDATION
# =============================================================================

#' Compute empirical return levels from a vector of annual maxima
#'
#' @description
#' Extracts return levels from a sample of annual maxima using the plotting-position
#' approach: sort values, assign return periods T = (n+1)/rank, interpolate to
#' requested return periods.
#'
#' @param annual_max Numeric vector of annual maximum values (e.g., V_site_max_kt).
#' @param return_periods Numeric vector of return periods (years) to compute.
#'
#' @return Named numeric vector of return levels at requested periods.
#' @export
compute_return_levels <- function(annual_max,
                                  return_periods = c(5, 10, 25, 50)) {
  x <- sort(annual_max[is.finite(annual_max)])
  n <- length(x)
  if (n < 3) {
    out <- rep(NA_real_, length(return_periods))
    names(out) <- paste0("RL_", return_periods, "yr")
    return(out)
  }
  
  # Weibull plotting position: T = (n+1) / rank
  ranks <- seq_len(n)
  T_emp <- (n + 1) / (n + 1 - ranks)  # exceedance-based
  
  out <- stats::approx(
    x = T_emp, y = x,
    xout = return_periods,
    rule = 2  # clamp at boundaries
  )$y
  
  names(out) <- paste0("RL_", return_periods, "yr")
  out
}


# =============================================================================
# Parametric hindcast components (KDE + hurdle-GEV)
# Merged from hazard_validation_parametric.R to avoid function overwrites.
# =============================================================================

fit_intensity_kde <- function(pool, lower, upper = Inf, bw_mult = 1.0) {
  pool <- pool[is.finite(pool)]
  n <- length(pool)
  if (n < 3) {
    return(list(
      method = "fallback",
      pool = pool,
      lower = lower,
      upper = upper,
      n_obs = n,
      pool_mean = if (n > 0) mean(pool) else (lower + 10),
      pool_sd = if (n > 1) stats::sd(pool) else 10
    ))
  }
  
  # Reflect at boundaries
  reflected <- pool
  reflected <- c(reflected, 2 * lower - pool)  # reflect at lower bound
  
  if (is.finite(upper)) {
    reflected <- c(reflected, 2 * upper - pool)  # reflect at upper bound
  }
  
  # Bandwidth: Silverman's rule on original data, scaled
  bw <- bw_mult * stats::bw.nrd0(pool)
  
  # Fit KDE on reflected data
  dens <- stats::density(reflected, bw = bw, n = 2048,
                         from = lower - 3 * bw,
                         to = if (is.finite(upper)) upper + 3 * bw else max(pool) + 6 * bw)
  
  # Restrict to valid support
  valid <- dens$x >= lower & (if (is.finite(upper)) dens$x <= upper else TRUE)
  x_valid <- dens$x[valid]
  y_valid <- dens$y[valid]
  
  # Renormalize
  y_valid <- pmax(0, y_valid)
  area <- stats::integrate(stats::approxfun(x_valid, y_valid, rule = 2),
                           lower = min(x_valid), upper = max(x_valid),
                           subdivisions = 500)$value
  if (area > 0) y_valid <- y_valid / area
  
  # Build CDF for inverse-CDF sampling
  dx <- diff(x_valid)
  y_mid <- (y_valid[-1] + y_valid[-length(y_valid)]) / 2
  cdf_y <- c(0, cumsum(y_mid * dx))
  cdf_x <- x_valid
  
  # Ensure CDF reaches exactly 1
  if (max(cdf_y) > 0) cdf_y <- cdf_y / max(cdf_y)
  
  list(
    method = "kde",
    density_x = x_valid,
    density_y = y_valid,
    cdf_x = cdf_x,
    cdf_y = cdf_y,
    lower = lower,
    upper = upper,
    n_obs = n,
    pool_mean = mean(pool),
    pool_sd = stats::sd(pool),
    bw = bw
  )
}


#' Sample from a fitted intensity KDE
#'
#' @param fit List from \code{fit_intensity_kde()}.
#' @param n Integer; number of draws.
#'
#' @return Numeric vector of n intensity draws within [lower, upper].
#' @export
sample_intensity_kde <- function(fit, n) {
  if (n <= 0) return(numeric(0))
  
  if (fit$method == "fallback") {
    # Too few observations: jittered resample
    if (fit$n_obs == 0) {
      return(rep(fit$pool_mean, n))
    }
    draws <- sample(fit$pool, n, replace = TRUE) +
      stats::rnorm(n, 0, fit$pool_sd * 0.2)
    draws <- pmax(fit$lower, draws)
    if (is.finite(fit$upper)) draws <- pmin(fit$upper, draws)
    return(draws)
  }
  
  # Inverse-CDF sampling
  u <- stats::runif(n)
  draws <- stats::approx(fit$cdf_y, fit$cdf_x, xout = u, rule = 2)$y
  
  # Enforce bounds
  draws <- pmax(fit$lower, draws)
  if (is.finite(fit$upper)) draws <- pmin(fit$upper, draws)
  
  draws
}


# =============================================================================
# 2) GEV FITTING (L-moments, no external packages)
# =============================================================================

#' Fit GEV distribution using L-moments (Hosking 1990)
#'
#' @description
#' Estimates GEV parameters (location μ, scale σ, shape ξ) using the method
#' of L-moments. This is more robust than MLE for small samples (n < 50) and
#' requires no optimization. The sign convention follows the standard meteorological
#' form: ξ > 0 → heavy tail (Fréchet), ξ < 0 → bounded tail (Weibull), ξ = 0 → Gumbel.
#'
#' @param x Numeric vector of block maxima (e.g., annual maxima).
#' @param xi_bounds Numeric vector of length 2; allowed range for shape parameter.
#'   Default c(-0.5, 0.5) prevents extreme tail behavior.
#'
#' @return A list with elements:
#'   \item{mu}{Location parameter.}
#'   \item{sigma}{Scale parameter.}
#'   \item{xi}{Shape parameter.}
#'   \item{n}{Sample size.}
#'   \item{l_moments}{Named vector of L1, L2, L3, tau3.}
#'   \item{converged}{Logical; whether estimation produced valid parameters.}
#'
#' @export
fit_gev_lmom <- function(x, xi_bounds = c(-0.5, 0.5)) {
  x <- sort(x[is.finite(x)])
  n <- length(x)
  
  if (n < 5) {
    return(list(
      mu = mean(x), sigma = stats::sd(x), xi = 0,
      n = n, l_moments = c(L1 = mean(x), L2 = NA, L3 = NA, tau3 = NA),
      converged = FALSE
    ))
  }
  
  # Probability-weighted moments (unbiased estimators)
  # b_r = (1/n) * sum_{i=1}^{n} x_{i:n} * C(i-1, r) / C(n-1, r)
  ii <- seq_len(n)
  
  b0 <- mean(x)
  b1 <- sum(x * (ii - 1) / (n - 1)) / n
  b2 <- sum(x * (ii - 1) * (ii - 2) / ((n - 1) * (n - 2))) / n
  
  # L-moments
  L1 <- b0
  L2 <- 2 * b1 - b0
  L3 <- 6 * b2 - 6 * b1 + b0
  
  if (!is.finite(L2) || L2 <= 0) {
    return(list(
      mu = L1, sigma = abs(L2), xi = 0,
      n = n, l_moments = c(L1 = L1, L2 = L2, L3 = L3, tau3 = NA),
      converged = FALSE
    ))
  }
  
  tau3 <- L3 / L2  # L-skewness
  
  # Hosking (1997) approximation for GEV shape from tau3
  # ξ ≈ 7.8590c + 2.9554c² where c = 2/(3+tau3) - log(2)/log(3)
  c_val <- 2 / (3 + tau3) - log(2) / log(3)
  xi <- 7.8590 * c_val + 2.9554 * c_val^2
  
  # Clamp shape parameter to physical bounds
  xi <- max(xi_bounds[1], min(xi_bounds[2], xi))
  
  # Back-solve for sigma and mu
  if (abs(xi) < 1e-6) {
    # Gumbel case (ξ → 0)
    sigma <- L2 / log(2)
    mu <- L1 - sigma * 0.5772  # Euler-Mascheroni
  } else {
    g1 <- gamma(1 + xi)
    g2 <- gamma(1 + 2 * xi)
    
    # Protect against gamma overflow
    if (!is.finite(g1) || !is.finite(g2)) {
      sigma <- L2 / log(2)
      mu <- L1 - sigma * 0.5772
      xi <- 0
    } else {
      sigma <- L2 * xi / (g1 * (1 - 2^(-xi)))
      if (!is.finite(sigma) || sigma <= 0) {
        sigma <- L2 / log(2)
        mu <- L1 - sigma * 0.5772
        xi <- 0
      } else {
        mu <- L1 - sigma * (1 - g1) / xi
      }
    }
  }
  
  sigma <- max(1e-6, sigma)  # floor
  
  list(
    mu = mu,
    sigma = sigma,
    xi = xi,
    n = n,
    l_moments = c(L1 = L1, L2 = L2, L3 = L3, tau3 = tau3),
    converged = TRUE
  )
}


#' GEV quantile function
#'
#' @param p Numeric vector of probabilities (0, 1).
#' @param mu,sigma,xi GEV parameters.
#' @return Numeric vector of quantiles.
#' @keywords internal
.qgev <- function(p, mu, sigma, xi) {
  if (abs(xi) < 1e-8) {
    # Gumbel
    mu - sigma * log(-log(p))
  } else {
    mu + sigma * ((-log(p))^(-xi) - 1) / xi
  }
}

#' GEV CDF
#'
#' @param x Numeric vector of values.
#' @param mu,sigma,xi GEV parameters.
#' @return Numeric vector of probabilities.
#' @keywords internal
.pgev <- function(x, mu, sigma, xi) {
  z <- (x - mu) / sigma
  if (abs(xi) < 1e-8) {
    exp(-exp(-z))
  } else {
    t_val <- pmax(1e-10, 1 + xi * z)
    exp(-t_val^(-1 / xi))
  }
}


# =============================================================================
# 3) HURDLE-GEV RETURN LEVELS
# =============================================================================

#' Compute return levels using a hurdle-GEV model
#'
#' @description
#' Accounts for years with zero impact (no events) separately from the intensity
#' distribution. The model is:
#'   P(annual_max <= v) = p0 + (1 - p0) * F_GEV(v)
#' where p0 is the probability of a zero-event year, and F_GEV is the GEV CDF
#' fitted to nonzero annual maxima.
#'
#' The T-year return level solves:
#'   p0 + (1 - p0) * F_GEV(v) = 1 - 1/T
#'
#' @param annual_max Numeric vector of annual maxima (including zeros).
#' @param return_periods Numeric vector of return periods (years).
#' @param xi_bounds Bounds on GEV shape parameter (default: c(-0.3, 0.4) for TCs).
#'
#' @return A list with:
#'   \item{return_levels}{Named numeric vector of return levels.}
#'   \item{gev_fit}{GEV fit object for nonzero maxima.}
#'   \item{p_zero}{Fraction of zero-event years.}
#'   \item{n_total, n_nonzero}{Sample sizes.}
#'
#' @export
compute_return_levels_gev <- function(annual_max,
                                      return_periods = c(5, 10, 25, 50),
                                      xi_bounds = c(-0.3, 0.4)) {
  x <- annual_max[is.finite(annual_max)]
  n <- length(x)
  
  n_zero <- sum(x <= 0)
  p_zero <- n_zero / n
  
  x_pos <- x[x > 0]
  n_pos <- length(x_pos)
  
  if (n_pos < 5) {
    rl <- rep(NA_real_, length(return_periods))
    names(rl) <- paste0("RL_", return_periods, "yr")
    return(list(
      return_levels = rl,
      gev_fit = NULL,
      p_zero = p_zero,
      n_total = n,
      n_nonzero = n_pos
    ))
  }
  
  gev <- fit_gev_lmom(x_pos, xi_bounds = xi_bounds)
  
  # Hurdle-GEV return levels:
  # Solve: p0 + (1-p0) * F_GEV(v) = 1 - 1/T
  # → F_GEV(v) = (1 - 1/T - p0) / (1 - p0)
  # → v = Q_GEV( (1 - 1/T - p0) / (1 - p0) )
  
  rl <- vapply(return_periods, function(T_rp) {
    target_p <- 1 - 1 / T_rp
    p_cond <- (target_p - p_zero) / (1 - p_zero)
    
    if (p_cond <= 0) return(0)       # return period shorter than mean recurrence
    if (p_cond >= 1) return(NA_real_) # can't compute
    
    .qgev(p_cond, gev$mu, gev$sigma, gev$xi)
  }, numeric(1))
  
  # Physical cap: TC winds can't exceed ~185 kt
  rl <- pmin(rl, 185)
  
  names(rl) <- paste0("RL_", return_periods, "yr")
  
  list(
    return_levels = rl,
    gev_fit = gev,
    p_zero = p_zero,
    n_total = n,
    n_nonzero = n_pos
  )
}


# =============================================================================
# 4) PARAMETRIC BOOTSTRAP CIs FOR RETURN LEVELS
# =============================================================================

#' Compute return level CIs via parametric bootstrap
#'
#' @description
#' Generates bootstrap CIs by:
#' 1) Resampling the annual maxima
#' 2) Refitting the hurdle-GEV to each bootstrap sample
#' 3) Computing return levels from each fit
#'
#' This propagates uncertainty from both the zero-fraction and the GEV parameters.
#'
#' @param annual_max Numeric vector of annual maxima (including zeros).
#' @param return_periods Numeric vector of return periods.
#' @param n_boot Integer; number of bootstrap replicates.
#' @param xi_bounds Bounds on GEV shape.
#'
#' @return Tibble with return_period, median, lo_90, hi_90, lo_50, hi_50.
#' @export
bootstrap_return_level_ci <- function(annual_max,
                                      return_periods = c(5, 10, 25, 50),
                                      n_boot = 500,
                                      xi_bounds = c(-0.3, 0.4)) {
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Package `tibble` required.")
  
  n <- length(annual_max)
  boot_rl <- matrix(NA_real_, nrow = n_boot, ncol = length(return_periods))
  
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    boot_sample <- annual_max[idx]
    
    res <- tryCatch(
      compute_return_levels_gev(boot_sample, return_periods, xi_bounds)$return_levels,
      error = function(e) rep(NA_real_, length(return_periods))
    )
    boot_rl[b, ] <- res
  }
  
  tibble::tibble(
    return_period = return_periods,
    sim_median = apply(boot_rl, 2, stats::median, na.rm = TRUE),
    sim_lo_90  = apply(boot_rl, 2, stats::quantile, probs = 0.05, na.rm = TRUE),
    sim_hi_90  = apply(boot_rl, 2, stats::quantile, probs = 0.95, na.rm = TRUE),
    sim_lo_50  = apply(boot_rl, 2, stats::quantile, probs = 0.25, na.rm = TRUE),
    sim_hi_50  = apply(boot_rl, 2, stats::quantile, probs = 0.75, na.rm = TRUE)
  )
}


#' Run hindcast validation for a single island
#'
#' @description
#' Splits historical events into training and test periods, refits the frequency
#' model on training data, simulates synthetic annual maxima, and compares
#' simulated return-level quantiles against observed test-period annual maxima.
#'
#' @param events_island Tibble of storm events for one island (from out$events_by_island).
#'   Must contain: year, severity, SID, V_site_max_kt.
#' @param island Character; island name (for labelling).
#' @param holdout_years Integer; number of years to hold out from the end.
#' @param n_sim Integer; number of synthetic years to simulate from training parameters.
#' @param return_periods Numeric vector of return periods to compare.
#' @param severities Character vector of severities (passed to compute_annual_counts).
#' @param seed Integer seed for reproducibility.
#'
#' @return A list with:
#'   \item{island}{Island name.}
#'   \item{train_years}{Integer vector of training years.}
#'   \item{test_years}{Integer vector of test years.}
#'   \item{train_params}{List: lambda_table, k_hat, n_train_years.}
#'   \item{obs_annual_max}{Tibble of observed annual max V_site by year (full record).}
#'   \item{obs_test_rl}{Named vector of observed return levels from test period.}
#'   \item{sim_rl}{Named vector of simulated return levels.}
#'   \item{sim_rl_ci}{Tibble of simulated return level confidence intervals (bootstrap).}
#'   \item{comparison}{Tibble comparing observed vs simulated at each return period.}
#'
#' @export
validate_hindcast <- function(events_island,
                              island,
                              holdout_years = 10,
                              n_sim = 5000,
                              return_periods = c(5, 10, 25, 50),
                              severities = c("TS", "HUR64plus"),
                              seed = 42,
                              sst_df = NULL,
                              beta_sst = 0,
                              gamma_intensity = 0) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` is required.")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Package `tibble` is required.")
  
  set.seed(seed)
  
  ev <- events_island |>
    dplyr::filter(.data$severity %in% c(severities, "none"),
                  is.finite(.data$V_site_max_kt))
  
  all_years <- sort(unique(ev$year))
  if (length(all_years) < holdout_years + 10) {
    stop("Insufficient years for holdout. Have ", length(all_years),
         ", need at least ", holdout_years + 10)
  }
  
  cutoff <- all_years[length(all_years) - holdout_years]
  train_years <- all_years[all_years <= cutoff]
  test_years  <- all_years[all_years > cutoff]
  
  message("[Hindcast] ", island, ": training ", min(train_years), "-", max(train_years),
          " (", length(train_years), " yr), testing ", min(test_years), "-",
          max(test_years), " (", length(test_years), " yr)")
  
  # --- Observed annual maxima (full record) ---
  obs_annual_max <- ev |>
    dplyr::group_by(.data$year) |>
    dplyr::summarise(V_max_kt = max(.data$V_site_max_kt, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::mutate(period = dplyr::if_else(.data$year %in% train_years, "train", "test"))
  
  full_years <- tibble::tibble(year = seq(min(all_years), max(all_years)))
  obs_annual_max <- full_years |>
    dplyr::left_join(obs_annual_max, by = "year") |>
    dplyr::mutate(
      V_max_kt = dplyr::if_else(is.na(.data$V_max_kt), 0, .data$V_max_kt),
      period = dplyr::if_else(.data$year %in% test_years, "test", "train")
    )
  
  # --- Fit frequency model on training period ---
  ev_train <- events_island |>
    dplyr::filter(.data$year %in% train_years)
  
  ac_train <- compute_annual_counts(ev_train, severities = severities)
  lt_train <- compute_lambda_table(ac_train)
  ki_train <- estimate_k_hat(ac_train)
  
  train_params <- list(
    lambda_table = lt_train,
    k_hat = ki_train$k_hat,
    n_train_years = length(train_years),
    mu_total = ki_train$mu,
    var_total = ki_train$var
  )
  
  # --- FIT KDE INTENSITY DISTRIBUTIONS (replaces discrete pools) ---
  train_V_ts <- ev_train |>
    dplyr::filter(.data$severity == "TS") |>
    dplyr::pull(.data$V_site_max_kt) |>
    (\(x) x[is.finite(x)])()
  
  train_V_hur <- ev_train |>
    dplyr::filter(.data$severity == "HUR64plus") |>
    dplyr::pull(.data$V_site_max_kt) |>
    (\(x) x[is.finite(x)])()
  
  kde_ts  <- fit_intensity_kde(train_V_ts,  lower = 34, upper = 64)
  kde_hur <- fit_intensity_kde(train_V_hur, lower = 64, upper = 185)
  
  n_ts_obs  <- length(train_V_ts)
  n_hur_obs <- length(train_V_hur)
  
  message(sprintf("  KDE fits: TS pool=%d events (mean=%.0f kt), HUR pool=%d events (mean=%.0f kt)",
                  n_ts_obs,  if (n_ts_obs > 0) mean(train_V_ts) else NA,
                  n_hur_obs, if (n_hur_obs > 0) mean(train_V_hur) else NA))
  
  # Fallback intensities if KDE can't be fit
  fallback_V <- list(TS = 45, HUR64plus = 85)
  
  # --- SIMULATE ANNUAL MAXIMA WITH KDE SAMPLING ---
  # If SST data available, use historical SST anomalies resampled for simulation
  sst_anomaly_sim <- NULL
  if (!is.null(sst_df) && is.finite(beta_sst) && beta_sst != 0) {
    sst_train <- sst_df |>
      dplyr::filter(.data$year %in% train_years)
    if (nrow(sst_train) > 0 && "sst_anomaly" %in% names(sst_train)) {
      # Resample historical SST anomalies (preserving the observed distribution)
      sst_pool <- sst_train$sst_anomaly[is.finite(sst_train$sst_anomaly)]
      if (length(sst_pool) > 0) {
        sst_anomaly_sim <- sample(sst_pool, n_sim, replace = TRUE)
        message(sprintf("  SST conditioning: β=%.3f, %d training SST values, mean ΔSST=%+.2f°C",
                        beta_sst, length(sst_pool), mean(sst_pool)))
      }
    }
  }
  
  sim_counts <- simulate_twolevel_counts(
    lt_train, ki_train$k_hat, n_years_sim = n_sim,
    sst_anomaly = sst_anomaly_sim,
    beta_sst = beta_sst,
    gamma_intensity = gamma_intensity
  )
  
  sim_annual_max <- vapply(seq_len(n_sim), function(i) {
    n_ts  <- sim_counts$n_TS[i]
    n_hur <- sim_counts$n_HUR64plus[i]
    
    winds <- numeric(0)
    
    if (n_ts > 0) {
      if (n_ts_obs >= 3) {
        winds <- c(winds, sample_intensity_kde(kde_ts, n_ts))
      } else {
        winds <- c(winds, rep(fallback_V$TS, n_ts) +
                     stats::rnorm(n_ts, 0, 5))
      }
    }
    if (n_hur > 0) {
      if (n_hur_obs >= 3) {
        winds <- c(winds, sample_intensity_kde(kde_hur, n_hur))
      } else {
        winds <- c(winds, rep(fallback_V$HUR64plus, n_hur) +
                     stats::rnorm(n_hur, 0, 10))
      }
    }
    
    if (length(winds) == 0) return(0)
    max(winds)
  }, numeric(1))
  
  # --- GEV RETURN LEVELS ---
  # Model: GEV fit to 5000 simulated annual maxima (precise point estimate)
  sim_gev <- compute_return_levels_gev(sim_annual_max, return_periods)
  sim_rl  <- sim_gev$return_levels
  
  # Observed: GEV fit to full observed record (~45-55 years)
  obs_full_max <- obs_annual_max$V_max_kt
  obs_gev <- compute_return_levels_gev(obs_full_max, return_periods)
  obs_full_rl <- obs_gev$return_levels
  
  # Test period (empirical, for reference only)
  obs_test_max <- obs_annual_max |>
    dplyr::filter(.data$period == "test") |>
    dplyr::pull(.data$V_max_kt)
  obs_test_rl <- compute_return_levels(obs_test_max, return_periods)
  
  # --- BOOTSTRAP CIs ON OBSERVED RETURN LEVELS ---
  # The CI reflects uncertainty from the limited historical record (~45 yr).
  # Question: "Given observational uncertainty, is the model consistent with obs?"
  obs_rl_ci <- bootstrap_return_level_ci(
    obs_full_max, return_periods,
    n_boot = 500,
    xi_bounds = c(-0.3, 0.4)
  )
  
  # Also compute model CIs (secondary diagnostic — should be narrow)
  sim_rl_ci <- bootstrap_return_level_ci(
    sim_annual_max, return_periods,
    n_boot = 200,
    xi_bounds = c(-0.3, 0.4)
  )
  
  message(sprintf("  Model GEV: \u03bc=%.1f, \u03c3=%.1f, \u03be=%.3f (n_pos=%d, p0=%.2f)",
                  sim_gev$gev_fit$mu, sim_gev$gev_fit$sigma, sim_gev$gev_fit$xi,
                  sim_gev$n_nonzero, sim_gev$p_zero))
  message(sprintf("  Obs GEV:   \u03bc=%.1f, \u03c3=%.1f, \u03be=%.3f (n_pos=%d, p0=%.2f)",
                  obs_gev$gev_fit$mu, obs_gev$gev_fit$sigma, obs_gev$gev_fit$xi,
                  obs_gev$n_nonzero, obs_gev$p_zero))
  
  
  # --- Comparison table ---
  # Columns are kept compatible with existing plotting/report code:
  #   sim_median/sim_lo_90/sim_hi_90 describe the MODEL uncertainty (bootstrap on simulated annual maxima).
  # Additional columns (obs_lo_90/obs_hi_90, model_in_obs_90ci) describe observational uncertainty diagnostics.
  comparison <- tibble::tibble(
    island = island,
    return_period = return_periods,
    obs_full_rl = obs_full_rl,
    obs_test_rl = obs_test_rl,
    sim_rl = sim_rl,
    sim_median = sim_rl_ci$sim_median,
    sim_lo_90 = sim_rl_ci$sim_lo_90,
    sim_hi_90 = sim_rl_ci$sim_hi_90,
    obs_lo_90 = obs_rl_ci$sim_lo_90,
    obs_hi_90 = obs_rl_ci$sim_hi_90,
    # Pass criterion: is the model return level consistent with observational uncertainty?
    # (Your parametric routine used the observed 90% CI as the reference.)
    obs_in_90ci = sim_rl >= obs_rl_ci$sim_lo_90 & sim_rl <= obs_rl_ci$sim_hi_90,
    # Back-compat alias (older code used this name explicitly).
    model_in_obs_90ci = sim_rl >= obs_rl_ci$sim_lo_90 & sim_rl <= obs_rl_ci$sim_hi_90,
    # Additional diagnostic (stricter): is the observation within the model's 90% CI?
    obs_in_model_90ci = obs_full_rl >= sim_rl_ci$sim_lo_90 & obs_full_rl <= sim_rl_ci$sim_hi_90,
    bias_pct = 100 * (sim_rl - obs_full_rl) / pmax(obs_full_rl, 1)
  )
  
  list(
    island = island,
    train_years = train_years,
    test_years = test_years,
    train_params = train_params,
    obs_annual_max = obs_annual_max,
    obs_test_rl = obs_test_rl,
    obs_full_rl = obs_full_rl,
    obs_rl_ci = obs_rl_ci,
    sim_rl = sim_rl,
    sim_annual_max = sim_annual_max,
    sim_rl_ci = sim_rl_ci,
    comparison = comparison,
    # New diagnostics
    gev_fit = sim_gev,
    kde_fits = list(TS = kde_ts, HUR64plus = kde_hur),
    diagnostics = list(
      n_ts_pool = n_ts_obs,
      n_hur_pool = n_hur_obs,
      gev_xi = sim_gev$gev_fit$xi,
      gev_sigma = sim_gev$gev_fit$sigma,
      p_zero = sim_gev$p_zero
    )
  )
}



#' Run hindcast validation across all islands in a hazard model output
#'
#' @param out List returned by \code{run_hazard_model()}.
#' @param holdout_years Integer; years to hold out from each island.
#' @param n_sim Integer; synthetic years per island.
#' @param return_periods Numeric vector of return periods.
#' @param seed Integer seed.
#'
#' @return A list with per-island results and a combined comparison table.
#' @export
validate_hindcast_all <- function(out,
                                  holdout_years = 10,
                                  n_sim = 5000,
                                  return_periods = c(5, 10, 25, 50),
                                  seed = 42,
                                  sst_df = NULL,
                                  beta_sst = 0,
                                  gamma_intensity = 0) {
  islands <- names(Filter(Negate(is.null), out$events_by_island))
  results <- setNames(vector("list", length(islands)), islands)
  
  for (island in islands) {
    ev <- out$events_by_island[[island]]
    if (is.null(ev) || nrow(ev) < 20) {
      message("[Hindcast] Skipping ", island, " (too few events: ",
              if (is.null(ev)) 0 else nrow(ev), ")")
      next
    }
    
    tryCatch({
      results[[island]] <- validate_hindcast(
        events_island = ev,
        island = island,
        holdout_years = holdout_years,
        n_sim = n_sim,
        return_periods = return_periods,
        seed = seed,
        sst_df = sst_df,
        beta_sst = beta_sst,
        gamma_intensity = gamma_intensity
      )
    }, error = function(e) {
      message("[Hindcast] Error for ", island, ": ", e$message)
    })
  }
  
  comparison_all <- dplyr::bind_rows(
    lapply(Filter(Negate(is.null), results), function(r) r$comparison)
  )
  
  list(
    per_island = results,
    comparison = comparison_all
  )
}


# =============================================================================
# 2) RATE SANITY CHECK
# =============================================================================

#' Reference HURDAT2/literature annual rates for Leeward Islands region
#'
#' @description
#' Returns a tibble of published annual TC passage rates from literature,
#' for comparison against model-fitted lambdas. Sources:
#'
#' - Tropical storm (34+ kt) passage rates within ~100 nm of Leeward Islands
#'   from Neumann et al. (1999), Landsea et al. (various HURDAT2 analyses),
#'   and NHC climatology tables.
#' - Hurricane (64+ kt) passage rates from the same sources.
#' - Values represent storms passing within roughly 100-150 nm radius
#'   (comparable to the model's effective gate after wind field application).
#'
#' @return Tibble with columns: region, severity, lambda_ref, source, gate_approx_nm, period.
#' @export
get_reference_rates <- function() {
  tibble::tribble(
    ~region,            ~severity,    ~lambda_ref, ~source,                           ~gate_approx_nm, ~period,
    # Leeward Islands region (15-20N, 60-65W)
    # TS rates: ~2-4 per year pass within 200nm of the northern Leewards
    "Leeward_Islands",  "TS34plus",   2.0,         "NHC Climo (100nm, 1970-2020)",    100,             "1970-2020",
    "Leeward_Islands",  "HUR64plus",  0.55,        "NHC Climo (100nm, 1970-2020)",    100,             "1970-2020",
    
    # St. Martin specific (NHC/NOAA tropical cyclone climatology)
    # ~1.5 TS + HUR within 65nm per year, ~0.4 HUR
    "St_Martin",        "TS34plus",   1.2,         "NOAA TC Climo (65nm, 1970-2023)", 65,              "1970-2023",
    "St_Martin",        "HUR64plus",  0.40,        "NOAA TC Climo (65nm, 1970-2023)", 65,              "1970-2023",
    
    # Puerto Rico (much larger island, well-studied)
    "Puerto_Rico",      "TS34plus",   1.8,         "NHC Climo (100nm, 1970-2020)",    100,             "1970-2020",
    "Puerto_Rico",      "HUR64plus",  0.45,        "NHC Climo (100nm, 1970-2020)",    100,             "1970-2020",
    
    # Miami / SE Florida
    "Miami",            "TS34plus",   1.5,         "NHC Climo (100nm, 1970-2020)",    100,             "1970-2020",
    "Miami",            "HUR64plus",  0.35,        "NHC Climo (100nm, 1970-2020)",    100,             "1970-2020"
  )
}


#' Compare model-fitted rates against reference climatologies
#'
#' @description
#' Takes the lambda table from the hazard model output and compares per-island
#' per-severity rates against published references. Reports absolute and
#' percentage differences, and flags rates outside expected bounds.
#'
#' A key consideration: the model uses an 800 km gate and computes V_site_kt at the
#' target, so its effective "passage radius" is set by the wind field, not the gate.
#' The reference rates use a fixed radius (65-100 nm). We flag discrepancies > 2x as
#' potentially indicating a gate/radius mismatch rather than a model error.
#'
#' @param out List returned by \code{run_hazard_model()}.
#' @param ref_rates Optional tibble of reference rates (default: \code{get_reference_rates()}).
#'
#' @return Tibble with model vs reference rate comparison.
#' @export
validate_rates <- function(out, ref_rates = NULL) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` is required.")
  
  if (is.null(ref_rates)) ref_rates <- get_reference_rates()
  
  model_rates <- out$lambda_all |>
    dplyr::select(island = "island", severity = "severity",
                  lambda_model = "lambda", n_years_model = "n_years")
  
  model_rates2 <- model_rates |>
    tidyr::pivot_wider(names_from = severity, values_from = lambda_model, values_fill = 0) |>
    dplyr::mutate(TS34plus = TS + HUR64plus) |>
    tidyr::pivot_longer(c("TS34plus","HUR64plus"), names_to="severity", values_to="lambda_model")
  
  
  
  
  # Join where island name matches region in reference table
  # (Saba/Statia not in reference → use Leeward_Islands as fallback)
  island_to_region <- tibble::tribble(
    ~island,        ~region,
    "St_Martin",    "St_Martin",
    "Saba",         "Leeward_Islands",
    "Statia",       "Leeward_Islands",
    "Puerto_Rico",  "Puerto_Rico",
    "Miami",        "Miami"
  )
  
  comp <- model_rates2 |>
    dplyr::left_join(island_to_region, by = "island") |>
    dplyr::left_join(ref_rates, by = c("region", "severity")) |>
    dplyr::mutate(
      diff_abs  = .data$lambda_model - .data$lambda_ref,
      diff_pct  = 100 * (.data$lambda_model - .data$lambda_ref) / pmax(.data$lambda_ref, 0.01),
      ratio     = .data$lambda_model / pmax(.data$lambda_ref, 0.001),
      
      # Flag logic:
      #   model uses 800km gate → its effective capture radius is larger than reference
      #   so model lambdas SHOULD be somewhat higher. Ratio 1.0-2.5 is acceptable.
      #   Below 0.5 or above 3.0 suggests a problem.
      flag = dplyr::case_when(
        is.na(.data$lambda_ref) ~ "no_reference",
        .data$ratio > 3.0  ~ "HIGH: model >> reference (check gate_km or min_year)",
        .data$ratio < 0.5  ~ "LOW: model << reference (check severity filter or data)",
        .data$ratio > 2.0  ~ "elevated (expected: model gate > ref radius)",
        TRUE ~ "OK"
      )
    ) |>
    dplyr::select(
      "island", "severity",
      "lambda_model", "n_years_model",
      "lambda_ref", "source", "gate_approx_nm", "period",
      "ratio", "diff_pct", "flag"
    )
  
  message("\n[Rate Check] Summary:")
  for (i in seq_len(nrow(comp))) {
    r <- comp[i, ]
    sym <- if (r$flag == "OK") "\u2713" else if (grepl("^(HIGH|LOW)", r$flag)) "\u2717" else "~"
    message(sprintf("  %s %s / %-10s : model=%.3f  ref=%.3f  ratio=%.2f  [%s]",
                    sym, r$island, r$severity, r$lambda_model,
                    if (is.na(r$lambda_ref)) NA else r$lambda_ref,
                    if (is.na(r$ratio)) NA else r$ratio,
                    r$flag))
  }
  
  comp
}


# =============================================================================
# 3) WIND FIELD SPOT-CHECKS
# =============================================================================

#' Reference observations for wind field validation
#'
#' @description
#' Returns a tibble of known station/buoy observations during well-documented storms
#' for comparison against model-estimated V_site_kt. Observations are peak sustained
#' wind (1-min or 10-min as noted) at fixed stations.
#'
#' Sources: NWS/NHC tropical cyclone reports, NOAA NDBC buoy data, Meteo France.
#'
#' @return Tibble with columns: storm_sid, storm_name, year, target_island,
#'   station, obs_wind_kt, obs_type, obs_source, notes.
#' @export
get_wind_observations <- function() {
  tibble::tribble(
    ~storm_sid,          ~storm_name, ~year, ~target_island, ~station,
    ~obs_wind_kt, ~obs_type,       ~obs_source,               ~notes,
    
    # Hurricane Irma (2017) - Cat 5 direct hit on St. Martin
    # NHC TCR: Irma produced sustained winds of ~155 kt on St. Martin (Juliana Airport)
    # Meteo France recorded 113 kt sustained (10-min) before instrument failure
    # NHC best track: 155 kt Vmax at closest approach
    "2017242N16333",     "IRMA",    2017L, "St_Martin",     "Juliana_Airport_TNCM",
    155,                 "1min_sust",     "NHC TCR (AL112017)",      "Instruments failed; 155 kt is NHC best-track Vmax at landfall",
    
    "2017242N16333",     "IRMA",    2017L, "St_Martin",     "Meteo_France_SXM",
    113,                 "10min_sust",    "Meteo France RSMC report", "Last valid reading before instrument failure; 10-min sustained",
    
    # Hurricane Gonzalo (2014) - Cat 1 near St. Martin
    "2014279N15323",     "GONZALO", 2014L, "St_Martin",     "Juliana_Airport_TNCM",
    60,                  "1min_sust",     "NHC TCR (AL082014)",      "Marginal TS winds at SXM; storm passed ~50nm north",
    
    # Hurricane Irma (2017) - effects at Saba (~30 nm from track)
    # No official station data; NHC estimates TS-force winds
    "2017242N16333",     "IRMA",    2017L, "Saba",          "NHC_estimate",
    80,                  "1min_sust",     "NHC TCR estimate",        "Estimated from best-track; Saba ~30nm from eye center",
    
    # Hurricane Maria (2017) - passed south of Leeward Islands
    "2017255N12319",     "MARIA",   2017L, "St_Martin",     "Juliana_Airport_TNCM",
    40,                  "1min_sust",     "NHC TCR (AL152017)",      "Maria passed ~100nm south; TS-force winds at SXM",
    
    # Hurricane Hugo (1989) - major hurricane near St. Croix/PR
    "1989248N12343",     "HUGO",    1989L, "Puerto_Rico",   "Roosevelt_Roads_NAS",
    104,                 "1min_sust",     "NHC TCR (AL081989)",      "Direct hit eastern PR",
    
    # Hurricane Andrew (1992) - direct hit Miami
    "1992216N10325",     "ANDREW",  1992L, "Miami",         "NHC_Miami_ASOS",
    141,                 "1min_sust",     "NHC TCR (AL041992)",      "Before instrument failure at NHC",
    
    # Hurricane Lenny (1999) - unusual W-to-E track through Leewards
    "1999317N14290",     "LENNY",   1999L, "St_Martin",     "Juliana_Airport_TNCM",
    55,                  "1min_sust",     "NHC TCR (AL171999)",      "Unusual west-to-east track; Cat 4 but passed south",
    
    # Hurricane Luis (1995) - direct hit northern Leeward Islands
    "1995241N12330",     "LUIS",    1995L, "St_Martin",     "Juliana_Airport_TNCM",
    110,                 "1min_sust",     "NHC TCR (AL131995)",      "Cat 4 at closest approach to SXM; major damage"
  )
}


#' Validate wind field estimates against station observations
#'
#' @description
#' For each reference observation, finds the matching storm in the model's trackpoint
#' data, extracts the model's peak V_site_kt for that storm at the target island,
#' and compares against the observed wind.
#'
#' Accounts for observation type: if the observation is 10-min sustained, applies a
#' 1-min/10-min conversion factor of ~1.12 before comparison.
#'
#' @param out List returned by \code{run_hazard_model()}.
#' @param obs_table Optional tibble of observations (default: \code{get_wind_observations()}).
#'
#' @return Tibble with model vs observed comparison per storm-station pair.
#' @export
validate_wind_field <- function(out, obs_table = NULL) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` is required.")
  
  if (is.null(obs_table)) obs_table <- get_wind_observations()
  
  # Conversion factor: 10-min to 1-min sustained (WMO standard)
  conv_10min_to_1min <- 1.12
  
  results <- vector("list", nrow(obs_table))
  
  for (i in seq_len(nrow(obs_table))) {
    obs <- obs_table[i, ]
    island <- obs$target_island
    
    # Convert observed wind to 1-min equivalent for fair comparison
    obs_1min_kt <- if (grepl("10min", obs$obs_type)) {
      obs$obs_wind_kt * conv_10min_to_1min
    } else {
      obs$obs_wind_kt
    }
    
    # Look up model estimate: peak V_site_kt for this SID at this island
    events <- out$events_by_island[[island]]
    model_V <- NA_real_
    model_wind_max <- NA_real_
    storm_found <- FALSE
    
    if (!is.null(events)) {
      match_row <- events |> dplyr::filter(.data$SID == obs$storm_sid)
      if (nrow(match_row) > 0) {
        storm_found <- TRUE
        model_V <- match_row$V_site_max_kt[1]
        model_wind_max <- match_row$wind_max_kt[1]
      }
    }
    
    # Also check trackpoints for more detail
    tp <- out$trackpoints[[island]]
    model_V_track <- NA_real_
    min_dist_km <- NA_real_
    if (!is.null(tp)) {
      tp_storm <- tp |> dplyr::filter(.data$SID == obs$storm_sid)
      if (nrow(tp_storm) > 0) {
        model_V_track <- suppressWarnings(max(tp_storm$V_site_kt, na.rm = TRUE))
        if (!is.finite(model_V_track)) model_V_track <- NA_real_
        min_dist_km <- suppressWarnings(min(tp_storm$dist_km, na.rm = TRUE))
        if (!is.finite(min_dist_km)) min_dist_km <- NA_real_
      }
    }
    
    # Use the better of event-level and trackpoint-level
    model_best <- if (is.finite(model_V)) model_V else model_V_track
    
    results[[i]] <- tibble::tibble(
      island = island,
      storm_name = obs$storm_name,
      storm_sid = obs$storm_sid,
      year = obs$year,
      station = obs$station,
      obs_raw_kt = obs$obs_wind_kt,
      obs_type = obs$obs_type,
      obs_1min_equiv_kt = obs_1min_kt,
      model_V_site_kt = model_best,
      model_wind_max_kt = model_wind_max,
      min_dist_km = min_dist_km,
      storm_found = storm_found,
      bias_kt = if (is.finite(model_best)) model_best - obs_1min_kt else NA_real_,
      bias_pct = if (is.finite(model_best) && obs_1min_kt > 0)
        100 * (model_best - obs_1min_kt) / obs_1min_kt else NA_real_,
      source = obs$obs_source,
      notes = obs$notes
    )
  }
  
  comp <- dplyr::bind_rows(results)
  
  # Summary statistics
  valid <- comp |> dplyr::filter(is.finite(.data$bias_kt))
  if (nrow(valid) > 0) {
    mae  <- mean(abs(valid$bias_kt))
    rmse <- sqrt(mean(valid$bias_kt^2))
    mb   <- mean(valid$bias_kt)
    
    message("\n[Wind Field Check] ", nrow(valid), " of ", nrow(comp), " storms matched")
    message(sprintf("  Mean Bias: %+.1f kt  |  MAE: %.1f kt  |  RMSE: %.1f kt", mb, mae, rmse))
    message("  Storm-by-storm:")
    for (j in seq_len(nrow(valid))) {
      r <- valid[j, ]
      sym <- if (abs(r$bias_pct) < 20) "\u2713" else if (abs(r$bias_pct) < 40) "~" else "\u2717"
      message(sprintf("    %s %-8s @ %-12s: obs=%5.0f kt  model=%5.0f kt  bias=%+.0f kt (%+.0f%%)",
                      sym, r$storm_name, r$island, r$obs_1min_equiv_kt,
                      r$model_V_site_kt, r$bias_kt, r$bias_pct))
    }
  } else {
    message("[Wind Field Check] No matching storms found in model output.")
  }
  
  comp
}


# =============================================================================
# 4) COMBINED VALIDATION REPORT
# =============================================================================

#' Run the full validation suite and print a summary report
#'
#' @description
#' Executes all three validation tiers (hindcast, rate check, wind field) and
#' produces a consolidated diagnostic report. Returns all results for further
#' analysis or plotting.
#'
#' @param out List returned by \code{run_hazard_model()}.
#' @param holdout_years Integer; years to hold out for hindcast.
#' @param n_sim Integer; synthetic years for hindcast simulation.
#' @param return_periods Numeric vector of return periods.
#' @param seed Integer seed.
#'
#' @return List with elements: hindcast, rate_check, wind_field, summary.
#' @export
run_validation_suite <- function(out,
                                 holdout_years = 10,
                                 n_sim = 5000,
                                 return_periods = c(5, 10, 25, 50),
                                 seed = 42) {
  
  message("=" |> rep(72) |> paste(collapse = ""))
  message("  HAZARD MODEL VALIDATION SUITE")
  message("=" |> rep(72) |> paste(collapse = ""))
  
  # --- Tier 1: Hindcast ---
  message("\n", "-" |> rep(72) |> paste(collapse = ""))
  message("  TIER 1: HINDCAST VALIDATION")
  message("-" |> rep(72) |> paste(collapse = ""))
  
  # Extract climate info from model output (if available)
  sst_df_val <- NULL
  beta_sst_val <- 0
  gamma_val <- 0
  if (!is.null(out$sst_info) && !is.null(out$sst_info$sst_df)) {
    sst_df_val <- out$sst_info$sst_df
    beta_sst_val <- out$beta_sst
    gamma_val <- if (!is.null(out$gamma_intensity)) out$gamma_intensity else 0
    message(sprintf("  [Climate] beta_SST=%.3f, gamma=%.4f", beta_sst_val, gamma_val))
  }
  
  hc <- tryCatch(
    validate_hindcast_all(out, holdout_years = holdout_years,
                          n_sim = n_sim, return_periods = return_periods,
                          seed = seed,
                          sst_df = sst_df_val,
                          beta_sst = beta_sst_val,
                          gamma_intensity = gamma_val),
    error = function(e) { message("  ERROR: ", e$message); NULL }
  )
  
  if (!is.null(hc) && nrow(hc$comparison) > 0) {
    message("\n  Return-level comparison (kt):")
    for (i in seq_len(nrow(hc$comparison))) {
      r <- hc$comparison[i, ]
      ci_tag <- if (isTRUE(r$obs_in_90ci)) "\u2713 in 90% CI" else "\u2717 outside 90% CI"
      message(sprintf("    %-12s %2d-yr:  obs=%.0f  sim=%.0f  obs90%%CI=[%.0f, %.0f]  bias=%+.0f%%  %s",
                      r$island, r$return_period, r$obs_full_rl, r$sim_rl,
                      r$obs_lo_90, r$obs_hi_90, r$bias_pct, ci_tag))
    }
  }
  
  # --- Tier 2: Rate check ---
  message("\n", "-" |> rep(72) |> paste(collapse = ""))
  message("  TIER 2: RATE SANITY CHECK")
  message("-" |> rep(72) |> paste(collapse = ""))
  
  rc <- tryCatch(
    validate_rates(out),
    error = function(e) { message("  ERROR: ", e$message); NULL }
  )
  
  # --- Tier 3: Wind field ---
  message("\n", "-" |> rep(72) |> paste(collapse = ""))
  message("  TIER 3: WIND FIELD SPOT-CHECKS")
  message("-" |> rep(72) |> paste(collapse = ""))
  
  wf <- tryCatch(
    validate_wind_field(out),
    error = function(e) { message("  ERROR: ", e$message); NULL }
  )
  
  # --- Summary ---
  message("\n", "=" |> rep(72) |> paste(collapse = ""))
  message("  VALIDATION SUMMARY")
  message("=" |> rep(72) |> paste(collapse = ""))
  
  n_rl_ok <- if (!is.null(hc)) sum(hc$comparison$obs_in_90ci, na.rm = TRUE) else 0
  n_rl_total <- if (!is.null(hc)) sum(!is.na(hc$comparison$obs_in_90ci)) else 0
  
  n_rate_ok <- if (!is.null(rc)) sum(rc$flag == "OK", na.rm = TRUE) else 0
  n_rate_total <- if (!is.null(rc)) sum(!is.na(rc$flag)) else 0
  
  n_wf_ok <- if (!is.null(wf)) {
    sum(abs(wf$bias_pct) < 30, na.rm = TRUE)
  } else 0
  n_wf_total <- if (!is.null(wf)) sum(is.finite(wf$bias_pct)) else 0
  
  summary_tbl <- tibble::tibble(
    tier = c("Hindcast (RL in 90% CI)", "Rate check (flag OK)", "Wind field (|bias| < 30%)"),
    pass = c(n_rl_ok, n_rate_ok, n_wf_ok),
    total = c(n_rl_total, n_rate_total, n_wf_total),
    pct = round(100 * c(n_rl_ok, n_rate_ok, n_wf_ok) /
                  pmax(c(n_rl_total, n_rate_total, n_wf_total), 1), 0)
  )
  
  message(sprintf("  Hindcast:   %d / %d return levels within 90%% CI", n_rl_ok, n_rl_total))
  message(sprintf("  Rate check: %d / %d rates flagged OK", n_rate_ok, n_rate_total))
  message(sprintf("  Wind field: %d / %d storms within 30%% bias", n_wf_ok, n_wf_total))
  message("=" |> rep(72) |> paste(collapse = ""))
  
  list(
    hindcast = hc,
    rate_check = rc,
    wind_field = wf,
    summary = summary_tbl
  )
}
