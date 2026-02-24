################################################################################
# hazard_validation.R
# Validation framework for the hurricane hazard model
#
# Three validation tiers:
#   1) Hindcast: hold-out last N years, compare simulated return levels vs observed
#   2) Rate sanity: compare model lambdas against published HURDAT2 climatologies
#   3) Wind field spot-checks: compare site wind against station observations
#
# User-facing API:
#   - make_validation_cfg(): typed config constructor (matches make_hazard_cfg())
#   - run_validation_suite(): all-in-one: run tiers + save plots/tables
#   - validate_hazard_model(): end-to-end: run model + validate
#   - plot_*() functions: standalone plotting (also called internally)
#
# All functions assume the hazard model source files are already loaded.
################################################################################


# Null-coalesce operator (available in base R >= 4.4; defined here for safety)
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}


# =============================================================================
# 0) VALIDATION CONFIGURATION
# =============================================================================

#' Create a validation configuration
#'
#' @description
#' Creates a typed configuration object for `run_validation_suite()`.
#' Follows the same pattern as `make_hazard_cfg()` and `make_sst_cfg()`:
#' common parameters up front, expert-only knobs in `advanced`.
#'
#' @param holdout_years Integer; number of years to hold out from the end of
#'   the historical record for train/test split (default: 10).
#' @param n_sim Integer; number of synthetic years to simulate for hindcast
#'   comparison (default: 5000).
#' @param return_periods Numeric vector of return periods (years) to compare
#'   (default: 5, 10, 25, 50).
#' @param seed Integer; random seed for reproducibility.
#' @param out_dir Character; output directory for saved plots and tables.
#' @param save_plots Logical; whether to save standard validation figures.
#' @param save_tables Logical; whether to save CSV + markdown tables.
#' @param advanced Optional named list of expert parameters. Most users should
#'   leave this as `NULL`. Supported names:
#'   \describe{
#'     \item{`xi_bounds`}{Numeric vector of length 2; allowed range for GEV
#'       shape parameter (default: `c(-0.3, 0.4)`).}
#'     \item{`n_boot`}{Integer; bootstrap replicates for return level CIs
#'       (default: 500).}
#'     \item{`base_size`}{Numeric; base font size for ggplot themes
#'       (default: 11).}
#'   }
#'
#' @return A list with class `c("validation_cfg", "list")`.
#' @export
#'
#' @examples
#' # Defaults — suitable for most users
#' val_cfg <- make_validation_cfg()
#'
#' # Custom holdout and output location
#' val_cfg <- make_validation_cfg(holdout_years = 15, out_dir = "results/val")
#'
#' # Expert tuning
#' val_cfg <- make_validation_cfg(
#'   n_sim = 10000,
#'   advanced = list(xi_bounds = c(-0.4, 0.5), n_boot = 1000, base_size = 13)
#' )
make_validation_cfg <- function(holdout_years  = 10L,
                                n_sim          = 5000L,
                                return_periods = c(5, 10, 25, 50),
                                seed           = 42L,
                                out_dir        = "output/validation",
                                save_plots     = TRUE,
                                save_tables    = TRUE,
                                advanced       = NULL) {

  defaults <- list(
    xi_bounds = c(-0.3, 0.4),
    n_boot    = 500L,
    base_size = 11
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

  # Input validation
  holdout_years  <- as.integer(holdout_years)
  n_sim          <- as.integer(n_sim)
  return_periods <- as.numeric(return_periods)
  seed           <- as.integer(seed)

  if (holdout_years < 1L) stop("holdout_years must be >= 1.", call. = FALSE)
  if (n_sim < 100L) stop("n_sim must be >= 100.", call. = FALSE)
  if (length(return_periods) == 0) stop("return_periods must have at least one value.", call. = FALSE)
  if (any(return_periods <= 1)) stop("return_periods must all be > 1.", call. = FALSE)

  cfg <- list(
    holdout_years  = holdout_years,
    n_sim          = n_sim,
    return_periods = return_periods,
    seed           = seed,
    out_dir        = as.character(out_dir),
    save_plots     = isTRUE(save_plots),
    save_tables    = isTRUE(save_tables),
    advanced       = advanced
  )
  class(cfg) <- c("validation_cfg", "list")
  cfg
}


#' @export
print.validation_cfg <- function(x, ...) {
  cat("Validation configuration\n")
  cat(sprintf("  Holdout       : %d years\n", x$holdout_years))
  cat(sprintf("  Simulation    : %s synthetic years\n",
              format(x$n_sim, big.mark = ",", scientific = FALSE, trim = TRUE)))
  cat(sprintf("  Return periods: %s yr\n", paste(x$return_periods, collapse = ", ")))
  cat(sprintf("  Seed          : %d\n", x$seed))
  cat(sprintf("  Output dir    : %s\n", x$out_dir))
  cat(sprintf("  Save plots    : %s\n", if (x$save_plots) "yes" else "no"))
  cat(sprintf("  Save tables   : %s\n", if (x$save_tables) "yes" else "no"))
  cat(sprintf("  GEV xi bounds : [%.2f, %.2f]\n",
              x$advanced$xi_bounds[1], x$advanced$xi_bounds[2]))
  invisible(x)
}


# =============================================================================
# 1) RETURN LEVEL COMPUTATION
# =============================================================================

#' Compute empirical return levels from a vector of annual maxima
#'
#' @description
#' Extracts return levels from a sample of annual maxima using the plotting-position
#' approach: sort values, assign return periods T = (n+1)/rank, interpolate to
#' requested return periods.
#'
#' @param annual_max Numeric vector of annual maximum values (e.g., peak_wind_kt).
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
# 2) KDE INTENSITY FITTING (internal)
# =============================================================================

.fit_intensity_kde <- function(pool, lower, upper = Inf, bw_mult = 1.0) {
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

  # Deduplicate CDF (approx() warns on tied x-values)
  dup <- duplicated(cdf_y)
  cdf_x <- cdf_x[!dup]
  cdf_y <- cdf_y[!dup]

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
#' @param fit List from `.fit_intensity_kde()`.
#' @param n Integer; number of draws.
#'
#' @return Numeric vector of n intensity draws within the interval
#'   `fit$lower` to `fit$upper`.
#' @keywords internal
.sample_intensity_kde <- function(fit, n) {
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
# 3) GEV FITTING (L-moments, no external packages)
# =============================================================================

#' Fit GEV distribution using L-moments (Hosking 1990)
#'
#' @description
#' Estimates GEV parameters (location \eqn{\mu}, scale \eqn{\sigma}, shape
#' \eqn{\xi}) using the method of L-moments. This is more robust than MLE for
#' small samples (n < 50) and requires no optimization. The sign convention
#' follows the standard meteorological form: \eqn{\xi > 0} = heavy tail
#' (Frechet), \eqn{\xi < 0} = bounded tail (Weibull), \eqn{\xi = 0} = Gumbel.
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
  c_val <- 2 / (3 + tau3) - log(2) / log(3)
  xi <- 7.8590 * c_val + 2.9554 * c_val^2

  # Clamp shape parameter to physical bounds
  xi <- max(xi_bounds[1], min(xi_bounds[2], xi))

  # Back-solve for sigma and mu
  if (abs(xi) < 1e-6) {
    # Gumbel case
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
#' @keywords internal
.qgev <- function(p, mu, sigma, xi) {
  if (abs(xi) < 1e-8) {
    mu - sigma * log(-log(p))
  } else {
    mu + sigma * ((-log(p))^(-xi) - 1) / xi
  }
}

#' GEV CDF
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
# 4) HURDLE-GEV RETURN LEVELS
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
# 5) PARAMETRIC BOOTSTRAP CIs FOR RETURN LEVELS
# =============================================================================

#' Compute return level CIs via parametric bootstrap
#'
#' @description
#' Generates bootstrap CIs by:
#' 1) Resampling the annual maxima
#' 2) Refitting the hurdle-GEV to each bootstrap sample
#' 3) Computing return levels from each fit
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


# =============================================================================
# 6) HINDCAST VALIDATION (internal workers)
# =============================================================================

#' Run hindcast validation for a single location
#' @keywords internal
.validate_hindcast <- function(events_island,
                               location,
                               holdout_years = 10,
                               n_sim = 5000,
                               return_periods = c(5, 10, 25, 50),
                               severities = c("TS", "HUR64plus"),
                               seed = 42,
                               sst_df = NULL,
                               beta_sst = 0,
                               gamma_intensity = 0,
                               xi_bounds = c(-0.3, 0.4),
                               n_boot = 500) {

  set.seed(seed)

  ev <- events_island |>
    dplyr::filter(.data$storm_class %in% c(severities, "none"),
                  is.finite(.data$peak_wind_kt))

  all_years <- sort(unique(ev$year))
  if (length(all_years) < holdout_years + 10) {
    stop("Insufficient years for holdout. Have ", length(all_years),
         ", need at least ", holdout_years + 10)
  }

  cutoff <- all_years[length(all_years) - holdout_years]
  train_years <- all_years[all_years <= cutoff]
  test_years  <- all_years[all_years > cutoff]

  message("[Hindcast] ", location, ": training ", min(train_years), "-", max(train_years),
          " (", length(train_years), " yr), testing ", min(test_years), "-",
          max(test_years), " (", length(test_years), " yr)")

  # --- Observed annual maxima (full record) ---
  obs_annual_max <- ev |>
    dplyr::group_by(.data$year) |>
    dplyr::summarise(V_max_kt = max(.data$peak_wind_kt, na.rm = TRUE),
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

  # --- FIT KDE INTENSITY DISTRIBUTIONS ---
  train_V_ts <- ev_train |>
    dplyr::filter(.data$storm_class == "TS") |>
    dplyr::pull(.data$peak_wind_kt) |>
    (\(x) x[is.finite(x)])()

  train_V_hur <- ev_train |>
    dplyr::filter(.data$storm_class == "HUR64plus") |>
    dplyr::pull(.data$peak_wind_kt) |>
    (\(x) x[is.finite(x)])()

  kde_ts  <- .fit_intensity_kde(train_V_ts,  lower = 34, upper = 64)
  # Tighter bandwidth for small hurricane samples
  hur_bw_mult <- if (length(train_V_hur) < 15) 0.7
  else if (length(train_V_hur) < 30) 0.85
  else 1.0
  kde_hur <- .fit_intensity_kde(train_V_hur, lower = 64, upper = 185, bw_mult = hur_bw_mult)

  n_ts_obs  <- length(train_V_ts)
  n_hur_obs <- length(train_V_hur)

  message(sprintf("  KDE fits: TS pool=%d events (mean=%.0f kt), HUR pool=%d events (mean=%.0f kt)",
                  n_ts_obs,  if (n_ts_obs > 0) mean(train_V_ts) else NA,
                  n_hur_obs, if (n_hur_obs > 0) mean(train_V_hur) else NA))

  # Fallback intensities if KDE can't be fit
  fallback_V <- list(TS = 45, HUR64plus = 85)

  # --- SIMULATE ANNUAL MAXIMA WITH KDE SAMPLING ---
  sst_anomaly_sim <- NULL
  if (!is.null(sst_df) && is.finite(beta_sst) && beta_sst != 0) {
    sst_train <- sst_df |>
      dplyr::filter(.data$year %in% train_years)
    if (nrow(sst_train) > 0 && "sst_anomaly" %in% names(sst_train)) {
      sst_pool <- sst_train$sst_anomaly[is.finite(sst_train$sst_anomaly)]
      if (length(sst_pool) > 0) {
        sst_anomaly_sim <- sample(sst_pool, n_sim, replace = TRUE)
        message(sprintf("  SST conditioning: \u03b2=%.3f, %d training SST values, mean \u0394SST=%+.2f\u00b0C",
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
    n_ts  <- sim_counts$n_ts[i]
    n_hur <- sim_counts$n_hur[i]

    winds <- numeric(0)

    if (n_ts > 0) {
      if (n_ts_obs >= 3) {
        winds <- c(winds, .sample_intensity_kde(kde_ts, n_ts))
      } else {
        winds <- c(winds, rep(fallback_V$TS, n_ts) +
                     stats::rnorm(n_ts, 0, 5))
      }
    }
    if (n_hur > 0) {
      if (n_hur_obs >= 3) {
        winds <- c(winds, .sample_intensity_kde(kde_hur, n_hur))
      } else {
        winds <- c(winds, rep(fallback_V$HUR64plus, n_hur) +
                     stats::rnorm(n_hur, 0, 10))
      }
    }

    if (length(winds) == 0) return(0)
    max(winds)
  }, numeric(1))

  # --- GEV RETURN LEVELS ---
  sim_gev <- compute_return_levels_gev(sim_annual_max, return_periods, xi_bounds)
  sim_rl  <- sim_gev$return_levels

  obs_full_max <- obs_annual_max$V_max_kt
  obs_gev <- compute_return_levels_gev(obs_full_max, return_periods, xi_bounds)
  obs_full_rl <- obs_gev$return_levels

  # Test period (empirical, for reference only)
  obs_test_max <- obs_annual_max |>
    dplyr::filter(.data$period == "test") |>
    dplyr::pull(.data$V_max_kt)
  obs_test_rl <- compute_return_levels(obs_test_max, return_periods)

  # --- BOOTSTRAP CIs ON OBSERVED RETURN LEVELS ---
  obs_rl_ci <- bootstrap_return_level_ci(
    obs_full_max, return_periods,
    n_boot = n_boot,
    xi_bounds = xi_bounds
  )

  # Model CIs (secondary diagnostic — should be narrow)
  sim_rl_ci <- bootstrap_return_level_ci(
    sim_annual_max, return_periods,
    n_boot = min(200L, n_boot),
    xi_bounds = xi_bounds
  )

  message(sprintf("  Model GEV: \u03bc=%.1f, \u03c3=%.1f, \u03be=%.3f (n_pos=%d, p0=%.2f)",
                  sim_gev$gev_fit$mu, sim_gev$gev_fit$sigma, sim_gev$gev_fit$xi,
                  sim_gev$n_nonzero, sim_gev$p_zero))
  message(sprintf("  Obs GEV:   \u03bc=%.1f, \u03c3=%.1f, \u03be=%.3f (n_pos=%d, p0=%.2f)",
                  obs_gev$gev_fit$mu, obs_gev$gev_fit$sigma, obs_gev$gev_fit$xi,
                  obs_gev$n_nonzero, obs_gev$p_zero))


  # --- Comparison table ---
  comparison <- tibble::tibble(
    location = location,
    return_period = return_periods,
    obs_full_rl = obs_full_rl,
    obs_test_rl = obs_test_rl,
    sim_rl = sim_rl,
    sim_median = sim_rl_ci$sim_median,
    sim_lo_90 = sim_rl_ci$sim_lo_90,
    sim_hi_90 = sim_rl_ci$sim_hi_90,
    obs_lo_90 = obs_rl_ci$sim_lo_90,
    obs_hi_90 = obs_rl_ci$sim_hi_90,
    obs_in_90ci = sim_rl >= obs_rl_ci$sim_lo_90 & sim_rl <= obs_rl_ci$sim_hi_90,
    model_in_obs_90ci = sim_rl >= obs_rl_ci$sim_lo_90 & sim_rl <= obs_rl_ci$sim_hi_90,
    obs_in_model_90ci = obs_full_rl >= sim_rl_ci$sim_lo_90 & obs_full_rl <= sim_rl_ci$sim_hi_90,
    bias_pct = 100 * (sim_rl - obs_full_rl) / pmax(obs_full_rl, 1)
  )

  list(
    location = location,
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


#' Run hindcast validation across all locations in a hazard model output
#' @keywords internal
.validate_hindcast_all <- function(out,
                                   holdout_years = 10,
                                   n_sim = 5000,
                                   return_periods = c(5, 10, 25, 50),
                                   seed = 42,
                                   sst_df = NULL,
                                   beta_sst = 0,
                                   gamma_intensity = 0,
                                   xi_bounds = c(-0.3, 0.4),
                                   n_boot = 500) {
  if (is.null(out$events)) stop("out$events is required.", call. = FALSE)
  locations <- sort(unique(out$events$location))
  results <- setNames(vector("list", length(locations)), locations)

  for (location in locations) {
    ev <- out$events |>
      dplyr::filter(.data$location == .env$location)
    if (is.null(ev) || nrow(ev) < 20) {
      message("[Hindcast] Skipping ", location, " (too few events: ",
              if (is.null(ev)) 0 else nrow(ev), ")")
      next
    }

    tryCatch({
      results[[location]] <- .validate_hindcast(
        events_island = ev,
        location = location,
        holdout_years = holdout_years,
        n_sim = n_sim,
        return_periods = return_periods,
        seed = seed,
        sst_df = sst_df,
        beta_sst = beta_sst,
        gamma_intensity = gamma_intensity,
        xi_bounds = xi_bounds,
        n_boot = n_boot
      )
    }, error = function(e) {
      message("[Hindcast] Error for ", location, ": ", e$message)
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
# 7) RATE SANITY CHECK
# =============================================================================

#' Reference HURDAT2/literature annual rates for Leeward Islands region
#'
#' @description
#' Returns a tibble of published annual TC passage rates from literature,
#' for comparison against model-fitted lambdas.
#'
#' @return Tibble with columns: region, storm_class, lambda_ref, source, gate_approx_nm, period.
#' @export
get_reference_rates <- function() {
  tibble::tribble(
    ~region,            ~storm_class,    ~lambda_ref, ~source,                           ~gate_approx_nm, ~period,
    "Leeward_Islands",  "TS34plus",   2.0,         "NHC Climo (100nm, 1970-2020)",    100,             "1970-2020",
    "Leeward_Islands",  "HUR64plus",  0.55,        "NHC Climo (100nm, 1970-2020)",    100,             "1970-2020",
    "St_Martin",        "TS34plus",   1.2,         "NOAA TC Climo (65nm, 1970-2023)", 65,              "1970-2023",
    "St_Martin",        "HUR64plus",  0.40,        "NOAA TC Climo (65nm, 1970-2023)", 65,              "1970-2023",
    "Puerto_Rico",      "TS34plus",   1.8,         "NHC Climo (100nm, 1970-2020)",    100,             "1970-2020",
    "Puerto_Rico",      "HUR64plus",  0.45,        "NHC Climo (100nm, 1970-2020)",    100,             "1970-2020",
    "Miami",            "TS34plus",   1.5,         "NHC Climo (100nm, 1970-2020)",    100,             "1970-2020",
    "Miami",            "HUR64plus",  0.35,        "NHC Climo (100nm, 1970-2020)",    100,             "1970-2020"
  )
}


#' Compare model-fitted rates against reference climatologies
#'
#' @description
#' Takes the lambda table from the hazard model output and compares per-location
#' per-storm_class rates against published references.
#'
#' @param out List returned by `run_hazard_model()`.
#' @param ref_rates Optional tibble of reference rates (default: `get_reference_rates()`).
#'
#' @return Tibble with model vs reference rate comparison.
#' @export
validate_rates <- function(out, ref_rates = NULL) {

  if (is.null(ref_rates)) ref_rates <- get_reference_rates()

  model_rates <- out$rates |>
    dplyr::select(location = "location", storm_class = "storm_class",
                  lambda_model = "lambda", n_years_model = "n_years")

  model_rates2 <- model_rates |>
    tidyr::pivot_wider(names_from = storm_class, values_from = lambda_model, values_fill = 0) |>
    dplyr::mutate(TS34plus = TS + HUR64plus) |>
    tidyr::pivot_longer(c("TS34plus","HUR64plus"), names_to="storm_class", values_to="lambda_model")

  island_to_region <- tibble::tribble(
    ~location,        ~region,
    "St_Martin",    "St_Martin",
    "Saba",         "Leeward_Islands",
    "Statia",       "Leeward_Islands",
    "Puerto_Rico",  "Puerto_Rico",
    "Miami",        "Miami"
  )

  comp <- model_rates2 |>
    dplyr::left_join(island_to_region, by = "location") |>
    dplyr::left_join(ref_rates, by = c("region", "storm_class")) |>
    dplyr::mutate(
      ratio     = .data$lambda_model / pmax(.data$lambda_ref, 0.001),

      expected_ratio = dplyr::case_when(
        .data$storm_class == "TS34plus"  & .data$gate_approx_nm >= 100 ~ 0.55,
        .data$storm_class == "TS34plus"  & .data$gate_approx_nm < 100  ~ 0.75,
        .data$storm_class == "HUR64plus" & .data$gate_approx_nm >= 100 ~ 0.30,
        .data$storm_class == "HUR64plus" & .data$gate_approx_nm < 100  ~ 0.45,
        TRUE ~ 0.50
      ),

      adj_ratio = .data$ratio / .data$expected_ratio,

      flag = dplyr::case_when(
        is.na(.data$lambda_ref)   ~ "no_reference",
        .data$adj_ratio > 2.5     ~ "HIGH: model >> expected",
        .data$adj_ratio < 0.4     ~ "LOW: model << expected",
        .data$adj_ratio > 1.8     ~ "elevated",
        .data$adj_ratio < 0.6     ~ "slightly_low",
        TRUE ~ "OK"
      )
    ) |>
    dplyr::select(
      "location", "storm_class",
      "lambda_model", "n_years_model",
      "lambda_ref", "source", "gate_approx_nm", "period",
      "expected_ratio", "ratio", "adj_ratio", "flag"
    )

  message("\n[Rate Check] Summary:")
  for (i in seq_len(nrow(comp))) {
    r <- comp[i, ]
    sym <- if (r$flag == "OK") "\u2713" else if (grepl("^(HIGH|LOW)", r$flag)) "\u2717" else "~"
    message(sprintf("  %s %s / %-10s : model=%.3f  ref=%.3f  raw_ratio=%.2f  exp_ratio=%.2f  adj_ratio=%.2f  [%s]",
                    sym, r$location, r$storm_class, r$lambda_model,
                    if (is.na(r$lambda_ref)) NA else r$lambda_ref,
                    if (is.na(r$ratio)) NA else r$ratio,
                    if (is.na(r$expected_ratio)) NA else r$expected_ratio,
                    if (is.na(r$adj_ratio)) NA else r$adj_ratio,
                    r$flag))
  }

  comp
}


# =============================================================================
# 8) WIND FIELD SPOT-CHECKS
# =============================================================================

#' Reference observations for wind field validation
#'
#' @description
#' Returns a tibble of known station/buoy observations during well-documented storms
#' for comparison against model-estimated site wind.
#'
#' @return Tibble with columns: storm_sid, storm_name, year, target_island,
#'   station, obs_wind_kt, obs_type, obs_source, notes.
#' @export
get_wind_observations <- function() {
  tibble::tribble(
    ~storm_sid,          ~storm_name, ~year, ~target_island, ~station,
    ~obs_wind_kt, ~obs_type,       ~obs_source,               ~notes,

    "2017242N16333",     "IRMA",    2017L, "St_Martin",     "Juliana_Airport_TNCM",
    155,                 "1min_sust",     "NHC TCR (AL112017)",      "Instruments failed; 155 kt is NHC best-track Vmax at landfall",

    "2017242N16333",     "IRMA",    2017L, "St_Martin",     "Meteo_France_SXM",
    113,                 "10min_sust",    "Meteo France RSMC report", "Last valid reading before instrument failure; 10-min sustained",

    "2014279N15323",     "GONZALO", 2014L, "St_Martin",     "Juliana_Airport_TNCM",
    60,                  "1min_sust",     "NHC TCR (AL082014)",      "Marginal TS winds at SXM; storm passed ~50nm north",

    "2017242N16333",     "IRMA",    2017L, "Saba",          "NHC_estimate",
    80,                  "1min_sust",     "NHC TCR estimate",        "Estimated from best-track; Saba ~30nm from eye center",

    "2017255N12319",     "MARIA",   2017L, "St_Martin",     "Juliana_Airport_TNCM",
    40,                  "1min_sust",     "NHC TCR (AL152017)",      "Maria passed ~100nm south; TS-force winds at SXM",

    "1989248N12343",     "HUGO",    1989L, "Puerto_Rico",   "Roosevelt_Roads_NAS",
    104,                 "1min_sust",     "NHC TCR (AL081989)",      "Direct hit eastern PR",

    "1992216N10325",     "ANDREW",  1992L, "Miami",         "NHC_Miami_ASOS",
    141,                 "1min_sust",     "NHC TCR (AL041992)",      "Before instrument failure at NHC",

    "1999317N14290",     "LENNY",   1999L, "St_Martin",     "Juliana_Airport_TNCM",
    55,                  "1min_sust",     "NHC TCR (AL171999)",      "Unusual west-to-east track; Cat 4 but passed south",

    "1995241N12330",     "LUIS",    1995L, "St_Martin",     "Juliana_Airport_TNCM",
    110,                 "1min_sust",     "NHC TCR (AL131995)",      "Cat 4 at closest approach to SXM; major damage"
  )
}


#' Validate wind field estimates against station observations
#'
#' @description
#' For each reference observation, finds the matching storm in the model's trackpoint
#' data, extracts the model's peak site wind for that storm at the target location,
#' and compares against the observed wind.
#'
#' @param out List returned by `run_hazard_model()`.
#' @param obs_table Optional tibble of observations (default: `get_wind_observations()`).
#'
#' @return Tibble with model vs observed comparison per storm-station pair.
#' @export
validate_wind_field <- function(out, obs_table = NULL) {

  if (is.null(obs_table)) obs_table <- get_wind_observations()

  # Conversion factor: 10-min to 1-min sustained (WMO standard)
  conv_10min_to_1min <- 1.12

  results <- vector("list", nrow(obs_table))

  for (i in seq_len(nrow(obs_table))) {
    obs <- obs_table[i, ]
    location <- obs$target_island

    # Convert observed wind to 1-min equivalent for fair comparison
    obs_1min_kt <- if (grepl("10min", obs$obs_type)) {
      obs$obs_wind_kt * conv_10min_to_1min
    } else {
      obs$obs_wind_kt
    }

    # Look up model estimate
    events <- out$events |>
      dplyr::filter(.data$location == .env$location)
    model_V <- NA_real_
    model_wind_max <- NA_real_
    storm_found <- FALSE

    if (!is.null(events)) {
      match_row <- events |> dplyr::filter(.data$storm_id == obs$storm_sid)
      if (nrow(match_row) > 0) {
        storm_found <- TRUE
        model_V <- match_row$peak_wind_kt[1]
        model_wind_max <- match_row$storm_intensity_kt[1]
      }
    }

    # Also check trackpoints for more detail
    tp <- out$trackpoints[[location]]
    model_V_track <- NA_real_
    min_dist_km <- NA_real_
    if (!is.null(tp)) {
      tp_storm <- tp |> dplyr::filter(.data$SID == obs$storm_sid)
      if (nrow(tp_storm) > 0) {
        model_V_track <- suppressWarnings(max(tp_storm$site_wind_kt, na.rm = TRUE))
        if (!is.finite(model_V_track)) model_V_track <- NA_real_
        min_dist_km <- suppressWarnings(min(tp_storm$dist_km, na.rm = TRUE))
        if (!is.finite(min_dist_km)) min_dist_km <- NA_real_
      }
    }

    model_best <- if (is.finite(model_V)) model_V else model_V_track

    results[[i]] <- tibble::tibble(
      location = location,
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
                      sym, r$storm_name, r$location, r$obs_1min_equiv_kt,
                      r$model_V_site_kt, r$bias_kt, r$bias_pct))
    }
  } else {
    message("[Wind Field Check] No matching storms found in model output.")
  }

  comp
}


# =============================================================================
# 9) COMBINED VALIDATION SUITE (all-in-one)
# =============================================================================

#' Run the full validation suite
#'
#' @description
#' Executes all three validation tiers (hindcast, rate check, wind field),
#' produces a consolidated diagnostic report, and optionally saves all plots
#' and tables to the configured output directory.
#'
#' @param out List returned by `run_hazard_model()`.
#' @param cfg Validation configuration from `make_validation_cfg()`. If `NULL`,
#'   default configuration is used.
#'
#' @return List with elements: `hindcast`, `rate_check`, `wind_field`, `summary`,
#'   and `artifacts` (paths of saved plots/tables, if any).
#' @export
run_validation_suite <- function(out, cfg = make_validation_cfg()) {
  if (!inherits(cfg, "validation_cfg")) {
    stop("cfg must be created by make_validation_cfg().", call. = FALSE)
  }

  holdout_years  <- cfg$holdout_years
  n_sim          <- cfg$n_sim
  return_periods <- cfg$return_periods
  seed           <- cfg$seed
  xi_bounds      <- cfg$advanced$xi_bounds
  n_boot         <- cfg$advanced$n_boot
  base_size      <- cfg$advanced$base_size
  out_dir        <- cfg$out_dir

  message("=" |> rep(72) |> paste(collapse = ""))
  message("  HAZARD MODEL VALIDATION SUITE")
  message("=" |> rep(72) |> paste(collapse = ""))

  # --- Extract climate info from model output (if available) ---
  sst_df_val <- NULL
  beta_sst_val <- 0
  gamma_val <- 0
  if (!is.null(out$fit) && nrow(out$fit) > 0) {
    beta_sst_val <- out$fit$beta_sst[1]
    gamma_val <- out$fit$gamma_intensity[1]
    message(sprintf("  [Climate] beta_SST=%.3f, gamma=%.4f", beta_sst_val, gamma_val))
  }

  # --- Tier 1: Hindcast ---
  message("\n", "-" |> rep(72) |> paste(collapse = ""))
  message("  TIER 1: HINDCAST VALIDATION")
  message("-" |> rep(72) |> paste(collapse = ""))

  hc <- tryCatch(
    .validate_hindcast_all(out,
                           holdout_years = holdout_years,
                           n_sim = n_sim,
                           return_periods = return_periods,
                           seed = seed,
                           sst_df = sst_df_val,
                           beta_sst = beta_sst_val,
                           gamma_intensity = gamma_val,
                           xi_bounds = xi_bounds,
                           n_boot = n_boot),
    error = function(e) { message("  ERROR: ", e$message); NULL }
  )

  if (!is.null(hc) && nrow(hc$comparison) > 0) {
    message("\n  Return-level comparison (kt):")
    for (i in seq_len(nrow(hc$comparison))) {
      r <- hc$comparison[i, ]
      ci_tag <- if (isTRUE(r$obs_in_90ci)) "\u2713 in 90% CI" else "\u2717 outside 90% CI"
      message(sprintf("    %-12s %2d-yr:  obs=%.0f  sim=%.0f  obs90%%CI=[%.0f, %.0f]  bias=%+.0f%%  %s",
                      r$location, r$return_period, r$obs_full_rl, r$sim_rl,
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

  val <- list(
    hindcast = hc,
    rate_check = rc,
    wind_field = wf,
    summary = summary_tbl
  )

  # --- Save artifacts ---
  artifacts <- list(plots = list(), tables = list())

  if (isTRUE(cfg$save_plots)) {
    .validate_dir_create(out_dir)
    artifacts$plots$hindcast        <- plot_hindcast_validation(val, out_dir = out_dir, base_size = base_size)
    artifacts$plots$rate_check      <- plot_rate_validation(val, out_dir = out_dir, base_size = base_size)
    artifacts$plots$wind_field      <- plot_wind_field_validation(val, out = out, out_dir = out_dir, base_size = base_size)
    artifacts$plots$bias_diagnostics <- plot_bias_diagnostics(val, out_dir = out_dir, base_size = base_size)
    artifacts$plots$qq_plots        <- plot_qq_validation(val, out_dir = out_dir, base_size = base_size)
    artifacts$plots$cdf_comparison  <- plot_cdf_comparison(val, out_dir = out_dir, base_size = base_size)
  }

  if (isTRUE(cfg$save_tables)) {
    .validate_dir_create(out_dir)
    tables <- list(
      `Hindcast Return Levels` = val$hindcast$comparison,
      `Rate Comparison`        = val$rate_check,
      `Wind Field Spot-Checks` = val$wind_field,
      `Summary`                = val$summary
    )

    if (!is.null(val$hindcast$comparison)) {
      artifacts$tables$hindcast_csv <- file.path(out_dir, "hindcast_return_levels.csv")
      .validate_write_csv(val$hindcast$comparison, artifacts$tables$hindcast_csv)
    }
    if (!is.null(val$rate_check)) {
      artifacts$tables$rate_check_csv <- file.path(out_dir, "rate_check.csv")
      .validate_write_csv(val$rate_check, artifacts$tables$rate_check_csv)
    }
    if (!is.null(val$wind_field)) {
      artifacts$tables$wind_field_csv <- file.path(out_dir, "wind_field.csv")
      .validate_write_csv(val$wind_field, artifacts$tables$wind_field_csv)
    }
    if (!is.null(val$summary)) {
      artifacts$tables$summary_csv <- file.path(out_dir, "validation_summary.csv")
      .validate_write_csv(val$summary, artifacts$tables$summary_csv)
    }

    artifacts$tables$tables_md <- file.path(out_dir, "validation_tables.md")
    .validate_write_md_tables(tables, artifacts$tables$tables_md)
  }

  val$artifacts <- artifacts
  val
}


# =============================================================================
# 10) END-TO-END CONVENIENCE WRAPPER
# =============================================================================

#' Validate the hazard model end-to-end
#'
#' @description
#' Convenience wrapper that runs the hazard model, executes the validation suite,
#' and saves all standard figures and tables to the configured directory.
#'
#' @param cfg Hazard configuration from `make_hazard_cfg()`.
#' @param targets Target locations tibble/data.frame.
#' @param validation_cfg Validation configuration from `make_validation_cfg()`.
#' @param severities Character vector of storm classes to simulate.
#' @param sst_cfg Optional SST configuration from `make_sst_cfg()`.
#'
#' @return A list with elements `out`, `val`, and `artifacts` (saved paths).
#' @export
validate_hazard_model <- function(cfg,
                                  targets,
                                  validation_cfg = make_validation_cfg(),
                                  severities = c("TS", "HUR64plus"),
                                  sst_cfg = NULL) {
  if (!inherits(validation_cfg, "validation_cfg")) {
    stop("validation_cfg must be created by make_validation_cfg().", call. = FALSE)
  }

  out <- run_hazard_model(
    cfg = cfg,
    targets = targets,
    severities = severities,
    sst_cfg = sst_cfg
  )

  val <- run_validation_suite(out = out, cfg = validation_cfg)

  list(out = out, val = val, artifacts = val$artifacts)
}


# =============================================================================
# 11) PRIVATE HELPERS (plotting infrastructure)
# =============================================================================

.validate_theme <- function(base_size = 11) {
  passaat_theme(base_size = base_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = base_size + 1))
}

.resolve_plot_cfg <- function(cfg, out_dir, base_size) {
  out_dir <- out_dir %||% (if (!is.null(cfg)) cfg$out_dir else "output/validation")
  base_size <- base_size %||% (if (!is.null(cfg)) cfg$advanced$base_size else 11)
  list(
    out_dir = out_dir,
    base_size = base_size,
    theme = .validate_theme(base_size = base_size)
  )
}

.validate_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

.validate_write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
  invisible(path)
}

.validate_write_md_tables <- function(tables, path) {
  if (!requireNamespace("knitr", quietly = TRUE)) return(invisible(NULL))
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  for (nm in names(tables)) {
    x <- tables[[nm]]
    if (is.null(x)) next
    writeLines(paste0("## ", nm, "\n"), con)
    writeLines(knitr::kable(x, format = "pipe"), con)
    writeLines("\n\n", con)
  }
  invisible(path)
}

.validate_save_plot <- function(p, path, width, height, dpi = 150) {
  if (is.null(p)) return(invisible(NULL))
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  ggplot2::ggsave(filename = path, plot = p, width = width, height = height, dpi = dpi)
  invisible(path)
}


# =============================================================================
# 12) PLOT FUNCTIONS
# =============================================================================

#' Plot hindcast validation figures
#'
#' @description
#' Creates the main hindcast validation figures (return-level comparison and
#' per-location annual-max distribution plots).
#'
#' @param val Output from `run_validation_suite()`.
#' @param cfg Optional `validation_cfg` object. When provided, `out_dir` and
#'   `base_size` are read from the config (explicit arguments override).
#' @param out_dir Directory to save plots in. Overrides `cfg$out_dir`.
#' @param base_size Base font size for plots. Overrides `cfg$advanced$base_size`.
#'
#' @return Named character vector of saved plot paths (invisibly) or `NULL`.
#' @export
plot_hindcast_validation <- function(val,
                                     cfg = NULL,
                                     out_dir = NULL,
                                     base_size = NULL) {
  plot_cfg <- .resolve_plot_cfg(cfg = cfg, out_dir = out_dir, base_size = base_size)
  out_dir <- plot_cfg$out_dir
  ggtheme <- plot_cfg$theme

  if (is.null(val$hindcast) || is.null(val$hindcast$comparison) || nrow(val$hindcast$comparison) == 0) {
    return(invisible(NULL))
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))

  .validate_dir_create(out_dir)

  comp <- val$hindcast$comparison
  p_rl <- ggplot2::ggplot(comp, ggplot2::aes(x = factor(return_period))) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = obs_lo_90, ymax = obs_hi_90),
      width = 0.3, color = "red", linewidth = 0.8, alpha = 0.7
    ) +
    ggplot2::geom_point(ggplot2::aes(y = obs_full_rl), size = 3, color = "red", shape = 17) +
    ggplot2::geom_point(ggplot2::aes(y = sim_rl), size = 3, color = "steelblue", shape = 16) +
    ggplot2::geom_hline(yintercept = 64, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    ggplot2::facet_wrap(~ location, scales = "free_y", ncol = 3) +
    ggplot2::labs(
      x = "Return period (years)",
      y = "Return level \u2014 peak site wind (kt)",
      title = "Hindcast Validation: Simulated vs Observed Return Levels",
      subtitle = "Blue dot = model; Red triangle [90% CI] = observed (full record); dashed = 64 kt HUR threshold"
    ) +
    ggtheme

  paths <- character(0)
  paths["hindcast_return_levels"] <- file.path(out_dir, "hindcast_return_levels.png")
  .validate_save_plot(p_rl, paths[["hindcast_return_levels"]], width = 12, height = 7, dpi = 150)

  # Per-location distribution plots
  if (!is.null(val$hindcast$per_island)) {
    for (isl in names(val$hindcast$per_island)) {
      hc_isl <- val$hindcast$per_island[[isl]]
      if (is.null(hc_isl)) next

      obs_df <- hc_isl$obs_annual_max
      obs_df <- obs_df[is.finite(obs_df$V_max_kt) & obs_df$V_max_kt > 0, , drop = FALSE]

      sim_df <- tibble::tibble(V_max_kt = hc_isl$sim_annual_max)
      sim_df <- sim_df[is.finite(sim_df$V_max_kt) & sim_df$V_max_kt > 0, , drop = FALSE]


      p_dist <- ggplot2::ggplot() +
        ggplot2::geom_histogram(
          data = sim_df,
          ggplot2::aes(x = V_max_kt, y = ggplot2::after_stat(density)),
          fill = "steelblue", alpha = 0.4, bins = 40
        ) +
        ggplot2::geom_density(data = sim_df, ggplot2::aes(x = V_max_kt), color = "steelblue", linewidth = 0.8) +
        ggplot2::geom_rug(data = obs_df, ggplot2::aes(x = V_max_kt, color = period), linewidth = 0.8, alpha = 0.8) +
        ggplot2::scale_color_manual(values = c(train = "grey40", test = "red")) +
        ggplot2::geom_vline(xintercept = 64, linetype = "dashed", color = "grey50") +
        ggplot2::labs(
          x = "Annual maximum site wind (kt)",
          y = "Density",
          title = paste("Annual Max Distribution:", isl),
          subtitle = "Blue = simulated; rug ticks = observed (grey=train, red=test)",
          color = "Observed"
        ) +
        ggtheme

      nm <- paste0("hindcast_dist_", tolower(isl))
      paths[[nm]] <- file.path(out_dir, paste0(nm, ".png"))
      .validate_save_plot(p_dist, paths[[nm]], width = 8, height = 5, dpi = 150)
    }
  }

  invisible(paths)
}

#' Plot rate-check validation figure
#'
#' @param val Output from `run_validation_suite()`.
#' @param cfg Optional `validation_cfg` object.
#' @param out_dir Directory to save plots in. Overrides `cfg$out_dir`.
#' @param base_size Base font size for plots. Overrides `cfg$advanced$base_size`.
#'
#' @return Path of saved plot (invisibly) or `NULL`.
#' @export
plot_rate_validation <- function(val,
                                 cfg = NULL,
                                 out_dir = NULL,
                                 base_size = NULL) {
  plot_cfg <- .resolve_plot_cfg(cfg = cfg, out_dir = out_dir, base_size = base_size)
  out_dir <- plot_cfg$out_dir
  ggtheme <- plot_cfg$theme

  if (is.null(val$rate_check) || nrow(val$rate_check) == 0) return(invisible(NULL))
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))

  .validate_dir_create(out_dir)

  rc <- dplyr::filter(val$rate_check, !is.na(lambda_ref))
  if (nrow(rc) == 0) return(invisible(NULL))

  p_rate <- ggplot2::ggplot(rc, ggplot2::aes(x = lambda_ref, y = lambda_model)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_abline(slope = 2, intercept = 0, linetype = "dotted", color = "grey70") +
    ggplot2::geom_abline(slope = 0.5, intercept = 0, linetype = "dotted", color = "grey70") +
    ggplot2::geom_point(ggplot2::aes(color = flag, shape = storm_class), size = 3.5) +
    ggplot2::geom_text(ggplot2::aes(label = location), nudge_y = 0.08, size = 3, check_overlap = TRUE) +
    ggplot2::scale_color_manual(values = c(
      "OK" = "forestgreen",
      "elevated (expected: model gate > ref radius)" = "orange",
      "HIGH: model >> reference (check search_radius_km or start_year)" = "red",
      "LOW: model << reference (check storm_class filter or data)" = "purple"
    )) +
    ggplot2::coord_equal(
      xlim = c(0, max(c(rc$lambda_model, rc$lambda_ref), na.rm = TRUE) + 0.3),
      ylim = c(0, max(c(rc$lambda_model, rc$lambda_ref), na.rm = TRUE) + 0.3)
    ) +
    ggplot2::labs(
      x = "Reference \u03bb (published climatology)",
      y = "Model \u03bb (fitted)",
      title = "Rate Sanity Check: Model vs Published Annual Rates",
      subtitle = "Dashed = 1:1; dotted = 0.5x and 2x bounds",
      color = "Flag", shape = "storm_class"
    ) +
    ggtheme

  path <- file.path(out_dir, "rate_comparison.png")
  .validate_save_plot(p_rate, path, width = 8, height = 6, dpi = 150)
  invisible(path)
}

#' Plot wind-field validation figures
#'
#' @param val Output from `run_validation_suite()`.
#' @param out Model output list from `run_hazard_model()`. Only needed for the
#'   Irma distance-profile plot.
#' @param cfg Optional `validation_cfg` object.
#' @param out_dir Directory to save plots in. Overrides `cfg$out_dir`.
#' @param base_size Base font size for plots. Overrides `cfg$advanced$base_size`.
#'
#' @return Named character vector of saved plot paths (invisibly) or `NULL`.
#' @export
plot_wind_field_validation <- function(val,
                                       out = NULL,
                                       cfg = NULL,
                                       out_dir = NULL,
                                       base_size = NULL) {
  plot_cfg <- .resolve_plot_cfg(cfg = cfg, out_dir = out_dir, base_size = base_size)
  out_dir <- plot_cfg$out_dir
  ggtheme <- plot_cfg$theme

  if (is.null(val$wind_field) || !any(is.finite(val$wind_field$model_V_site_kt))) {
    return(invisible(NULL))
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))

  .validate_dir_create(out_dir)

  wf <- dplyr::filter(val$wind_field, is.finite(model_V_site_kt))
  if (nrow(wf) == 0) return(invisible(NULL))

  p_wf <- ggplot2::ggplot(wf, ggplot2::aes(x = obs_1min_equiv_kt, y = model_V_site_kt)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_abline(slope = 1, intercept = 20, linetype = "dotted", color = "grey70") +
    ggplot2::geom_abline(slope = 1, intercept = -20, linetype = "dotted", color = "grey70") +
    ggplot2::geom_point(ggplot2::aes(color = location, shape = location), size = 3.5) +
    ggplot2::geom_text(ggplot2::aes(label = storm_name), nudge_y = 4, size = 3, check_overlap = TRUE) +
    ggplot2::coord_equal(
      xlim = c(0, max(c(wf$obs_1min_equiv_kt, wf$model_V_site_kt), na.rm = TRUE) + 20),
      ylim = c(0, max(c(wf$obs_1min_equiv_kt, wf$model_V_site_kt), na.rm = TRUE) + 20)
    ) +
    ggplot2::labs(
      x = "Observed wind (1-min equiv, kt)",
      y = "Model V_site_kt (kt)",
      title = "Wind Field Validation: Model vs Station Observations",
      subtitle = "Dashed = 1:1; dotted = \u00b120 kt bounds",
      color = "location", shape = "location"
    ) +
    ggtheme

  paths <- character(0)
  paths[["wind_field_scatter"]] <- file.path(out_dir, "wind_field_scatter.png")
  .validate_save_plot(p_wf, paths[["wind_field_scatter"]], width = 8, height = 6, dpi = 150)

  # Optional per-storm distance profile for Irma at St. Martin
  if (!is.null(out) && !is.null(out$trackpoints$St_Martin)) {
    irma_tp <- dplyr::filter(out$trackpoints$St_Martin, SID == "2017242N16333")
    wind_col <- if ("site_wind_kt" %in% names(irma_tp)) {
      "site_wind_kt"
    } else if ("V_site_kt" %in% names(irma_tp)) {
      "V_site_kt"
    } else {
      NULL
    }
    if (!is.null(wind_col)) {
      irma_tp <- dplyr::filter(irma_tp, is.finite(.data[[wind_col]]), is.finite(dist_km))
    } else {
      irma_tp <- irma_tp[0, , drop = FALSE]
    }

    if (nrow(irma_tp) > 0) {
      p_irma <- ggplot2::ggplot(irma_tp, ggplot2::aes(x = dist_km, y = .data[[wind_col]])) +
        ggplot2::geom_line(color = "steelblue", linewidth = 0.8) +
        ggplot2::geom_point(color = "steelblue", size = 1.5) +
        ggplot2::geom_hline(yintercept = c(34, 64), linetype = "dashed", color = c("orange", "red")) +
        ggplot2::geom_hline(yintercept = 155, linetype = "solid", color = "red", alpha = 0.5) +
        ggplot2::annotate(
          "text",
          x = max(irma_tp$dist_km) * 0.7,
          y = 155,
          label = "NHC best-track Vmax = 155 kt",
          size = 3,
          color = "red"
        ) +
        ggplot2::labs(
          x = "Distance from St. Martin (km)",
          y = "Model site wind (kt)",
          title = "Hurricane Irma (2017): Wind Profile at St. Martin",
          subtitle = "Modeled site wind vs distance at each 6-hourly track point"
        ) +
        ggtheme

      paths[["irma_stmartin_profile"]] <- file.path(out_dir, "irma_stmartin_profile.png")
      .validate_save_plot(p_irma, paths[["irma_stmartin_profile"]], width = 8, height = 5, dpi = 150)
    }
  }

  invisible(paths)
}


#' Plot bias decomposition: frequency vs intensity contributions
#'
#' @param val Output from `run_validation_suite()`.
#' @param cfg Optional `validation_cfg` object.
#' @param out_dir Directory to save plots in. Overrides `cfg$out_dir`.
#' @param base_size Base font size for plots. Overrides `cfg$advanced$base_size`.
#'
#' @return Named character vector of saved plot paths (invisibly) or `NULL`.
#' @export
plot_bias_diagnostics <- function(val,
                                  cfg = NULL,
                                  out_dir = NULL,
                                  base_size = NULL) {
  plot_cfg <- .resolve_plot_cfg(cfg = cfg, out_dir = out_dir, base_size = base_size)
  out_dir <- plot_cfg$out_dir
  ggtheme <- plot_cfg$theme

  if (is.null(val$hindcast) || is.null(val$hindcast$per_island)) return(invisible(NULL))
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))

  .validate_dir_create(out_dir)
  paths <- character(0)

  # --- Collect per-location bias decomposition ---
  bias_rows <- list()
  for (isl in names(val$hindcast$per_island)) {
    hc <- val$hindcast$per_island[[isl]]
    if (is.null(hc)) next

    obs_am <- hc$obs_annual_max$V_max_kt
    sim_am <- hc$sim_annual_max

    obs_freq <- mean(obs_am > 0)
    sim_freq <- mean(sim_am > 0)
    obs_int <- mean(obs_am[obs_am > 0])
    sim_int <- mean(sim_am[sim_am > 0])
    obs_mean <- mean(obs_am)
    sim_mean <- mean(sim_am)

    # Decomposition:  E[max] = P(any event) * E[max | event]
    freq_contrib <- (sim_freq - obs_freq) * obs_int
    int_contrib  <- obs_freq * (sim_int - obs_int)
    interact     <- (sim_freq - obs_freq) * (sim_int - obs_int)

    bias_rows[[isl]] <- tibble::tibble(
      location = isl,
      obs_event_rate  = obs_freq,
      sim_event_rate  = sim_freq,
      freq_bias_pct   = 100 * (sim_freq - obs_freq) / pmax(obs_freq, 0.01),
      obs_mean_int_kt = obs_int,
      sim_mean_int_kt = sim_int,
      int_bias_pct    = 100 * (sim_int - obs_int) / pmax(obs_int, 1),
      total_bias_kt   = sim_mean - obs_mean,
      freq_contrib_kt = freq_contrib,
      int_contrib_kt  = int_contrib,
      interact_kt     = interact
    )
  }

  if (length(bias_rows) == 0) return(invisible(NULL))
  bias_df <- dplyr::bind_rows(bias_rows)

  # --- Panel 1: Stacked bar chart of bias contributions ---
  bias_long <- bias_df |>
    dplyr::select("location", Frequency = "freq_contrib_kt",
                  Intensity = "int_contrib_kt", Interaction = "interact_kt") |>
    tidyr::pivot_longer(-"location", names_to = "source", values_to = "bias_kt")

  p_decomp <- ggplot2::ggplot(bias_long, ggplot2::aes(x = location, y = bias_kt, fill = source)) +
    ggplot2::geom_col(position = "stack", width = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = c(
      Frequency = "#E69F00", Intensity = "#56B4E9", Interaction = "#999999"
    )) +
    ggplot2::labs(
      x = NULL, y = "Bias contribution (kt)",
      title = "Bias Decomposition: Frequency vs Intensity",
      subtitle = "Positive = model overestimates; Negative = model underestimates",
      fill = "Source"
    ) +
    ggtheme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

  paths[["bias_decomposition"]] <- file.path(out_dir, "bias_decomposition.png")
  .validate_save_plot(p_decomp, paths[["bias_decomposition"]], width = 8, height = 5, dpi = 150)

  # --- Panel 2: Frequency & intensity scatter ---
  p_freq <- ggplot2::ggplot(bias_df, ggplot2::aes(x = obs_event_rate, y = sim_event_rate)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(size = 3.5, color = "steelblue") +
    ggplot2::geom_text(ggplot2::aes(label = location), nudge_y = 0.02, size = 3, check_overlap = TRUE) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = "Observed P(any event)", y = "Simulated P(any event)",
      title = "Frequency Bias: Event Occurrence Rate",
      subtitle = "Dashed = 1:1 line"
    ) +
    ggtheme

  p_int <- ggplot2::ggplot(bias_df, ggplot2::aes(x = obs_mean_int_kt, y = sim_mean_int_kt)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(size = 3.5, color = "#D55E00") +
    ggplot2::geom_text(ggplot2::aes(label = location), nudge_y = 1.5, size = 3, check_overlap = TRUE) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = "Observed E[max | event] (kt)", y = "Simulated E[max | event] (kt)",
      title = "Intensity Bias: Conditional Mean Wind",
      subtitle = "Dashed = 1:1 line"
    ) +
    ggtheme

  paths[["freq_scatter"]] <- file.path(out_dir, "bias_frequency_scatter.png")
  .validate_save_plot(p_freq, paths[["freq_scatter"]], width = 6, height = 5, dpi = 150)
  paths[["int_scatter"]] <- file.path(out_dir, "bias_intensity_scatter.png")
  .validate_save_plot(p_int, paths[["int_scatter"]], width = 6, height = 5, dpi = 150)

  # --- Save decomposition table ---
  paths[["bias_table"]] <- file.path(out_dir, "bias_decomposition.csv")
  .validate_write_csv(bias_df, paths[["bias_table"]])

  invisible(paths)
}


#' Plot QQ plots comparing simulated vs observed annual maxima
#'
#' @param val Output from `run_validation_suite()`.
#' @param cfg Optional `validation_cfg` object.
#' @param out_dir Directory to save plots in. Overrides `cfg$out_dir`.
#' @param base_size Base font size for plots. Overrides `cfg$advanced$base_size`.
#'
#' @return Path of saved plot (invisibly) or `NULL`.
#' @export
plot_qq_validation <- function(val,
                               cfg = NULL,
                               out_dir = NULL,
                               base_size = NULL) {
  plot_cfg <- .resolve_plot_cfg(cfg = cfg, out_dir = out_dir, base_size = base_size)
  out_dir <- plot_cfg$out_dir
  ggtheme <- plot_cfg$theme

  if (is.null(val$hindcast) || is.null(val$hindcast$per_island)) return(invisible(NULL))
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))

  .validate_dir_create(out_dir)

  qq_rows <- list()
  for (isl in names(val$hindcast$per_island)) {
    hc <- val$hindcast$per_island[[isl]]
    if (is.null(hc)) next

    obs_am <- sort(hc$obs_annual_max$V_max_kt)
    n_obs  <- length(obs_am)
    if (n_obs < 5) next

    probs <- (seq_len(n_obs) - 0.5) / n_obs
    sim_q <- stats::quantile(hc$sim_annual_max, probs = probs, na.rm = TRUE)

    qq_rows[[isl]] <- tibble::tibble(
      location = isl,
      obs_quantile = obs_am,
      sim_quantile = as.numeric(sim_q),
      prob = probs
    )
  }

  if (length(qq_rows) == 0) return(invisible(NULL))
  qq_df <- dplyr::bind_rows(qq_rows)

  p_qq <- ggplot2::ggplot(qq_df, ggplot2::aes(x = obs_quantile, y = sim_quantile)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(ggplot2::aes(color = prob), size = 2, alpha = 0.8) +
    ggplot2::scale_color_viridis_c(name = "Quantile", option = "C", limits = c(0, 1)) +
    ggplot2::facet_wrap(~ location, scales = "free", ncol = 3) +
    ggplot2::labs(
      x = "Observed annual max (kt)",
      y = "Simulated annual max (kt)",
      title = "QQ Plot: Simulated vs Observed Annual Maxima",
      subtitle = "Points above dashed line = model overprediction; color = quantile level"
    ) +
    ggtheme

  path <- file.path(out_dir, "qq_annual_max.png")
  .validate_save_plot(p_qq, path, width = 12, height = 7, dpi = 150)
  invisible(path)
}


#' Plot CDF comparison of simulated vs observed annual maxima
#'
#' @param val Output from `run_validation_suite()`.
#' @param cfg Optional `validation_cfg` object.
#' @param out_dir Directory to save plots in. Overrides `cfg$out_dir`.
#' @param base_size Base font size for plots. Overrides `cfg$advanced$base_size`.
#'
#' @return Path of saved plot (invisibly) or `NULL`.
#' @export
plot_cdf_comparison <- function(val,
                                cfg = NULL,
                                out_dir = NULL,
                                base_size = NULL) {
  plot_cfg <- .resolve_plot_cfg(cfg = cfg, out_dir = out_dir, base_size = base_size)
  out_dir <- plot_cfg$out_dir
  ggtheme <- plot_cfg$theme

  if (is.null(val$hindcast) || is.null(val$hindcast$per_island)) return(invisible(NULL))
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))

  .validate_dir_create(out_dir)

  cdf_rows <- list()
  gev_rows <- list()
  for (isl in names(val$hindcast$per_island)) {
    hc <- val$hindcast$per_island[[isl]]
    if (is.null(hc)) next

    obs_am <- sort(hc$obs_annual_max$V_max_kt)
    sim_am <- sort(hc$sim_annual_max)
    n_obs  <- length(obs_am)
    n_sim  <- length(sim_am)
    if (n_obs < 5) next

    cdf_rows[[paste0(isl, "_obs")]] <- tibble::tibble(
      location = isl,
      source = "Observed",
      wind_kt = obs_am,
      ecdf = seq_len(n_obs) / n_obs
    )
    idx_sim <- seq(1, n_sim, by = max(1, n_sim %/% 500))
    cdf_rows[[paste0(isl, "_sim")]] <- tibble::tibble(
      location = isl,
      source = "Simulated",
      wind_kt = sim_am[idx_sim],
      ecdf = idx_sim / n_sim
    )

    if (!is.null(hc$gev_fit) && !is.null(hc$gev_fit$gev_fit)) {
      gev <- hc$gev_fit$gev_fit
      p0  <- hc$gev_fit$p_zero
      x_grid <- seq(0, max(obs_am, sim_am[idx_sim]) * 1.1, length.out = 200)
      cdf_gev <- vapply(x_grid, function(v) {
        if (v <= 0) return(p0)
        p0 + (1 - p0) * .pgev(v, gev$mu, gev$sigma, gev$xi)
      }, numeric(1))

      gev_rows[[isl]] <- tibble::tibble(
        location = isl,
        wind_kt = x_grid,
        cdf_gev = cdf_gev
      )
    }
  }

  if (length(cdf_rows) == 0) return(invisible(NULL))
  cdf_df <- dplyr::bind_rows(cdf_rows)
  gev_df <- if (length(gev_rows) > 0) dplyr::bind_rows(gev_rows) else NULL

  p_cdf <- ggplot2::ggplot(cdf_df, ggplot2::aes(x = wind_kt, y = ecdf, color = source)) +
    ggplot2::geom_step(linewidth = 0.7, alpha = 0.85) +
    ggplot2::scale_color_manual(values = c(Observed = "red", Simulated = "steelblue")) +
    ggplot2::geom_vline(xintercept = 64, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    ggplot2::facet_wrap(~ location, scales = "free_x", ncol = 3) +
    ggplot2::labs(
      x = "Annual maximum site wind (kt)",
      y = "Cumulative probability",
      title = "CDF Comparison: Simulated vs Observed Annual Maxima",
      subtitle = "Simulated CDF right of observed = overprediction; dashed = 64 kt HUR threshold",
      color = "Source"
    ) +
    ggtheme

  if (!is.null(gev_df)) {
    p_cdf <- p_cdf +
      ggplot2::geom_line(
        data = gev_df,
        ggplot2::aes(x = wind_kt, y = cdf_gev),
        color = "steelblue", linetype = "dotted", linewidth = 0.6,
        inherit.aes = FALSE
      )
  }

  paths <- character(0)
  paths[["cdf_comparison"]] <- file.path(out_dir, "cdf_comparison.png")
  .validate_save_plot(p_cdf, paths[["cdf_comparison"]], width = 12, height = 7, dpi = 150)

  # --- Exceedance probability plot (1-CDF, log scale) ---
  p_exceed <- ggplot2::ggplot(cdf_df, ggplot2::aes(x = wind_kt, y = 1 - ecdf, color = source)) +
    ggplot2::geom_step(linewidth = 0.7, alpha = 0.85) +
    ggplot2::scale_y_log10(
      limits = c(0.001, 1),
      labels = scales::label_percent(accuracy = 0.1)
    ) +
    ggplot2::scale_color_manual(values = c(Observed = "red", Simulated = "steelblue")) +
    ggplot2::geom_vline(xintercept = 64, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    ggplot2::facet_wrap(~ location, scales = "free_x", ncol = 3) +
    ggplot2::labs(
      x = "Annual maximum site wind (kt)",
      y = "Exceedance probability (log scale)",
      title = "Exceedance Probability: Simulated vs Observed",
      subtitle = "Higher exceedance at a given wind speed = model overpredicts that intensity",
      color = "Source"
    ) +
    ggtheme

  paths[["exceedance_plot"]] <- file.path(out_dir, "exceedance_comparison.png")
  .validate_save_plot(p_exceed, paths[["exceedance_plot"]], width = 12, height = 7, dpi = 150)

  invisible(paths)
}
