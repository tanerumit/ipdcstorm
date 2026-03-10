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
#' @param n_sim Integer or `NULL`; number of synthetic years to simulate for
#'   hindcast comparison. If `NULL` (default), `run_validation_suite()`
#'   inherits the simulation length from the hazard model output.
#' @param return_periods Numeric vector of return periods (years) to compare
#'   (default: 5, 10, 25, 50).
#' @param conf_level Numeric; confidence level for return-level CIs
#'   (default: 0.90). Must be in `(0, 1)`. This controls the width of the
#'   confidence interval around GEV-based return level estimates, computed
#'   via the delta method or parametric bootstrap of the hurdle-GEV model.
#'   Common choices: 0.90 (standard), 0.95 (stricter), 0.80 (exploratory).
#'   **Note**: this quantifies estimation uncertainty in return levels given
#'   the observed record length --- it is not a prediction interval for
#'   future storm intensities.
#' @param seed Integer; random seed for reproducibility.
#' @param out_dir Character; output directory for saved plots and tables.
#' @param save_plots Logical; whether to save standard validation figures.
#' @param save_tables Logical; whether to save CSV + markdown tables.
#' @param advanced Optional named list of expert parameters. Most users should
#'   leave this as `NULL`. Supported names:
#'   \describe{
#'     \item{`xi_bounds`}{Numeric vector of length 2; allowed range for GEV
#'       shape parameter (default: `c(-0.3, 0.4)`).}
#'     \item{`base_size`}{Numeric; base font size for ggplot themes
#'       (default: 11).}
#'     \item{`hindcast_use_raw_rates`}{Logical; whether Tier 1A hindcast uses
#'       raw fitted rates instead of adjusted rates (default: `TRUE`).}
#'   }
#'
#' @return A list with class `c("validation_cfg", "list")`.
#' @export
#'
#' @examples
#' # Defaults - suitable for most users
#' val_cfg <- make_validation_cfg()
#'
#' # Custom holdout and output location
#' val_cfg <- make_validation_cfg(holdout_years = 15, out_dir = "results/val")
#'
#' # Expert tuning
#' val_cfg <- make_validation_cfg(
#'   n_sim = 10000,
#'   advanced = list(xi_bounds = c(-0.4, 0.5), base_size = 13)
#' )
make_validation_cfg <- function(holdout_years  = 10L,
                                n_sim          = NULL,
                                return_periods = c(5, 10, 25, 50),
                                conf_level     = 0.90,
                                seed           = 42L,
                                out_dir        = "output/validation",
                                save_plots     = TRUE,
                                save_tables    = TRUE,
                                advanced       = NULL) {

  defaults <- list(
    xi_bounds = c(-0.3, 0.4),
    base_size = 11,
    hindcast_use_raw_rates = TRUE
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
  if (!is.null(n_sim)) {
    n_sim <- as.integer(n_sim)
  }
  return_periods <- as.numeric(return_periods)
  seed           <- as.integer(seed)

  if (holdout_years < 1L) stop("holdout_years must be >= 1.", call. = FALSE)
  if (!is.null(n_sim) && n_sim < 100L) stop("n_sim must be >= 100.", call. = FALSE)
  if (length(return_periods) == 0) stop("return_periods must have at least one value.", call. = FALSE)
  if (any(return_periods <= 1)) stop("return_periods must all be > 1.", call. = FALSE)

  conf_level <- as.numeric(conf_level)
  if (length(conf_level) != 1 || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be a single number in (0, 1).", call. = FALSE)
  }
  if (conf_level < 0.50 || conf_level > 0.99) {
    warning("conf_level = ", conf_level, " is outside the typical range [0.50, 0.99]. ",
            "Are you sure?", call. = FALSE)
  }

  cfg <- list(
    holdout_years  = holdout_years,
    n_sim          = n_sim,
    return_periods = return_periods,
    conf_level     = conf_level,
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
  sim_label <- if (is.null(x$n_sim)) {
    "inherit from hazard model output"
  } else {
    sprintf(
      "%s synthetic years",
      format(x$n_sim, big.mark = ",", scientific = FALSE, trim = TRUE)
    )
  }
  cat(sprintf("  Simulation    : %s\n", sim_label))
  cat(sprintf("  Return periods: %s yr\n", paste(x$return_periods, collapse = ", ")))
  cat(sprintf("  Conf. level   : %.0f%%\n", (x$conf_level %||% 0.90) * 100))
  cat(sprintf("  Seed          : %d\n", x$seed))
  cat(sprintf("  Output dir    : %s\n", x$out_dir))
  cat(sprintf("  Save plots    : %s\n", if (x$save_plots) "yes" else "no"))
  cat(sprintf("  Save tables   : %s\n", if (x$save_tables) "yes" else "no"))
  cat(sprintf("  GEV xi bounds : [%.2f, %.2f]\n",
              x$advanced$xi_bounds[1], x$advanced$xi_bounds[2]))
  invisible(x)
}

#' Resolve validation simulation length with inheritance
#' @keywords internal
.resolve_validation_n_sim <- function(cfg, out) {
  if (!is.null(cfg$n_sim)) {
    n_sim <- as.integer(cfg$n_sim)
    source <- "validation_cfg$n_sim"
  } else {
    out_config <- out$config %||% out$cfg
    if (is.null(out_config)) {
      stop(
        "Validation n_sim is NULL and model output does not contain out$config or out$cfg.",
        call. = FALSE
      )
    }
    n_sim <- out_config$n_sim
    if (is.null(n_sim)) {
      stop(
        "Validation n_sim is NULL and model output config does not contain n_sim.",
        call. = FALSE
      )
    }
    n_sim <- as.integer(n_sim)
    source <- "model output config$n_sim"
  }
  if (!is.finite(n_sim) || length(n_sim) != 1L || n_sim < 100L) {
    stop("Effective n_sim must be a single integer >= 100.", call. = FALSE)
  }
  list(n_sim = n_sim, source = source)
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
#' @param xi_bounds Bounds on GEV shape parameter (default: c(-0.3, 0.4) for TSs).
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

  # Physical cap: TS winds can't exceed ~185 kt
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
#' @param conf_level Numeric; confidence level for CIs (default: 0.90).
#'
#' @return Tibble with columns: `return_period`, `sim_median`, `sim_ci_lo`,
#'   `sim_ci_hi`, `sim_lo_50`, `sim_hi_50`. The `ci_lo`/`ci_hi` columns
#'   correspond to the specified `conf_level`.
#' @export
bootstrap_return_level_ci <- function(annual_max,
                                      return_periods = c(5, 10, 25, 50),
                                      n_boot = 500,
                                      xi_bounds = c(-0.3, 0.4),
                                      conf_level = 0.90) {

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

  alpha <- (1 - conf_level) / 2
  tibble::tibble(
    return_period = return_periods,
    sim_median = apply(boot_rl, 2, stats::median, na.rm = TRUE),
    sim_ci_lo  = apply(boot_rl, 2, stats::quantile, probs = alpha, na.rm = TRUE),
    sim_ci_hi  = apply(boot_rl, 2, stats::quantile, probs = 1 - alpha, na.rm = TRUE),
    sim_lo_50  = apply(boot_rl, 2, stats::quantile, probs = 0.25, na.rm = TRUE),
    sim_hi_50  = apply(boot_rl, 2, stats::quantile, probs = 0.75, na.rm = TRUE)
  )
}

# Internal GEV negative log-likelihood (par = c(mu, log_sigma, xi)).
.gev_nll_transformed <- function(par, x, xi_bounds = c(-0.3, 0.4)) {
  mu <- par[1]
  sigma <- exp(par[2])
  xi <- par[3]
  if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0 || !is.finite(xi)) return(Inf)
  if (xi < xi_bounds[1] || xi > xi_bounds[2]) return(Inf)

  z <- (x - mu) / sigma
  if (abs(xi) < 1e-8) {
    return(length(x) * log(sigma) + sum(z + exp(-z)))
  }

  t_val <- 1 + xi * z
  if (any(!is.finite(t_val)) || any(t_val <= 0)) return(Inf)
  length(x) * log(sigma) + (1 + 1 / xi) * sum(log(t_val)) + sum(t_val^(-1 / xi))
}

# Internal GEV random generator from inverse-CDF.
.rgev <- function(n, mu, sigma, xi) {
  u <- stats::runif(n, min = 1e-12, max = 1 - 1e-12)
  .qgev(u, mu = mu, sigma = sigma, xi = xi)
}

# Internal RL gradient wrt c(mu, sigma, xi) for a fixed conditional probability.
.gev_rl_gradient <- function(p_cond, sigma, xi) {
  y <- -log(p_cond)
  lgy <- log(y)
  if (abs(xi) < 1e-8) {
    a <- -lgy
    dq_dxi <- sigma * 0.5 * (lgy^2)
  } else {
    y_pow <- y^(-xi)
    a <- (y_pow - 1) / xi
    dq_dxi <- sigma * ((-lgy * y_pow * xi) - (y_pow - 1)) / (xi^2)
  }
  c(1, a, dq_dxi)
}

# Internal observed RL CI using delta-method first, deterministic bootstrap fallback.
.compute_obs_return_level_ci <- function(annual_max,
                                         return_periods = c(5, 10, 25, 50),
                                         xi_bounds = c(-0.3, 0.4),
                                         conf_level = 0.90,
                                         n_boot_fallback = 300L,
                                         min_n_pos_ci = 10L,
                                         seed = NULL,
                                         force_bootstrap = FALSE) {
  if (!is.null(seed) && is.finite(seed)) set.seed(as.integer(seed))

  base_fit <- compute_return_levels_gev(annual_max, return_periods, xi_bounds = xi_bounds)
  out <- tibble::tibble(
    return_period = return_periods,
    sim_median = as.numeric(base_fit$return_levels),
    sim_ci_lo = NA_real_,
    sim_ci_hi = NA_real_,
    sim_lo_50 = NA_real_,
    sim_hi_50 = NA_real_,
    ci_method = "unavailable"
  )

  n_total <- base_fit$n_total
  n_pos <- base_fit$n_nonzero
  p_zero <- base_fit$p_zero
  gev <- base_fit$gev_fit

  if (!is.finite(n_pos) || n_pos < as.integer(min_n_pos_ci) || is.null(gev)) {
    return(out)
  }
  if (!all(is.finite(c(gev$mu, gev$sigma, gev$xi))) || gev$sigma <= 0) {
    return(out)
  }

  vcov_theta <- NULL
  start_par <- c(gev$mu, log(gev$sigma), gev$xi)
  x_pos <- annual_max[is.finite(annual_max) & annual_max > 0]
  mle <- NULL
  if (!isTRUE(force_bootstrap)) {
    mle <- tryCatch(
      stats::optim(
        par = start_par,
        fn = .gev_nll_transformed,
        x = x_pos,
        xi_bounds = xi_bounds,
        method = "BFGS",
        control = list(maxit = 500, reltol = 1e-10)
      ),
      error = function(e) NULL
    )
    if (!is.null(mle) && is.finite(mle$value)) {
      hess <- tryCatch(
        stats::optimHess(par = mle$par, fn = .gev_nll_transformed, x = x_pos, xi_bounds = xi_bounds),
        error = function(e) NULL
      )
      if (!is.null(hess) && all(is.finite(hess))) {
        inv_hess <- tryCatch(solve(hess), error = function(e) NULL)
        if (!is.null(inv_hess) && all(is.finite(inv_hess))) {
          sigma_hat <- exp(mle$par[2])
          jac <- diag(c(1, sigma_hat, 1))
          vcov_theta <- jac %*% inv_hess %*% jac
          eig <- tryCatch(eigen(vcov_theta, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NULL)
          if (is.null(eig) || any(!is.finite(eig)) || any(eig <= 0)) vcov_theta <- NULL
        }
      }
    }
  }

  z <- stats::qnorm((1 + conf_level) / 2)
  delta_ok <- FALSE
  if (!is.null(vcov_theta)) {
    fit_theta <- c(gev$mu, gev$sigma, gev$xi)
    for (i in seq_along(return_periods)) {
      T_rp <- return_periods[i]
      target_p <- 1 - 1 / T_rp
      p_cond <- (target_p - p_zero) / (1 - p_zero)
      if (!is.finite(p_cond) || p_cond <= 0 || p_cond >= 1) next
      q_hat <- out$sim_median[i]
      if (!is.finite(q_hat)) next
      grad <- .gev_rl_gradient(p_cond, fit_theta[2], fit_theta[3])
      var_q <- as.numeric(t(grad) %*% vcov_theta %*% grad)
      if (!is.finite(var_q) || var_q < 0) next
      se_q <- sqrt(var_q)
      lo <- q_hat - z * se_q
      hi <- q_hat + z * se_q
      if (is.finite(lo) && is.finite(hi)) {
        lo <- max(0, min(185, lo))
        hi <- max(0, min(185, hi))
        out$sim_ci_lo[i] <- min(lo, hi)
        out$sim_ci_hi[i] <- max(lo, hi)
      }
    }
    if (all(is.finite(out$sim_ci_lo) & is.finite(out$sim_ci_hi))) {
      out$ci_method <- "delta"
      delta_ok <- TRUE
    }
  }

  if (delta_ok) return(out)

  n_boot <- max(200L, min(500L, as.integer(n_boot_fallback)))
  if (!is.finite(n_boot) || n_boot < 1L) return(out)

  boot_rl <- matrix(NA_real_, nrow = n_boot, ncol = length(return_periods))
  for (b in seq_len(n_boot)) {
    is_zero <- stats::runif(n_total) < p_zero
    x_boot <- numeric(n_total)
    n_pos_boot <- sum(!is_zero)
    if (n_pos_boot > 0) {
      x_boot[!is_zero] <- .rgev(n_pos_boot, gev$mu, gev$sigma, gev$xi)
      x_boot[!is_zero] <- pmax(0, pmin(185, x_boot[!is_zero]))
    }
    fit_b <- tryCatch(
      compute_return_levels_gev(x_boot, return_periods = return_periods, xi_bounds = xi_bounds)$return_levels,
      error = function(e) rep(NA_real_, length(return_periods))
    )
    boot_rl[b, ] <- fit_b
  }

  for (j in seq_along(return_periods)) {
    rlj <- boot_rl[, j]
    ok <- is.finite(rlj)
    if (sum(ok) < max(30L, ceiling(0.25 * n_boot))) next
    alpha_fb <- (1 - conf_level) / 2
    qq <- stats::quantile(rlj[ok], probs = c(alpha_fb, 0.25, 0.5, 0.75, 1 - alpha_fb), na.rm = TRUE, names = FALSE, type = 7)
    out$sim_ci_lo[j] <- max(0, min(185, qq[1]))
    out$sim_ci_hi[j] <- max(0, min(185, qq[5]))
    out$sim_lo_50[j] <- max(0, min(185, qq[2]))
    out$sim_median[j] <- max(0, min(185, qq[3]))
    out$sim_hi_50[j] <- max(0, min(185, qq[4]))
  }
  if (all(is.finite(out$sim_ci_lo) & is.finite(out$sim_ci_hi))) {
    out$ci_method <- "bootstrap"
  }
  out
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
                               conf_level = 0.90,
                               storm_classes = c("TS", "HUR"),
                               seed = 42,
                               sst_df = NULL,
                               beta_sst = 0,
                               gamma_intensity = 0,
                               use_raw_rates = TRUE,
                               xi_bounds = c(-0.3, 0.4),
                               n_boot = 500) {

  set.seed(seed)
  storm_classes <- .normalize_storm_classes(storm_classes = storm_classes)

  ev <- events_island |>
    dplyr::filter(.data$storm_class %in% c(storm_classes, "none"),
                  is.finite(.data$peak_wind_kt))

  # Tier 1A alignment with Tier 1B "modern TS+ annual-max" regime:
  # restrict to modern era and evaluate annual maxima over TS+ storms with explicit zero years.
  min_year <- 1970L
  ev <- ev |>
    dplyr::filter(is.finite(.data$year), .data$year >= min_year)

  if (nrow(ev) == 0) {
    message("[Hindcast] Skipping ", location, ": no events at/after ", min_year, ".")
    return(list(
      location = location,
      skipped = TRUE,
      skip_reason = "no_modern_years",
      min_year = min_year
    ))
  }

  year_max <- max(as.integer(ev$year), na.rm = TRUE)
  all_years <- seq(min_year, year_max)

  # TS+ events for observed annual-max construction (years with no TS+ -> 0 via join below)
  ev_tsplus <- ev |>
    dplyr::filter(.data$storm_class %in% storm_classes)

  min_train_years <- 10L
  min_holdout_years <- 3L

  n_years <- length(all_years)
  max_holdout <- n_years - min_train_years
  holdout_years_site <- min(as.integer(holdout_years), max_holdout)

  if (!is.finite(holdout_years_site) || holdout_years_site < min_holdout_years) {
    message(
      "[Hindcast] Skipping ", location, ": insufficient years for holdout. Have ",
      n_years, ", need at least ", (min_train_years + min_holdout_years),
      " (min_train_years=", min_train_years, ", min_holdout_years=", min_holdout_years, ")."
    )
    return(list(
      location = location,
      skipped = TRUE,
      skip_reason = "insufficient_years_for_holdout",
      n_years_available = n_years,
      min_train_years = min_train_years,
      min_holdout_years = min_holdout_years,
      requested_holdout_years = as.integer(holdout_years),
      used_holdout_years = as.integer(holdout_years_site)
    ))
  }

  cutoff <- all_years[n_years - holdout_years_site]
  train_years <- all_years[all_years <= cutoff]
  test_years  <- all_years[all_years > cutoff]


  message("[Hindcast] ", location, ": training ", min(train_years), "-", max(train_years),
          " (", length(train_years), " yr), testing ", min(test_years), "-",
          max(test_years), " (", length(test_years), " yr)")

  # --- Observed annual maxima (full record) ---
  obs_annual_max <- ev_tsplus |>
    dplyr::group_by(.data$year) |>
    dplyr::summarise(V_max_kt = max(.data$peak_wind_kt, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::mutate(period = dplyr::if_else(.data$year %in% train_years, "train", "test"))

  full_years <- tibble::tibble(year = all_years)
  obs_annual_max <- full_years |>
    dplyr::left_join(obs_annual_max, by = "year") |>
    dplyr::mutate(
      V_max_kt = dplyr::if_else(is.na(.data$V_max_kt), 0, .data$V_max_kt),
      period = dplyr::if_else(.data$year %in% test_years, "test", "train")
    )

  # --- Fit frequency model on training period ---
  ev_train <- events_island |>
    dplyr::filter(.data$year %in% train_years)

  ac_train <- compute_annual_counts(ev_train, storm_classes = storm_classes)
  lt_train <- compute_lambda_table(ac_train)
  rate_check_train <- .build_rate_check_table(list(
    rates = lt_train |>
      dplyr::mutate(location = location) |>
      dplyr::relocate("location", .before = "storm_class")
  ))
  lambda_scalers_train <- .lambda_scalers_from_rate_check(rate_check_train)
  lt_train_adj <- .apply_lambda_scalers_to_lambda_table(
    lambda_table = lt_train,
    location = location,
    lambda_scalers = lambda_scalers_train
  )

  lt_train_for_sim <- if (isTRUE(use_raw_rates)) lt_train else lt_train_adj
  ki_train <- estimate_k_hat(ac_train)

  train_params <- list(
    lambda_table = lt_train,
    lambda_table_adjusted = lt_train_adj,
    lambda_table_for_sim = lt_train_for_sim,
    use_raw_rates = isTRUE(use_raw_rates),
    rate_check = rate_check_train,
    lambda_scalers = lambda_scalers_train,
    k_hat = ki_train$k_hat,
    n_train_years = length(train_years),
    mu_total = ki_train$mu,
    var_total = ki_train$var
  )

  message("  Hindcast frequency mode: ", if (isTRUE(use_raw_rates)) "RAW (point)" else "SCALED (reference)")

  message("  Rate calibration (training):")
  for (i in seq_len(nrow(rate_check_train))) {
    rr <- rate_check_train[i, ]
    message(sprintf(
      "    %-10s raw=%.3f  ref=%.3f  scale=%.3f  adj=%.3f  [%s]",
      rr$storm_class,
      rr$lambda_model_raw,
      if (is.na(rr$lambda_ref)) NA else rr$lambda_ref,
      rr$lambda_scale,
      rr$lambda_adj,
      rr$scale_status
    ))
  }

  # --- FIT KDE INTENSITY DISTRIBUTIONS ---
  train_V_ts <- ev_train |>
    dplyr::filter(.data$storm_class == "TS") |>
    dplyr::pull(.data$peak_wind_kt) |>
    (\(x) x[is.finite(x)])()

  train_V_hur <- ev_train |>
    dplyr::filter(.data$storm_class == "HUR") |>
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
  fallback_V <- list(TS = 45, HUR = 85)

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
    lt_train_for_sim, ki_train$k_hat, n_years_sim = n_sim,
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
        winds <- c(winds, rep(fallback_V$HUR, n_hur) +
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

  # --- OBSERVED RL CIs (delta first, deterministic bootstrap fallback) ---
  obs_rl_ci <- .compute_obs_return_level_ci(
    annual_max = obs_full_max,
    return_periods = return_periods,
    xi_bounds = xi_bounds,
    conf_level = conf_level,
    n_boot_fallback = if (isTRUE(as.integer(n_boot) > 0L)) n_boot else 300L,
    min_n_pos_ci = 10L,
    seed = seed + 1000L
  )

  if (isTRUE(as.integer(n_boot) > 0L)) {
    # Model CIs (secondary diagnostic - should be narrow)
    sim_rl_ci <- bootstrap_return_level_ci(
      sim_annual_max, return_periods,
      n_boot = min(200L, n_boot),
      xi_bounds = xi_bounds,
      conf_level = conf_level
    )
  } else {
    sim_rl_ci <- tibble::tibble(
      return_period = return_periods,
      sim_median = sim_rl,
      sim_ci_lo = NA_real_,
      sim_ci_hi = NA_real_,
      sim_lo_50 = NA_real_,
      sim_hi_50 = NA_real_
    )
  }

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
    sim_ci_lo = sim_rl_ci$sim_ci_lo,
    sim_ci_hi = sim_rl_ci$sim_ci_hi,
    obs_ci_lo = obs_rl_ci$sim_ci_lo,
    obs_ci_hi = obs_rl_ci$sim_ci_hi,
    obs_ci_method = obs_rl_ci$ci_method,
    obs_in_ci = dplyr::if_else(
      is.finite(obs_rl_ci$sim_ci_lo) & is.finite(obs_rl_ci$sim_ci_hi),
      sim_rl >= obs_rl_ci$sim_ci_lo & sim_rl <= obs_rl_ci$sim_ci_hi,
      NA
    ),
    model_in_obs_ci = dplyr::if_else(
      is.finite(obs_rl_ci$sim_ci_lo) & is.finite(obs_rl_ci$sim_ci_hi),
      sim_rl >= obs_rl_ci$sim_ci_lo & sim_rl <= obs_rl_ci$sim_ci_hi,
      NA
    ),
    obs_in_model_ci = dplyr::if_else(
      is.finite(sim_rl_ci$sim_ci_lo) & is.finite(sim_rl_ci$sim_ci_hi),
      obs_full_rl >= sim_rl_ci$sim_ci_lo & obs_full_rl <= sim_rl_ci$sim_ci_hi,
      NA
    ),
    bias_pct = 100 * (sim_rl - obs_full_rl) / pmax(obs_full_rl, 1)
  )
  attr(comparison, "conf_level") <- conf_level

  list(
    location = location,
    train_years = train_years,
    test_years = test_years,
    train_params = train_params,
    obs_annual_max = obs_annual_max,
    obs_test_rl = obs_test_rl,
    obs_full_rl = obs_full_rl,
    obs_rl_ci = obs_rl_ci,
    obs_gev = obs_gev,
    sim_rl = sim_rl,
    sim_annual_max = sim_annual_max,
    sim_rl_ci = sim_rl_ci,
    comparison = comparison,
    gev_fit = sim_gev,
    kde_fits = list(TS = kde_ts, HUR = kde_hur),
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
                                   conf_level = 0.90,
                                   seed = 42,
                                   sst_df = NULL,
                                   beta_sst = 0,
                                   gamma_intensity = 0,
                                   use_raw_rates = TRUE,
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
        conf_level = conf_level,
        seed = seed,
        sst_df = sst_df,
        beta_sst = beta_sst,
        gamma_intensity = gamma_intensity,
        use_raw_rates = use_raw_rates,
        xi_bounds = xi_bounds,
        n_boot = n_boot
      )

      if (isTRUE(results[[location]]$skipped)) next

    }, error = function(e) {
      message("[Hindcast] Error for ", location, ": ", e$message)
    })
  }

  ok_results <- Filter(function(r) {
    !is.null(r) && !isTRUE(r$skipped) && !is.null(r$comparison)
  }, results)

  comparison_all <- dplyr::bind_rows(lapply(ok_results, function(r) r$comparison))

  list(
    per_island = results,
    comparison = comparison_all
  )
}


# =============================================================================
# 7) RATE CHECK
# =============================================================================

#' Reference HURDAT2/literature annual rates for Leeward Islands region
#'
#' @description
#' Returns a tibble of published annual TS passage rates from literature,
#' for comparison against model-fitted lambdas.
#'
#' @return Tibble with columns: region, storm_class, lambda_ref, source, gate_approx_nm, period.
#' @export
get_reference_rates <- function() {
  tibble::tribble(
    ~region,            ~storm_class,    ~lambda_ref, ~source,                           ~gate_approx_nm, ~period,
    "Leeward_Islands",  "TS",   2.0,         "NHC Climo (100nm, 1970-2020)",    100,             "1970-2020",
    "Leeward_Islands",  "HUR",  0.55,        "NHC Climo (100nm, 1970-2020)",    100,             "1970-2020",
    "St_Martin",        "TS",   1.2,         "NOAA TC Climo (65nm, 1970-2023)", 65,              "1970-2023",
    "St_Martin",        "HUR",  0.40,        "NOAA TC Climo (65nm, 1970-2023)", 65,              "1970-2023",
    "Puerto_Rico",      "TS",   1.8,         "NHC Climo (100nm, 1970-2020)",    100,             "1970-2020",
    "Puerto_Rico",      "HUR",  0.45,        "NHC Climo (100nm, 1970-2020)",    100,             "1970-2020",
    "Miami",            "TS",   1.5,         "NHC Climo (100nm, 1970-2020)",    100,             "1970-2020",
    "Miami",            "HUR",  0.35,        "NHC Climo (100nm, 1970-2020)",    100,             "1970-2020"
  )
}


.build_rate_check_table <- function(out, ref_rates = NULL, lambda_scaling_mode = NULL, ...) {
  if (is.null(ref_rates)) ref_rates <- get_reference_rates()
  if (is.null(out$rates)) stop("out$rates is required.", call. = FALSE)
  if (nrow(out$rates) == 0) {
    return(tibble::tibble(
      location = character(0),
      storm_class = character(0),
      lambda_model = numeric(0),
      lambda_model_raw = numeric(0),
      n_years_model = integer(0),
      lambda_ref = numeric(0),
      lambda_target = numeric(0),
      lambda_scale = numeric(0),
      lambda_adj = numeric(0),
      scale_status = character(0),
      scale_clamped = logical(0),
      source = character(0),
      gate_approx_nm = numeric(0),
      period = character(0),
      expected_ratio = numeric(0),
      raw_ratio = numeric(0),
      ratio = numeric(0),
      adj_ratio = numeric(0),
      flag = character(0)
    ))
  }

  model_rates <- out$rates |>
    dplyr::select(
      location = "location",
      storm_class = "storm_class",
      lambda_model_raw = "lambda",
      n_years_model = "n_years"
    )

  model_rates2 <- model_rates |>
    tidyr::pivot_wider(names_from = storm_class, values_from = lambda_model_raw, values_fill = 0)
  if (!("TS" %in% names(model_rates2))) model_rates2$TS <- 0
  if (!("HUR" %in% names(model_rates2))) model_rates2$HUR <- 0

  island_to_region <- tibble::tribble(
    ~location,        ~region,
    "St_Martin",      "St_Martin",
    "Saba",           "Leeward_Islands",
    "Statia",         "Leeward_Islands",
    "Puerto_Rico",    "Puerto_Rico",
    "Miami",          "Miami"
  )

  ref_wide <- ref_rates |>
    dplyr::mutate(storm_class = as.character(.data$storm_class)) |>
    tidyr::pivot_wider(
      names_from = storm_class,
      values_from = c(lambda_ref, source, gate_approx_nm, period),
      names_sep = "__"
    )

  comp_base <- model_rates2 |>
    dplyr::left_join(island_to_region, by = "location") |>
    dplyr::left_join(ref_wide, by = "region") |>
    dplyr::mutate(
      expected_ratio_total = dplyr::if_else(.data$gate_approx_nm__TS >= 100, 0.55, 0.75),
      expected_ratio_hur = dplyr::if_else(.data$gate_approx_nm__HUR >= 100, 0.30, 0.45),
      lambda_ref_ts = pmax(.data$lambda_ref__TS - .data$lambda_ref__HUR, 0),
      lambda_target_total = .data$lambda_ref__TS * .data$expected_ratio_total,
      lambda_target_hur = .data$lambda_ref__HUR * .data$expected_ratio_hur,
      lambda_target_ts = pmax(.data$lambda_target_total - .data$lambda_target_hur, 0),
      expected_ratio_ts = dplyr::if_else(
        is.finite(.data$lambda_ref_ts) & .data$lambda_ref_ts > 0,
        .data$lambda_target_ts / .data$lambda_ref_ts,
        .data$expected_ratio_total
      )
    ) |>
    dplyr::transmute(
      location = .data$location,
      region = .data$region,
      storm_class = "TS",
      lambda_model_raw = .data$TS,
      lambda_model_total_raw = .data$TS + .data$HUR,
      lambda_ref = .data$lambda_ref_ts,
      lambda_ref_total = .data$lambda_ref__TS,
      lambda_target_total = .data$lambda_target_total,
      source = .data$source__TS,
      gate_approx_nm = .data$gate_approx_nm__TS,
      period = .data$period__TS,
      expected_ratio = .data$expected_ratio_ts,
      expected_ratio_total = .data$expected_ratio_total,
      n_years_model = .data$n_years_model
    ) |>
    dplyr::bind_rows(
      model_rates2 |>
        dplyr::left_join(island_to_region, by = "location") |>
        dplyr::left_join(ref_wide, by = "region") |>
        dplyr::mutate(
          expected_ratio = dplyr::if_else(.data$gate_approx_nm__HUR >= 100, 0.30, 0.45)
        ) |>
        dplyr::transmute(
          location = .data$location,
          region = .data$region,
          storm_class = "HUR",
          lambda_model_raw = .data$HUR,
          lambda_model_total_raw = .data$TS + .data$HUR,
          lambda_ref = .data$lambda_ref__HUR,
          lambda_ref_total = .data$lambda_ref__TS,
          lambda_target_total = .data$lambda_ref__TS * dplyr::if_else(.data$gate_approx_nm__TS >= 100, 0.55, 0.75),
          source = .data$source__HUR,
          gate_approx_nm = .data$gate_approx_nm__HUR,
          period = .data$period__HUR,
          expected_ratio = .data$expected_ratio,
          expected_ratio_total = dplyr::if_else(.data$gate_approx_nm__TS >= 100, 0.55, 0.75),
          n_years_model = .data$n_years_model
        )
    )

  lambda_scalers <- out$lambda_scalers
  if (is.null(lambda_scalers) || nrow(lambda_scalers) == 0) {
    lambda_scalers <- .lambda_scalers_from_rate_check(comp_base)
  }

  comp_base |>
    dplyr::left_join(
      dplyr::select(
        tibble::as_tibble(lambda_scalers),
        dplyr::any_of(c(
          "location", "storm_class", "lambda_scaling_mode",
          "lambda_target", "lambda_scale", "lambda_adj",
          "scale_status", "scale_clamped"
        ))
      ),
      by = c("location", "storm_class")
    ) |>
    dplyr::mutate(
      lambda_scale = dplyr::if_else(
        is.finite(.data$lambda_scale) & .data$lambda_scale > 0,
        .data$lambda_scale,
        1
      ),
      lambda_target = dplyr::if_else(
        is.finite(.data$lambda_target),
        .data$lambda_target,
        dplyr::if_else(
          is.finite(.data$lambda_ref),
          .data$lambda_ref * .data$expected_ratio,
          .data$lambda_model_raw
        )
      ),
      lambda_adj = dplyr::if_else(
        is.finite(.data$lambda_adj),
        .data$lambda_adj,
        .data$lambda_model_raw * .data$lambda_scale
      ),
      lambda_model = .data$lambda_model_raw,
      raw_ratio = dplyr::if_else(
        is.finite(.data$lambda_ref),
        .data$lambda_model_raw / pmax(.data$lambda_ref, 0.001),
        NA_real_
      ),
      ratio = dplyr::if_else(
        is.finite(.data$lambda_ref),
        .data$lambda_adj / pmax(.data$lambda_ref, 0.001),
        NA_real_
      ),
      adj_ratio = dplyr::if_else(
        is.finite(.data$ratio) & is.finite(.data$expected_ratio) & .data$expected_ratio > 0,
        .data$ratio / .data$expected_ratio,
        NA_real_
      ),
      flag = dplyr::case_when(
        is.na(.data$lambda_ref)   ~ "no_reference",
        .data$adj_ratio > 2.5     ~ "HIGH: model >> expected",
        .data$adj_ratio < 0.4     ~ "LOW: model << expected",
        .data$adj_ratio > 1.8     ~ "elevated",
        .data$adj_ratio < 0.6     ~ "slightly_low",
        TRUE ~ "OK"
      )
    ) |>
    dplyr::mutate(
      lambda_scaling_mode = if (is.null(lambda_scaling_mode)) NA_character_ else as.character(lambda_scaling_mode)
    ) |>
    dplyr::select(
      "location", "storm_class",
      "lambda_model", "lambda_model_raw", "n_years_model",
      "lambda_ref", "lambda_target", "lambda_scale", "lambda_adj",
      "scale_status", "scale_clamped",
      "source", "gate_approx_nm", "period",
      "expected_ratio", "raw_ratio", "ratio", "adj_ratio", "flag"
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
  comp <- .build_rate_check_table(out = out, ref_rates = ref_rates)

  message("\n[Rate Check] Summary:")
  for (i in seq_len(nrow(comp))) {
    r <- comp[i, ]
    sym <- if (r$flag == "OK") "\u2713" else if (grepl("^(HIGH|LOW)", r$flag)) "\u2717" else "~"
    message(sprintf(
      "  %s %s / %-10s : raw=%.3f  ref=%.3f  scale=%.3f  adj=%.3f  raw_ratio=%.2f  exp_ratio=%.2f  adj_ratio=%.2f  [%s; %s]",
      sym, r$location, r$storm_class, r$lambda_model_raw,
      if (is.na(r$lambda_ref)) NA else r$lambda_ref,
      r$lambda_scale,
      r$lambda_adj,
      if (is.na(r$raw_ratio)) NA else r$raw_ratio,
      if (is.na(r$expected_ratio)) NA else r$expected_ratio,
      if (is.na(r$adj_ratio)) NA else r$adj_ratio,
      r$flag,
      r$scale_status
    ))
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
#' for comparison against model-estimated site wind. Each observation includes a
#' quality tier (`obs_quality`) that determines the acceptable bias threshold:
#' - **A** (Direct): Station measurement, instrument survived (\eqn{\pm}15%).
#' - **B** (Converted): 10-min reading with averaging-period conversion (\eqn{\pm}25%).
#' - **C** (Estimated): NHC best-track or indirect estimate (\eqn{\pm}35%).
#'
#' @return Tibble with columns: storm_sid, storm_name, year, target_island,
#'   station, obs_wind_kt, obs_type, obs_quality, obs_source, notes.
#' @export
get_wind_observations <- function() {
  tibble::tribble(
    ~storm_sid,          ~storm_name, ~year, ~target_island, ~station,
    ~obs_wind_kt, ~obs_type,       ~obs_quality, ~obs_source,               ~notes,

    "2017242N16333",     "IRMA",    2017L, "St_Martin",     "Meteo_France_SXM",
    113,                 "10min_sust",    "B",  "Meteo France RSMC report", "Last valid reading before instrument failure; 10-min sustained",

    "2014279N15323",     "GONZALO", 2014L, "St_Martin",     "Juliana_Airport_TNCM",
    60,                  "1min_sust",     "A",  "NHC TCR (AL082014)",      "Marginal TS winds at SXM; storm passed ~50nm north",

    "2017242N16333",     "IRMA",    2017L, "Saba",          "NHC_estimate",
    80,                  "1min_sust",     "C",  "NHC TCR estimate",        "Estimated from best-track; Saba ~30nm from eye center",

    "2017255N12319",     "MARIA",   2017L, "St_Martin",     "Juliana_Airport_TNCM",
    40,                  "1min_sust",     "A",  "NHC TCR (AL152017)",      "Maria passed ~100nm south; TS-force winds at SXM",

    "1989248N12343",     "HUGO",    1989L, "Puerto_Rico",   "Roosevelt_Roads_NAS",
    104,                 "1min_sust",     "A",  "NHC TCR (AL081989)",      "Direct hit eastern PR",

    "1992216N10325",     "ANDREW",  1992L, "Miami",         "NHC_Miami_ASOS",
    141,                 "1min_sust",     "A",  "NHC TCR (AL041992)",      "Before instrument failure at NHC",

    "1999317N14290",     "LENNY",   1999L, "St_Martin",     "Juliana_Airport_TNCM",
    55,                  "1min_sust",     "A",  "NHC TCR (AL171999)",      "Unusual west-to-east track; Cat 4 but passed south",

    "1995241N12330",     "LUIS",    1995L, "St_Martin",     "Juliana_Airport_TNCM",
    110,                 "1min_sust",     "A",  "NHC TCR (AL131995)",      "Cat 4 at closest approach to SXM; major damage"
  )
}


#' Observation quality bias thresholds for wind field validation
#'
#' @description
#' Returns acceptable absolute bias percentages per observation quality tier.
#' Quality tiers reflect measurement uncertainty:
#' - **A** (Direct): Station measurement, instrument survived. Threshold: 15%.
#' - **B** (Converted): 10-min reading with averaging-period conversion. Threshold: 25%.
#' - **C** (Estimated): NHC best-track or indirect estimate. Threshold: 35%.
#'
#' @return Named numeric vector: quality tier -> acceptable |bias| (%).
#' @keywords internal
.wf_quality_thresholds <- function() {
  c(A = 15, B = 25, C = 35)
}


#' Validate wind field estimates against station observations
#'
#' @description
#' For each reference observation, finds the matching storm in the model's trackpoint
#' data, extracts the model's peak site wind for that storm at the target location,
#' and compares against the observed wind. Each observation carries a quality tier
#' (A/B/C) that determines the acceptable bias threshold.
#'
#' @param out List returned by `run_hazard_model()`.
#' @param obs_table Optional tibble of observations (default: `get_wind_observations()`).
#'
#' @return Tibble with model vs observed comparison per storm-station pair,
#'   including `obs_quality`, `bias_threshold_pct`, and `bias_ok` columns.
#' @export
validate_wind_field <- function(out, obs_table = NULL) {

  if (is.null(obs_table)) obs_table <- get_wind_observations()

  # Ensure obs_quality column exists (backward compat with user-supplied tables)
  if (!"obs_quality" %in% names(obs_table)) {
    obs_table$obs_quality <- "B"
  }

  quality_thresh <- .wf_quality_thresholds()

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

    # Quality-dependent threshold
    q <- toupper(as.character(obs$obs_quality))
    if (!q %in% names(quality_thresh)) q <- "B"
    thresh <- as.numeric(quality_thresh[q])

    bias_pct_val <- if (is.finite(model_best) && obs_1min_kt > 0) {
      100 * (model_best - obs_1min_kt) / obs_1min_kt
    } else {
      NA_real_
    }

    results[[i]] <- tibble::tibble(
      location = location,
      storm_name = obs$storm_name,
      storm_sid = obs$storm_sid,
      year = obs$year,
      station = obs$station,
      obs_raw_kt = obs$obs_wind_kt,
      obs_type = obs$obs_type,
      obs_quality = q,
      obs_1min_equiv_kt = obs_1min_kt,
      model_V_site_kt = model_best,
      model_wind_max_kt = model_wind_max,
      min_dist_km = min_dist_km,
      storm_found = storm_found,
      bias_kt = if (is.finite(model_best)) model_best - obs_1min_kt else NA_real_,
      bias_pct = bias_pct_val,
      bias_threshold_pct = thresh,
      bias_ok = if (is.finite(bias_pct_val)) abs(bias_pct_val) <= thresh else NA,
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
      sym <- if (isTRUE(r$bias_ok)) "\u2713" else "\u2717"
      message(sprintf("    %s %-8s @ %-12s: obs=%5.0f kt  model=%5.0f kt  bias=%+.0f kt (%+.0f%%)  [%s \u00b1%d%%]",
                      sym, r$storm_name, r$location, r$obs_1min_equiv_kt,
                      r$model_V_site_kt, r$bias_kt, r$bias_pct,
                      r$obs_quality, r$bias_threshold_pct))
    }
  } else {
    message("[Wind Field Check] No matching storms found in model output.")
  }

  comp
}


# Internal helper for scalar quantiles used in compact validation summaries.
.safe_quantile_scalar <- function(x, prob) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE, type = 7))
}

.validation_summary_from_hindcast <- function(hc) {
  if (is.null(hc) || is.null(hc$per_island)) {
    return(tibble::tibble(
      location = character(0),
      storm_class = character(0),
      delta_top1_p50 = numeric(0),
      delta_overall_p99 = numeric(0)
    ))
  }

  rows <- list()
  idx <- 0L
  for (loc in names(hc$per_island)) {
    hloc <- hc$per_island[[loc]]
    if (is.null(hloc) || isTRUE(hloc$skipped)) next
    obs <- hloc$obs_annual_max$V_max_kt
    sim <- hloc$sim_annual_max
    obs <- obs[is.finite(obs)]
    sim <- sim[is.finite(sim)]

    obs_p99 <- .safe_quantile_scalar(obs, 0.99)
    sim_p99 <- .safe_quantile_scalar(sim, 0.99)
    obs_top <- obs[is.finite(obs_p99) & obs >= obs_p99]
    sim_top <- sim[is.finite(sim_p99) & sim >= sim_p99]

    idx <- idx + 1L
    rows[[idx]] <- tibble::tibble(
      location = loc,
      storm_class = "TS",
      delta_top1_p50 = .safe_quantile_scalar(sim_top, 0.50) - .safe_quantile_scalar(obs_top, 0.50),
      delta_overall_p99 = sim_p99 - obs_p99
    )
  }

  if (!length(rows)) {
    return(tibble::tibble(
      location = character(0),
      storm_class = character(0),
      delta_top1_p50 = numeric(0),
      delta_overall_p99 = numeric(0)
    ))
  }
  dplyr::bind_rows(rows)
}

#' Run the full validation suite
#'
#' @description
#' Executes all three validation tiers (hindcast, rate check, wind field),
#' produces a consolidated diagnostic report, and optionally saves all plots
#' and tables to the configured output directory.
#'
# =============================================================================
# CONSOLE FORMATTING HELPERS
# =============================================================================

.val_header <- function(title, width = 72) {
  message(paste(rep("\u2550", width), collapse = ""))
  message("  ", title)
}

.val_header_close <- function(width = 72) {
  message(paste(rep("\u2550", width), collapse = ""))
}

.val_section <- function(title, width = 72) {
  pad <- width - nchar(title) - 4
  message("\n\u2500\u2500 ", title, " ", paste(rep("\u2500", max(1, pad)), collapse = ""))
}

.val_table_row <- function(..., widths = NULL) {
  args <- list(...)
  if (is.null(widths)) {
    message("  ", paste(args, collapse = "  "))
  } else {
    parts <- mapply(function(a, w) {
      formatC(a, width = w, flag = "-")
    }, args, widths, SIMPLIFY = TRUE)
    message("  ", paste(parts, collapse = " "))
  }
}

.val_table_sep <- function(width = 68) {
  message("  ", paste(rep("\u2500", width), collapse = ""))
}

#' RUN VALIDATION SUITE
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

  n_sim_info <- .resolve_validation_n_sim(cfg = cfg, out = out)
  holdout_years  <- cfg$holdout_years
  n_sim          <- n_sim_info$n_sim
  return_periods <- cfg$return_periods
  conf_level     <- cfg$conf_level %||% 0.90
  seed           <- cfg$seed
  xi_bounds      <- cfg$advanced$xi_bounds
  base_size      <- cfg$advanced$base_size
  use_raw_rates  <- isTRUE(cfg$advanced$hindcast_use_raw_rates)
  out_dir        <- cfg$out_dir

  ci_pct <- sprintf("%.0f%%", conf_level * 100)

  # â”€â”€ Header â”€â”€
  .val_header("HAZARD MODEL VALIDATION SUITE")
  message(sprintf("  CI: %s | Holdout: %d yr | Sim: %s yr (%s)",
                  ci_pct, holdout_years,
                  format(n_sim, big.mark = ",", scientific = FALSE, trim = TRUE),
                  n_sim_info$source))
  .val_header_close()

  # --- Extract climate info ---
  sst_df_val <- NULL
  beta_sst_val <- 0
  gamma_val <- 0
  if (!is.null(out$fit) && nrow(out$fit) > 0) {
    beta_sst_val <- out$fit$beta_sst[1]
    gamma_val <- out$fit$gamma_intensity[1]
  }

  # â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
  # TIER 1A: HINDCAST
  # â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
  .val_section("HINDCAST VALIDATION (return levels)")

  hc <- tryCatch(
    suppressMessages(.validate_hindcast_all(out,
                                            holdout_years = holdout_years,
                                            n_sim = n_sim,
                                            return_periods = return_periods,
                                            conf_level = conf_level,
                                            seed = seed,
                                            sst_df = sst_df_val,
                                            beta_sst = beta_sst_val,
                                            gamma_intensity = gamma_val,
                                            use_raw_rates = use_raw_rates,
                                            xi_bounds = xi_bounds,
                                            n_boot = 0L)),
    error = function(e) { message("  ERROR: ", e$message); NULL }
  )

  if (!is.null(hc) && !is.null(hc$per_island) && length(hc$per_island) > 0) {
    # Training splits summary
    message("")
    for (isl in names(hc$per_island)) {
      hi <- hc$per_island[[isl]]
      if (is.null(hi) || isTRUE(hi$skipped)) next
      train_rng <- range(hi$train_years)
      n_train <- length(hi$train_years)
      n_test  <- length(hi$test_years %||% numeric(0))
      end_yr  <- if (n_test > 0) max(hi$test_years) else max(hi$train_years)
      # GEV summary
      gev_obs <- hi$obs_gev
      gev_sim <- hi$gev_fit
      obs_str <- if (!is.null(gev_obs) && !is.null(gev_obs$gev_fit))
        sprintf("\u03be=%.2f n=%d", gev_obs$gev_fit$xi, gev_obs$n_nonzero) else "NA"
      sim_str <- if (!is.null(gev_sim) && !is.null(gev_sim$gev_fit))
        sprintf("\u03be=%.2f n=%d", gev_sim$gev_fit$xi, gev_sim$n_nonzero) else "NA"
      message(sprintf("  %-14s %d\u2013%d  obs GEV: %s | sim GEV: %s",
                      isl, train_rng[1], end_yr, obs_str, sim_str))
    }

    # Return-level table
    if (nrow(hc$comparison) > 0) {
      message("")
      w <- c(14, 4, 7, 7, 14, 6, 4)
      .val_table_row("Location", "RP", "Obs kt", "Sim kt",
                     paste0(ci_pct, " CI"), "Bias", "",
                     widths = w)
      .val_table_sep()

      for (i in seq_len(nrow(hc$comparison))) {
        r <- hc$comparison[i, ]
        ci_str <- if (is.finite(r$obs_ci_lo) && is.finite(r$obs_ci_hi)) {
          sprintf("[%3.0f, %3.0f]", r$obs_ci_lo, r$obs_ci_hi)
        } else {
          "  [NA,  NA]"
        }
        status <- if (!is.finite(r$obs_ci_lo) || !is.finite(r$obs_ci_hi)) {
          "\u2014"
        } else if (isTRUE(r$obs_in_ci)) {
          "\u2713"
        } else {
          "\u2717"
        }
        .val_table_row(
          r$location,
          sprintf("%d", r$return_period),
          sprintf("%.0f", r$obs_full_rl),
          sprintf("%.0f", r$sim_rl),
          ci_str,
          sprintf("%+.0f%%", r$bias_pct),
          status,
          widths = w
        )
      }
    }
  }

  # â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
  # TIER 2: RATE CHECK
  # â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
  rc <- tryCatch(
    suppressMessages(validate_rates(out)),
    error = function(e) NULL
  )

  # â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
  # TIER 3: WIND FIELD
  # â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
  .val_section("WIND FIELD SPOT-CHECKS")

  wf <- tryCatch(
    suppressMessages(validate_wind_field(out)),
    error = function(e) { message("  ERROR: ", e$message); NULL }
  )

  if (!is.null(wf)) {
    valid_wf <- wf[is.finite(wf$bias_kt), , drop = FALSE]
    if (nrow(valid_wf) > 0) {
      mae  <- mean(abs(valid_wf$bias_kt))
      rmse <- sqrt(mean(valid_wf$bias_kt^2))
      mb   <- mean(valid_wf$bias_kt)
      message(sprintf("  Mean bias: %+.1f kt | MAE: %.1f kt | RMSE: %.1f kt (%d storms)",
                      mb, mae, rmse, nrow(valid_wf)))
      message("")
      w_wf <- c(9, 14, 6, 6, 7, 6, 4)
      .val_table_row("Storm", "Location", "Obs", "Model", "Bias", "Tol", "",
                     widths = w_wf)
      .val_table_sep(55)
      for (j in seq_len(nrow(valid_wf))) {
        r <- valid_wf[j, ]
        sym <- if (isTRUE(r$bias_ok)) "\u2713" else "\u2717"
        tol_str <- sprintf("%s\u00b1%d%%", r$obs_quality, r$bias_threshold_pct)
        .val_table_row(
          r$storm_name,
          r$location,
          sprintf("%.0f", r$obs_1min_equiv_kt),
          sprintf("%.0f", r$model_V_site_kt),
          sprintf("%+.0f%%", r$bias_pct),
          tol_str,
          sym,
          widths = w_wf
        )
      }
      message("  Tol: A=direct station (\u00b115%) | B=converted/10-min (\u00b125%) | C=estimated (\u00b135%)")
    } else {
      message("  No matching storms found.")
    }
  }

  # â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
  # SUMMARY
  # â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
  hc_comp_scored <- NULL
  if (!is.null(hc) && !is.null(hc$comparison) && nrow(hc$comparison) > 0) {
    hc_comp_scored <- hc$comparison
  }

  n_rl_ok <- if (!is.null(hc_comp_scored)) sum(hc_comp_scored$obs_in_ci, na.rm = TRUE) else 0
  n_rl_total <- if (!is.null(hc_comp_scored)) sum(!is.na(hc_comp_scored$obs_in_ci)) else 0

  n_rate_ok <- if (!is.null(rc)) sum(rc$flag == "OK", na.rm = TRUE) else 0
  n_rate_total <- if (!is.null(rc)) sum(!is.na(rc$flag)) else 0

  n_wf_ok <- if (!is.null(wf)) {
    sum(wf$bias_ok, na.rm = TRUE)
  } else 0
  n_wf_total <- if (!is.null(wf)) sum(!is.na(wf$bias_ok)) else 0

  summary_tbl <- .validation_summary_from_hindcast(hc)

  rl_sym  <- if (n_rl_total > 0 && n_rl_ok == n_rl_total) "\u2713" else if (n_rl_total == 0) "\u2014" else "\u2717"
  rc_sym  <- if (n_rate_total > 0 && n_rate_ok == n_rate_total) "\u2713" else if (n_rate_total == 0) "\u2014" else "\u2717"
  wf_sym  <- if (n_wf_total > 0 && n_wf_ok == n_wf_total) "\u2713" else if (n_wf_total == 0) "\u2014" else "\u2717"

  message("")
  .val_header("SUMMARY")
  message(sprintf("  Hindcast         : %2d / %2d within %s CI     %s",
                  n_rl_ok, n_rl_total, ci_pct, rl_sym))
  message(sprintf("  Rate check       : %2d / %2d flagged OK        %s",
                  n_rate_ok, n_rate_total, rc_sym))
  message(sprintf("  Wind field       : %2d / %2d within tolerance  %s",
                  n_wf_ok, n_wf_total, wf_sym))
  .val_header_close()

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
    hindcast_plot <- plot_hindcast_validation(val, cfg = cfg, out_dir = out_dir, base_size = base_size)
    if (!is.null(hindcast_plot) && length(hindcast_plot) > 0) {
      artifacts$plots$hindcast_return_levels <- unname(hindcast_plot[[1]])
    }

    bias_plot <- plot_bias_diagnostics(val, cfg = cfg, out_dir = out_dir, base_size = base_size)
    if (!is.null(bias_plot) && length(bias_plot) > 0) {
      artifacts$plots$bias_decomposition <- unname(bias_plot[[1]])
    }

    wind_plot <- plot_wind_field_validation(val, cfg = cfg, out_dir = out_dir, base_size = base_size)
    if (!is.null(wind_plot) && length(wind_plot) > 0) {
      artifacts$plots$wind_field_scatter <- unname(wind_plot[[1]])
    }

    qq_plot <- plot_qq_validation(val, cfg = cfg, out_dir = out_dir, base_size = base_size)
    if (!is.null(qq_plot) && length(qq_plot) > 0) {
      artifacts$plots$qq_annual_max <- unname(qq_plot[[1]])
    }

    exceedance_plot <- plot_cdf_comparison(val, cfg = cfg, out_dir = out_dir, base_size = base_size)
    if (!is.null(exceedance_plot) && length(exceedance_plot) > 0) {
      artifacts$plots$exceedance_comparison <- unname(exceedance_plot[[1]])
    }
  }

  if (isTRUE(cfg$save_tables)) {
    .validate_dir_create(out_dir)
    artifacts$tables$hindcast_csv <- file.path(out_dir, "hindcast_return_levels.csv")
    .validate_write_csv(
      if (!is.null(val$hindcast) && !is.null(val$hindcast$comparison)) val$hindcast$comparison else data.frame(),
      artifacts$tables$hindcast_csv
    )

    artifacts$tables$rate_check_csv <- file.path(out_dir, "rate_check.csv")
    .validate_write_csv(
      if (!is.null(val$rate_check)) val$rate_check else data.frame(),
      artifacts$tables$rate_check_csv
    )

    artifacts$tables$wind_field_csv <- file.path(out_dir, "wind_field.csv")
    .validate_write_csv(
      if (!is.null(val$wind_field)) val$wind_field else data.frame(),
      artifacts$tables$wind_field_csv
    )

    artifacts$tables$summary_csv <- file.path(out_dir, "validation_summary.csv")
    .validate_write_csv(
      if (!is.null(val$summary)) val$summary else data.frame(),
      artifacts$tables$summary_csv
    )
  }

  val$artifacts <- artifacts
  val
}


# =============================================================================
# 10) END-TO-END CONVENIENCE WRAPPER
# =============================================================================

#' Run hazard model and validate in one step
#'
#' @description
#' End-to-end convenience wrapper: runs `run_hazard_model()` then
#' `run_validation_suite()`. This is the recommended entry point for
#' standard validation workflows.
#'
#' The validation uses a three-tier framework:
#' - **Tier 1 (Hindcast)**: Hold out the last N years of the historical
#'   record, simulate synthetic annual maxima from the training period, and
#'   compare return levels. The simulated return level should fall within the
#'   observed return level's confidence interval (computed via delta-method
#'   or parametric bootstrap of the hurdle-GEV model).
#' - **Tier 2 (Rate check)**: Compare model storm rates (lambda) against
#'   published HURDAT2/IBTrACS climatologies for the target region.
#' - **Tier 3 (Wind field)**: Spot-check site-level winds from the Holland
#'   wind profile model against historical station observations.
#'
#' @param cfg A `hazard_cfg` object from `make_hazard_cfg()`. Controls
#'   basin, year range, search radius, and core model parameters.
#' @param targets Data frame of target locations with columns:
#'   - `location` (character): Site name (e.g., "Saba")
#'   - `lat` / `lon` (numeric): Coordinates in decimal degrees
#'   - `search_radius_km` (numeric): Radius for storm proximity filtering
#' @param validation_cfg A `validation_cfg` object from
#'   `make_validation_cfg()`. Controls holdout period, simulation size,
#'   return periods, confidence level, and output paths. Default uses 90%
#'   CI with 10-year holdout and inherited synthetic years from the hazard model output.
#' @param storm_classes Character vector of storm classes to include
#'   (default: `c("TS", "HUR")`). Must match IBTrACS classification codes.
#'   Use `"HUR"` alone for hurricane-only analysis.
#' @param climate_cfg Optional climate configuration from `make_climate_cfg()`
#'   for hazard-model climate conditioning. When `NULL`, climate conditioning
#'   is disabled.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{`out`}{Full hazard model output from `run_hazard_model()`,
#'       containing `events`, `fit`, `annual_max`, and `config` components.}
#'     \item{`val`}{Validation results from `run_validation_suite()`:
#'       - `hindcast`: Per-island return level comparison table with CI
#'         membership flags. Access via `val$hindcast$comparison`.
#'       - `rate_check`: Rate sanity check results with OK/FLAG status.
#'       - `wind_field`: Wind field spot-check results with bias estimates.
#'       - `summary`: Tibble with compact hindcast tail-diagnostic columns.}
#'     \item{`artifacts`}{Named list of saved file paths:
#'       - `plots`: Named character vector of PNG paths.
#'       - `tables`: Named character vector of CSV/MD paths.}
#'   }
#'
#' @export
#'
#' @examples
#' # --- Basic validation with defaults (90% CI) ---
#' targets_df <- data.frame(
#'   location         = c("Saba", "St. Eustatius", "St. Martin"),
#'   lat              = c(17.63, 17.49, 18.04),
#'   lon              = c(-63.23, -62.97, -63.07),
#'   search_radius_km = c(200, 200, 200)
#' )
#'
#' result <- validate_hazard_model(
#'   cfg     = make_hazard_cfg(simulation_years = 2000),
#'   targets = targets_df
#' )
#'
#' # Inspect results
#' result$val$summary                # Compact hindcast summary diagnostics
#' result$val$hindcast$comparison    # Detailed return-level table
#' result$artifacts$plots            # Paths to saved figures
#'
#' # --- Stricter validation with 95% CI ---
#' result_95 <- validate_hazard_model(
#'   cfg            = make_hazard_cfg(simulation_years = 4000),
#'   targets        = targets_df,
#'   validation_cfg = make_validation_cfg(
#'     conf_level     = 0.95,
#'     n_sim          = 10000,
#'     return_periods = c(10, 25, 50, 100)
#'   )
#' )
#'
#' # --- With climate conditioning ---
#' result_climate <- validate_hazard_model(
#'   cfg     = make_hazard_cfg(),
#'   targets = targets_df,
#'   climate_cfg = make_climate_cfg()
#' )
validate_hazard_model <- function(cfg,
                                  targets,
                                  validation_cfg = make_validation_cfg(),
                                  storm_classes = c("TS", "HUR"),
                                  climate_cfg = NULL) {
  if (!inherits(validation_cfg, "validation_cfg")) {
    stop("validation_cfg must be created by make_validation_cfg().", call. = FALSE)
  }
  storm_classes <- .normalize_storm_classes(storm_classes = storm_classes)

  out <- run_hazard_model(
    cfg = cfg,
    targets = targets,
    storm_classes = storm_classes,
    climate_cfg = climate_cfg
  )

  val <- run_validation_suite(out = out, cfg = validation_cfg)

  list(out = out, val = val, artifacts = val$artifacts)
}


# =============================================================================
# 11) PRIVATE HELPERS (plotting infrastructure)
# =============================================================================

.validate_theme <- function(base_size = 11) {
  plot_theme(base_size = base_size) +
    theme(plot.title = element_text(face = "bold", size = base_size + 1))
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

.validate_save_plot <- function(p, path, width, height, dpi = 150) {
  if (is.null(p)) return(invisible(NULL))
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  ggsave(filename = path, plot = p, width = width, height = height, dpi = dpi)
  invisible(path)
}


# =============================================================================
# 12) PLOT FUNCTIONS
# =============================================================================

#' Plot hindcast validation figures
#'
#' @description
#' Creates the retained hindcast return-level comparison figure.
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
  p_rl <- ggplot(comp, aes(x = factor(return_period))) +
    geom_point(aes(y = obs_full_rl), size = 3, color = "red", shape = 17) +
    geom_point(aes(y = sim_rl), size = 3, color = "steelblue", shape = 16) +
    geom_hline(yintercept = 64, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    facet_wrap(~ location, scales = "free_y", ncol = 3) +
    labs(
      x = "Return period (years)",
      y = "Return level \u2014 peak site wind (kt)",
      title = "Hindcast Validation: Simulated vs Observed Return Levels",
      subtitle = "Blue dot = model; Red triangle = observed (full record); dashed = 64 kt HUR threshold"
    ) +
    ggtheme
  if (any(is.finite(comp$obs_ci_lo) & is.finite(comp$obs_ci_hi))) {
    ci_pct <- attr(comp, "conf_level") %||% 0.90
    ci_label <- sprintf("%.0f%%", ci_pct * 100)
    p_rl <- p_rl +
      geom_errorbar(
        aes(ymin = obs_ci_lo, ymax = obs_ci_hi),
        width = 0.3, color = "red", linewidth = 0.8, alpha = 0.7
      ) +
      labs(
        subtitle = paste0("Blue dot = model; Red triangle [", ci_label, " CI] = observed (full record); dashed = 64 kt HUR threshold")
      )
  }

  paths <- character(0)
  paths["hindcast_return_levels"] <- file.path(out_dir, "hindcast_return_levels.png")
  .validate_save_plot(p_rl, paths[["hindcast_return_levels"]], width = 12, height = 7, dpi = 150)

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

  y_col <- if ("lambda_adj" %in% names(rc)) "lambda_adj" else "lambda_model"
  y_vals <- rc[[y_col]]

  p_rate <- ggplot(rc, aes(x = lambda_ref, y = .data[[y_col]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_abline(slope = 2, intercept = 0, linetype = "dotted", color = "grey70") +
    geom_abline(slope = 0.5, intercept = 0, linetype = "dotted", color = "grey70") +
    geom_point(aes(color = flag, shape = storm_class), size = 3.5) +
    geom_text(aes(label = location), nudge_y = 0.08, size = 3, check_overlap = TRUE) +
    scale_color_manual(values = c(
      "OK" = "forestgreen",
      "elevated" = "orange",
      "HIGH: model >> expected" = "red",
      "LOW: model << expected" = "purple",
      "slightly_low" = "goldenrod3",
      "no_reference" = "grey60"
    )) +
    coord_equal(
      xlim = c(0, max(c(y_vals, rc$lambda_ref), na.rm = TRUE) + 0.3),
      ylim = c(0, max(c(y_vals, rc$lambda_ref), na.rm = TRUE) + 0.3)
    ) +
    labs(
      x = "Reference \u03bb (published climatology)",
      y = if (identical(y_col, "lambda_adj")) "Adjusted model \u03bb" else "Model \u03bb (fitted)",
      title = "Rate Sanity Check: Adjusted Model vs Published Annual Rates",
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
#' @param out Model output list from `run_hazard_model()`.
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

  p_wf <- ggplot(wf, aes(x = obs_1min_equiv_kt, y = model_V_site_kt)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_abline(slope = 1, intercept = 20, linetype = "dotted", color = "grey70") +
    geom_abline(slope = 1, intercept = -20, linetype = "dotted", color = "grey70") +
    geom_point(aes(color = location, shape = location), size = 3.5) +
    geom_text(aes(label = storm_name), nudge_y = 4, size = 3, check_overlap = TRUE) +
    coord_equal(
      xlim = c(0, max(c(wf$obs_1min_equiv_kt, wf$model_V_site_kt), na.rm = TRUE) + 20),
      ylim = c(0, max(c(wf$obs_1min_equiv_kt, wf$model_V_site_kt), na.rm = TRUE) + 20)
    ) +
    labs(
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

  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))

  .validate_dir_create(out_dir)

  if (is.null(val$hindcast) || is.null(val$hindcast$per_island)) return(invisible(NULL))
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

  p_decomp <- ggplot(bias_long, aes(x = location, y = bias_kt, fill = source)) +
    geom_col(position = "stack", width = 0.6) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    scale_fill_manual(values = c(
      Frequency = "#E69F00", Intensity = "#56B4E9", Interaction = "#999999"
    )) +
    labs(
      x = NULL, y = "Bias contribution (kt)",
      title = "Bias Decomposition: Frequency vs Intensity",
      subtitle = "Positive = model overestimates; Negative = model underestimates",
      fill = "Source"
    ) +
    ggtheme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  paths[["bias_decomposition"]] <- file.path(out_dir, "bias_decomposition.png")
  .validate_save_plot(p_decomp, paths[["bias_decomposition"]], width = 8, height = 5, dpi = 150)

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

  p_qq <- ggplot(qq_df, aes(x = obs_quantile, y = sim_quantile)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = prob), size = 2, alpha = 0.8) +
    scale_color_viridis_c(name = "Quantile", option = "C", limits = c(0, 1)) +
    facet_wrap(~ location, scales = "free", ncol = 3) +
    labs(
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

  p_cdf <- ggplot(cdf_df, aes(x = wind_kt, y = ecdf, color = source)) +
    geom_step(linewidth = 0.7, alpha = 0.85) +
    scale_color_manual(values = c(Observed = "red", Simulated = "steelblue")) +
    geom_vline(xintercept = 64, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    facet_wrap(~ location, scales = "free_x", ncol = 3) +
    labs(
      x = "Annual maximum site wind (kt)",
      y = "Cumulative probability",
      title = "CDF Comparison: Simulated vs Observed Annual Maxima",
      subtitle = "Simulated CDF right of observed = overprediction; dashed = 64 kt HUR threshold",
      color = "Source"
    ) +
    ggtheme

  if (!is.null(gev_df)) {
    p_cdf <- p_cdf +
      geom_line(
        data = gev_df,
        aes(x = wind_kt, y = cdf_gev),
        color = "steelblue", linetype = "dotted", linewidth = 0.6,
        inherit.aes = FALSE
      )
  }

  paths <- character(0)
  paths[["exceedance_comparison"]] <- file.path(out_dir, "exceedance_comparison.png")
  .validate_save_plot(p_cdf, paths[["exceedance_comparison"]], width = 12, height = 7, dpi = 150)

  invisible(paths)
}
