# =============================================================================
# Script overview: model orchestration
# - make_hazard_cfg(): user-facing hazard config constructor with presets.
# - print.hazard_cfg(): concise human-readable config summary.
# - run_hazard_model(): end-to-end hazard workflow across multiple targets.
# =============================================================================

# =============================================================================
# 0) Console output helpers (internal)
# =============================================================================

#' Print a styled section header
#' @keywords internal
.cli_h <- function(title, width = 60) {
  pad <- max(2, width - nchar(title) - 2)
  message(sprintf("\n\u2500\u2500 %s %s", title, strrep("\u2500", pad)))
}

#' Print an indented info line
#' @keywords internal
.cli_info <- function(...) {
  message(sprintf("   %s", paste0(...)))
}

#' Print a success line with checkmark
#' @keywords internal
.cli_ok <- function(...) {
  message(sprintf("   \u2714 %s", paste0(...)))
}

#' Format a number with comma separators
#' @keywords internal
.fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

#' Resolve an IBTrACS path in dev and installed package contexts
#' @keywords internal
.resolve_ibtracs_path <- function(path) {
  path <- as.character(path)[1]
  if (is.na(path) || !nzchar(path)) {
    stop("cfg$data_path must be a non-empty character path.", call. = FALSE)
  }
  if (file.exists(path)) {
    return(path)
  }

  base_name <- basename(path)
  candidates <- c(
    file.path("inst", "extdata", "ibtracs", base_name),
    file.path("inst", "extdata", base_name),
    system.file("extdata", "ibtracs", base_name, package = "ipdcstorm"),
    system.file("extdata", base_name, package = "ipdcstorm")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0) {
    return(existing[[1]])
  }

  stop("File not found: ", path, call. = FALSE)
}

#' Normalize user-supplied targets to the required schema
#' @keywords internal
.normalize_hazard_targets <- function(targets) {
  if (!is.data.frame(targets)) {
    stop("targets must be a data frame with columns name, lat, lon.", call. = FALSE)
  }
  if (!("name" %in% names(targets)) && "location" %in% names(targets)) {
    targets <- dplyr::rename(targets, name = "location")
  }

  required_cols <- c("name", "lat", "lon")
  missing_cols <- setdiff(required_cols, names(targets))
  if (length(missing_cols) > 0) {
    stop(
      "targets must contain columns: ",
      paste(required_cols, collapse = ", "),
      ". Missing: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  targets
}

# =============================================================================
# 1) Hazard configuration
# =============================================================================

#' Create a hazard model configuration
#'
#' @description
#' Creates a user-facing configuration object for `run_hazard_model()`.
#' Most users only need to set the IBTrACS path and simulation horizon.
#'
#' The `preset = "default"` setting applies standard threshold and cap values.
#' Expert tuning knobs can be provided via `advanced`, but are optional.
#'
#' @param data_path Character; path to IBTrACS CSV input data.
#' @param search_radius_km Numeric; maximum distance from each target used to
#'   include track points.
#' @param historical_start_year Integer; first historical year used to fit rates.
#' @param simulation_years Integer; number of synthetic years to simulate.
#' @param preset Character; currently only `"default"`.
#' @param advanced Optional named list of expert parameters. Most users should
#'   leave this as `NULL`. Supported names:
#'   `ts_threshold_kt`, `strong_storm_threshold_kt`, `hurricane_threshold_kt`,
#'   `r34_cap_nm`, `r50_cap_nm`, `r64_cap_nm`,
#'   `lambda_scaling_mode` (`"target"` or `"down_only"`).
#'
#' @return A list with class `c("hazard_cfg", "list")`.
#' @export
make_hazard_cfg <- function(data_path = "data/ibtracs/ibtracs.NA.list.v04r01.csv",
                            search_radius_km = 800,
                            historical_start_year = 1970L,
                            simulation_years = 1000L,
                            preset = "default",
                            advanced = NULL) {

  preset <- match.arg(preset, choices = c("default"))

  defaults <- list(
    ts_threshold_kt = 34,
    strong_storm_threshold_kt = 50,
    hurricane_threshold_kt = 64,
    r34_cap_nm = 600,
    r50_cap_nm = 400,
    r64_cap_nm = 250,
    lambda_scaling_mode = "target"
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

  cfg <- list(
    preset = preset,
    data_path = as.character(data_path),
    search_radius_km = as.numeric(search_radius_km),
    start_year = as.integer(historical_start_year),
    n_sim = as.integer(simulation_years),
    advanced = advanced,
    resampling_method = "stratified",
    copula_min_n = 30L,
    copula_k = 1L,
    copula_robust_scale = TRUE
  )
  class(cfg) <- c("hazard_cfg", "list")
  cfg$advanced$lambda_scaling_mode <- .normalize_lambda_scaling_mode(cfg$advanced$lambda_scaling_mode)
  cfg
}

#' @export
print.hazard_cfg <- function(x, ...) {
  cat(sprintf("Hazard configuration (preset: \"%s\")\n", x$preset))
  cat(sprintf("  IBTrACS data  : %s\n", x$data_path))
  cat(sprintf("  Study period  : %d - present\n", x$start_year))
  cat(sprintf("  Search radius : %s km\n", format(round(x$search_radius_km, 2), trim = TRUE)))
  cat(sprintf("  Simulation    : %s synthetic years\n", format(x$n_sim, big.mark = ",", scientific = FALSE, trim = TRUE)))
  cat(sprintf(
    "  Thresholds    : WMO standard (TS >= %s kt, Hurricane >= %s kt) [preset]\n",
    format(x$advanced$ts_threshold_kt, trim = TRUE),
    format(x$advanced$hurricane_threshold_kt, trim = TRUE)
  ))
  cat(sprintf("  Lambda scaling: %s\n", x$advanced$lambda_scaling_mode))
  invisible(x)
}

# =============================================================================
# 8) Orchestrator
# =============================================================================

#' Run hazard model for multiple target locations
#'
#' @description
#' End-to-end workflow:
#' 1) Read IBTrACS and filter to basin/season
#' 2) Gate track points within cfg$search_radius_km of each target
#' 3) Compute site winds and storm event summaries
#' 4) Classify storm classes and compute annual rates
#' 5) Simulate annual counts with optional stationary climate-shift adjustments
#'
#' @param cfg Hazard configuration from `make_hazard_cfg()`.
#' @param targets Data frame/tibble with columns: name, lat, lon.
#' @param storm_classes Character vector of storm classes to model (default TS and HUR).
#' @param climate Optional climate input from `make_climate_input()`. When
#'   `NULL`, the run is stationary baseline with `delta_sst = 0`,
#'   `beta_sst = 0`, `gamma = 0`, and no Level 3 perturbation.
#' @param seed Optional integer seed. When `NULL`, `run_hazard_model()`
#'   generates one, sets it at entry, and records it in `run_metadata$seed`.
#' @param verbose Logical; if TRUE (default), print progress to console.
#'
#' @return A list with elements:
#'   \code{sim}, \code{events}, \code{trackpoints}, \code{rates}, \code{fit},
#'   and \code{cfg}. The same configuration is also returned in \code{config}.
#'
#' @export
run_hazard_model <- function(cfg, targets,
                             storm_classes = c("TS", "HUR"),
                             climate = NULL,
                             seed = NULL,
                             verbose = TRUE) {
  if (!inherits(cfg, "hazard_cfg")) {
    stop("cfg must be created by make_hazard_cfg().", call. = FALSE)
  }
  storm_classes <- .normalize_storm_classes(storm_classes = storm_classes)
  targets <- .normalize_hazard_targets(targets)

  ts_threshold_kt <- as.numeric(cfg$advanced$ts_threshold_kt)
  hurricane_threshold_kt <- as.numeric(cfg$advanced$hurricane_threshold_kt)
  lambda_scaling_mode <- .normalize_lambda_scaling_mode(cfg$advanced$lambda_scaling_mode)
  data_path <- .resolve_ibtracs_path(cfg$data_path)
  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
  } else {
    if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed)) {
      stop("seed must be NULL or a single finite numeric value.", call. = FALSE)
    }
    seed <- as.integer(seed)
  }
  set.seed(seed)

  if (!is.null(climate) && !inherits(climate, "climate_input")) {
    stop("climate must be NULL or created by make_climate_input().", call. = FALSE)
  }

  # --- 1) Load IBTrACS (suppress sub-function messages) --------------------
  if (verbose) .cli_h("Loading data")

  ib_sub <- read_ibtracs_clean(
    ibtracs_csv = data_path,
    basin = "NA",
    season = NULL,
    keep_all = TRUE,
    verbose = FALSE
  )

  yr_range <- range(as.integer(ib_sub$SEASON), na.rm = TRUE)
  if (verbose) {
    .cli_info(sprintf("Input file       : %s", basename(data_path)))
    .cli_info(sprintf(
      "Raw NA input     : %s track points (%d\u2013%d)",
      .fmt_n(nrow(ib_sub)), yr_range[1], yr_range[2]
    ))
    .cli_info(sprintf("Model start year : %d", cfg$start_year))
  }

  # --- 2) Filter storms per location --------------------------------------
  if (verbose) .cli_h("Target/location filtering")

  trackpoints_list    <- setNames(vector("list", nrow(targets)), targets$name)
  events_list         <- setNames(vector("list", nrow(targets)), targets$name)
  annual_counts_list  <- setNames(vector("list", nrow(targets)), targets$name)
  rates_list          <- setNames(vector("list", nrow(targets)), targets$name)
  sim_list            <- setNames(vector("list", nrow(targets)), targets$name)
  fit_list            <- setNames(vector("list", nrow(targets)), targets$name)

  # Collect per-location summary for table display
  loc_summary <- list()

  for (i in seq_len(nrow(targets))) {
    loc <- targets[i, ]
    loc_name <- as.character(loc$name)

    dat_loc <- ib_sub |>
      dplyr::mutate(dist_km = dist_to_target(.data$lat, .data$lon, loc$lat, loc$lon)) |>
      dplyr::filter(.data$dist_km <= cfg$search_radius_km)

    if (nrow(dat_loc) == 0) {
      warning("No trackpoints within search_radius_km for location: ", loc_name, " (search_radius_km=", cfg$search_radius_km, "). Skipping.")
      next
    }

    dat_loc <- compute_site_winds_full(dat_loc, loc$lat, loc$lon)
    trackpoints_list[[loc_name]] <- dat_loc |>
      dplyr::rename(
        site_wind_kt = "V_site_kt",
        site_wind_symmetric_kt = "V_site_symmetric_kt",
        storm_wind_kt = "Vmax_kt",
        bearing_deg = "bearing_to_target"
      ) |>
      dplyr::select(-dplyr::any_of("R34_missing"))

    ev <- make_storm_events(dat_loc) |>
      dplyr::mutate(
        location = loc_name,
        storm_class = classify_severity(
          .data$peak_wind_kt,
          ts_threshold_kt = ts_threshold_kt,
          hurricane_threshold_kt = hurricane_threshold_kt
        )
      ) |>
      dplyr::relocate("location", .before = "storm_id")

    # Collect summary stats for table
    n_ts  <- sum(ev$peak_wind_kt >= ts_threshold_kt, na.rm = TRUE)
    n_hur <- sum(ev$peak_wind_kt >= hurricane_threshold_kt, na.rm = TRUE)
    loc_summary[[loc_name]] <- list(n_ts = n_ts, n_hur = n_hur)

    ev <- ev |>
      dplyr::filter(.data$year >= cfg$start_year) |>
      dplyr::filter(is.finite(.data$peak_wind_kt)) |>
      dplyr::filter(.data$storm_class != "unknown")

    events_list[[loc_name]] <- ev

    ac <- compute_annual_counts(ev, storm_classes = storm_classes)
    lt <- compute_lambda_table(ac)
    kinfo <- estimate_k_hat(ac)

    annual_counts_list[[loc_name]] <- ac |>
      dplyr::mutate(location = loc_name) |>
      dplyr::relocate("location", .before = "year")

    rates_list[[loc_name]] <- lt |>
      dplyr::mutate(location = loc_name) |>
      dplyr::relocate("location", .before = "storm_class")

    fit_list[[loc_name]] <- tibble::tibble(
      location = loc_name,
      k_hat = kinfo$k_hat,
      mean_annual_total = kinfo$mu,
      var_annual_total = kinfo$var
    )
  }

  # Print location summary table
  if (verbose && length(loc_summary) > 0) {
    max_name_len <- max(nchar(names(loc_summary)))
    col_w <- max(max_name_len, 10)
    .cli_info(sprintf("%-*s  %5s  %5s", col_w, "Location", "TS", "HUR"))
    for (nm in names(loc_summary)) {
      s <- loc_summary[[nm]]
      .cli_info(sprintf("%-*s  %5d  %5d", col_w, nm, s$n_ts, s$n_hur))
    }
  }

  events_all        <- dplyr::bind_rows(Filter(Negate(is.null), events_list))
  annual_counts_all <- dplyr::bind_rows(Filter(Negate(is.null), annual_counts_list))
  rates_all         <- dplyr::bind_rows(Filter(Negate(is.null), rates_list))
  fit_all           <- dplyr::bind_rows(Filter(Negate(is.null), fit_list))
  if (verbose && nrow(events_all) > 0) {
    sample_years <- range(events_all$year, na.rm = TRUE)
    .cli_info(sprintf(
      "Model sample     : %s target-event records (%d\u2013%d) across %d locations",
      .fmt_n(nrow(events_all)),
      sample_years[1],
      sample_years[2],
      dplyr::n_distinct(events_all$location)
    ))
  }
  if (verbose) .cli_h("Rate calibration")
  rate_check_seed <- .build_rate_check_table(
    list(rates = rates_all),
    lambda_scaling_mode = lambda_scaling_mode
  )
  lambda_scalers <- .lambda_scalers_from_rate_check(
    rate_check_seed,
    scaling_mode = lambda_scaling_mode
  )
  lambda_scaler_id <- .lambda_scaler_id(lambda_scalers)

  if (verbose && nrow(lambda_scalers) > 0) {
    n_clamped <- sum(lambda_scalers$scale_clamped, na.rm = TRUE)
    n_missing_ref <- sum(lambda_scalers$scale_status == "no_ref", na.rm = TRUE)
    .cli_info(sprintf(
      "Adjustments      : %d location/class adjustments (mode=%s, clamped=%d, missing_ref=%d)",
      nrow(lambda_scalers), lambda_scaling_mode, n_clamped, n_missing_ref
    ))
  }
  if (verbose && identical(lambda_scaling_mode, "target") && any(lambda_scalers$was_upscaled, na.rm = TRUE)) {
    .cli_info(
      "Some location/class rates were increased to match reference rates; ",
      "set advanced$lambda_scaling_mode='down_only' to prevent upscaling."
    )
  }

  # =========================================================================
  # CLIMATE ADJUSTMENTS
  # =========================================================================
  climate_info <- NULL
  delta_sst <- 0
  beta_sst <- 0
  gamma_intensity <- 0
  perturb_cfg <- NULL

  if (!is.null(climate)) {
    climate_info <- climate
    delta_sst <- climate$shift$delta_sst
    beta_sst <- climate$response$beta_sst
    gamma_intensity <- climate$response$gamma
    perturb_cfg <- climate$response$perturb
  }

  if (verbose) {
    .cli_h("Climate")
    if (!is.null(climate_info)) {
      bl_range <- climate$shift$baseline_years
      .cli_info(sprintf("delta_sst         : %.3f", delta_sst))
      .cli_info(sprintf("Baseline years    : %d-%d", min(bl_range), max(bl_range)))
      .cli_info(sprintf("beta_sst          : %.3f", beta_sst))
      .cli_info(sprintf("gamma             : %.3f", gamma_intensity))
      .cli_info(sprintf("Storm perturbation: %s", if (is.null(perturb_cfg)) "disabled" else "enabled"))
    } else {
      .cli_info("delta_sst         : 0.000")
      .cli_info("beta_sst          : 0.000")
      .cli_info("gamma             : 0.000")
      .cli_info("Storm perturbation: disabled")
    }
    .cli_info(sprintf("seed              : %s", as.character(seed)))
  }

  # =========================================================================
  # Simulate annual counts (with climate modifications)
  # =========================================================================
  if (verbose) {
    .cli_h("Simulation")
    .cli_info(sprintf("Synthetic years  : %s", .fmt_n(cfg$n_sim)))
  }

  sim_list <- setNames(vector("list", nrow(targets)), targets$name)

  for (loc_name in names(Filter(Negate(is.null), fit_list))) {
    lt <- rates_list[[loc_name]] |>
      dplyr::select(-"location")
    lt_sim <- .apply_lambda_scalers_to_lambda_table(
      lambda_table = lt,
      location = loc_name,
      lambda_scalers = lambda_scalers
    )
    k <- fit_list[[loc_name]]$k_hat

    p_hur_island <- compute_p_hur_base(lt_sim)

    sim <- simulate_twolevel_counts(
      lambda_table = lt_sim,
      k_hat = k,
      n_years_sim = cfg$n_sim,
      delta_sst = delta_sst,
      beta_sst = beta_sst,
      gamma_intensity = gamma_intensity,
      p_hur_base = p_hur_island
    )

    sim_list[[loc_name]] <- sim |>
      dplyr::mutate(location = loc_name) |>
      dplyr::relocate("location", .before = "sim_year")
  }

  sim_all <- dplyr::bind_rows(Filter(Negate(is.null), sim_list))

  # Print simulation summary
  if (verbose) {
    if (beta_sst != 0 || delta_sst != 0) {
      rs <- exp(beta_sst * delta_sst)
      .cli_info(sprintf("Rate scaling      : %.3fx", rs))
    }
    .cli_ok("Done")
  }

  fit_all <- fit_all |>
    dplyr::mutate(
      beta_sst = beta_sst,
      gamma_intensity = gamma_intensity,
      p_hurricane_base = NA_real_
    )
  attr(fit_all, "perturb") <- perturb_cfg
  attr(fit_all, "cc_params") <- perturb_cfg
  attr(fit_all, "delta_sst") <- delta_sst

  cfg_out <- cfg
  cfg_out$advanced$lambda_scaling_mode <- lambda_scaling_mode
  cfg_out$data_path <- data_path
  cfg_out$climate <- climate_info

  data_file <- basename(data_path)
  data_rows <- nrow(ib_sub)
  data_id <- paste0(data_file, "|rows=", format(data_rows, scientific = FALSE, trim = TRUE))
  param_fields <- c(
    cfg$preset,
    data_path,
    cfg$search_radius_km,
    cfg$start_year,
    cfg$n_sim,
    cfg$advanced$ts_threshold_kt,
    cfg$advanced$hurricane_threshold_kt,
    cfg$advanced$r34_cap_nm,
    cfg$advanced$r50_cap_nm,
    cfg$advanced$r64_cap_nm,
    lambda_scaling_mode
  )
  parameter_id <- .checksum_id_from_text(param_fields, prefix = "params")
  run_metadata <- list(
    seed = seed,
    ibtracs_file = data_file,
    ibtracs_rows = as.integer(data_rows),
    ibtracs_data_id = data_id,
    parameter_id = parameter_id,
    lambda_scaling_mode = lambda_scaling_mode
  )
  if (verbose) {
    meta_parts <- c(
      sprintf("seed=%s", as.character(run_metadata$seed)),
      sprintf("lambda_mode=%s", run_metadata$lambda_scaling_mode)
    )
    .cli_h("Run metadata")
    .cli_info(paste(meta_parts, collapse = " "))
  }

  list(
    sim = sim_all,
    events = events_all,
    trackpoints = trackpoints_list,
    rates = rates_all,
    lambda_scalers = lambda_scalers,
    lambda_scaler_id = lambda_scaler_id,
    fit = fit_all,
    cfg = cfg_out,
    config = cfg_out,
    run_metadata = run_metadata
  )
}
