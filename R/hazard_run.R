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

#' Run the stochastic tropical-cyclone hazard model for one or more target locations
#'
#' @description
#' Runs the end-to-end site-level hurricane hazard workflow using cleaned IBTrACS
#' North Atlantic track data and a user-supplied set of target locations.
#'
#' The function reads and filters the IBTrACS archive, selects storm track points
#' within a specified search radius of each target, computes location-specific
#' wind exposure metrics from the passing storms, aggregates those track points
#' into storm-event summaries, estimates historical annual occurrence rates by
#' storm class, calibrates dispersion for annual counts, optionally applies
#' climate-change adjustments, and then simulates synthetic annual storm counts
#' for each location over \code{cfg$n_sim} years.
#'
#' This is the main user-level entry point for the hazard model. It returns both
#' the historical intermediate products used for calibration
#' (\code{events}, \code{trackpoints}, \code{rates}, \code{fit}) and the final
#' stochastic simulation output (\code{sim}), along with resolved configuration
#' and run metadata for reproducibility.
#'
#' @details
#' The workflow consists of the following major steps:
#'
#' \enumerate{
#'   \item \strong{IBTrACS data loading and filtering.}
#'   The function resolves the IBTrACS input path from \code{cfg$data_path},
#'   reads the cleaned North Atlantic archive with
#'   \code{\link{read_ibtracs_clean}()}, and keeps all available fields needed
#'   by downstream wind-field and event-processing functions.
#'
#'   \item \strong{Target-based spatial gating.}
#'   For each target location, all track points are assigned a great-circle
#'   distance to the site using \code{\link{dist_to_target}()}. Only points
#'   within \code{cfg$search_radius_km} are retained for that location. This is
#'   a site-centred filtering step: the same storm may be retained for one
#'   island and excluded for another, depending on distance.
#'
#'   \item \strong{Site wind estimation.}
#'   The retained track points are passed to
#'   \code{\link{compute_site_winds_full}()}, which estimates the wind at the
#'   target site from storm intensity, geometry, and wind-radii information.
#'   The returned track-point table includes location-specific wind metrics such
#'   as asymmetric site wind, symmetric site wind, storm maximum wind, and
#'   bearing to target.
#'
#'   \item \strong{Storm-event construction and classification.}
#'   Filtered track points are aggregated into per-storm event summaries using
#'   \code{\link{make_storm_events}()}. Each event is then classified into the
#'   requested storm classes with \code{\link{classify_severity}()} using the
#'   thresholds stored in \code{cfg$advanced}, currently tropical-storm and
#'   hurricane thresholds in knots.
#'
#'   \item \strong{Historical annual rates and dispersion calibration.}
#'   The event catalogue is reduced to annual counts with
#'   \code{\link{compute_annual_counts}()}. Mean annual rates by storm class are
#'   estimated with \code{\link{compute_lambda_table}()}, and over-dispersion of
#'   annual total activity is estimated with \code{\link{estimate_k_hat}()}.
#'   These provide the frequency inputs for the stochastic simulation.
#'
#'   \item \strong{Rate-scaling adjustment.}
#'   Before simulation, the function computes location/class-specific lambda
#'   scalers from the calibrated rate table using the configured
#'   \code{lambda_scaling_mode}. This can downscale or, depending on mode,
#'   upscale location/class rates relative to the reference calibration target.
#'   The applied scalers and their identifier are returned.
#'
#'   \item \strong{Optional climate adjustments.}
#'   If \code{climate} inherits from class \code{"climate_cfg"}, the climate
#'   inputs are resolved with \code{resolve_climate_inputs()}. This produces a
#'   scalar \code{delta_sst}, a frequency sensitivity \code{beta_sst}, a
#'   hurricane-fraction sensitivity \code{gamma}, and optional storm
#'   perturbation settings. These affect the stochastic simulation but do not
#'   rewrite the historical event catalogue itself.
#'
#'   \item \strong{Simulation of synthetic annual counts.}
#'   For each location, the rate table is adjusted, the island-specific base
#'   hurricane fraction is computed, and
#'   \code{\link{simulate_twolevel_counts}()} is used to generate
#'   \code{cfg$n_sim} synthetic years of annual storm activity. The simulation
#'   produces counts such as total storms, tropical storms, and hurricanes by
#'   year and location.
#' }
#'
#' Notes:
#'
#' \itemize{
#'   \item \strong{Historical sample versus synthetic output.}
#'   The historical component is built from IBTrACS North Atlantic storms gated
#'   around each site. The synthetic component does not resample full tracks;
#'   it simulates annual storm counts conditional on the calibrated local rate
#'   structure.
#'
#'   \item \strong{Storm classes.}
#'   The model currently uses normalized storm-class labels such as
#'   \code{"TS"} and \code{"HUR"}. These are assigned from site-level peak wind
#'   at the target location, not basin-wide best-track status labels.
#'
#'   \item \strong{Climate interpretation.}
#'   When \code{climate = NULL}, the climate module is fully disabled:
#'   \code{delta_sst = 0}, \code{beta_sst = 0}, \code{gamma = 0}, and storm
#'   perturbation is off. When a baseline climate configuration is supplied
#'   through \code{make_climate_cfg(...)} with \code{delta_sst = 0}, the
#'   climate resolver still runs and records climate metadata, but the climate
#'   shift itself is zero.
#'
#'   \item \strong{Reproducibility metadata.}
#'   The returned object includes \code{run_metadata} with the seed, IBTrACS
#'   file identifier, row count, derived data identifier, parameter hash, and
#'   lambda-scaling mode.
#' }
#'
#' @param cfg A \code{hazard_cfg} object created by
#'   \code{\link{make_hazard_cfg}()}. This controls the IBTrACS source file,
#'   target search radius, simulation length, start year for the calibration
#'   sample, and advanced storm/wind/rate settings.
#'
#' @param targets A data frame or tibble of target locations with required
#'   columns:
#'   \describe{
#'     \item{\code{name}}{Character location name used as the location key in
#'       output tables and lists.}
#'     \item{\code{lat}}{Numeric latitude in decimal degrees.}
#'     \item{\code{lon}}{Numeric longitude in decimal degrees.}
#'   }
#'   Each row defines one site for which a separate hazard calibration and
#'   stochastic simulation will be performed.
#'
#' @param storm_classes Character vector of storm classes to model. These are
#'   normalized internally and currently intended to include \code{"TS"} and
#'   \code{"HUR"}. Only the requested classes are carried into the annual-count
#'   and rate workflow.
#'
#' @param climate Optional climate configuration object produced by
#'   \code{\link{make_climate_cfg}()} and inheriting from class
#'   \code{"climate_cfg"}.
#'
#'   Use a baseline climate run with \code{delta_sst = 0} when you want the
#'   climate resolver and its metadata to remain active but represent present-day
#'   conditions. Use \code{NULL} only when you explicitly want the climate
#'   module disabled.
#'
#'   When \code{NULL}, the run uses:
#'   \itemize{
#'     \item \code{delta_sst = 0}
#'     \item \code{beta_sst = 0}
#'     \item \code{gamma = 0}
#'     \item storm perturbation disabled
#'   }
#'
#' @param seed Optional integer scalar random seed. If \code{NULL}, a seed is
#'   generated internally, set at function entry, and recorded in
#'   \code{run_metadata$seed}. All stochastic simulation in the run inherits
#'   from this seed.
#'
#' @param verbose Logical; if \code{TRUE} (default), print progress, summary
#'   statistics, climate settings, and run metadata to the console.
#'
#' @return
#' A named list with the following elements:
#' \describe{
#'   \item{\code{sim}}{Tibble of synthetic annual storm counts for all simulated
#'     years and locations. Contains one row per simulation year and location.}
#'   \item{\code{events}}{Tibble of historical storm-event summaries across all
#'     locations after filtering, event construction, and classification.}
#'   \item{\code{trackpoints}}{Named list of tibbles, one per location,
#'     containing filtered IBTrACS track points and computed site-wind
#'     diagnostics.}
#'   \item{\code{rates}}{Tibble of calibrated mean annual rates
#'     (\code{lambda}) by location and storm class.}
#'   \item{\code{lambda_scalers}}{Tibble of location/class-specific lambda
#'     scaling factors used to adjust the calibrated rate table before
#'     simulation.}
#'   \item{\code{lambda_scaler_id}}{Character identifier summarizing the applied
#'     lambda-scaler configuration.}
#'   \item{\code{fit}}{Tibble of fitted dispersion and summary parameters by
#'     location, with climate-related attributes attached.}
#'   \item{\code{cfg}}{Resolved hazard configuration used for the run, including
#'     normalized data path and climate metadata.}
#'   \item{\code{config}}{Duplicate of \code{cfg}, returned for compatibility.}
#'   \item{\code{run_metadata}}{List with reproducibility metadata including
#'     \code{seed}, \code{ibtracs_file}, \code{ibtracs_rows},
#'     \code{ibtracs_data_id}, \code{parameter_id}, and
#'     \code{lambda_scaling_mode}.}
#' }
#'
#' @seealso
#' \code{\link{make_hazard_cfg}},
#' \code{\link{make_climate_cfg}},
#' \code{\link{read_ibtracs_clean}},
#' \code{\link{compute_site_winds_full}},
#' \code{\link{make_storm_events}},
#' \code{\link{compute_annual_counts}},
#' \code{\link{compute_lambda_table}},
#' \code{\link{estimate_k_hat}},
#' \code{\link{simulate_twolevel_counts}}
#'
#' @examples
#' \dontrun{
#' targets <- tibble::tribble(
#'   ~name,        ~lat,     ~lon,
#'   "Saba",        17.63,  -63.23,
#'   "St_Martin",   18.07,  -63.05
#' )
#'
#' cfg <- make_hazard_cfg(
#'   data_path = "inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv",
#'   start_year = 1980,
#'   n_sim = 1000
#' )
#'
#' # Historical/stationary run with climate module disabled
#' out <- run_hazard_model(
#'   cfg = cfg,
#'   targets = targets,
#'   storm_classes = c("TS", "HUR"),
#'   climate = NULL,
#'   seed = 123
#' )
#'
#' out$sim
#' out$rates
#' out$run_metadata
#'
#' # Baseline climate run with climate resolver active and delta_sst = 0
#' climate_cfg <- make_climate_cfg(
#'   scenario = "baseline",
#'   delta_sst = 0
#' )
#'
#' out_climate <- run_hazard_model(
#'   cfg = cfg,
#'   targets = targets,
#'   climate = climate_cfg,
#'   seed = 123
#' )
#' }
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
  if (!is.null(climate) && !inherits(climate, "climate_cfg")) {
    stop("climate must be NULL or inherit from \"climate_cfg\".", call. = FALSE)
  }

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
  climate_cfg <- NULL
  climate_resolved <- NULL
  delta_sst <- 0
  beta_sst <- 0
  gamma_intensity <- 0
  perturb_cfg <- NULL
  if (inherits(climate, "climate_cfg")) {
    climate_cfg <- climate
    future_period <- seq.int(
      from = as.integer(climate_cfg$start_year),
      length.out = 30L
    )
    climate_resolved <- resolve_climate_inputs(
      climate_cfg = climate_cfg,
      annual_counts = annual_counts_all,
      lambda_table = rates_all |>
        dplyr::group_by(.data$storm_class) |>
        dplyr::summarise(lambda = sum(.data$lambda), .groups = "drop"),
      min_year = cfg$start_year,
      future_period = future_period,
      verbose = verbose
    )
    delta_sst <- climate_resolved$delta_sst
    beta_sst <- climate_resolved$beta_sst
    gamma_intensity <- climate_resolved$gamma
    perturb_cfg <- climate_resolved$perturb
    climate_info <- list(
      enabled = climate_resolved$enabled,
      scenario = climate_resolved$scenario,
      source = climate_resolved$source,
      delta_sst = delta_sst,
      baseline_years = climate_resolved$baseline_years,
      future_period = climate_resolved$future_period,
      sensitivity_mode = climate_resolved$sensitivity_mode,
      k_beta = climate_resolved$k_beta,
      k_gamma = climate_resolved$k_gamma,
      beta_0 = climate_resolved$beta_0,
      gamma_0 = climate_resolved$gamma_0,
      beta_sst = beta_sst,
      gamma = gamma_intensity,
      p_hur_base = climate_resolved$p_hur_base,
      perturb = perturb_cfg,
      perturb_state = climate_resolved$perturb_state,
      climate_mode = climate_resolved$climate_mode,
      config = climate_cfg
    )
  }

  if (verbose) {
    .fmt_scenario_label <- function(x) {
      if (!is.character(x) || length(x) != 1L || !nzchar(x)) {
        return("Stationary")
      }
      if (grepl("^ssp", x, ignore.case = TRUE)) {
        return(toupper(x))
      }
      x
    }
    .fmt_perturb_status <- function(perturb_state, perturb_cfg) {
      if (identical(perturb_state, "default")) {
        return("enabled (defaults)")
      }
      if (identical(perturb_state, "custom")) {
        return("enabled (custom)")
      }
      if (is.null(perturb_cfg)) {
        return("disabled")
      }
      "enabled"
    }
    .cli_h("Climate")
    if (!is.null(climate_info)) {
      bl_range <- climate_info$baseline_years
      scenario_label <- .fmt_scenario_label(climate_info$scenario)
      perturb_state <- climate_info$perturb_state
      rate_pct <- 100 * (exp(beta_sst) - 1)
      intensity_pct <- 100 * gamma_intensity
      .cli_info(sprintf(
        "Climate mode      : %s",
        switch(
          climate_info$climate_mode,
          baseline = "baseline climate run (delta_sst = 0)",
          future = sprintf("future climate run (delta_sst = %+.2f C)", delta_sst),
          "climate run"
        )
      ))
      .cli_info(sprintf("Climate scenario  : %s", scenario_label))
      .cli_info(sprintf("SST baseline     : %d-%d", min(bl_range), max(bl_range)))
      .cli_info(sprintf("Sensitivity mode : %s", climate_info$sensitivity_mode))
      .cli_info(sprintf(
        "Rate effect      : beta_0=%.2f -> beta=%.2f (%+.0f%% per +1\u00B0C)",
        climate_info$beta_0, beta_sst, rate_pct
      ))
      .cli_info(sprintf(
        "Intensity effect : gamma_0=%.3f -> gamma=%.3f (%+.0f%% per +1\u00B0C)",
        climate_info$gamma_0, gamma_intensity, intensity_pct
      ))
      .cli_info(sprintf(
        "Storm perturbation: %s",
        .fmt_perturb_status(perturb_state, perturb_cfg)
      ))
    } else {
      .cli_info("Climate mode      : off (no climate module)")
      .cli_info("Climate scenario  : Off")
      .cli_info("Rate effect      : +0% activity (\u03B2 = 0.00)")
      .cli_info("Intensity effect : +0% hurricane fraction (\u03B3 = 0.000)")
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
  attr(fit_all, "climate_mode") <- if (!is.null(climate_info)) climate_info$climate_mode else "off"
  attr(fit_all, "climate_scenario") <- if (!is.null(climate_info)) climate_info$scenario else "off"
  attr(fit_all, "climate_source") <- if (!is.null(climate_info) && !is.null(climate_info$source)) climate_info$source else "off"
  attr(fit_all, "climate_cfg") <- climate_cfg

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
