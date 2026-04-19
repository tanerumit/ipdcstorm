#' Save daily hazard output to per-location CSV files
#'
#' @description
#' Writes a named list of daily hazard tables to disk, creating one CSV per
#' location. Filenames include both the scenario label and the location name.
#'
#' @param daily Named list of per-location tibbles, typically returned by
#'   \code{\link{generate_daily_hazard_impact_spatial}()}.
#' @param scenario Character scalar scenario label to embed in each filename.
#' @param out_dir Character scalar output directory. Created if it does not
#'   already exist.
#'
#' @return Named character vector of file paths, one per location.
#'
#' @examples
#' daily <- list(
#'   Saba = tibble::tibble(
#'     sim_year = 1L,
#'     doy = 1L,
#'     wind_kt = 0,
#'     surge_m = NA_real_,
#'     event_id = NA_character_,
#'     pressure_hpa = NA_real_,
#'     pressure_deficit_hpa = NA_real_,
#'     rmw_km = NA_real_,
#'     damage_intensity = 0,
#'     damage_rate = 0,
#'     cum_damage = 0
#'   )
#' )
#' save_daily_hazard_csvs(daily, scenario = "baseline", out_dir = tempdir())
#' @export
save_daily_hazard_csvs <- function(daily,
                                   scenario,
                                   out_dir = file.path("output", "raw")) {

  if (!is.list(daily) || is.data.frame(daily) || is.null(names(daily)) ||
      any(is.na(names(daily))) || any(!nzchar(names(daily)))) {
    stop("`daily` must be a named list of per-location tables.", call. = FALSE)
  }
  if (!is.character(scenario) || length(scenario) != 1L ||
      is.na(scenario) || !nzchar(scenario)) {
    stop("`scenario` must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(out_dir) || length(out_dir) != 1L ||
      is.na(out_dir) || !nzchar(out_dir)) {
    stop("`out_dir` must be a single non-empty character string.", call. = FALSE)
  }

  required_cols <- c(
    "sim_year", "doy", "wind_kt", "surge_m", "event_id",
    "pressure_hpa", "pressure_deficit_hpa", "rmw_km",
    "damage_intensity", "damage_rate", "cum_damage"
  )

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(out_dir)) {
    stop("Failed to create `out_dir`.", call. = FALSE)
  }

  sanitize <- function(x) {
    x <- gsub("[^A-Za-z0-9_]+", "_", x)
    x <- gsub("_+", "_", x)
    gsub("^_|_$", "", x)
  }

  paths <- stats::setNames(character(length(daily)), names(daily))
  scenario_id <- sanitize(scenario)

  for (loc in names(daily)) {
    tbl <- tibble::as_tibble(daily[[loc]])
    missing_cols <- setdiff(required_cols, names(tbl))
    if (length(missing_cols) > 0L) {
      stop(
        "`daily[['", loc, "']]` is missing required columns: ",
        paste(missing_cols, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    loc_id <- sanitize(loc)
    path <- file.path(out_dir, sprintf("daily_%s__%s.csv", scenario_id, loc_id))
    utils::write.csv(tbl[, required_cols, drop = FALSE], file = path, row.names = FALSE)
    paths[[loc]] <- path
  }

  paths
}
