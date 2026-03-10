
# =============================================================================
# Script overview: model orchestration
# - run_hazard_model(): end-to-end hazard workflow across multiple targets.
# =============================================================================

# =============================================================================
# 8) Orchestrator
# =============================================================================

#' Run hazard model for multiple targets (islands/locations)
#'
#' @description
#' End-to-end workflow:
#' 1) Read IBTrACS and filter to basin/season
#' 2) Gate track points within cfg$gate_km of each target
#' 3) Compute site winds (V_site_kt) and storm event summaries
#' 4) Classify severities and compute annual rates
#' 5) Simulate annual counts with a shared activity factor (two-level model)
#'
#' @param cfg List with required fields:
#'   ibtracs_path, gate_km, thr_ts, thr_hur, min_year, n_years_sim.
#' @param targets Data frame/tibble with columns: name, lat, lon.
#' @param per_target_cfg Named list of per-target threshold overrides (optional).
#' @param severities Character vector of event severities to model (default TS and HUR).
#'
#' @return A list with elements:
#'   trackpoints (list by island),
#'   events_by_island (list),
#'   events_all, annual_counts_all, lambda_all, sim_all, kinfo_all.
#'
#' @export
run_hazard_model <- function(cfg, targets, per_target_cfg = list(), severities = c("TC", "HUR")) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` is required.")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Package `tibble` is required.")
  
  ib_sub <- read_ibtracs_clean(
    ibtracs_csv = cfg$ibtracs_path,
    basin = "NA",
    season = NULL,
    keep_all = TRUE,
    verbose = TRUE
  )
  
  trackpoints_list    <- setNames(vector("list", nrow(targets)), targets$name)
  events_list         <- setNames(vector("list", nrow(targets)), targets$name)
  annual_counts_list  <- setNames(vector("list", nrow(targets)), targets$name)
  lambda_list         <- setNames(vector("list", nrow(targets)), targets$name)
  sim_list            <- setNames(vector("list", nrow(targets)), targets$name)
  kinfo_list          <- setNames(vector("list", nrow(targets)), targets$name)
  
  for (i in seq_len(nrow(targets))) {
    loc <- targets[i, ]
    island <- as.character(loc$name)
    
    thr <- list(
      thr_ts  = cfg$thr_ts,
      thr_hur = cfg$thr_hur,
      thr_port  = NA_real_,
      thr_infra = NA_real_
    )
    if (!is.null(per_target_cfg[[island]])) {
      thr[names(per_target_cfg[[island]])] <- per_target_cfg[[island]]
    }
    
    dat_loc <- ib_sub |>
      dplyr::mutate(dist_km = dist_to_target(.data$lat, .data$lon, loc$lat, loc$lon)) |>
      dplyr::filter(.data$dist_km <= cfg$gate_km)
    
    if (nrow(dat_loc) == 0) {
      warning("No trackpoints within gate_km for island: ", island, " (gate_km=", cfg$gate_km, "). Skipping.")
      next
    }
    
    dat_loc <- compute_site_winds_full(dat_loc, loc$lat, loc$lon)
    trackpoints_list[[island]] <- dat_loc
    
    ev <- make_storm_events(dat_loc) |>
      dplyr::mutate(island = island) |>
      dplyr::relocate("island", .before = "SID") |>
      dplyr::mutate(
        severity = classify_severity(.data$V_site_max_kt, thr_ts = cfg$thr_ts, thr_hur = cfg$thr_hur)
      )
    
    message("[", island, "] storms=", nrow(ev),
            " finite V_site_max=", sum(is.finite(ev$V_site_max_kt)),
            " any_Vsite_known=", sum(isTRUE(ev$any_Vsite_known)),
            " >=34=", sum(ev$V_site_max_kt >= cfg$thr_ts, na.rm=TRUE),
            " >=64=", sum(ev$V_site_max_kt >= cfg$thr_hur, na.rm=TRUE))
    
    
    ev <- ev |>
      dplyr::filter(.data$year >= cfg$min_year) |>
      dplyr::filter(is.finite(.data$V_site_max_kt)) |>
      dplyr::filter(.data$severity != "unknown")
    
    events_list[[island]] <- ev
    
    ac <- compute_annual_counts(ev, severities = severities)
    lt <- compute_lambda_table(ac)
    kinfo <- estimate_k_hat(ac)
    sim <- simulate_twolevel_counts(lt, kinfo$k_hat, n_years_sim = cfg$n_years_sim)
    
    annual_counts_list[[island]] <- ac |>
      dplyr::mutate(island = island) |>
      dplyr::relocate("island", .before = "year")
    
    lambda_list[[island]] <- lt |>
      dplyr::mutate(island = island) |>
      dplyr::relocate("island", .before = "severity")
    
    sim_list[[island]] <- sim |>
      dplyr::mutate(island = island) |>
      dplyr::relocate("island", .before = "sim_year")
    
    kinfo_list[[island]] <- tibble::tibble(
      island = island,
      k_hat = kinfo$k_hat,
      annual_mean = kinfo$mu,
      annual_var = kinfo$var
    )
  }
  
  events_all        <- dplyr::bind_rows(Filter(Negate(is.null), events_list))
  annual_counts_all <- dplyr::bind_rows(Filter(Negate(is.null), annual_counts_list))
  lambda_all        <- dplyr::bind_rows(Filter(Negate(is.null), lambda_list))
  sim_all           <- dplyr::bind_rows(Filter(Negate(is.null), sim_list))
  kinfo_all         <- dplyr::bind_rows(Filter(Negate(is.null), kinfo_list))
  
  list(
    trackpoints = trackpoints_list,
    events_by_island = events_list,
    events_all = events_all,
    annual_counts_all = annual_counts_all,
    lambda_all = lambda_all,
    sim_all = sim_all,
    kinfo_all = kinfo_all,
    cfg = cfg
  )
}

