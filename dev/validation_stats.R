# dev/validation_stats.R
# Focused validation statistics for hindcast diagnosis:
# - annual maxima series
# - train/test p0 and quantiles
# - top years and drivers (SID, climo radii flag)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# ---- helpers ----

# Blocked CV split: contiguous blocks of fixed size (years).
# Returns a list of folds, each fold is list(train=<int>, test=<int>, fold_id=<int>).
.split_years_blocked_cv <- function(years, block_size = 5L) {
  years <- sort(unique(years[is.finite(years)]))
  if (!length(years)) return(list())

  # make contiguous blocks
  blocks <- split(years, ceiling(seq_along(years) / block_size))

  folds <- vector("list", length(blocks))
  for (i in seq_along(blocks)) {
    test_years <- as.integer(blocks[[i]])
    train_years <- as.integer(setdiff(years, test_years))
    folds[[i]] <- list(train = train_years, test = test_years, fold_id = i)
  }
  folds
}

.safe_quantile <- function(x, probs) {
  x <- x[is.finite(x)]
  if (!length(x)) return(rep(NA_real_, length(probs)))
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7))
}

.safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

# Default holdout split: last `holdout_years` years are test.
.split_years_last_block <- function(years, holdout_years) {
  years <- sort(unique(years[is.finite(years)]))
  if (length(years) == 0) return(list(train = integer(), test = integer()))
  if (length(years) <= holdout_years) {
    return(list(train = integer(), test = years))
  }
  n_train <- length(years) - holdout_years
  list(train = years[seq_len(n_train)], test = years[(n_train + 1L):length(years)])
}

# ---- main computations ----

# Build annual-max series from out$trackpoints for one site.
# Requires compute_site_winds_full() to be available in the session (load_all()).
.compute_site_annual_max <- function(dat_loc, lat, lon, site_name,
                                     min_year = 1970,
                                     storm_vmax_min = 34) {
  tmp <- compute_site_winds_full(df = dat_loc, target_lat = lat, target_lon = lon)

  # Require these columns; fail fast with informative error
  needed <- c("SID", "iso_time", "V_site_kt", "Vmax_kt", "dist_km")
  miss <- setdiff(needed, names(tmp))
  if (length(miss) > 0) {
    stop("compute_site_winds_full() output missing columns: ", paste(miss, collapse = ", "))
  }

  tmp <- tmp %>%
    mutate(
      year = as.integer(format(.data$iso_time, "%Y")),
      # ensure these flags exist even if compute_site_winds_full changes
      R34_is_climo = if ("R34_is_climo" %in% names(tmp)) .data$R34_is_climo else NA,
      R34_eff_km   = if ("R34_eff_km" %in% names(tmp)) .data$R34_eff_km else NA_real_,
      RMW_used_km  = if ("RMW_used_km" %in% names(tmp)) .data$RMW_used_km else NA_real_
    )

  # Restrict to modern era for hindcast-style diagnostics (IBTrACS reliability).
  if (is.finite(min_year)) {
    tmp <- tmp %>% filter(is.finite(.data$year), .data$year >= min_year)
  }

  # Build annual maxima over TS+ storms only (storm intensity threshold),
  # so annual-max series and p0 are consistent with TS34plus definition.
  if (is.finite(storm_vmax_min)) {
    tmp <- tmp %>% filter(is.finite(.data$Vmax_kt), .data$Vmax_kt >= storm_vmax_min)
  }

  # Annual maxima at the site (one value per year), plus the “driver row” metadata
  ann <- tmp %>%
    filter(is.finite(.data$year), is.finite(.data$V_site_kt)) %>%
    group_by(.data$year) %>%
    slice_max(order_by = .data$V_site_kt, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      site = site_name,
      year = .data$year,
      SID = .data$SID,
      iso_time = .data$iso_time,
      ann_max_kt = .data$V_site_kt,
      Vmax_kt = .data$Vmax_kt,
      dist_km = .data$dist_km,
      RMW_used_km = .data$RMW_used_km,
      R34_eff_km = .data$R34_eff_km,
      R34_is_climo = .data$R34_is_climo
    )

  ann
}

# Summarise annual-max series for a given train/test year split.
.summarise_split <- function(ann, years, threshold_kt = 34) {
  if (!length(years)) {
    return(tibble(
      n_years = 0L,
      n_pos_years = 0L,
      p0 = NA_real_,
      q50 = NA_real_, q75 = NA_real_, q90 = NA_real_, q95 = NA_real_
    ))
  }
  x <- ann %>%
    filter(.data$year %in% years) %>%
    pull(.data$ann_max_kt)

  n_years <- length(years)
  n_pos <- sum(is.finite(x) & x >= threshold_kt, na.rm = TRUE)
  p0 <- 1 - (n_pos / n_years)

  qs <- .safe_quantile(x, c(0.50, 0.75, 0.90, 0.95))
  tibble(
    n_years = as.integer(n_years),
    n_pos_years = as.integer(n_pos),
    p0 = as.numeric(p0),
    q50 = qs[1], q75 = qs[2], q90 = qs[3], q95 = qs[4]
  )
}

# Top K annual maxima in the test years (with drivers)
.top_k_test <- function(ann, test_years, k = 5) {
  ann %>%
    filter(.data$year %in% test_years) %>%
    arrange(desc(.data$ann_max_kt)) %>%
    slice_head(n = k)
}

# ---- entrypoint ----

run_validation_stats <- function(out, targets, sites = NULL,
                                 holdout_years = 5,
                                 threshold_kt = 34,
                                 top_k = 5,
                                 min_year = 1970,
                                 storm_vmax_min = 34,
                                 split_mode = c("last_block", "blocked_cv"),
                                 block_size = 5) {

  split_mode <- match.arg(split_mode)

  if (is.null(out$trackpoints) || !is.list(out$trackpoints)) {
    stop("`out$trackpoints` must be a named list of per-site trackpoint tables.")
  }
  if (is.null(targets) || !all(c("name", "lat", "lon") %in% names(targets))) {
    stop("`targets` must contain columns: name, lat, lon.")
  }

  if (is.null(sites)) sites <- as.character(targets$name)
  sites <- as.character(sites)

  all_ann <- list()
  out_tbl <- list()
  top_tbl <- list()

  for (site in sites) {
    if (is.null(out$trackpoints[[site]]) || nrow(out$trackpoints[[site]]) == 0) {
      message("[skip] ", site, ": no trackpoints")
      next
    }

    lat <- targets$lat[targets$name == site][1]
    lon <- targets$lon[targets$name == site][1]

    ann <- .compute_site_annual_max(out$trackpoints[[site]], lat, lon, site,
                                    min_year = min_year, storm_vmax_min = storm_vmax_min)


    all_years <- sort(unique(ann$year))

    folds <- if (split_mode == "last_block") {
      list(c(.split_years_last_block(all_years, holdout_years = holdout_years), fold_id = 1L))
    } else {
      .split_years_blocked_cv(all_years, block_size = as.integer(block_size))
    }

    for (fold in folds) {
      train_sum <- .summarise_split(ann, fold$train, threshold_kt)
      test_sum  <- .summarise_split(ann, fold$test, threshold_kt)

      summ <- tibble(
        site = site,
        split_mode = split_mode,
        fold_id = as.integer(fold$fold_id),
        holdout_years = if (split_mode == "last_block") holdout_years else NA_integer_,
        block_size = if (split_mode == "blocked_cv") as.integer(block_size) else NA_integer_,
        min_year = min_year,
        storm_vmax_min = storm_vmax_min,
        threshold_kt = threshold_kt,
        train_years = paste(fold$train, collapse = ","),
        test_years  = paste(fold$test, collapse = ",")
      ) %>%
        bind_cols(
          train_sum %>% rename_with(~paste0("train_", .x)),
          test_sum  %>% rename_with(~paste0("test_", .x))
        )

      topk <- .top_k_test(ann, fold$test, k = top_k) %>%
        mutate(
          r_norm_used = ifelse(is.finite(.data$RMW_used_km) & .data$RMW_used_km > 0,
                               .data$dist_km / .data$RMW_used_km, NA_real_),
          v_ratio = ifelse(is.finite(.data$Vmax_kt) & .data$Vmax_kt > 0,
                           .data$ann_max_kt / .data$Vmax_kt, NA_real_),
          split_mode = split_mode,
          fold_id = as.integer(fold$fold_id)
        )

      out_tbl[[length(out_tbl) + 1L]] <- summ
      top_tbl[[length(top_tbl) + 1L]] <- topk
    }

    if (is.null(all_ann[[site]])) all_ann[[site]] <- ann
  }

  summary_tbl <- bind_rows(out_tbl)

  if (split_mode == "blocked_cv") {
    message("\n(Note) blocked_cv: one row per site × fold. Use dplyr to summarise across folds.")
  }

  top_tbl_all <- bind_rows(top_tbl)
  ann_all <- bind_rows(all_ann)

  # print concise outputs
  message("\n==============================")
  message("VALIDATION STATS (annual maxima)")
  message("==============================")
  print(summary_tbl, n = nrow(summary_tbl))

  message("\n------------------------------")
  message("Top ", top_k, " test-year annual maxima (drivers)")
  message("------------------------------")
  print(top_tbl_all, n = nrow(top_tbl_all))

  invisible(list(summary = summary_tbl, top_test = top_tbl_all, annual_max = ann_all))
}

# ---- example usage (uncomment and edit to run) ----
# pkgload::load_all(".")
# cfg <- make_hazard_cfg(data_path = "inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv")
# targets <- tibble::tribble(
#   ~name, ~lat, ~lon,
#   "St_Martin", 18.0708, -63.0501,
#   "Saba", 17.6350, -63.2300,
#   "Statia", 17.4890, -62.9740,
#   "Puerto_Rico", 18.2208, -66.5901,
#   "Miami", 25.7617, -80.1918
# )
# out <- run_hazard_model(cfg = cfg, targets = targets, verbose = FALSE)
# res <- run_validation_stats(out, targets, sites = c("Miami","Saba","Statia"), holdout_years = 5, threshold_kt = 34, top_k = 5)
