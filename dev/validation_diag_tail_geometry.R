# Focused validation diagnostics for tail-geometry checks in the TS+ regime.
.run_validation_diagnostics_tail_geometry <- function(out_trackpoints,
                                                      targets,
                                                      sites = c("Saba", "Statia", "St_Martin"),
                                                      tail_prob = 0.95,
                                                      top_k = 10,
                                                      tail_mode = c("quantile", "top_k", "both"),
                                                      storm_filter = c("site_wind", "storm_vmax")) {
  fmt_num <- function(x, digits = 2) {
    if (length(x) == 0 || !is.finite(x[1])) return("NA")
    format(round(x[1], digits), trim = TRUE, nsmall = digits)
  }

  safe_quantile <- function(x, prob) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE, type = 7))
  }

  safe_median <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    stats::median(x, na.rm = TRUE)
  }

  if (is.null(out_trackpoints) || !is.list(out_trackpoints)) {
    stop("out_trackpoints must be a named list like out$trackpoints.", call. = FALSE)
  }
  if (is.null(targets) || !all(c("name", "lat", "lon") %in% names(targets))) {
    stop("targets must contain columns `name`, `lat`, and `lon`.", call. = FALSE)
  }
  tail_mode <- match.arg(tail_mode)
  storm_filter <- match.arg(storm_filter)
  tail_prob <- as.numeric(tail_prob)
  top_k <- as.integer(top_k)
  if (!is.finite(tail_prob) || tail_prob <= 0 || tail_prob >= 1) {
    stop("tail_prob must be a finite number strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.finite(top_k) || top_k < 1L) {
    stop("top_k must be an integer >= 1.", call. = FALSE)
  }

  target_names <- as.character(targets$name)
  req_sites <- as.character(sites)
  n_sites <- length(req_sites)

  rows <- vector("list", n_sites)

  message("\n", "-" |> rep(72) |> paste(collapse = ""))
  message("  TAIL GEOMETRY DIAGNOSTICS (TS+ tail)")
  message("-" |> rep(72) |> paste(collapse = ""))

  for (i in seq_len(n_sites)) {
    site_name <- req_sites[i]
    site_idx <- match(site_name, target_names)

    row_i <- tibble::tibble(
      site = site_name,
      n_pos = 0L,
      top_n = 0L,
      R34_p50 = NA_real_,
      R34_p90 = NA_real_,
      R34_unique = 0L,
      dist_p50_top = NA_real_,
      dist_over_R34_p50_top = NA_real_,
      dist_over_R34mean_p50_top = NA_real_,
      frac_rnorm_lt1 = NA_real_,
      rnorm_p10 = NA_real_,
      rnorm_p50 = NA_real_,
      rnorm_p90 = NA_real_
    )

    if (!is.finite(site_idx)) {
      message(sprintf("  %s: skipped (site not found in `targets`).", site_name))
      rows[[i]] <- row_i
      next
    }
    if (is.null(out_trackpoints[[site_name]]) || nrow(out_trackpoints[[site_name]]) == 0) {
      message(sprintf("  %s: skipped (no per-site trackpoints in `out_trackpoints`).", site_name))
      rows[[i]] <- row_i
      next
    }

    dat_loc <- out_trackpoints[[site_name]]
    lat <- targets$lat[site_idx]
    lon <- targets$lon[site_idx]

    tmp <- tryCatch(
      compute_site_winds_full(df = dat_loc, target_lat = lat, target_lon = lon),
      error = function(e) e
    )
    if (inherits(tmp, "error")) {
      message(sprintf("  %s: skipped (compute_site_winds_full failed: %s).", site_name, tmp$message))
      rows[[i]] <- row_i
      next
    }

    keep <- if (storm_filter == "storm_vmax") {
      is.finite(tmp$Vmax_kt) &
        (tmp$Vmax_kt >= 34) &
        is.finite(tmp$R34_km) &
        (tmp$R34_km > 0) &
        is.finite(tmp$dist_km) &
        (tmp$dist_km >= 0)
    } else {
      is.finite(tmp$V_site_kt) &
        (tmp$V_site_kt >= 34) &
        is.finite(tmp$R34_km) &
        (tmp$R34_km > 0) &
        is.finite(tmp$dist_km) &
        (tmp$dist_km >= 0)
    }

    tmp_pos <- tmp[keep, , drop = FALSE]
    row_i$n_pos <- nrow(tmp_pos)

    if (nrow(tmp_pos) == 0) {
      message(sprintf("  %s: skipped (no TS+ rows with finite positive R34_km).", site_name))
      rows[[i]] <- row_i
      next
    }

    row_i$R34_p50 <- safe_median(tmp_pos$R34_km)
    row_i$R34_p90 <- safe_quantile(tmp_pos$R34_km, 0.90)

    r34_round <- round(tmp_pos$R34_km[is.finite(tmp_pos$R34_km)], 2)
    row_i$R34_unique <- length(unique(r34_round))

    if (nrow(tmp_pos) < 5) {
      message(sprintf("  %s: insufficient TS+ rows (n=%d); skipping tail stats.", site_name, nrow(tmp_pos)))
      rows[[i]] <- row_i
      next
    }

    tail_q <- rep(FALSE, nrow(tmp_pos))
    q_tail <- safe_quantile(tmp_pos$V_site_kt, tail_prob)
    if (is.finite(q_tail)) {
      tail_q <- is.finite(tmp_pos$V_site_kt) & (tmp_pos$V_site_kt >= q_tail)
    }

    tail_k <- rep(FALSE, nrow(tmp_pos))
    k_use <- min(top_k, nrow(tmp_pos))
    if (k_use > 0) {
      tmp_rank <- tmp_pos |>
        dplyr::mutate(.tail_row_id = seq_len(nrow(tmp_pos)))
      top_rows <- dplyr::slice_max(
        tmp_rank,
        order_by = .data$V_site_kt,
        n = k_use,
        with_ties = TRUE
      )
      tail_k[top_rows$.tail_row_id] <- TRUE
    }

    top <- switch(
      tail_mode,
      quantile = tail_q,
      top_k = tail_k,
      both = tail_q | tail_k
    )

    row_i$top_n <- sum(top, na.rm = TRUE)

    r_norm <- rep(NA_real_, nrow(tmp_pos))
    ok_rnorm <- is.finite(tmp_pos$dist_km) &
      is.finite(tmp_pos$RMW_km) &
      (tmp_pos$RMW_km > 0)
    r_norm[ok_rnorm] <- tmp_pos$dist_km[ok_rnorm] / tmp_pos$RMW_km[ok_rnorm]

    if (any(top, na.rm = TRUE)) {
      row_i$dist_p50_top <- safe_median(tmp_pos$dist_km[top])

      ratio_dir <- rep(NA_real_, sum(top))
      r34_top <- tmp_pos$R34_km[top]
      dist_top <- tmp_pos$dist_km[top]
      ok_dir <- is.finite(dist_top) & is.finite(r34_top) & (r34_top > 0)
      ratio_dir[ok_dir] <- dist_top[ok_dir] / r34_top[ok_dir]
      row_i$dist_over_R34_p50_top <- safe_median(ratio_dir)

      if ("R34_mean_km" %in% names(tmp_pos)) {
        ratio_mean <- rep(NA_real_, sum(top))
        r34_mean_top <- tmp_pos$R34_mean_km[top]
        ok_mean <- is.finite(dist_top) & is.finite(r34_mean_top) & (r34_mean_top > 0)
        ratio_mean[ok_mean] <- dist_top[ok_mean] / r34_mean_top[ok_mean]
        row_i$dist_over_R34mean_p50_top <- safe_median(ratio_mean)
      }

      r_norm_top <- r_norm[top]
      row_i$frac_rnorm_lt1 <- mean(r_norm_top < 1, na.rm = TRUE)
      row_i$rnorm_p10 <- safe_quantile(r_norm_top, 0.10)
      row_i$rnorm_p50 <- safe_quantile(r_norm_top, 0.50)
      row_i$rnorm_p90 <- safe_quantile(r_norm_top, 0.90)
    }

    message(sprintf(
      "  %s: n_pos=%d; top_n=%d; R34_p50=%s km; R34_p90=%s km; R34_unique=%d; dist_p50_top=%s km; dist/R34_p50_top=%s; dist/R34_mean_p50_top=%s; frac_rnorm_lt1=%s; rnorm_p10=%s; rnorm_p50=%s; rnorm_p90=%s",
      site_name,
      row_i$n_pos,
      row_i$top_n,
      fmt_num(row_i$R34_p50),
      fmt_num(row_i$R34_p90),
      row_i$R34_unique,
      fmt_num(row_i$dist_p50_top),
      fmt_num(row_i$dist_over_R34_p50_top, digits = 3),
      fmt_num(row_i$dist_over_R34mean_p50_top, digits = 3),
      fmt_num(row_i$frac_rnorm_lt1, digits = 3),
      fmt_num(row_i$rnorm_p10, digits = 3),
      fmt_num(row_i$rnorm_p50, digits = 3),
      fmt_num(row_i$rnorm_p90, digits = 3)
    ))

    if (!all(c("SID", "iso_time") %in% names(tmp_pos))) {
      message(sprintf("    %s event-max: skipped (missing SID or iso_time).", site_name))
    } else {
      evt <- tmp_pos |>
        dplyr::group_by(.data$SID) |>
        dplyr::slice_max(order_by = .data$V_site_kt, n = 1, with_ties = FALSE) |>
        dplyr::ungroup()
      evt_top <- dplyr::slice_max(evt, order_by = .data$V_site_kt, n = 10, with_ties = FALSE)

      cols_evt <- c("SID", "iso_time", "V_site_kt", "dist_km", "RMW_km", "R34_km", "Vmax_kt", "quadrant")
      for (nm in cols_evt) {
        if (!(nm %in% names(evt_top))) evt_top[[nm]] <- NA
      }
      evt_top <- evt_top |>
        dplyr::select(dplyr::all_of(cols_evt))

      message(sprintf("    %s event-max top 10:", site_name))
      print(evt_top, n = nrow(evt_top))
    }

    rows[[i]] <- row_i
  }

  out <- dplyr::bind_rows(rows)
  print(out, n = nrow(out))
  invisible(out)
}
