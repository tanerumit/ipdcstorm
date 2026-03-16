  # =============================================================================
  # Script overview: IBTrACS ingestion
  # - read_ibtracs_clean(): read/clean IBTrACS CSV and standardize USA best-track fields.
  # =============================================================================

  #' Read and clean IBTrACS CSV (USA fields + pressure/structure when present)
  #'
  #' @description
  #' Reads an IBTrACS CSV and returns a cleaned track-point tibble using USA best-track
  #' variables. Adds central pressure `pres_hpa` (USA_PRES with WMO_PRES fallback)
  #' and optional structure fields (USA_POCI/USA_ROCI/USA_RMW) when present.
  #'
  #' @param ibtracs_csv Path to an IBTrACS CSV file.
  #' @param basin Character vector or NULL; filter BASIN if provided (e.g. "NA").
  #' @param season Integer vector or NULL; filter SEASON if provided.
  #' @param keep_all Logical; keep original columns if TRUE.
  #' @param verbose Logical; print read diagnostics if TRUE.
  #'
  #' @return Tibble of cleaned track points.
  #'
  #' @export
  read_ibtracs_clean <- function(ibtracs_csv,
                                 basin = "NA",
                                 season = NULL,
                                 keep_all = FALSE,
                                 verbose = TRUE) {
    
    if (!is.character(ibtracs_csv) || length(ibtracs_csv) != 1L || !nzchar(ibtracs_csv)) {
      stop("`ibtracs_csv` must be a non-empty character scalar path.")
    }
    if (!file.exists(ibtracs_csv)) stop("File not found: ", ibtracs_csv)
    
    if (!requireNamespace("readr", quietly = TRUE)) stop("Package `readr` is required.")
    if (!requireNamespace("lubridate", quietly = TRUE)) stop("Package `lubridate` is required.")
    
    col_num_or_na <- function(dat, col) {
      if (col %in% names(dat)) .to_num_quiet(dat[[col]]) else rep(NA_real_, nrow(dat))
    }
    col_chr_or_na <- function(dat, col) {
      if (col %in% names(dat)) as.character(dat[[col]]) else rep(NA_character_, nrow(dat))
    }
    
    if (isTRUE(verbose)) message("[IBTrACS] Reading CSV: ", ibtracs_csv)
    
    df <- readr::read_csv(
      file = ibtracs_csv,
      na = c("", " ", ".", "-"),   # do NOT include "NA" here
      show_col_types = FALSE,
      progress = FALSE,
      guess_max = 200000
    )
    
    if ("BASIN" %in% names(df))    df[["BASIN"]]    <- toupper(trimws(as.character(df[["BASIN"]])))
    if ("SUBBASIN" %in% names(df)) df[["SUBBASIN"]] <- toupper(trimws(as.character(df[["SUBBASIN"]])))
    
    if (!is.null(basin)) basin <- toupper(trimws(as.character(basin)))
    
    if (isTRUE(verbose) && "BASIN" %in% names(df)) {
      message("[IBTrACS] Unique BASIN (first 20): ",
              paste(utils::head(sort(unique(df$BASIN)), 20), collapse = ", "))
    }

    cleaned_required_cols <- c(
      "SID", "SEASON", "BASIN", "iso_time",
      "lat", "lon", "wind_kt",
      "r34_ne_nm", "r34_se_nm", "r34_sw_nm", "r34_nw_nm"
    )
    if (all(cleaned_required_cols %in% names(df))) {
      if (!is.null(basin))  df <- dplyr::filter(df, .data$BASIN %in% basin)
      if (!is.null(season)) df <- dplyr::filter(df, .data$SEASON %in% season)

      if (nrow(df) == 0) {
        present <- sort(unique(toupper(trimws(as.character(df$BASIN)))))
        present <- present[!is.na(present)]
        stop("After filtering, 0 rows remain. Requested basin=",
             paste(basin, collapse = ", "),
             ". File BASIN values include: ",
             paste(utils::head(present, 20), collapse = ", "))
      }

      df2 <- df |>
        dplyr::mutate(
          iso_time = lubridate::ymd_hms(.data$iso_time, tz = "UTC", quiet = TRUE),
          lat = .to_num_quiet(.data$lat),
          lon = .to_num_quiet(.data$lon),
          wind_kt = .to_num_quiet(.data$wind_kt),
          r34_ne_nm = .to_num_quiet(.data$r34_ne_nm),
          r34_se_nm = .to_num_quiet(.data$r34_se_nm),
          r34_sw_nm = .to_num_quiet(.data$r34_sw_nm),
          r34_nw_nm = .to_num_quiet(.data$r34_nw_nm),
          pres_hpa = if ("pres_hpa" %in% names(df)) .to_num_quiet(.data$pres_hpa) else NA_real_,
          pres_source = if ("pres_source" %in% names(df)) as.character(.data$pres_source) else NA_character_,
          usa_pres_hpa = if ("usa_pres_hpa" %in% names(df)) .to_num_quiet(.data$usa_pres_hpa) else NA_real_,
          wmo_pres_hpa = if ("wmo_pres_hpa" %in% names(df)) .to_num_quiet(.data$wmo_pres_hpa) else NA_real_,
          poci_hpa = if ("poci_hpa" %in% names(df)) .to_num_quiet(.data$poci_hpa) else NA_real_,
          roci_km = if ("roci_km" %in% names(df)) .to_num_quiet(.data$roci_km) else NA_real_,
          rmw_km = if ("rmw_km" %in% names(df)) .to_num_quiet(.data$rmw_km) else NA_real_,
          storm_name = if ("storm_name" %in% names(df)) as.character(.data$storm_name) else NA_character_,
          storm_status = if ("storm_status" %in% names(df)) as.character(.data$storm_status) else NA_character_,
          r50_ne_nm = if ("r50_ne_nm" %in% names(df)) .to_num_quiet(.data$r50_ne_nm) else NA_real_,
          r50_se_nm = if ("r50_se_nm" %in% names(df)) .to_num_quiet(.data$r50_se_nm) else NA_real_,
          r50_sw_nm = if ("r50_sw_nm" %in% names(df)) .to_num_quiet(.data$r50_sw_nm) else NA_real_,
          r50_nw_nm = if ("r50_nw_nm" %in% names(df)) .to_num_quiet(.data$r50_nw_nm) else NA_real_,
          r64_ne_nm = if ("r64_ne_nm" %in% names(df)) .to_num_quiet(.data$r64_ne_nm) else NA_real_,
          r64_se_nm = if ("r64_se_nm" %in% names(df)) .to_num_quiet(.data$r64_se_nm) else NA_real_,
          r64_sw_nm = if ("r64_sw_nm" %in% names(df)) .to_num_quiet(.data$r64_sw_nm) else NA_real_,
          r64_nw_nm = if ("r64_nw_nm" %in% names(df)) .to_num_quiet(.data$r64_nw_nm) else NA_real_,
          storm_speed_kt = if ("storm_speed_kt" %in% names(df)) .to_num_quiet(.data$storm_speed_kt) else NA_real_,
          storm_dir_deg = if ("storm_dir_deg" %in% names(df)) .to_num_quiet(.data$storm_dir_deg) else NA_real_
        ) |>
        dplyr::arrange(.data$SID, .data$iso_time)

      out <- if (isTRUE(keep_all)) {
        df2
      } else {
        df2 |>
          dplyr::select(
            SID, SEASON, BASIN,
            iso_time, storm_name, storm_status,
            lat, lon, wind_kt,
            pres_hpa, pres_source, usa_pres_hpa, wmo_pres_hpa,
            poci_hpa, roci_km, rmw_km,
            r34_ne_nm, r34_se_nm, r34_sw_nm, r34_nw_nm,
            r50_ne_nm, r50_se_nm, r50_sw_nm, r50_nw_nm,
            r64_ne_nm, r64_se_nm, r64_sw_nm, r64_nw_nm,
            storm_speed_kt, storm_dir_deg
          )
      }

      out <- tibble::as_tibble(out)

      if (isTRUE(verbose)) {
        message("[IBTrACS] Rows: ", nrow(out), " | Pc available: ", sum(is.finite(out$pres_hpa)))
        message("[IBTrACS] Non-missing coords: lat=", sum(is.finite(out$lat)),
                " lon=", sum(is.finite(out$lon)))
      }

      return(out)
    }
    
    required_cols <- c("SID", "SEASON", "BASIN", "ISO_TIME",
                       "USA_LAT", "USA_LON", "USA_WIND",
                       "USA_R34_NE", "USA_R34_SE", "USA_R34_SW", "USA_R34_NW")
    missing_cols <- setdiff(required_cols, names(df))
    if (length(missing_cols) > 0) {
      stop("IBTrACS file is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    
    if (!is.null(basin))  df <- dplyr::filter(df, .data$BASIN %in% basin)
    if (!is.null(season)) df <- dplyr::filter(df, .data$SEASON %in% season)
    
    if (nrow(df) == 0) {
      present <- sort(unique(toupper(trimws(readr::read_csv(
        ibtracs_csv, show_col_types = FALSE, progress = FALSE
      )[["BASIN"]]))))
      present <- present[!is.na(present)]
      stop("After filtering, 0 rows remain. Requested basin=",
           paste(basin, collapse = ", "),
           ". File BASIN values include: ",
           paste(utils::head(present, 20), collapse = ", "))
    }
    
    # Coords: prefer USA coords, fallback to LAT/LON
    lat_raw <- col_num_or_na(df, "USA_LAT")
    lon_raw <- col_num_or_na(df, "USA_LON")
    lat_wmo <- col_num_or_na(df, "LAT")
    lon_wmo <- col_num_or_na(df, "LON")
    
    lat_clean <- ifelse(is.finite(lat_raw), lat_raw, lat_wmo)
    lon_clean <- ifelse(is.finite(lon_raw), lon_raw, lon_wmo)
    
    # Optional columns
    usa_pres <- col_num_or_na(df, "USA_PRES")
    wmo_pres <- col_num_or_na(df, "WMO_PRES")
    poci     <- col_num_or_na(df, "USA_POCI")
    roci     <- col_num_or_na(df, "USA_ROCI")
    rmw_nm   <- col_num_or_na(df, "USA_RMW")
    rmw      <- rmw_nm * 1.852
    rmw[!is.finite(rmw) | rmw <= 5 | rmw >= 150] <- NA_real_
    
    usa_status <- col_chr_or_na(df, "USA_STATUS")
    name_chr   <- col_chr_or_na(df, "NAME")
    
    storm_speed_kt <- col_num_or_na(df, "STORM_SPEED")
    storm_dir_deg  <- col_num_or_na(df, "STORM_DIR")
    
    df2 <- df |>
      dplyr::mutate(
        iso_time = lubridate::ymd_hms(.data$ISO_TIME, tz = "UTC", quiet = TRUE),
        
        lat = lat_clean,
        lon = lon_clean,
        wind_kt = .to_num_quiet(.data$USA_WIND),
        
        r34_ne_nm = .to_radii_nm(df[["USA_R34_NE"]], cap_nm = 600),
        r34_se_nm = .to_radii_nm(df[["USA_R34_SE"]], cap_nm = 600),
        r34_sw_nm = .to_radii_nm(df[["USA_R34_SW"]], cap_nm = 600),
        r34_nw_nm = .to_radii_nm(df[["USA_R34_NW"]], cap_nm = 600),
        
        r50_ne_nm = .to_radii_nm(col_chr_or_na(df, "USA_R50_NE"), cap_nm = 400),
        r50_se_nm = .to_radii_nm(col_chr_or_na(df, "USA_R50_SE"), cap_nm = 400),
        r50_sw_nm = .to_radii_nm(col_chr_or_na(df, "USA_R50_SW"), cap_nm = 400),
        r50_nw_nm = .to_radii_nm(col_chr_or_na(df, "USA_R50_NW"), cap_nm = 400),
        
        r64_ne_nm = .to_radii_nm(col_chr_or_na(df, "USA_R64_NE"), cap_nm = 250),
        r64_se_nm = .to_radii_nm(col_chr_or_na(df, "USA_R64_SE"), cap_nm = 250),
        r64_sw_nm = .to_radii_nm(col_chr_or_na(df, "USA_R64_SW"), cap_nm = 250),
        r64_nw_nm = .to_radii_nm(col_chr_or_na(df, "USA_R64_NW"), cap_nm = 250),
        
        usa_pres_hpa = usa_pres,
        wmo_pres_hpa = wmo_pres,
        
        pres_hpa = dplyr::case_when(
          is.finite(.data$usa_pres_hpa) ~ .data$usa_pres_hpa,
          is.finite(.data$wmo_pres_hpa) ~ .data$wmo_pres_hpa,
          TRUE ~ NA_real_
        ),
        pres_source = dplyr::case_when(
          is.finite(.data$usa_pres_hpa) ~ "USA",
          is.finite(.data$wmo_pres_hpa) ~ "WMO",
          TRUE ~ NA_character_
        ),
        
        poci_hpa = poci,
        roci_km  = roci,
        rmw_km   = rmw,
        
        storm_status = usa_status,
        storm_name   = name_chr,
        
        storm_speed_kt = storm_speed_kt,
        storm_dir_deg  = storm_dir_deg
      ) |>
      dplyr::arrange(.data$SID, .data$iso_time)
    
    out <- if (isTRUE(keep_all)) {
      df2
    } else {
      df2 |>
        dplyr::select(
          SID, SEASON, BASIN,
          iso_time, storm_name, storm_status,
          lat, lon, wind_kt,
          pres_hpa, pres_source, usa_pres_hpa, wmo_pres_hpa,
          poci_hpa, roci_km, rmw_km,
          r34_ne_nm, r34_se_nm, r34_sw_nm, r34_nw_nm,
          r50_ne_nm, r50_se_nm, r50_sw_nm, r50_nw_nm,
          r64_ne_nm, r64_se_nm, r64_sw_nm, r64_nw_nm,
          storm_speed_kt, storm_dir_deg
        )
    }
    
    out <- tibble::as_tibble(out)
    
    if (isTRUE(verbose)) {
      message("[IBTrACS] Rows: ", nrow(out), " | Pc available: ", sum(is.finite(out$pres_hpa)))
      message("[IBTrACS] Non-missing coords: lat=", sum(is.finite(out$lat)),
              " lon=", sum(is.finite(out$lon)))
    }
    
    out
  }
  
