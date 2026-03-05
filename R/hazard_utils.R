  ################################################################################
  # Script overview: shared utilities
  # - .to_num_quiet(): robust numeric parsing for IBTrACS strings.
  # - .max_na(): NA-safe max for numeric vectors.
  # - .min_na(): NA-safe min for numeric vectors.
  # - .mean_na(): NA-safe mean for numeric vectors.
  # - .to_radii_nm(): parse wind radii fields with caps and sentinels.
  #
  # Conventions used:
  # - Internal helpers: .name() + @keywords internal + not exported
  # - User-facing functions: exported with @export
  # - Package deps referenced with :: and requireNamespace() where needed
  ################################################################################
  
  # =============================================================================
  # 0) Internal utilities (type checks, NA-safe summaries)
  # =============================================================================
  
  #' Quiet numeric parsing for messy IBTrACS string fields
  #'
  #' @description
  #' Converts character inputs to numeric robustly by stripping common text
  #' decorations (e.g., "deg", "degrees_north") and treating known placeholders
  #' as missing.
  #'
  #' @param x A vector (character/numeric) to parse as numeric.
  #'
  #' @return Numeric vector with non-parsable values set to NA.
  #'
  #' @keywords internal
  .to_num_quiet <- function(x) {
    x <- stringr::str_trim(as.character(x))
    x[x %in% c("", "NA", "N/A", "NULL", "null", ".", "-")] <- NA_character_
    x <- stringr::str_replace_all(
      x,
      stringr::regex("\\b(degrees?_north|degrees?_south|degrees?_east|degrees?_west)\\b", ignore_case = TRUE),
      ""
    )
    x <- stringr::str_replace_all(x, stringr::regex("\\b(deg|degrees?)\\b", ignore_case = TRUE), "")
    x <- stringr::str_trim(x)
    
    out <- suppressWarnings(readr::parse_number(x, na = c("", "NA", "N/A", "NULL", "null")))
    out[out %in% c(-999, -99, 999)] <- NA_real_
    out
  }

  #' NA-safe max for numeric vectors
  #' @param x Numeric vector.
  #' @return Scalar numeric; NA if no finite values.
  #' @keywords internal
  .max_na <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else max(x)
  }
  
  #' NA-safe min for numeric vectors
  #' @param x Numeric vector.
  #' @return Scalar numeric; NA if no finite values.
  #' @keywords internal
  .min_na <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else min(x)
  }
  
  #' NA-safe mean for numeric vectors
  #' @param x Numeric vector.
  #' @return Scalar numeric; NA if no finite values.
  #' @keywords internal
  .mean_na <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else mean(x)
  }

  #' Parse IBTrACS wind radii fields (nautical miles) with caps and sentinels
  #'
  #' @param x A vector (character/numeric) of radii values.
  #' @param cap_nm Numeric scalar; values greater than this are treated as NA.
  #'
  #' @return Numeric vector of radii (nm) with invalid/sentinel values set to NA.
  #'
  #' @keywords internal
.to_radii_nm <- function(x, cap_nm) {
    x <- stringr::str_trim(as.character(x))
    x[x %in% c("", "NA", "N/A", "NULL", "null", ".", "-")] <- NA_character_
    
    out <- suppressWarnings(readr::parse_number(x, na = c("", "NA", "N/A", "NULL", "null")))
    out[out %in% c(-999, -99, 0, 99, 999, 9999, 99999)] <- NA_real_
    out[is.finite(out) & out > cap_nm] <- NA_real_
  out
}

#' Normalize lambda scaling mode
#'
#' @param mode Character scalar mode.
#'
#' @return Character scalar in {"target", "down_only"}.
#'
#' @keywords internal
.normalize_lambda_scaling_mode <- function(mode) {
  mode_chr <- if (is.null(mode) || length(mode) == 0L || !is.character(mode)) {
    "target"
  } else {
    as.character(mode[1])
  }
  if (!isTRUE(nzchar(mode_chr))) {
    mode_chr <- "target"
  }
  if (!(mode_chr %in% c("target", "down_only"))) {
    stop("lambda scaling mode must be one of: target, down_only.", call. = FALSE)
  }
  mode_chr
}

#' Build multiplicative lambda scalers from a rate-check table
#'
#' @param rate_tbl Tibble with location, storm_class, lambda_model or
#'   lambda_model_raw, lambda_ref, and optional expected_ratio.
#' @param scale_min,scale_max Numeric scalar clamp bounds for lambda_scale.
#' @param scaling_mode Character scalar. `"target"` preserves historical behavior;
#'   `"down_only"` prevents lambda upscaling above modeled rates.
#'
#' @return Tibble keyed by location and storm_class with raw lambda, target
#'   lambda, lambda_scale, adjusted lambda, and scaler status.
#'
#' @keywords internal
.lambda_scalers_from_rate_check <- function(rate_tbl,
                                            scale_min = 0.25,
                                            scale_max = 4,
                                            scaling_mode = "target") {
  rt <- tibble::as_tibble(rate_tbl)
  required <- c("location", "storm_class", "lambda_ref")
    missing_cols <- setdiff(required, names(rt))
    if (length(missing_cols) > 0) {
      stop(
        "rate_tbl is missing required columns: ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }

    lambda_col <- if ("lambda_model_raw" %in% names(rt)) {
      "lambda_model_raw"
    } else if ("lambda_model" %in% names(rt)) {
      "lambda_model"
    } else {
      stop(
        "rate_tbl must contain either lambda_model_raw or lambda_model.",
        call. = FALSE
      )
    }

    if (!("expected_ratio" %in% names(rt))) {
      rt$expected_ratio <- 1
    }

  scale_min <- as.numeric(scale_min)
  scale_max <- as.numeric(scale_max)
  if (!is.finite(scale_min) || !is.finite(scale_max) || scale_min <= 0 || scale_max < scale_min) {
    stop("scale_min and scale_max must be finite with 0 < scale_min <= scale_max.", call. = FALSE)
  }
  scaling_mode <- .normalize_lambda_scaling_mode(scaling_mode)

  if (any(is.finite(as.numeric(rt[[lambda_col]])) & as.numeric(rt[[lambda_col]]) < 0, na.rm = TRUE)) {
    stop("lambda_model_raw must be >= 0 for all finite values.", call. = FALSE)
  }

  out <- rt |>
    dplyr::mutate(
      lambda_model_raw = as.numeric(.data[[lambda_col]]),
      expected_ratio = dplyr::if_else(
          is.finite(.data$expected_ratio) & .data$expected_ratio > 0,
          as.numeric(.data$expected_ratio),
          1
        ),
        lambda_target = dplyr::if_else(
          is.finite(.data$lambda_ref),
          as.numeric(.data$lambda_ref) * .data$expected_ratio,
          .data$lambda_model_raw
        ),
        scale_raw_base = dplyr::case_when(
          !is.finite(.data$lambda_ref) ~ 1,
          !is.finite(.data$lambda_model_raw) ~ 1,
          .data$lambda_model_raw <= 0 & .data$lambda_target > 0 ~ scale_max,
          .data$lambda_model_raw <= 0 ~ 1,
          TRUE ~ .data$lambda_target / .data$lambda_model_raw
        ),
        scale_raw = if (identical(scaling_mode, "down_only")) {
          pmin(1, .data$scale_raw_base)
        } else {
          .data$scale_raw_base
        },
        scale_raw = dplyr::if_else(
          is.finite(.data$scale_raw) & .data$scale_raw > 0,
          .data$scale_raw,
          1
        ),
        lambda_scale = pmax(scale_min, pmin(scale_max, .data$scale_raw)),
        scale_status = dplyr::case_when(
          !is.finite(.data$lambda_ref) ~ "no_ref",
          !is.finite(.data$lambda_model_raw) ~ "no_model",
          .data$lambda_model_raw <= 0 & .data$lambda_target > 0 ~ "zero_model",
          .data$scale_raw < scale_min ~ "clamped_low",
          .data$scale_raw > scale_max ~ "clamped_high",
          TRUE ~ "ok"
        ),
        scale_clamped = .data$scale_status %in% c("clamped_low", "clamped_high", "zero_model"),
        lambda_adj = .data$lambda_model_raw * .data$lambda_scale,
        lambda_scaling_mode = scaling_mode,
        lambda_scale_applied = .data$lambda_scale,
        was_upscaled = .data$lambda_scale_applied > (1 + 1e-12)
      ) |>
      dplyr::select(
        "location", "storm_class",
        "lambda_model_raw", "lambda_ref", "expected_ratio",
        "lambda_target", "lambda_scale", "lambda_adj",
        "scale_status", "scale_clamped",
        "lambda_scaling_mode", "lambda_scale_applied", "was_upscaled"
      )

  if (any(!is.finite(out$lambda_adj), na.rm = TRUE) || any(out$lambda_adj < 0, na.rm = TRUE)) {
    stop("lambda_adj must be finite and >= 0 for all rows.", call. = FALSE)
  }

  n_clamped <- sum(out$scale_clamped, na.rm = TRUE)
  if (n_clamped > 0) {
    message(sprintf("[lambda] Clamp applied to %d site/class rows.", n_clamped))
  }
  if (identical(scaling_mode, "down_only") && any(out$was_upscaled, na.rm = TRUE)) {
    stop("down_only scaling mode produced upscaling, which is not allowed.", call. = FALSE)
  }

  out
}

  #' Apply site/class lambda scalers to a TS/HUR lambda table
  #'
  #' @param lambda_table Tibble from compute_lambda_table().
  #' @param location Character scalar; site name.
  #' @param lambda_scalers Optional output from .lambda_scalers_from_rate_check().
  #'
  #' @return Tibble with adjusted lambda values. Existing columns are preserved;
  #'   lambda_raw, lambda_scale, and lambda_adj are added for traceability.
  #'
  #' @keywords internal
  .apply_lambda_scalers_to_lambda_table <- function(lambda_table,
                                                    location,
                                                    lambda_scalers = NULL) {
    lt <- tibble::as_tibble(lambda_table)
    if (!all(c("storm_class", "lambda") %in% names(lt))) {
      stop("lambda_table must contain storm_class and lambda.", call. = FALSE)
    }

    location <- as.character(location)[1]
    lt$storm_class <- as.character(lt$storm_class)
    lt$lambda <- as.numeric(lt$lambda)
    lt$lambda_raw <- lt$lambda
    lt$lambda_scale <- 1
    lt$lambda_adj <- lt$lambda

    if (is.null(lambda_scalers) || nrow(lambda_scalers) == 0) {
      return(lt)
    }

    scaler_tbl <- tibble::as_tibble(lambda_scalers) |>
      dplyr::filter(.data$location == .env$location)

    if (nrow(scaler_tbl) == 0) {
      return(lt)
    }

    scale_total <- scaler_tbl |>
      dplyr::filter(.data$storm_class == "TS34plus") |>
      dplyr::pull("lambda_scale")
    scale_hur <- scaler_tbl |>
      dplyr::filter(.data$storm_class == "HUR") |>
      dplyr::pull("lambda_scale")

    scale_total <- if (length(scale_total) > 0 && is.finite(scale_total[1])) scale_total[1] else 1
    scale_hur <- if (length(scale_hur) > 0 && is.finite(scale_hur[1])) scale_hur[1] else 1

    idx_ts <- which(lt$storm_class == "TS")
    idx_hur <- which(lt$storm_class == "HUR")

    lambda_ts_raw <- if (length(idx_ts) > 0) sum(lt$lambda_raw[idx_ts], na.rm = TRUE) else 0
    lambda_hur_raw <- if (length(idx_hur) > 0) sum(lt$lambda_raw[idx_hur], na.rm = TRUE) else 0
    lambda_total_raw <- lambda_ts_raw + lambda_hur_raw

    lambda_total_adj <- lambda_total_raw * scale_total
    lambda_hur_adj <- lambda_hur_raw * scale_hur

    consistency_clamped <- FALSE
    if (lambda_hur_adj > lambda_total_adj) {
      lambda_hur_adj <- lambda_total_adj
      consistency_clamped <- TRUE
    }

    lambda_ts_adj <- max(0, lambda_total_adj - lambda_hur_adj)

    if (length(idx_ts) > 0) {
      lt$lambda_adj[idx_ts] <- lambda_ts_adj
      lt$lambda_scale[idx_ts] <- if (lambda_ts_raw > 0) lambda_ts_adj / lambda_ts_raw else 1
      lt$lambda[idx_ts] <- lambda_ts_adj
    }
    if (length(idx_hur) > 0) {
      lt$lambda_adj[idx_hur] <- lambda_hur_adj
      lt$lambda_scale[idx_hur] <- if (lambda_hur_raw > 0) lambda_hur_adj / lambda_hur_raw else 1
      lt$lambda[idx_hur] <- lambda_hur_adj
    }

    attr(lt, "lambda_scaler_context") <- list(
      location = location,
      scale_total = scale_total,
      scale_hur = scale_hur,
      consistency_clamped = consistency_clamped
    )

    lt
  }

  #' Deterministic identifier for a lambda scaler table
  #'
  #' @param lambda_scalers Tibble from .lambda_scalers_from_rate_check().
  #'
  #' @return Character scalar ID stable for identical scaler content.
  #'
  #' @keywords internal
.lambda_scaler_id <- function(lambda_scalers) {
    if (is.null(lambda_scalers) || nrow(lambda_scalers) == 0) {
      return("lambda-scalers-none")
    }

    scaler_tbl <- tibble::as_tibble(lambda_scalers) |>
      dplyr::mutate(
        location = as.character(.data$location),
        storm_class = as.character(.data$storm_class),
        lambda_scale = as.numeric(.data$lambda_scale),
        scale_status = as.character(.data$scale_status)
      ) |>
      dplyr::arrange(.data$location, .data$storm_class)

    keys <- paste(
      scaler_tbl$location,
      scaler_tbl$storm_class,
      sprintf("%.6f", scaler_tbl$lambda_scale),
      scaler_tbl$scale_status,
      sep = "|"
    )
    code_points <- utf8ToInt(paste(keys, collapse = ";"))
    weights <- (seq_along(code_points) %% 251) + 1
    checksum <- sum(as.numeric(code_points) * weights)
    checksum <- checksum %% 4294967295

  sprintf("lambda-scalers-%08x", as.integer(checksum))
}

#' Deterministic checksum ID from text
#'
#' @param text Character vector.
#' @param prefix Character scalar prefix for the ID.
#'
#' @return Character scalar deterministic ID.
#'
#' @keywords internal
.checksum_id_from_text <- function(text, prefix = "id") {
  text <- paste(as.character(text), collapse = "|")
  code_points <- utf8ToInt(enc2utf8(text))
  if (length(code_points) == 0L) {
    return(sprintf("%s-%08x", prefix, as.integer(0)))
  }
  weights <- (seq_along(code_points) %% 251) + 1
  checksum <- sum(as.numeric(code_points) * weights)
  checksum <- checksum %% 4294967295
  sprintf("%s-%08x", prefix, as.integer(checksum))
}
