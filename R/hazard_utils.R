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
