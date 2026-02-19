#' @keywords internal
"_PACKAGE"

#' Package-wide imports and NSE bindings.
#'
#' @name ipdcstorm-imports
#' @keywords internal
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @import purrr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats dpois median na.omit qnorm quantile rnorm setNames
#' @importFrom utils modifyList globalVariables
NULL

utils::globalVariables(
  c(
    "BASIN", "SEASON", "SID", "annual_rate", "count", "doy", "dur_days",
    "end_date", "event_class", "event_id", "exceedance_prob", "expected",
    "extreme", "iso_time", "lat", "location", "lon", "max_wind",
    "max_wind_kt", "mean_wind", "month", "n", "n_days", "n_events",
    "poci_hpa", "pres_hpa", "pres_source", "r34_ne_nm", "r34_nw_nm",
    "r34_se_nm", "r34_sw_nm", "r50_ne_nm", "r50_nw_nm", "r50_se_nm",
    "r50_sw_nm", "r64_ne_nm", "r64_nw_nm", "r64_se_nm", "r64_sw_nm",
    "return_period", "rmw_km", "roci_km", "sim_year", "start_date",
    "start_doy", "start_month", "stat", "storm_name", "storm_status",
    "upper", "usa_pres_hpa", "wind_kt", "wmo_pres_hpa"
  )
)
