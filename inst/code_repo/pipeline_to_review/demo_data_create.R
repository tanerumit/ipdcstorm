# Run once to create inst/extdata/ibtracs_demo.csv
library(readr)
library(dplyr)
library(geosphere)

# --- Config ---
full_ibtracs  <- "inst/extdata/ibtracs/ibtracs.NA.list.v04r01.csv"
output_path   <- "inst/extdata/ibtracs_demo.csv"
target_coords <- c(lon = -63.05, lat = 17.63)  # St. Martin
radius_km     <- 800
start_year    <- 1970L

# --- Read & filter ---
df <- read_csv(full_ibtracs, show_col_types = FALSE, guess_max = 5000)

required_cols <- c("SID", "SEASON", "BASIN", "ISO_TIME",
                   "USA_LAT", "USA_LON", "USA_WIND",
                   "USA_R34_NE", "USA_R34_SE", "USA_R34_SW", "USA_R34_NW")

df <- df |>
  mutate(
    SEASON  = as.integer(SEASON),
    USA_LAT = as.numeric(USA_LAT),
    USA_LON = as.numeric(USA_LON)
  ) |>
  filter(SEASON >= start_year, !is.na(USA_LAT), !is.na(USA_LON))

# Find storms that pass within radius
df <- df |>
  mutate(dist_km = distHaversine(
    cbind(USA_LON, USA_LAT),
    matrix(target_coords[c("lon", "lat")], nrow = nrow(df), ncol = 2, byrow = TRUE)
  ) / 1000)

nearby_sids <- df |>
  filter(dist_km <= radius_km) |>
  pull(SID) |>
  unique()

# Keep ALL track points for those storms (needed for wind field computation)
demo <- df |>
  filter(SID %in% nearby_sids) |>
  select(all_of(required_cols)) |>
  arrange(SID, ISO_TIME)

# --- Write ---
write_csv(demo, output_path)



# Verify
#demo <- read_csv(output_path, show_col_types = FALSE)
#cat("Columns:", paste(names(demo), collapse = ", "), "\n")
#cat("Storms:", n_distinct(demo$SID), "| Rows:", nrow(demo), "\n")
#cat("Year range:", range(as.integer(demo$SEASON), na.rm = TRUE), "\n")


read_ibtracs_clean(ibtracs_csv = "C:/Users/taner/AppData/Local/R/win-library/4.5/ipdcstorm/extdata/ibtracs_demo.csv")

path <- system.file("extdata", "ibtracs_demo.csv", package = "ipdcstorm")
df <- readr::read_csv(path, na = c("", " ", ".", "-"), show_col_types = FALSE, n_max = 5)
names(df)

