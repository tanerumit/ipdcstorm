################################################################################
# Hazard model + temporal downscaling + damage forcing
# Legacy wrapper: source split hazard_* files
################################################################################

this_file <- sys.frame(1)$ofile
this_dir <- if (!is.null(this_file)) dirname(this_file) else getwd()
base_dir <- normalizePath(file.path(this_dir, ".."), mustWork = FALSE)

source(file.path(base_dir, "hazard_utils.R"))
source(file.path(base_dir, "hazard_ibtracs.R"))
source(file.path(base_dir, "hazard_core.R"))
source(file.path(base_dir, "hazard_downscale.R"))
source(file.path(base_dir, "hazard_run.R"))
