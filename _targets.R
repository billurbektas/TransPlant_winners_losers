library(targets)

source("config/analysis.R")
source("config/canonical_analysis.R")
active_pipeline_files <- canonical_pipeline_source_files()
if (any(!file.exists(active_pipeline_files))) {
  stop("One or more canonical upstream pipeline modules are missing.")
}
purrr::walk(active_pipeline_files, source)
purrr::walk(list.files("R/functions", pattern = "[.]R$", full.names = TRUE), source)
source("R/targets/data_targets.R")
source("R/targets/analysis_targets.R")
source("R/targets/publication_targets.R")

if (!nzchar(Sys.which("quarto"))) {
  bundled_quarto <- "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto"
  if (file.exists(bundled_quarto)) Sys.setenv(QUARTO_PATH = bundled_quarto)
}

workers <- as.integer(Sys.getenv("TRANSPLANT_CHELSA_WORKERS", "2"))
controller <- if (requireNamespace("crew", quietly = TRUE)) {
  crew::crew_controller_local(workers = workers)
} else {
  NULL
}

tar_option_set(
  packages = c(
    "ape", "cowplot", "DescTools", "TNRS", "dplyr", "emmeans", "geometry", "ggalluvial",
    "ggplot2", "ggrepel", "lme4", "Matrix", "mclogit", "metafor", "nnet",
    "patchwork", "purrr", "readr", "sandwich", "scatterpie", "sf", "stringi",
    "stringr", "terra", "tibble", "tidyTRY", "tidyr"
  ),
  controller = controller,
  error = "stop",
  seed = 20260824L
)

c(
  canonical_data_targets(),
  canonical_analysis_targets(),
  canonical_publication_targets()
)
