library(targets)

source("config/analysis.R")
source("config/workflow.R")
targets::tar_source("R")

if (!nzchar(Sys.which("quarto"))) {
  bundled_quarto <- "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto"
  if (file.exists(bundled_quarto)) Sys.setenv(QUARTO_PATH = bundled_quarto)
}

workers <- as.integer(Sys.getenv("TRANSPLANT_WORKERS", "2"))
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
  source_targets(),
  community_targets(),
  covariate_targets(),
  response_targets(),
  resident_model_targets(),
  predictor_model_targets(),
  evaluation_targets(),
  publication_targets(),
  release_targets(),
  legacy_certification_targets()
)
