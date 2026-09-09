PROJECT_ROOT <- normalizePath(testthat::test_path("..", ".."))
project_path <- function(...) file.path(PROJECT_ROOT, ...)
source(project_path("config", "analysis.R"))
source(project_path("config", "canonical_analysis.R"))
active_pipeline_files <- canonical_pipeline_source_files(PROJECT_ROOT)
if (any(!file.exists(active_pipeline_files))) {
  stop("One or more canonical upstream pipeline modules are missing.")
}
purrr::walk(active_pipeline_files, source)
canonical_function_files <- list.files(
  project_path("R", "functions"), pattern = "[.]R$", full.names = TRUE
)
purrr::walk(canonical_function_files, source)
