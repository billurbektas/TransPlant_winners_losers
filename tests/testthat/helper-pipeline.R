PROJECT_ROOT <- normalizePath(testthat::test_path("..", ".."))
project_path <- function(...) file.path(PROJECT_ROOT, ...)
source(project_path("config", "analysis.R"))
purrr::walk(list.files(project_path("R", "pipeline"), pattern = "[.]R$", full.names = TRUE), source)
