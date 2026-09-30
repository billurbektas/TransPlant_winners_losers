test_that("active modules and functions carry documentation contracts", {
  files <- c(
    active_function_files(PROJECT_ROOT),
    project_path("R", "targets", "stages.R"),
    project_path("R", "targets", "definitions", c("data.R", "models.R", "publication.R")),
    project_path("config", c("analysis.R", "workflow.R"))
  )
  files <- unique(files[file.exists(files)])

  for (path in files) {
    lines <- readLines(path, warn = FALSE)
    expect_true(
      any(grepl("^# Module:", lines)),
      info = paste("missing module header:", path)
    )
    function_lines <- grep("<- function", lines, fixed = TRUE)
    function_lines <- function_lines[
      grepl("^[A-Za-z.][A-Za-z0-9._]*[[:space:]]*<-", lines[function_lines])
    ]
    for (line in function_lines) {
      previous <- line - 1L
      while (previous > 0L && !nzchar(trimws(lines[[previous]]))) {
        previous <- previous - 1L
      }
      expect_true(
        previous > 0L && startsWith(trimws(lines[[previous]]), "#' @return"),
        info = paste("missing Roxygen return contract:", path, lines[[line]])
      )
    }
  }
})

test_that("active function and target names follow the semantic naming policy", {
  files <- c(
    active_function_files(PROJECT_ROOT),
    list.files(project_path("R", "targets"), pattern = "[.]R$",
               recursive = TRUE, full.names = TRUE)
  )
  definitions <- unlist(lapply(unique(files), function(path) {
    expressions <- parse(path)
    vapply(expressions, function(expression) {
      is_function <-
        is.call(expression) && identical(expression[[1]], as.name("<-")) &&
        is.symbol(expression[[2]]) && is.call(expression[[3]]) &&
        identical(expression[[3]][[1]], as.name("function"))
      if (is_function) as.character(expression[[2]]) else NA_character_
    }, character(1))
  }), use.names = FALSE)
  definitions <- stats::na.omit(definitions)
  forbidden <- "(^|_)(upstream|downstream|canonical|step[0-9]*|mboth)(_|$)|_526($|_)"
  expect_false(any(grepl(forbidden, definitions, ignore.case = TRUE)))

  previous_directory <- setwd(PROJECT_ROOT)
  on.exit(setwd(previous_directory), add = TRUE)
  target_names <- targets::tar_manifest(callr_function = NULL)$name
  expect_false(any(grepl(forbidden, target_names, ignore.case = TRUE)))
})

test_that("the public scripts and reports are the consolidated interfaces", {
  expect_setequal(
    basename(list.files(project_path("scripts"), pattern = "[.]R$")),
    c("run_pipeline.R", "render_figure.R", "refresh_dependencies.R")
  )
  expect_setequal(
    basename(list.files(project_path("reports"), pattern = "[.]qmd$")),
    c(
      "data_foundation.qmd", "response_estimation.qmd",
      "resident_climate_model.qmd", "family_context_model.qmd",
      "predictor_models.qmd", "model_synthesis.qmd",
      "release_validation.qmd"
    )
  )
})

test_that("current reports and artifact names have no historical identifiers", {
  reports <- list.files(
    project_path("reports"), pattern = "[.]qmd$", full.names = TRUE
  )
  report_text <- unlist(lapply(reports, readLines, warn = FALSE), use.names = FALSE)
  expect_false(any(grepl("tar_read\\((legacy_|.*equivalence)", report_text)))
  expect_false(any(grepl("_526|M_both|Step (19|22)", report_text)))

  migration <- readr::read_csv(
    project_path("docs", "name_migration.csv"), show_col_types = FALSE
  )
  artifacts <- dplyr::filter(migration, .data$category == "artifact")
  expect_true(nrow(artifacts) > 0L)
  expect_false(any(grepl(
    "_526|downstream|M_both|M_family|canonical",
    artifacts$current_name, ignore.case = TRUE
  )))
  expect_false(anyDuplicated(migration[c("category", "former_name")]) > 0L)
})
