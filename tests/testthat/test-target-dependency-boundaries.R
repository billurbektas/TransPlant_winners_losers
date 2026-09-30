testthat::test_that("current scientific targets have no retired inputs", {
  previous_directory <- setwd(PROJECT_ROOT)
  on.exit(setwd(previous_directory), add = TRUE)
  manifest <- targets::tar_manifest(
    fields = c("name", "command", "pattern", "description"),
    callr_function = NULL
  )
  testthat::expect_false(anyDuplicated(manifest$name) > 0L)
  testthat::expect_false(any(grepl(
    "legacy|equivalence|step(19|22)|archive/|tests/fixtures|reports/reference_text",
    paste(manifest$name, manifest$command), ignore.case = TRUE
  )))
  testthat::expect_true(all(c(
    "community", "analysis_data", "resident_auc_data",
    "resident_sampling_variance_matrix", "resident_climate_model",
    "family_model_bundle", "phylogeny_model_pair", "release_validation",
    "report_analysis_workflow"
  ) %in% manifest$name))
  testthat::expect_identical(
    manifest$name[startsWith(manifest$name, "report_")],
    "report_analysis_workflow"
  )
})

testthat::test_that("functional status and phylogeny derive from reviewed inputs", {
  text <- paste(readLines(project_path("analysis_workflow.qmd"), warn = FALSE),
                collapse = "\n")
  testthat::expect_match(text, "functional_group_overrides", fixed = TRUE)
  testthat::expect_match(text, "resident_species_overrides", fixed = TRUE)
  testthat::expect_match(text, "build_phylogeny_from_placement_decisions", fixed = TRUE)
  testthat::expect_match(text, "phylogeny_placement_decisions", fixed = TRUE)
  testthat::expect_false(grepl("tar_read(legacy_", text, fixed = TRUE))
})

testthat::test_that("the source inventory contains the active review document", {
  files <- active_scientific_source_files(PROJECT_ROOT)
  testthat::expect_true(all(file.exists(files)))
  testthat::expect_true(all(active_function_files(PROJECT_ROOT) %in% files))
  testthat::expect_true(project_path("analysis_workflow.qmd") %in% files)
  testthat::expect_false(any(grepl("/archive/|/tests/fixtures/", files)))
})

testthat::test_that("each active top-level function is defined once", {
  files <- active_function_files(PROJECT_ROOT)
  definitions <- unlist(lapply(files, function(path) {
    expressions <- parse(path)
    names <- vapply(expressions, function(expression) {
      is_function_assignment <-
        is.call(expression) && identical(expression[[1]], as.name("<-")) &&
        is.symbol(expression[[2]]) && is.call(expression[[3]]) &&
        identical(expression[[3]][[1]], as.name("function"))
      if (is_function_assignment) as.character(expression[[2]]) else NA_character_
    }, character(1))
    stats::na.omit(names)
  }), use.names = FALSE)
  testthat::expect_identical(anyDuplicated(definitions), 0L)
})
