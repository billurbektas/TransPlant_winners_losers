test_that("ordinary thematic reports do not read historical assemblies", {
  report_paths <- project_path("reports", c(
    "data_foundation.qmd", "response_estimation.qmd",
    "resident_climate_model.qmd", "family_context_model.qmd",
    "predictor_models.qmd", "model_synthesis.qmd",
    "release_validation.qmd"
  ))
  ordinary_paths <- report_paths[!grepl("release_validation[.]qmd$", report_paths)]
  report_text <- unlist(lapply(ordinary_paths, readLines, warn = FALSE), use.names = FALSE)
  expect_false(any(grepl("legacy_comparison_candidate", report_text, fixed = TRUE)))
  expect_false(any(grepl("legacy_reference_validation", report_text, fixed = TRUE)))
})

test_that("the seven thematic reports retain auditable scientific sections", {
  reports <- list.files(
    project_path("reports"), pattern = "^[a-z].+[.]qmd$", full.names = TRUE
  )
  expect_setequal(basename(reports), c(
    "data_foundation.qmd", "response_estimation.qmd",
    "resident_climate_model.qmd", "family_context_model.qmd",
    "predictor_models.qmd", "model_synthesis.qmd",
    "release_validation.qmd"
  ))
  text <- vapply(reports, function(path) {
    paste(readLines(path, warn = FALSE), collapse = "\n")
  }, character(1))
  expect_true(all(grepl("PASS/FAIL", text, fixed = TRUE)))
  expect_true(all(grepl("Principal objects|Principal outputs", text)))
})

test_that("trait, functional, and phylogeny reports read branch-specific tables", {
  read_report <- function(name) {
    paste(readLines(project_path("reports", name), warn = FALSE), collapse = "\n")
  }
  trait <- read_report("predictor_models.qmd")
  functional <- trait
  phylogeny <- trait

  expect_match(trait, "trait_model_tables", fixed = TRUE)
  expect_false(grepl("tar_read(predictor_model_fit)", trait, fixed = TRUE))
  expect_match(functional, "functional_model_tables", fixed = TRUE)
  expect_match(functional, "woodiness_support_mask_audit", fixed = TRUE)
  expect_false(grepl("woodiness_figure_file", functional, fixed = TRUE))
  expect_false(grepl("tar_read(predictor_model_fit)", functional, fixed = TRUE))
  expect_match(phylogeny, "phylogeny_model_tables", fixed = TRUE)
  expect_false(grepl("tar_read(predictor_model_fit)", phylogeny, fixed = TRUE))
})

test_that("the obsolete standalone Woodiness figure is not generated", {
  source_text <- paste(
    readLines(
      project_path("R", "targets", "definitions", "publication.R"), warn = FALSE
    ),
    collapse = "\n"
  )
  expect_false(grepl("woodiness_figure,", source_text, fixed = TRUE))
  expect_false(grepl("woodiness_figure_file,", source_text, fixed = TRUE))
  expect_match(source_text, "woodiness_support_mask,", fixed = TRUE)
})

test_that("current exports are independent of the late legacy assembly", {
  source_text <- readLines(
    project_path("R", "targets", "definitions", "models.R"), warn = FALSE
  )
  late_assembly_line <- grep("legacy_comparison_candidate", source_text, fixed = TRUE)[[1]]
  primary_exports <- c(
    "resident_data_file,", "resident_climate_model_file,",
    "species_classification_file,", "family_model_file,",
    "trait_scaling_file,", "functional_species_file,",
    "phylogeny_file,", "woodiness_model_file,",
    "heterogeneity_synthesis_file,"
  )
  export_lines <- vapply(primary_exports, function(target_name) {
    grep(target_name, source_text, fixed = TRUE)[[1]]
  }, integer(1))
  expect_true(all(export_lines < late_assembly_line))

  manifest_start <- grep("release_artifact_files,", source_text,
                         fixed = TRUE)[[1]]
  manifest_end <- grep("release_artifact_manifest_file,", source_text,
                       fixed = TRUE)[[1]]
  manifest_block <- paste(source_text[manifest_start:manifest_end], collapse = "\n")
  expect_false(grepl("legacy_comparison_candidate", manifest_block, fixed = TRUE))
  expect_match(manifest_block, "vapply(unname(paths), sha256_file", fixed = TRUE)
})

test_that("functional assignments use the active source plus the explicit ledger", {
  source_text <- readLines(
    project_path("R", "targets", "definitions", "models.R"), warn = FALSE
  )
  start <- grep("functional_species,", source_text, fixed = TRUE)[[1]]
  end <- grep("functional_source_coverage_audit,", source_text, fixed = TRUE)[[1]]
  functional_block <- paste(source_text[start:(end - 1L)], collapse = "\n")
  expect_match(functional_block, "functional_groups_input", fixed = TRUE)
  expect_match(functional_block, "resident_species_overrides", fixed = TRUE)
  expect_false(grepl("legacy_trait_species", functional_block, fixed = TRUE))
})

test_that("the candidate phylogeny is rebuilt rather than adopted from history", {
  source_text <- paste(
    readLines(
      project_path("R", "targets", "definitions", "models.R"), warn = FALSE
    ),
    collapse = "\n"
  )
  expect_match(
    source_text,
    "resident_phylogeny_candidate, resident_phylogeny_rebuild$tree",
    fixed = TRUE
  )
  expect_match(
    source_text,
    "build_phylogeny_from_placement_decisions",
    fixed = TRUE
  )
  expect_match(source_text, "legacy_file_artifacts", fixed = TRUE)
  expect_match(
    source_text, "phylogeny_package_installation_files", fixed = TRUE
  )
})

test_that("current command-line entry points do not source historical runners", {
  scripts <- project_path(
    "scripts", c("run_pipeline.R", "refresh_dependencies.R")
  )
  text <- paste(
    unlist(lapply(scripts, readLines, warn = FALSE), use.names = FALSE),
    collapse = "\n"
  )
  expect_false(grepl("source(\"R/", text, fixed = TRUE))
  expect_match(text, "jinyizju/V.PhyloMaker2@", fixed = TRUE)
  expect_match(text, "tarchetypes", fixed = TRUE)
  expect_match(text, "quarto", fixed = TRUE)
})

test_that("the target graph loads the stage-organized R tree", {
  target_text <- paste(
    readLines(project_path("_targets.R"), warn = FALSE), collapse = "\n"
  )
  config_text <- paste(
    readLines(project_path("config", "workflow.R"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(target_text, "targets::tar_source(\"R\")", fixed = TRUE)
  expect_false(grepl("list.files(\"R/pipeline\"", target_text, fixed = TRUE))
  expect_false(dir.exists(project_path("R", "pipeline")))
  expect_false(dir.exists(project_path("R", "reporting")))
  expect_false(grepl("additional_analysis.R", config_text, fixed = TRUE))
  expect_false(grepl("frequentist_hierarchical_auc.R", config_text, fixed = TRUE))
  expect_false(grepl('"figures.R"', config_text, fixed = TRUE))
  expect_true(all(file.exists(c(
    project_path("R", "data", "analysis_dataset.R"),
    project_path("R", "responses", "auc_model_helpers.R"),
    project_path("R", "figures", "temporal.R"),
    project_path("R", "figures", "climate.R"),
    project_path("R", "figures", "consistency.R")
  ))))
  expect_false(grepl("biological_context_review.R", config_text, fixed = TRUE))
  expect_false(grepl("climate_response_model.R", config_text, fixed = TRUE))
  expect_false(grepl("hierarchical_trait_trial.R", config_text, fixed = TRUE))
})

test_that("active modules define each top-level function exactly once", {
  files <- unique(c(
    active_function_files(PROJECT_ROOT),
    list.files(
      project_path("R"), pattern = "[.]R$", recursive = TRUE,
      full.names = TRUE
    )
  ))
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
  expect_identical(anyDuplicated(definitions), 0L)
})

test_that("the data-foundation report inventories active source hashes", {
  report_text <- paste(
    readLines(project_path("reports", "data_foundation.qmd"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(report_text, "active_source_inventory", fixed = TRUE)

  files <- active_scientific_source_files(PROJECT_ROOT)
  expect_true(all(file.exists(files)))
  expect_true(all(active_function_files(PROJECT_ROOT) %in% files))
  expect_false(any(grepl("/archive/", files, fixed = TRUE)))
  expect_false(any(grepl("/_targets/", files, fixed = TRUE)))
})

test_that("current analysis targets do not read superseded output trees", {
  target_source <- paste(
    readLines(
      testthat::test_path("..", "..", "R", "targets", "definitions", "models.R"),
      warn = FALSE
    ),
    collapse = "\n"
  )
  expect_false(grepl("output/step_by_step_modeling", target_source, fixed = TRUE))
  expect_false(grepl("output/frequentist_hierarchical_auc", target_source,
                     fixed = TRUE))
  expect_match(target_source, "canonical_transition_reference.rds", fixed = TRUE)
})
