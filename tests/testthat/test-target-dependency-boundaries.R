test_that("ordinary reports do not read the late Step 22 assembly", {
  report_paths <- list.files(
    project_path("reports"), pattern = "^[0-9]{2}_.+[.]qmd$", full.names = TRUE
  )
  ordinary_paths <- report_paths[!grepl("16_reference_and_release[.]qmd$", report_paths)]
  report_text <- unlist(lapply(ordinary_paths, readLines, warn = FALSE), use.names = FALSE)
  expect_false(any(grepl("candidate_step22_results", report_text, fixed = TRUE)))
  expect_false(any(grepl("step22_reference_validation", report_text, fixed = TRUE)))
})

test_that("canonical reports retain the required scientific sections", {
  cleaning <- sprintf("%02d", 1:8)
  analysis <- sprintf("%02d", 9:15)

  report_text <- function(step) {
    path <- list.files(
      project_path("reports"), pattern = paste0("^", step, "_.*[.]qmd$"),
      full.names = TRUE
    )
    expect_length(path, 1L)
    paste(readLines(path, warn = FALSE), collapse = "\n")
  }

  purrr::walk(vapply(cleaning, report_text, character(1)), function(text) {
    expect_match(text, "## Purpose", fixed = TRUE)
    expect_match(text, "## Inputs", fixed = TRUE)
    expect_match(text, "## Decisions and exclusions", fixed = TRUE)
    expect_match(text, "## Before/after counts and validation checks", fixed = TRUE)
    expect_match(text, "## Problems flagged", fixed = TRUE)
    expect_match(text, "## Canonical outputs", fixed = TRUE)
    expect_match(text, "## PASS/FAIL", fixed = TRUE)
  })
  purrr::walk(vapply(analysis, report_text, character(1)), function(text) {
    expect_match(text, "## Question", fixed = TRUE)
    expect_match(text, "## Dataset and sample size", fixed = TRUE)
    expect_match(text, "## Exact model specification", fixed = TRUE)
    expect_match(text, "## Assumptions", fixed = TRUE)
    expect_match(text, "## Convergence and covariance diagnostics", fixed = TRUE)
    expect_match(text, "## Main numerical results", fixed = TRUE)
    expect_match(text, "## Interpretation boundary", fixed = TRUE)
    expect_match(text, "## Canonical objects", fixed = TRUE)
    expect_match(text, "## PASS/FAIL", fixed = TRUE)
  })
})

test_that("trait, functional, and phylogeny reports read branch-specific tables", {
  read_report <- function(name) {
    paste(readLines(project_path("reports", name), warn = FALSE), collapse = "\n")
  }
  trait <- read_report("12_trait_information.qmd")
  functional <- read_report("13_functional_information.qmd")
  phylogeny <- read_report("14_phylogeny.qmd")

  expect_match(trait, "trait_model_tables_526", fixed = TRUE)
  expect_false(grepl("tar_read(predictor_model_fit)", trait, fixed = TRUE))
  expect_match(functional, "functional_model_tables_526", fixed = TRUE)
  expect_match(functional, "woodiness_figure_file", fixed = TRUE)
  expect_false(grepl("tar_read(predictor_model_fit)", functional, fixed = TRUE))
  expect_match(phylogeny, "phylogeny_model_tables_526", fixed = TRUE)
  expect_false(grepl("tar_read(predictor_model_fit)", phylogeny, fixed = TRUE))
})

test_that("the Woodiness figure has a narrowly scoped plotting dependency", {
  source_text <- paste(
    readLines(
      project_path("R", "targets", "publication_targets.R"), warn = FALSE
    ),
    collapse = "\n"
  )
  start <- regexpr("woodiness_figure,", source_text, fixed = TRUE)[[1]]
  end <- regexpr("woodiness_figure_file,", source_text, fixed = TRUE)[[1]]
  expect_gt(start, 0L)
  expect_gt(end, start)
  block <- substring(source_text, start, end - 1L)
  expect_match(block, "woodiness_model_sequence_526", fixed = TRUE)
  expect_match(block, "woodiness_support_mask", fixed = TRUE)
  expect_false(grepl("candidate_step22_results", block, fixed = TRUE))
})

test_that("canonical exports are independent of the late validation object", {
  source_text <- readLines(
    project_path("R", "targets", "analysis_targets.R"), warn = FALSE
  )
  late_assembly_line <- grep("candidate_step22_results", source_text, fixed = TRUE)[[1]]
  primary_exports <- c(
    "canonical_downstream_data_file,", "canonical_mboth_model_file,",
    "canonical_species_classification_file,", "canonical_family_model_file,",
    "canonical_trait_scaling_file,", "canonical_functional_species_file,",
    "canonical_phylogeny_file,", "canonical_woodiness_model_file,",
    "canonical_heterogeneity_synthesis_file,"
  )
  export_lines <- vapply(primary_exports, function(target_name) {
    grep(target_name, source_text, fixed = TRUE)[[1]]
  }, integer(1))
  expect_true(all(export_lines < late_assembly_line))

  manifest_start <- grep("canonical_numerical_export_files,", source_text,
                         fixed = TRUE)[[1]]
  manifest_end <- grep("canonical_output_manifest_file,", source_text,
                       fixed = TRUE)[[1]]
  manifest_block <- paste(source_text[manifest_start:manifest_end], collapse = "\n")
  expect_false(grepl("candidate_step22_results", manifest_block, fixed = TRUE))
})

test_that("functional assignments use the active source plus the explicit ledger", {
  source_text <- readLines(
    project_path("R", "targets", "analysis_targets.R"), warn = FALSE
  )
  start <- grep("functional_species_526,", source_text, fixed = TRUE)[[1]]
  end <- grep("functional_source_coverage_audit,", source_text, fixed = TRUE)[[1]]
  functional_block <- paste(source_text[start:(end - 1L)], collapse = "\n")
  expect_match(functional_block, "canonical_functional_groups_input", fixed = TRUE)
  expect_match(functional_block, "downstream_species_overrides", fixed = TRUE)
  expect_false(grepl("trait_species_reference_526", functional_block, fixed = TRUE))
})

test_that("the candidate phylogeny is rebuilt rather than adopted from Step 19", {
  source_text <- paste(
    readLines(
      project_path("R", "targets", "analysis_targets.R"), warn = FALSE
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
})

test_that("canonical command-line entry points do not source historical runners", {
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

test_that("the target graph loads an explicit canonical module list", {
  target_text <- paste(
    readLines(project_path("_targets.R"), warn = FALSE), collapse = "\n"
  )
  config_text <- paste(
    readLines(project_path("config", "canonical_analysis.R"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(target_text, "canonical_pipeline_source_files()", fixed = TRUE)
  expect_false(grepl("list.files(\"R/pipeline\"", target_text, fixed = TRUE))
  expect_false(grepl("additional_analysis.R", config_text, fixed = TRUE))
  expect_false(grepl("frequentist_hierarchical_auc.R", config_text, fixed = TRUE))
  expect_false(grepl('"figures.R"', config_text, fixed = TRUE))
  expect_true(all(file.exists(project_path("R", "functions", c(
    "upstream_analysis.R", "auc_model_helpers.R", "temporal_figures.R",
    "climate_display_figure.R", "consistency_figures.R"
  )))))
  expect_false(grepl("biological_context_review.R", config_text, fixed = TRUE))
  expect_false(grepl("climate_response_model.R", config_text, fixed = TRUE))
  expect_false(grepl("hierarchical_trait_trial.R", config_text, fixed = TRUE))
})

test_that("active modules define each top-level function exactly once", {
  files <- c(
    canonical_pipeline_source_files(PROJECT_ROOT),
    list.files(
      project_path("R", "functions"), pattern = "[.]R$", full.names = TRUE
    )
  )
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

test_that("Step 01 inventories and reports active scientific source hashes", {
  report_text <- paste(
    readLines(project_path("reports", "01_inputs_and_scope.qmd"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(report_text, "canonical_scientific_source_inventory", fixed = TRUE)

  files <- canonical_scientific_source_files(PROJECT_ROOT)
  expect_true(all(file.exists(files)))
  expect_true(all(canonical_pipeline_source_files(PROJECT_ROOT) %in% files))
  expect_false(any(grepl("/archive/", files, fixed = TRUE)))
  expect_false(any(grepl("/_targets/", files, fixed = TRUE)))
})

test_that("canonical analysis targets do not read superseded output trees", {
  target_source <- paste(
    readLines(
      testthat::test_path("..", "..", "R", "targets", "analysis_targets.R"),
      warn = FALSE
    ),
    collapse = "\n"
  )
  expect_false(grepl("output/step_by_step_modeling", target_source, fixed = TRUE))
  expect_false(grepl("output/frequentist_hierarchical_auc", target_source,
                     fixed = TRUE))
  expect_match(target_source, "canonical_transition_reference.rds", fixed = TRUE)
})
