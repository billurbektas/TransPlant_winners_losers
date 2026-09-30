testthat::test_that("display chunks are assignment-free and tangling is exact", {
  path <- project_path("analysis_workflow.qmd")
  lines <- readLines(path, warn = FALSE)
  fence_starts <- grep("^```\\{r", lines)
  display_assignments <- character()
  intended_targets <- character()
  for (start in fence_starts) {
    relative_end <- which(lines[(start + 1L):length(lines)] == "```")[[1]]
    end <- start + relative_end
    chunk <- lines[(start + 1L):(end - 1L)]
    code <- chunk[!grepl("^#\\|", chunk)]
    expressions <- parse(text = code)
    assigned <- vapply(expressions, function(expression) {
      is.call(expression) &&
        as.character(expression[[1]])[[1]] %in% c("<-", "=", "->")
    }, logical(1))
    if (any(grepl("^#\\| eval: true$", chunk))) {
      if (any(assigned)) {
        display_assignments <- c(
          display_assignments,
          vapply(expressions[assigned], deparse1, character(1))
        )
      }
    } else if (any(assigned)) {
      intended_targets <- c(
        intended_targets,
        vapply(expressions[assigned], function(expression) {
          if (identical(as.character(expression[[1]])[[1]], "->")) {
            as.character(expression[[3]])
          } else {
            as.character(expression[[2]])
          }
        }, character(1))
      )
    }
  }
  testthat::expect_length(display_assignments, 0L)

  tangled <- tarchetypes::tar_tangle(path)
  names <- vapply(tangled, function(target) target$name, character(1))
  testthat::expect_setequal(names, intended_targets)
  testthat::expect_false(anyDuplicated(names) > 0L)
  testthat::expect_setequal(
    intersect(
      names,
      c("input_inventory", "release_validation", "figure01_winners_losers")
    ),
    character()
  )
})

testthat::test_that("the manifest contains the review workflow and no history", {
  previous_directory <- setwd(PROJECT_ROOT)
  on.exit(setwd(previous_directory), add = TRUE)
  manifest <- targets::tar_manifest(callr_function = NULL)
  testthat::expect_false(anyDuplicated(manifest$name) > 0L)
  testthat::expect_false(any(grepl(
    "legacy|equivalence|step(19|22)", manifest$name, ignore.case = TRUE
  )))
  major <- c(
    "community", "analysis_data", "resident_sampling_variance_matrix",
    "resident_climate_model", "ldmc_model_pair",
    "model_cross_validation_metrics", "figure06_predictor_synthesis",
    "release_validation", "report_analysis_workflow"
  )
  selected <- manifest[match(major, manifest$name), ]
  testthat::expect_false(anyNA(selected$name))
  testthat::expect_true(all(
    !is.na(selected$description) & nzchar(selected$description)
  ))
})

testthat::test_that("regrouped results keep their separate scientific exports", {
  previous_directory <- setwd(PROJECT_ROOT)
  on.exit(setwd(previous_directory), add = TRUE)
  manifest <- targets::tar_manifest(callr_function = NULL)
  removed_projections <- c(
    "resident_model_data", "climate_scaling", "species_classification",
    "family_auc_data", "family_sampling_variance_matrix",
    "predictor_model_fit", "predictor_model_comparisons", "predictor_CR1",
    "heterogeneity_components", "heterogeneity_synthesis",
    "climate_context_experiments", "family_random_coefficients",
    "resident_phylogeny_candidate", "phylogenetic_distance",
    "pcoa_eigen_summary", "pev_scaling"
  )
  testthat::expect_false(any(removed_projections %in% manifest$name))
  retained_exports <- c(
    "resident_data_file", "species_classification_file",
    "family_data_file", "family_sampling_variance_file",
    "predictor_fit_file", "predictor_comparison_file", "predictor_CR1_file",
    "heterogeneity_components_file", "heterogeneity_synthesis_file",
    "climate_context_experiments_file", "family_random_file",
    "phylogeny_file", "phylogeny_newick_file",
    "phylogenetic_distance_file", "pcoa_eigen_file", "pev_scaling_file"
  )
  testthat::expect_true(all(retained_exports %in% manifest$name))
  testthat::expect_match(
    manifest$command[match("family_sampling_variance_file", manifest$name)],
    "family_auc_bundle$sampling_variance_matrix", fixed = TRUE
  )
  testthat::expect_match(
    manifest$command[match("predictor_comparison_file", manifest$name)],
    "predictor_model_tables$comparisons", fixed = TRUE
  )
})

testthat::test_that("model specifications are literal in the review document", {
  text <- paste(
    readLines(project_path("analysis_workflow.qmd"), warn = FALSE),
    collapse = "\n"
  )
  required <- c(
    "resident_climate_formula <- ~ T_z * P_z * W_z",
    "~ 1 | species_id",
    "~ 0 + W_z | species_id",
    "~ 0 + T_z + P_z + TP_z | species_id",
    "model_random_covariance_structure <- \"GEN\"",
    "trait_full_formula <- ~ Trait_z * T_z * P_z * W_z",
    "graminoid_full_formula <- ~ GraminoidStatus * T_z * P_z * W_z",
    "woodiness_full_formula <- ~ Woodiness * T_z * P_z * W_z",
    "method = \"ML\"",
    "method = \"REML\"",
    "comparison_method = \"ML\"",
    "final_method = \"REML\"",
    "sparse = TRUE",
    "robust_cluster_column <- \"experiment_id\"",
    "robust_small_sample_adjustment <- TRUE",
    "sampling_variance_matrix = resident_sampling_variance_matrix"
  )
  for (value in required) {
    testthat::expect_match(text, value, fixed = TRUE)
  }
  testthat::expect_false(grepl(
    "additive_formula|additive_ML|additive_REML|M1_model", text
  ))
})

testthat::test_that("migration map gives every historical path one disposition", {
  history <- readr::read_csv(
    project_path("docs", "history_retirement_manifest.csv"),
    show_col_types = FALSE
  )
  migration <- readr::read_csv(
    project_path("docs", "workflow_migration_map.csv"),
    show_col_types = FALSE
  )
  allowed <- c("migrated", "merged", "renamed", "intentionally retired")
  testthat::expect_false(anyDuplicated(migration$old_component) > 0L)
  testthat::expect_true(all(migration$status %in% allowed))
  testthat::expect_setequal(
    history$path,
    intersect(history$path, migration$old_component)
  )
  testthat::expect_true(all(
    nzchar(migration$new_component) & nzchar(migration$note)
  ))
  replacements <- migration[
    migration$status %in% c("migrated", "merged", "renamed"),
  ]
  testthat::expect_false(any(startsWith(replacements$new_component, "none")))

  renamed <- migration[migration$status == "renamed", ]
  function_files <- active_function_files(PROJECT_ROOT)
  function_names <- unlist(lapply(function_files, function(path) {
    expressions <- parse(path)
    names <- vapply(expressions, function(expression) {
      is_definition <-
        is.call(expression) && identical(expression[[1]], as.name("<-")) &&
        is.symbol(expression[[2]]) && is.call(expression[[3]]) &&
        identical(expression[[3]][[1]], as.name("function"))
      if (is_definition) as.character(expression[[2]]) else NA_character_
    }, character(1))
    stats::na.omit(names)
  }), use.names = FALSE)
  previous_directory <- setwd(PROJECT_ROOT)
  on.exit(setwd(previous_directory), add = TRUE)
  target_names <- targets::tar_manifest(callr_function = NULL)$name
  resolves <- ifelse(
    renamed$component_type == "function",
    renamed$new_component %in% function_names,
    renamed$new_component %in% target_names
  )
  testthat::expect_true(all(resolves))

  nonretired_targets <- migration[
    migration$component_type == "target" &
      migration$status != "intentionally retired",
  ]
  replacement_targets <- sub(
    " .*", "", sub("^target:", "", nonretired_targets$new_component)
  )
  replacement_targets <- sub("\\$.*", "", replacement_targets)
  testthat::expect_true(all(replacement_targets %in% target_names))

  nonretired_files <- migration[
    migration$component_type == "file" &
      migration$status != "intentionally retired",
  ]
  file_replacements_resolve <- vapply(
    nonretired_files$new_component,
    function(replacement) {
      if (startsWith(replacement, "target:")) {
        target <- sub("\\$.*", "", sub(" .*", "", sub("^target:", "", replacement)))
        return(target %in% target_names)
      }
      path <- sub(" (plus|/) .*", "", replacement)
      file.exists(project_path(path))
    },
    logical(1)
  )
  testthat::expect_true(all(file_replacements_resolve))

  nonretired_functions <- migration[
    migration$component_type == "function" &
      migration$status != "intentionally retired",
  ]
  function_replacements_resolve <- vapply(
    nonretired_functions$new_component,
    function(replacement) {
      if (grepl("analysis_workflow.qmd", replacement, fixed = TRUE)) {
        return(file.exists(project_path("analysis_workflow.qmd")))
      }
      replacement %in% function_names
    },
    logical(1)
  )
  testthat::expect_true(all(function_replacements_resolve))

  factory_targets <- vapply(
    c(data_target_definitions(), model_target_definitions(),
      publication_target_definitions()),
    function(target) target$name,
    character(1)
  )
  testthat::expect_setequal(
    paste0("target:", factory_targets),
    intersect(paste0("target:", factory_targets), migration$old_component)
  )
})
