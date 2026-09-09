source(project_path("R", "functions", "figures.R"))
source(project_path("R", "functions", "provenance.R"))

test_that("functional coverage is computed independently by group", {
  data <- tibble::tibble(
    species_id = c("a", "b", "c", "d"),
    experiment_id = c("e1", "e2", "e1", "e3"),
    T_z = c(-1, 0, 2, 3), P_z = c(-2, 0, 1, 4), W_z = c(-1, 1, -2, 2),
    Woodiness = factor(c("non_woody", "non_woody", "woody", "woody"))
  )
  coverage <- functional_climate_coverage(data, "Woodiness", "Woodiness")
  expect_setequal(as.character(coverage$group), c("non_woody", "woody"))
  expect_equal(coverage$T_z_min[coverage$group == "woody"], 2)
  expect_equal(coverage$T_z_max[coverage$group == "non_woody"], 0)
})

test_that("the Woodiness figure fails closed without a group-specific mask", {
  function_text <- paste(
    deparse(body(build_woodiness_response_surface_figure)), collapse = "\n"
  )
  expect_match(function_text, "group-specific observed climate support", fixed = TRUE)
  expect_match(function_text, "canonical_predictor_formulas()$Woodiness_full", fixed = TRUE)
})

test_that("all eleven retained figure copies match their approved sources", {
  audit <- audit_retained_figure_copies(project_path())
  expect_identical(nrow(audit), 11L)
  expect_true(all(audit$source_matches_expected))
  expect_true(all(audit$copy_matches_source))
})

test_that("figure provenance exposes source and duplicate relationships", {
  path <- project_path("output", "manifests", "figure_provenance.csv")
  skip_if_not(file.exists(path))
  manifest <- readr::read_csv(path, show_col_types = FALSE)
  expect_true(all(c(
    "source_code_available", "duplicate_stem",
    "duplicate_format_relationships", "proposed_status"
  ) %in% names(manifest)))
  approved <- readr::read_csv(
    project_path("config", "retained_figures.csv"), show_col_types = FALSE
  )
  original_rows <- manifest$figure_path %in% approved$figure_path
  expect_identical(sum(original_rows), 11L)
  expect_true(all(
    manifest$proposed_status[original_rows] == "KEEP_FOR_MANUSCRIPT"
  ))
})

test_that("regeneratable upstream figures use canonical candidate paths", {
  manifest <- figure_manifest(PROJECT_ROOT)
  canonical_names <- c(
    "compare_rel.abund_vs.sqrt.rel.abund.pdf",
    "destination_climate_warming_gradient.pdf",
    "repeated_species_consistency.pdf",
    "temporal_category_flows.pdf",
    "temporal_transition_matrices.pdf",
    "transient_winning_losing_strength.pdf"
  )
  canonical <- manifest[basename(manifest$filename) %in% canonical_names, ]
  expect_setequal(basename(canonical$filename), canonical_names)
  expect_true(all(grepl("/figures/manuscript/", canonical$filename, fixed = TRUE)))
  expect_true(all(
    canonical$status %in% c(
      "CURRENT_CANONICAL_ANALYSIS_FIGURE",
      "implemented target; not generated"
    )
  ))
})
