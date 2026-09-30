source(project_path("R", "figures", "species_context.R"))
source(project_path("R", "infrastructure", "provenance.R"))

test_that("Figure 6 uses direct climate-only versus full support", {
  key <- predictor_synthesis_key()
  fit <- dplyr::bind_rows(
    tibble::tibble(
      predictor = key$source_predictor,
      model = "climate only", method = "ML", AIC = 100
    ),
    tibble::tibble(
      predictor = key$source_predictor,
      model = "full predictor-climate", method = "ML",
      AIC = 100 + c(8.2879, -4.0575, -3.4711, -6.4543, 10.3332, 11.7257)
    )
  )
  comparison <- tibble::tibble(
    predictor = key$source_predictor,
    comparison = "climate only vs full predictor-climate",
    likelihood_ratio = 1,
    df = 1,
    p_value = c(0.4621, 0.0101, 0.0125, 0.00414, 0.6845, 0.2485),
    valid_nested_comparison = TRUE,
    same_rows_and_sampling_variances = TRUE,
    method = "ML"
  )
  heterogeneity <- tibble::tibble(
    Predictor = key$source_predictor,
    `Overall heterogeneity resolved (%)` = 1,
    `Warming heterogeneity resolved (%)` = 2,
    `Context heterogeneity resolved (%)` = 3,
    `Total heterogeneity resolved (%)` = 4
  )
  fit_file <- withr::local_tempfile(fileext = ".csv")
  comparison_file <- withr::local_tempfile(fileext = ".csv")
  heterogeneity_file <- withr::local_tempfile(fileext = ".csv")
  readr::write_csv(fit, fit_file)
  readr::write_csv(comparison, comparison_file)
  readr::write_csv(heterogeneity, heterogeneity_file)

  prepared <- prepare_predictor_synthesis(
    fit_file, comparison_file, heterogeneity_file
  )
  expect_equal(
    prepared$delta_AIC_full_minus_climate,
    c(8.2879, -4.0575, -3.4711, -6.4543, 10.3332, 11.7257)
  )
  expect_equal(prepared$lrt_p_value, comparison$p_value)
  expect_false(any(grepl("additive|M[0-3]", names(prepared))))
  figure <- build_predictor_synthesis_figure(prepared)
  expect_s3_class(figure, "ggplot")
})

test_that("functional coverage is computed independently by group", {
  data <- tibble::tibble(
    species_id = c("a", "b", "c", "d"),
    experiment_id = c("e1", "e2", "e1", "e3"),
    T_z = c(-1, 0, 2, 3), P_z = c(-2, 0, 1, 4), W_z = c(-1, 1, -2, 2),
    Woodiness = factor(c("non_woody", "non_woody", "woody", "woody"))
  )
  coverage <- summarize_functional_climate_coverage(data, "Woodiness", "Woodiness")
  expect_setequal(as.character(coverage$group), c("non_woody", "woody"))
  expect_equal(coverage$T_z_min[coverage$group == "woody"], 2)
  expect_equal(coverage$T_z_max[coverage$group == "non_woody"], 0)
})

test_that("the Woodiness figure fails closed without a group-specific mask", {
  function_text <- paste(
    deparse(body(build_woodiness_response_surface_figure)), collapse = "\n"
  )
  expect_match(function_text, "group-specific observed climate support", fixed = TRUE)
  expect_match(function_text, "predictor_formulas()$Woodiness_full", fixed = TRUE)
})

test_that("all currently retained figures match their approved checksums", {
  approved <- readr::read_csv(
    project_path("config", "retained_figures.csv"), show_col_types = FALSE
  )
  audit <- audit_retained_figure_copies(project_path())
  expect_identical(nrow(audit), nrow(approved))
  expect_true(all(audit$source_matches_expected))
  expect_true(all(audit$copy_matches_source))
  release <- retained_figure_release_check(audit)
  expect_identical(nrow(release), nrow(approved))
  expect_true(all(release$pass))
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
  retained_rows <- manifest$figure_path %in% approved$figure_path
  expect_identical(sum(retained_rows), nrow(approved))
  expect_true(all(
    manifest$proposed_status[retained_rows] == "KEEP_FOR_MANUSCRIPT"
  ))
})

test_that("regeneratable upstream figures use canonical manuscript paths", {
  manifest <- figure_manifest(PROJECT_ROOT)
  canonical_names <- c(
    "compare_rel.abund_vs.sqrt.rel.abund.pdf",
    "appendix_destination_climate_warming_gradient.pdf",
    "appendix_destination_climate_warming_gradient.png",
    "fig02_repeated_species_consistency.pdf",
    "fig02_repeated_species_consistency.png",
    "fig03_temporal_category_flows.pdf",
    "fig03_temporal_category_flows.png",
    "fig04_species_minmax_classification.pdf",
    "fig04_species_minmax_classification.png",
    "fig05_climate_only_fixed_effect_prediction_surfaces.pdf",
    "fig05_climate_only_fixed_effect_prediction_surfaces.png",
    "fig06_predictor_synthesis_support_and_heterogeneity.pdf",
    "fig06_predictor_synthesis_support_and_heterogeneity.png",
    "fig07_predictor_response_atlas.pdf",
    "fig07_predictor_response_atlas.png",
    "appendix_family_minmax_classification.pdf",
    "appendix_family_minmax_classification.png",
    "appendix_model_predictive_performance.pdf",
    "appendix_model_predictive_performance.png",
    "appendix_pev2_pev3_prediction_atlas.pdf",
    "appendix_pev2_pev3_prediction_atlas.png",
    "appendix_phylogeny_pcoa_summary.pdf",
    "appendix_phylogeny_pcoa_summary.png",
    "appendix_temporal_winning_losing_strength_correlations.pdf",
    "appendix_temporal_winning_losing_strength_correlations.png",
    "appendix_temporal_transition_matrices.pdf",
    "appendix_temporal_transition_matrices.png"
  )
  canonical <- manifest[basename(manifest$filename) %in% canonical_names, ]
  expect_setequal(basename(canonical$filename), canonical_names)
  expect_true(all(grepl("/figures/manuscript/", canonical$filename, fixed = TRUE)))
  expect_true(all(
    canonical$status %in% c(
      "CURRENT_ANALYSIS_FIGURE",
      "FINAL_MANUSCRIPT",
      "implemented target; not generated"
    )
  ))
})

test_that("Figure 4 uses enlarged typography and manuscript labels", {
  expect_gte(formals(build_species_context_callout_figure)$classification_base_size, 14)
  expect_gte(formals(build_species_context_callout_figure)$classification_label_size, 4.8)
  expect_gte(formals(build_species_context_callout_figure)$inset_base_size, 14)
  body_text <- paste(
    deparse(body(build_species_context_callout_figure)), collapse = "\n"
  )
  expect_true(grepl(
    "classification_base_size[[:space:]]*[+][[:space:]]*2", body_text
  ))
  expect_true(grepl(
    "inset_base_size[[:space:]]*[+][[:space:]]*2", body_text
  ))
  expect_false(grepl("Delta T =", body_text, fixed = TRUE))
  expect_match(body_text, "Minimum fitted winning-losing strength", fixed = TRUE)
  expect_match(body_text, "Maximum fitted winning-losing strength", fixed = TRUE)
  expect_match(
    body_text,
    "Warming sensitivity\\n(shift with stronger warming)",
    fixed = TRUE
  )
  expect_match(body_text, "Destination precipitation (mm)", fixed = TRUE)
  expect_match(body_text, "Fitted winning/losing strength", fixed = TRUE)
})

test_that("Figure 5 climate-only renderer is fixed-effect-only and support-masked", {
  renderer <- readLines(
    project_path(
      "R", "figures", "render_climate_response.R"
    ),
    warn = FALSE
  )
  renderer_text <- paste(renderer, collapse = "\n")
  expect_false(grepl("targets::tar_make(", renderer_text, fixed = TRUE))
  expect_false(grepl("rma.mv(", renderer_text, fixed = TRUE))
  expect_match(renderer_text, "Fixed-effect predictions", fixed = TRUE)
  expect_match(renderer_text, "Canonical elevation labels", fixed = TRUE)
  expect_match(renderer_text, "Support-mask enforcement", fixed = TRUE)
  expect_match(renderer_text, "All experiments annotated", fixed = TRUE)
})

test_that("Figure 7 renderer is figure-only and support-masked", {
  renderer <- readLines(
    project_path("R", "figures", "render_predictor_atlas.R"),
    warn = FALSE
  )
  renderer_text <- paste(renderer, collapse = "\n")
  expect_false(grepl("targets::tar_make(", renderer_text, fixed = TRUE))
  expect_false(grepl("rma.mv(", renderer_text, fixed = TRUE))
  expect_match(renderer_text, "Support-mask enforcement", fixed = TRUE)
  expect_match(renderer_text, "group-specific observed three-dimensional climate hull", fixed = TRUE)
  expect_match(renderer_text, "Predicted winning-losing strength", fixed = TRUE)
  expect_match(renderer_text, "Destination temperature (°C)", fixed = TRUE)
  expect_match(renderer_text, "Destination precipitation (mm)", fixed = TRUE)
})

test_that("the destination-climate renderer drops inherited grouping", {
  source(project_path("R", "figures", "climate.R"))
  body_text <- paste(
    deparse(body(plot_destination_climate_warming_gradient)), collapse = "\n"
  )
  expect_match(body_text, "dplyr::ungroup(data)", fixed = TRUE)
})

test_that("the destination-climate appendix uses continuous axes without cutoffs", {
  body_text <- paste(
    deparse(body(plot_destination_climate_warming_gradient_continuous)),
    collapse = "\n"
  )
  expect_match(body_text, "dplyr::ungroup(data)", fixed = TRUE)
  expect_false(grepl("geom_rect", body_text, fixed = TRUE))
  expect_false(grepl("geom_vline", body_text, fixed = TRUE))
  expect_false(grepl("geom_hline", body_text, fixed = TRUE))
  expect_false(grepl("temperature_context", body_text, fixed = TRUE))
  expect_false(grepl("precipitation_context", body_text, fixed = TRUE))
  expect_match(body_text, 'legend.position = "bottom"', fixed = TRUE)
})
