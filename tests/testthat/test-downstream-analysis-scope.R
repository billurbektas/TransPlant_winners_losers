testthat::test_that("declared taxa remain upstream but are excluded downstream", {
  results_path <- project_path("output", "tables", "winner_loser_results.csv")
  testthat::skip_if_not(file.exists(results_path))
  exclusions <- resident_downstream_model_exclusions()
  testthat::expect_equal(length(exclusions), 6L)

  upstream <- readr::read_csv(results_path, show_col_types = FALSE) |>
    dplyr::mutate(horizon = as.character(.data$horizon)) |>
    dplyr::filter(
      .data$transform == "raw", .data$contrast == "to_origin",
      .data$pool == "resident", .data$horizon == "max"
    )
  taxonomy <- readr::read_csv(
    project_path("output", "tables", "analysis_dataset.csv"),
    show_col_types = FALSE
  ) |>
    dplyr::distinct(
      experiment_id = as.character(.data$experiment_id),
      source_species = as.character(.data$species),
      species_id = as.character(.data$accepted_name)
    )
  upstream <- upstream |>
    dplyr::transmute(
      experiment_id = as.character(.data$experiment_id),
      source_species = as.character(.data$species)
    ) |>
    dplyr::left_join(taxonomy, by = c("experiment_id", "source_species"))
  testthat::expect_setequal(
    intersect(unique(upstream$species_id), exclusions), exclusions
  )
  testthat::expect_equal(sum(upstream$species_id %in% exclusions), 21L)

  downstream <- readRDS(project_path(
    "tests", "fixtures", "canonical_transition_reference.rds"
  ))$downstream_auc_data
  testthat::expect_equal(nrow(downstream), 1604L)
  testthat::expect_equal(dplyr::n_distinct(downstream$species_id), 526L)
  testthat::expect_false(any(downstream$species_id %in% exclusions))
})

testthat::test_that("graminoid and woody status are independent binary fields", {
  analysis_path <- project_path("output", "tables", "analysis_dataset.csv")
  override_path <- project_path("config", "functional_group_overrides.csv")
  testthat::skip_if_not(file.exists(analysis_path) && file.exists(override_path))
  taxonomy <- readr::read_csv(analysis_path, show_col_types = FALSE) |>
    dplyr::distinct(.data$accepted_name, .data$family)
  overrides <- readr::read_csv(
    override_path, show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )
  downstream_overrides <- readr::read_csv(
    project_path("config", "downstream_species_overrides.csv"),
    show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )
  statuses <- derive_functional_groups(
    taxonomy, overrides, downstream_overrides
  )
  testthat::expect_setequal(
    unique(statuses$graminoid_status), c("graminoid", "non_graminoid")
  )
  testthat::expect_setequal(
    unique(statuses$woody_status), c("woody", "non_woody")
  )
  testthat::expect_true(all(
    statuses$graminoid_status[statuses$family %in%
      c("Poaceae", "Cyperaceae", "Juncaceae")] == "graminoid"
  ))
  testthat::expect_equal(
    statuses$woody_status[statuses$accepted_name == "Salix herbacea"],
    "woody"
  )
})

testthat::test_that("functional-information models test both statuses separately", {
  formulas <- canonical_predictor_formulas()
  testthat::expect_equal(
    formulas$Graminoid_full,
    stats::as.formula("~ GraminoidStatus * T_z * P_z * W_z")
  )
  testthat::expect_equal(
    formulas$Woodiness_full,
    stats::as.formula("~ Woodiness * T_z * P_z * W_z")
  )
  all_formulas <- vapply(formulas, function(x) paste(deparse(x), collapse = " "),
                         character(1))
  testthat::expect_false(any(grepl("FunctionalGroup", all_formulas, fixed = TRUE)))
})

testthat::test_that("updated phylogeny exactly represents the downstream set", {
  fixture <- readRDS(project_path(
    "tests", "fixtures", "canonical_transition_reference.rds"
  ))$phylogeny
  tree <- fixture$tree_rds
  scores <- fixture$pev_scores
  species_order <- fixture$species_order
  zero_distance <- fixture$zero_distance_audit
  placement <- readr::read_csv(
    project_path("config", "phylogeny_placement_decisions.csv"),
    show_col_types = FALSE
  )

  testthat::expect_equal(length(tree$tip.label), 526L)
  testthat::expect_setequal(
    tree$tip.label,
    gsub("_", " ", species_order$newick_tip_label, fixed = TRUE)
  )
  testthat::expect_false(any(species_order$species_id %in%
                               resident_downstream_model_exclusions()))
  testthat::expect_equal(
    sum(placement$match_type ==
          "infraspecific -> parent species placement"),
    10L
  )
  testthat::expect_equal(nrow(scores), 526L)
  testthat::expect_true(all(is.finite(as.matrix(scores[c("PEV1", "PEV2", "PEV3")]))))
  testthat::expect_equal(nrow(zero_distance), 5L)
  testthat::expect_true(all(
    zero_distance$classification == "intentional shared placement"
  ))
  testthat::expect_gt(scores$PEV1[scores$species_id == "Poa trivialis"], 0)
  testthat::expect_lt(
    scores$PEV2[scores$species_id == "Erigeron formosissimus"], 0
  )
  testthat::expect_gt(scores$PEV3[scores$species_id == "Euphrasia minima"], 0)
})
