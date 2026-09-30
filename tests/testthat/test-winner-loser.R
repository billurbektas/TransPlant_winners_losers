test_that("experiment pools are fixed from the complete experimental duration", {
  fixture <- tibble::tribble(
    ~ODT, ~species, ~relative_year, ~relative_cover,
    "origin_controls", "resident", 1, 0.2,
    "destination_controls", "colonizer", 1, 0.2,
    "warmed", "resident", 1, 0.2,
    "warmed", "colonizer", 1, 0.2,
    "warmed", "outsider", 1, 0.2,
    "warmed", "late", 5, 0.2
  )
  pool <- experiment_pool(fixture)
  expect_equal(pool$pool[pool$species == "resident"], "resident")
  expect_equal(pool$pool[pool$species == "colonizer"], "colonizer")
  expect_equal(pool$pool[pool$species == "outsider"], "outsider")
  expect_equal(pool$pool[pool$species == "late"], "outsider")
})

test_that("validation rejects changing pools across horizons", {
  classified <- tibble::tibble(
    experiment_id = "experiment", horizon = c("2", "4"), transform = "raw",
    species = "species", pool = c("outsider", "resident"),
    category = factor(c("non_significant_outsider", "non_significant_winner")),
    significant = FALSE, estimate = 0.1
  )
  pools <- dplyr::select(classified, experiment_id, horizon, species, pool)
  expect_error(
    validate_winner_loser_outputs(classified, pools),
    "Pool membership changes across horizons"
  )
})

test_that("validation rejects resident categories for nonresident pools", {
  classified <- tibble::tibble(
    experiment_id = "experiment", horizon = "2", transform = "raw",
    species = "species", pool = "outsider", category = factor("winner"),
    significant = TRUE, estimate = 0.1
  )
  pools <- dplyr::select(classified, experiment_id, horizon, species, pool)
  expect_error(
    validate_winner_loser_outputs(classified, pools),
    "classified as a resident response"
  )
})

test_that("legacy observation filter removes singleton species-treatment groups", {
  fixture <- tibble::tibble(
    ODT = c("warmed", "warmed", "origin_controls"),
    species = c("replicated", "replicated", "singleton"),
    relative_year = c(1, 2, 1),
    relative_cover = c(0.1, 0.2, 0.3)
  )
  filtered <- filter_minimum_plot_observations(fixture)
  expect_equal(unique(filtered$species), "replicated")
})

test_that("colonizers and outsiders are split by significance", {
  contrasts <- tibble::tibble(
    experiment_id = "experiment", horizon = "max", transform = "raw",
    contrast = "to_origin", species = c("sig", "nonsig", "sig_outsider", "nonsig_outsider"),
    estimate = c(0.2, 0.1, 0.3, 0.1), p.value = c(0.01, 0.2, 0.02, 0.4)
  )
  pools <- tibble::tibble(
    experiment_id = "experiment", horizon = "max",
    species = c("sig", "nonsig", "sig_outsider", "nonsig_outsider"),
    pool = c("colonizer", "colonizer", "outsider", "outsider")
  )
  result <- classify_winners_losers(contrasts, pools)
  expect_equal(
    as.character(result$category),
    c("colonizer", "non_significant_colonizer", "colonizer", "non_significant_outsider")
  )
  expect_equal(
    as.character(result$analysis_category),
    c("colonizer", "non_significant_colonizer", "colonizer", "non_significant_colonizer")
  )
})

test_that("non-significant residents retain contrast direction", {
  contrasts <- tibble::tibble(
    experiment_id = "experiment", horizon = "max", transform = "raw",
    contrast = "to_origin", species = c("negative", "positive"),
    estimate = c(-0.2, 0.2), p.value = c(0.3, 0.3)
  )
  pools <- tibble::tibble(
    experiment_id = "experiment", horizon = "max",
    species = c("negative", "positive"), pool = "resident"
  )
  result <- classify_winners_losers(contrasts, pools)
  expect_equal(
    as.character(result$category),
    c("non_significant_loser", "non_significant_winner")
  )
})

test_that("flow summaries conserve all transitions", {
  transitions <- tibble::tibble(
    later_horizon = c("4", "4", "6"), state_2 = c("non_significant_winner", "loser", "non_significant_winner"),
    state_later = c("non_significant_winner", "non_significant_loser", "winner"), same_state = c(TRUE, FALSE, FALSE)
  )
  summary <- transition_summary(transitions)
  expect_true(all(validate_flow_conservation(transitions, summary$counts)$pass))
})

test_that("pairwise temporal flows require both horizons", {
  classified <- tibble::tribble(
    ~transform, ~experiment_id, ~species, ~horizon, ~category,
    "raw", "experiment", "outsider_then_resident", "2", "non_significant_outsider",
    "raw", "experiment", "outsider_then_resident", "4", "non_significant_loser",
    "raw", "experiment", "later_only", "4", "winner",
    "raw", "experiment", "matched_resident", "2", "loser",
    "raw", "experiment", "matched_resident", "4", "loser"
  )
  result <- temporal_transitions(classified)
  horizon_4 <- dplyr::filter(result, .data$later_horizon == "4")
  expect_true("non_significant_outsider" %in% horizon_4$state_2)
  expect_false("later_only" %in% horizon_4$species)
  expect_false(any(horizon_4$state_2 == "not_yet_detected"))
})

test_that("linked temporal flow retains the first-two-year cohort", {
  classified <- tibble::tribble(
    ~transform, ~experiment_id, ~species, ~horizon, ~category,
    "raw", "experiment", "complete_outsider", "2", "non_significant_outsider",
    "raw", "experiment", "complete_outsider", "4", "non_significant_outsider",
    "raw", "experiment", "complete_outsider", "6", "colonizer",
    "raw", "experiment", "missing_four", "2", "non_significant_winner",
    "raw", "experiment", "missing_four", "6", "non_significant_winner"
  )
  result <- temporal_linked_transitions(classified)
  expect_equal(sort(result$species), c("complete_outsider", "missing_four"))
  expect_equal(result$state_2[result$species == "complete_outsider"], "non_significant_outsider")
  expect_equal(result$state_4[result$species == "missing_four"], "follow_up_unavailable")
  expect_equal(result$state_6[result$species == "missing_four"], "non_significant_winner")
})

test_that("linked temporal flow rejects changing pools", {
  classified <- tibble::tribble(
    ~transform, ~experiment_id, ~species, ~horizon, ~category, ~pool,
    "raw", "experiment", "species", "2", "non_significant_outsider", "outsider",
    "raw", "experiment", "species", "4", "non_significant_winner", "resident"
  )
  expect_error(
    temporal_linked_transitions(classified),
    "Temporal results contain pool changes"
  )
})

test_that("alluvial pair percentages exactly match transition-matrix rows", {
  linked <- tibble::tribble(
    ~experiment_id, ~species, ~state_2, ~state_4, ~state_6, ~pool_2, ~pool_4, ~pool_6,
    "e1", "a", "non_significant_colonizer", "non_significant_colonizer", "non_significant_colonizer", "outsider", "outsider", "outsider",
    "e1", "b", "non_significant_colonizer", "non_significant_colonizer", "colonizer", "outsider", "outsider", "outsider",
    "e1", "c", "non_significant_colonizer", "colonizer", "colonizer", "outsider", "outsider", "outsider",
    "e1", "d", "non_significant_colonizer", "follow_up_unavailable", "follow_up_unavailable", "outsider", "outsider", "outsider",
    "e1", "e", "non_significant_winner", "winner", "winner", "resident", "resident", "resident",
    "e1", "f", "non_significant_winner", "non_significant_winner", "winner", "resident", "resident", "resident"
  )
  flow <- temporal_pairwise_flow_summary(linked) |>
    dplyr::filter(.data$visible) |>
    dplyr::transmute(
      comparison = .data$comparison,
      from_category = .data$source_state,
      to_category = .data$destination_state,
      flow_n = .data$n,
      flow_proportion = .data$row_proportion
    )
  matrix <- build_temporal_transition_matrix(linked) |>
    dplyr::filter(
      .data$comparison %in% c("First 2 → First 4", "First 4 → First 6"),
      .data$n > 0
    ) |>
    dplyr::transmute(
      comparison = .data$comparison,
      from_category = .data$from_category,
      to_category = .data$to_category,
      matrix_n = .data$n, matrix_proportion = .data$row_proportion
    )
  comparison <- dplyr::full_join(
    flow, matrix, by = c("comparison", "from_category", "to_category")
  )
  expect_false(anyNA(comparison))
  expect_equal(comparison$flow_n, comparison$matrix_n)
  expect_equal(comparison$flow_proportion, comparison$matrix_proportion)
})

test_that("temporal analysis uses the merged analysis category", {
  classified <- tibble::tribble(
    ~transform, ~experiment_id, ~species, ~horizon, ~category, ~analysis_category, ~pool,
    "raw", "experiment", "species", "2", "non_significant_outsider", "non_significant_colonizer", "outsider",
    "raw", "experiment", "species", "4", "non_significant_outsider", "non_significant_colonizer", "outsider",
    "raw", "experiment", "species", "6", "colonizer", "colonizer", "outsider"
  )
  result <- temporal_linked_transitions(classified)
  expect_equal(result$state_2, "non_significant_colonizer")
  expect_equal(result$state_4, "non_significant_colonizer")
  expect_equal(result$state_6, "colonizer")
  expect_equal(result$pool_2, "outsider")
})

test_that("temporal composition preserves no-change and reports yearly sampling", {
  classified <- tibble::tribble(
    ~transform, ~experiment_id, ~species, ~horizon, ~category, ~analysis_category,
    "raw", "e1", "a", "2", "winner", "winner",
    "raw", "e1", "b", "2", "no_change", "no_change",
    "raw", "e1", "a", "4", "winner", "winner",
    "raw", "e1", "b", "4", "non_significant_loser", "non_significant_loser",
    "raw", "e1", "a", "6", "non_significant_winner", "non_significant_winner",
    "raw", "e1", "b", "6", NA_character_, NA_character_
  )
  result <- temporal_category_composition(classified)
  expect_true(all(c("2", "4", "6") %in% result$horizon))
  expect_true(all(result$category[result$n == 0] %in% unique(result$category)))
  expect_equal(
    result |>
      dplyr::summarise(total = sum(.data$proportion), .by = "horizon") |>
      dplyr::pull(.data$total),
    rep(1, 3)
  )
  expect_equal(
    result$n[result$horizon == "6" & result$category == "no_change"],
    1L
  )
  expect_equal(unique(result$n_species[result$horizon == "2"]), 2L)
  expect_equal(unique(result$n_experiments[result$horizon == "2"]), 1L)
  expect_equal(unique(result$total_species_experiment[result$horizon == "2"]), 2L)
})

test_that("category persistence reports exact same-state denominators", {
  linked <- tibble::tribble(
    ~experiment_id, ~species, ~state_2, ~state_4, ~state_6,
    "e1", "a", "winner", "winner", "winner",
    "e1", "b", "winner", "loser", "loser",
    "e2", "c", "loser", "follow_up_unavailable", "follow_up_unavailable"
  )
  result <- summarise_temporal_category_persistence(linked)
  winner_24 <- result$by_category |>
    dplyr::filter(.data$comparison == "First 2 → First 4", .data$from_category == "winner")
  expect_equal(winner_24$available_n, 2L)
  expect_equal(winner_24$same_n, 1L)
  expect_equal(winner_24$same_proportion, 0.5)
  all_three <- result$overall |>
    dplyr::filter(.data$comparison == "Same at First 2, 4, and 6")
  expect_equal(all_three$available_n, 2L)
  expect_equal(all_three$same_n, 1L)
})

test_that("transition matrix distinguishes significance and direction changes", {
  linked <- tibble::tribble(
    ~experiment_id, ~species, ~state_2, ~state_4, ~state_6, ~pool_2, ~pool_4, ~pool_6,
    "e1", "a", "non_significant_winner", "winner", "winner", "resident", "resident", "resident",
    "e1", "b", "non_significant_winner", "loser", "loser", "resident", "resident", "resident",
    "e1", "c", "colonizer", "non_significant_colonizer", "non_significant_colonizer", "outsider", "outsider", "outsider"
  )
  result <- build_temporal_transition_matrix(linked) |>
    dplyr::filter(.data$comparison == "First 2 → First 4", .data$n > 0)
  expect_equal(
    result$transition_type[result$to_category == "winner"],
    "significance only"
  )
  expect_equal(
    result$transition_type[result$to_category == "loser"],
    "loser-winner reversal"
  )
  expect_equal(
    result$transition_type[result$to_category == "non_significant_colonizer"],
    "significance only"
  )
})

test_that("maximum-duration AUC does not pass a NULL endpoint", {
  skip_if_not_installed("DescTools")
  fixture <- tibble::tibble(
    ODT = factor("warmed", levels = c("destination_controls", "origin_controls", "warmed")),
    stable_plot_id = "plot-1", species = "Species one", transform = "raw",
    relative_year = c(1, 2, 3), relative_cover = c(0, 1, 1)
  ) |>
    dplyr::select(-transform)
  result <- auc_by_plot(fixture, "max", "experiment")
  expect_equal(result$auc[result$transform == "raw"], 1.5)
})

test_that("contrast seeds are stable and group-specific", {
  first <- stable_contrast_seed("region::origin::destination", "2", "raw")
  expect_identical(first, stable_contrast_seed("region::origin::destination", "2", "raw"))
  expect_false(first == stable_contrast_seed("region::origin::destination", "3", "raw"))
  expect_true(first > 0L)
})

test_that("current predictor formulas retain the settled interaction structures", {
  formulas <- predictor_formulas()
  continuous_terms <- attr(stats::terms(formulas$continuous_full), "term.labels")
  phylogeny_terms <- attr(stats::terms(formulas$Phylogeny_full), "term.labels")

  expect_true("Trait_z:T_z:P_z:W_z" %in% continuous_terms)
  expect_true(all(c(
    "W_z:PEV1_z", "W_z:PEV2_z", "W_z:PEV3_z",
    "T_z:PEV1_z", "T_z:PEV2_z", "T_z:PEV3_z",
    "P_z:PEV1_z", "P_z:PEV2_z", "P_z:PEV3_z",
    "PEV1_z:TP_z", "PEV2_z:TP_z", "PEV3_z:TP_z"
  ) %in% phylogeny_terms))
})
