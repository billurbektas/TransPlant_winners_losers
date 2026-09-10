test_that("destination climate uses the five reviewed groups and review gate", {
  fixture <- tibble::tibble(
    experiment_id = paste0("e", 1:6),
    destination_T_JJA = c(12, 14, 18, 14, 14, 18),
    destination_P_JJA = c(400, 400, 400, 600, 900, 600),
    climate_summary_status = "complete"
  )
  result <- reviewed_five_climate_contexts(fixture)
  expect_equal(
    as.character(result$assigned$climate_context),
    c("cold_dry", "intermediate_dry", "warm_dry", "cold_intermediate", "cold_wet")
  )
  expect_equal(result$review$experiment_id, "e6")
  expect_error(
    validate_five_climate_contexts(result),
    "climate_context_review.csv"
  )
  expect_equal(
    validate_five_climate_contexts(reviewed_five_climate_contexts(fixture[1:5, ])),
    reviewed_five_climate_contexts(fixture[1:5, ])$assigned
  )
})

test_that("current climate snapshot has reviewed assignments and named moves", {
  data <- readr::read_csv(
    file.path(PROJECT_ROOT, "output", "tables", "analysis_dataset.csv"),
    show_col_types = FALSE
  )
  result <- reviewed_five_climate_contexts(data)
  expect_equal(nrow(result$review), 0L)
  expect_equal(
    dplyr::count(result$assigned, .data$climate_context)$n,
    c(15L, 9L, 11L, 2L, 5L)
  )
  named <- result$assigned |>
    dplyr::filter(.data$experiment_id %in% c(
      "DE_Grainau::HIGH::LOW", "DE_Susalps::EB::GW", "CH_Calanda::Pea::Cal"
    )) |>
    dplyr::select("experiment_id", "climate_context")
  expect_equal(
    as.character(named$climate_context[match(
      c("DE_Grainau::HIGH::LOW", "DE_Susalps::EB::GW", "CH_Calanda::Pea::Cal"),
      named$experiment_id
    )]),
    c("warm_dry", "warm_dry", "cold_dry")
  )
})

test_that("repeated-species consistency retains and orders significant identities only", {
  fixture <- tibble::tibble(
    species = rep(paste0("raw ", LETTERS[1:5]), each = 4),
    accepted_name = rep(paste0("Species ", LETTERS[1:5]), each = 4),
    experiment_id = rep(paste0("e", 1:4), times = 5),
    category = c(
      "loser", "loser", "loser", "non_significant_loser",
      "loser", "winner", "non_significant_loser", "non_significant_winner",
      "winner", "winner", "non_significant_winner", "no_change",
      "colonizer", "colonizer", "non_significant_colonizer", "no_change",
      rep("no_change", 4)
    )
  )
  result <- signed_species_consistency(fixture)
  expect_equal(
    as.character(result$consistency_group),
    c("Loser only", "Experiment-dependent responses", "Winner only", "Colonizer only")
  )
  expect_equal(result$species, paste0("Species ", LETTERS[1:4]))
  expect_false("Species E" %in% result$species)
  expect_true(all(result$significant_total > 0L))
  expect_equal(result$significant_total, c(3L, 2L, 2L, 2L))

  residents <- signed_species_consistency(fixture, include_colonizers = FALSE)
  expect_equal(
    as.character(residents$consistency_group),
    c("Loser only", "Experiment-dependent responses", "Winner only")
  )
  expect_equal(residents$species, paste0("Species ", LETTERS[1:3]))
  expect_true(all(residents$colonizer == 0L))
  expect_true(all(residents$non_significant_colonizer == 0L))
  expect_false(any(as.character(residents$consistency_group) == "Colonizer only"))
})

test_that("family consistency uses unique species-experiment combinations and percentages", {
  fixture <- tibble::tibble(
    species = c("raw a", "raw a", "raw b", "raw c", "raw d", "raw e"),
    accepted_name = c("A a", "A a", "B b", "C c", "D d", "E e"),
    experiment_id = c("e1", "e1", "e1", "e2", "e2", "e3"),
    family = c("Family one", "Family one", "Family one", "Family one", "Family two", "Family three"),
    category = c("loser", "loser", "winner", "non_significant_winner", "winner", "no_change")
  )
  result <- signed_family_consistency(fixture)
  family_one <- dplyr::filter(result, .data$family == "Family one")
  expect_equal(family_one$total_combinations, 3L)
  expect_equal(family_one$n_species, 3L)
  expect_equal(family_one$n_experiments, 2L)
  expect_equal(family_one$loser_percent, 100 / 3)
  expect_equal(family_one$winner_percent, 100 / 3)
  expect_equal(family_one$neutral_percent, 100 / 3)
  expect_equal(
    family_one$colonizer_high - family_one$loser_low,
    100
  )
  expect_equal(as.character(family_one$consistency_group), "Mixed significant responses")
  expect_false("Family three" %in% result$family)
  percentage_plot <- family_consistency_plot_data(result, metric = "percentage") |>
    dplyr::summarise(total_height = sum(abs(.data$value)), .by = "family")
  count_plot <- family_consistency_plot_data(result, metric = "count") |>
    dplyr::summarise(total_height = sum(abs(.data$value)), .by = "family") |>
    dplyr::left_join(
      dplyr::select(result, "family", "total_combinations"), by = "family"
    )
  expect_equal(percentage_plot$total_height, rep(100, nrow(percentage_plot)))
  expect_equal(count_plot$total_height, count_plot$total_combinations)
  expect_equal(
    family_consistency_plot_data(result, metric = "count") |>
      dplyr::distinct(.data$family, .data$consistency_group) |>
      dplyr::count(.data$family) |>
      dplyr::pull(.data$n),
    rep(1L, nrow(result))
  )
})

test_that("clean experiment climate map keeps one finite record per experiment", {
  fixture <- tibble::tibble(
    experiment_id = rep(c("Region::origin::destination", "Region::local::warm"), each = 2),
    destination_T_JJA = rep(c(10, 12), each = 2),
    destination_P_JJA = rep(c(300, 450), each = 2),
    diff_T_JJA = rep(c(2, 3), each = 2),
    diff_P_JJA = rep(c(-50, 25), each = 2),
    climate_summary_status = rep(c("complete", "partial"), each = 2),
    species = letters[1:4]
  )
  result <- destination_climate_experiment_map_data(fixture)
  expect_equal(nrow(result), 2L)
  expect_equal(
    result$experiment_label,
    c("Region\norigin -> destination", "Region\nlocal -> warm*")
  )
  invalid <- fixture
  invalid$destination_T_JJA[[2]] <- 11
  expect_error(
    destination_climate_experiment_map_data(invalid),
    "one unique climate-map record"
  )
})

test_that("climate density data retain raw significant-response counts", {
  fixture <- tibble::tibble(
    experiment_id = c(rep("e1", 4), rep("e2", 3), rep("e3", 2)),
    category = c(
      "loser", "loser", "loser", "winner",
      "winner", "winner", "colonizer",
      "colonizer", "loser"
    ),
    destination_T_JJA = c(rep(10, 4), rep(12, 3), rep(15, 2)),
    destination_P_JJA = c(rep(300, 4), rep(450, 3), rep(700, 2))
  )
  result <- category_density_data_xy(
    fixture, "destination_T_JJA", "destination_P_JJA", grid_size = 20L
  )
  expect_false("proportion" %in% names(result$counts))
  expect_equal(
    result$counts$n[
      result$counts$experiment_id == "e1" & result$counts$category == "loser"
    ],
    3L
  )
  expect_equal(
    result$counts |>
      dplyr::summarise(total = sum(.data$n), .by = "experiment_id") |>
      dplyr::arrange(.data$experiment_id) |>
      dplyr::pull(.data$total),
    c(4L, 3L, 2L)
  )
  expect_true(all(is.finite(result$density$density)))
})
