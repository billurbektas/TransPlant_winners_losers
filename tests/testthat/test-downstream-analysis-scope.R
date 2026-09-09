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

  downstream <- read_resident_auc_data(PROJECT_ROOT, horizons = "max")
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
  statuses <- derive_functional_groups(taxonomy, overrides)
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
  script <- paste(readLines(
    project_path("R", "29_step_17_species_information_heterogeneity.R"),
    warn = FALSE
  ), collapse = "\n")
  testthat::expect_match(
    script, "~ GraminoidStatus * T_z * P_z * W_z", fixed = TRUE
  )
  testthat::expect_match(
    script, "~ WoodyStatus * T_z * P_z * W_z", fixed = TRUE
  )
  testthat::expect_false(grepl(
    "~ FunctionalGroup * T_z * P_z * W_z", script, fixed = TRUE
  ))
})

testthat::test_that("updated phylogeny exactly represents the downstream set", {
  output_dir <- project_path(
    "output", "step_by_step_modeling", "resident_phylogeny_pcoa"
  )
  required <- file.path(
    output_dir,
    c("resident_phylogeny.rds", "tnrs_backbone_crosswalk.csv",
      "pcoa_axis_extreme_species.csv", "species_pev_scores.csv",
      "infraspecific_placement_audit.csv",
      "phylogeny_pcoa_audit_report.txt")
  )
  testthat::skip_if_not(all(file.exists(required)))
  tree <- readRDS(required[[1]])
  crosswalk <- readr::read_csv(required[[2]], show_col_types = FALSE)
  extremes <- readr::read_csv(required[[3]], show_col_types = FALSE)
  scores <- readr::read_csv(required[[4]], show_col_types = FALSE)
  infraspecific_audit <- readr::read_csv(required[[5]], show_col_types = FALSE)
  report <- paste(readLines(required[[6]], warn = FALSE), collapse = "\n")

  testthat::expect_equal(length(tree$tip.label), 526L)
  testthat::expect_false(any(
    tree$tip.label %in% resident_downstream_model_exclusions()
  ))
  testthat::expect_equal(
    sum(crosswalk$match_type ==
          "infraspecific -> parent species placement"),
    10L
  )
  testthat::expect_equal(nrow(extremes), 30L)
  testthat::expect_equal(
    unname(as.integer(table(extremes$PEV))), rep(10L, 3L)
  )
  testthat::expect_equal(nrow(infraspecific_audit), 16L)
  testthat::expect_true(
    "Parent species tip present in GBOTB.extended.LCVP?" %in%
      names(infraspecific_audit)
  )
  testthat::expect_match(
    report,
    paste0(
      "Absence from GBOTB.extended.LCVP indicates absence from the ",
      "phylogenetic backbone and does not imply that the taxon is unrecognized"
    ),
    fixed = TRUE
  )

  shared <- crosswalk |>
    dplyr::filter(.data$shared_parent_placement) |>
    dplyr::group_by(.data$`backbone matched name`) |>
    dplyr::filter(dplyr::n() > 1L) |>
    dplyr::ungroup()
  if (nrow(shared)) {
    coordinate_columns <- grep("^PEV", names(scores), value = TRUE)
    joined <- shared |>
      dplyr::select("species_id", parent = "backbone matched name") |>
      dplyr::left_join(scores, by = "species_id")
    maximum_within_parent_difference <- joined |>
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(coordinate_columns),
          ~ max(.x) - min(.x)
        ),
        .by = "parent"
      ) |>
      dplyr::select(-"parent") |>
      unlist(use.names = FALSE) |>
      max()
    testthat::expect_lt(maximum_within_parent_difference, 1e-10)
  }
})
