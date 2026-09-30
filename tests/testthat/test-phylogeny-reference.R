source(project_path("R", "models", "phylogeny.R"))
source(project_path("R", "infrastructure", "provenance.R"))

test_that("the reviewed phylogeny placement ledger is frozen and complete", {
  path <- project_path("config", "phylogeny_placement_decisions.csv")
  expect_identical(
    sha256_file(path),
    "adc146dffede99200b95b4e3bb43a8d316ed49edbaec25e00316a7d4b81e4e63"
  )
  decisions <- readr::read_csv(path, show_col_types = FALSE)
  expect_identical(nrow(decisions), 526L)
  expect_identical(dplyr::n_distinct(decisions$species_id), 526L)
  expect_false(any(decisions$match_type == "unresolved"))
})

test_that("the installed phylogeny package is the pinned Git revision", {
  provenance <- phylogeny_package_provenance()
  checks <- validate_phylogeny_package_provenance(provenance)
  expect_true(all(checks$pass))
})

test_that("the external-source ledger records the pinned phylogeny commit", {
  sources <- readr::read_csv(
    project_path("config", "external_sources.csv"), show_col_types = FALSE
  )
  phylogeny <- dplyr::filter(sources, .data$component == "V.PhyloMaker2")
  expect_equal(nrow(phylogeny), 1L)
  expect_identical(
    phylogeny$repository[[1]],
    "https://github.com/jinyizju/V.PhyloMaker2"
  )
  expect_identical(
    phylogeny$git_sha[[1]],
    "7af3fb5152f691af2e4ec9d5e2e467d1b50505e9"
  )
})

test_that("the current PEV orientation follows the documented convention", {
  path <- project_path(
    "output", "data", "phylogeny_predictor_scores.csv"
  )
  expect_true(file.exists(path), info = "Build pev_scores_file")
  scores <- readr::read_csv(path, show_col_types = FALSE)
  oriented <- orient_pev_axes(scores)
  expect_gt(oriented$PEV1[oriented$species_id == "Poa trivialis"], 0)
  expect_lt(oriented$PEV2[oriented$species_id == "Erigeron formosissimus"], 0)
  expect_gt(oriented$PEV3[oriented$species_id == "Euphrasia minima"], 0)
  expect_true(all(compare_pev_coordinates(oriented, scores)$pass))
})

test_that("the current zero-distance pairs are intentional", {
  distance_path <- project_path(
    "output", "data", "phylogeny_distances.rds"
  )
  audit_path <- project_path(
    "output", "audits", "zero_distance_pair_audit.csv"
  )
  expect_true(
    all(file.exists(c(distance_path, audit_path))),
    info = "Build canonical phylogenetic-distance audit files"
  )
  checks <- validate_zero_distance_placements(
    readRDS(distance_path),
    readr::read_csv(audit_path, show_col_types = FALSE)
  )
  expect_true(all(checks$pass))
})

test_that("the PEV modeling interface validates while retaining reconstructed values", {
  reference <- tibble::tibble(
    species_id = c(
      "Poa trivialis", "Erigeron formosissimus", "Euphrasia minima"
    ),
    PEV1 = c(1, 0, -1),
    PEV2 = c(2, -2, 0),
    PEV3 = c(-2, -1, 3)
  )
  recomputed <- reference
  recomputed$PEV1 <- recomputed$PEV1 + c(1e-14, 0, -1e-14)
  interface <- build_legacy_pev_validation_interface(recomputed, reference)
  expect_identical(interface$data, recomputed)
  expect_false(interface$audit$historical_values_used_for_modeling)
  expect_true(all(interface$audit$pass))

  drifted <- recomputed
  drifted$PEV1[[1]] <- drifted$PEV1[[1]] + 1e-4
  expect_error(
    build_legacy_pev_validation_interface(drifted, reference),
    "PEV"
  )
})
