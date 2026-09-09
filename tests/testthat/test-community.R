test_that("current adapter reproduces the reviewed snapshot", {
  input_path <- project_path("data", "TransPlantNetwork_060824.RData")
  skip_if_not(file.exists(input_path))
  input <- read_current_transplant_input(input_path)
  metadata <- site_metadata_from_regions(input$regional)
  decisions <- readr::read_csv(
    project_path("config", "species_mapping_decisions.csv"),
    col_types = readr::cols(.default = readr::col_character())
  )
  mapping <- resolve_local_species_mappings(input$cleanedspecies, decisions)
  community <- canonical_community(input$dat, mapping, metadata)
  sites <- build_sites(community, input$regional)
  experiments <- build_experiments(community, sites)

  expect_equal(nrow(sites), 59L)
  expect_equal(nrow(experiments), 44L)
  expect_equal(sum(experiments$included), 42L)
  expect_setequal(experiments$exclusion_reason[!experiments$included],
                  c("insufficient duration", "warmed observations have missing Year"))
  expect_true(all(experiments$has_origin_control))
  expect_true(all(experiments$has_destination_control))
  expect_false(any(grepl("::(19|20)[0-9]{2}$", community$stable_plot_id)))
})

test_that("the nine ambiguous botanical keys are explicitly reviewed", {
  input_path <- project_path("data", "TransPlantNetwork_060824.RData")
  skip_if_not(file.exists(input_path))
  input <- read_current_transplant_input(input_path)
  decisions <- readr::read_csv(project_path("config", "species_mapping_decisions.csv"), show_col_types = FALSE)
  review <- local_species_mapping_review(input$cleanedspecies, decisions)
  expect_equal(nrow(dplyr::distinct(review, region, site_id, local_species_code)), 9L)
  expect_false(any(is.na(review$decision)))
})

test_that("explicit taxonomy exclusions do not enter analysis community", {
  community <- tibble::tibble(
    species = c("Accepted species", "Carex", "Carex sp.", "Asteraceae", "Region_NID_1"),
    cover = 1
  )
  taxonomy <- tibble::tibble(
    name_original = c("Accepted species", "Carex", "Carex sp.", "Asteraceae", "Region_NID_1"),
    accepted_name_reviewed = c("Accepted species", "Carex", "Carex", "Asteraceae", NA_character_)
  )
  result <- filter_community_to_reviewed_taxa(community, taxonomy)
  expect_equal(result$species, "Accepted species")
})

test_that("taxonomic resolution separates species-, genus-, and family-level names", {
  taxonomy <- tibble::tibble(
    name_original = c("Carex nigra", "Carex", "Carex sp.", "Asteraceae"),
    accepted_name_reviewed = c("Carex nigra", "Carex", "Carex", "Asteraceae")
  )
  result <- taxonomic_resolution_review(taxonomy)
  expect_equal(
    result$identification_level,
    c("species_or_lower", "genus_only", "genus_only", "family_only")
  )
  expect_equal(result$include_species_level, c(TRUE, FALSE, FALSE, FALSE))
})
