# Settled scientific choices for the canonical 526-species analysis.

canonical_pipeline_source_files <- function(root = ".") {
  file.path(
    root, "R", "functions",
    c(
      "upstream_utils.R", "upstream_community.R", "upstream_climate.R",
      "upstream_traits.R", "upstream_winner_loser.R", "upstream_temporal.R"
    )
  )
}

canonical_scientific_source_files <- function(root = ".") {
  fixed <- file.path(
    root,
    c(
      "_targets.R", ".gitignore", "README.md", "renv.lock",
      "TransPlant_winners_losers.Rproj", "external/get_chelsa.R",
      "data/README.md", "data/source_manifest.csv"
    )
  )
  inventory_trees <- c(
    "R/functions" = "[.]R$",
    "R/targets" = "[.]R$",
    "config" = "[.](R|csv)$",
    "reports" = "[.]qmd$",
    "scripts" = "[.]R$",
    "tests" = "[.](R|rds|csv|md)$"
  )
  trees <- unlist(Map(
    function(directory, pattern) {
      path <- file.path(root, directory)
      if (!dir.exists(path)) return(character())
      files <- list.files(
        path, pattern = pattern, recursive = TRUE, full.names = TRUE,
        all.files = TRUE, no.. = TRUE
      )
      files[file.exists(files) & !dir.exists(files)]
    }, names(inventory_trees), unname(inventory_trees)
  ), use.names = FALSE)
  sort(unique(c(fixed[file.exists(fixed)], canonical_pipeline_source_files(root), trees)))
}

canonical_downstream_config <- function(root = here::here()) {
  list(
    step = 22L,
    response = "maximum-duration raw relative-cover Delta AUC (warmed minus origin)",
    exclusions = resident_downstream_model_exclusions(),
    expected = list(effects = 1604L, species = 526L, experiments = 42L),
    climate_formula = stats::as.formula("~ T_z * P_z * W_z"),
    species_random = list(
      stats::as.formula("~ 1 | species_id"),
      stats::as.formula("~ 0 + W_z | species_id"),
      stats::as.formula("~ 0 + T_z + P_z + TP_z | species_id")
    ),
    random_structure = "GEN",
    classification = list(
      minimum_experiments = 4L,
      required_context_rank = 4L,
      winner = "minimum fitted Delta AUC > 0",
      loser = "maximum fitted Delta AUC < 0",
      otherwise = "Context-dependent"
    ),
    trait_scaling = "natural log of unique-species median, standardized before row expansion",
    phylogeny_axes = c("PEV1", "PEV2", "PEV3"),
    pev_orientation = c(
      PEV1 = "Poa trivialis positive",
      PEV2 = "Erigeron formosissimus negative",
      PEV3 = "Euphrasia minima positive"
    ),
    reference_path = file.path(
      root, "output", "step_by_step_modeling", "downstream_526",
      "downstream_526_reference_results.rds"
    ),
    fixture_path = file.path(root, "tests", "fixtures", "step22_reference_results.rds"),
    canonical_output = list(
      data = file.path(root, "output", "data"),
      models = file.path(root, "output", "models"),
      tables = file.path(root, "output", "tables"),
      audits = file.path(root, "output", "audits"),
      manifests = file.path(root, "output", "manifests")
    )
  )
}
