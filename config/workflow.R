# Module: config/workflow
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Active workflow configuration and source inventory.

#' Active Function Files.
#'
#' @param root A project-relative or absolute filesystem path.
#' @return An R object described by the function name; its structure is validated before later stages use it.
active_function_files <- function(root = ".") {
  stage_directories <- c(
    "infrastructure", "data", "responses", "models", "evaluation",
    "figures", "validation"
  )
  sort(unlist(lapply(stage_directories, function(directory) {
    path <- file.path(root, "R", directory)
    list.files(path, pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
  }), use.names = FALSE))
}

#' Active Scientific Source Files.
#'
#' @param root A project-relative or absolute filesystem path.
#' @return An R object described by the function name; its structure is validated before later stages use it.
active_scientific_source_files <- function(root = ".") {
  fixed <- file.path(
    root,
    c(
      "_targets.R", ".gitignore", "README.md", "renv.lock",
      "TransPlant_winners_losers.Rproj", "external/get_chelsa.R",
      "data/README.md", "data/source_manifest.csv"
    )
  )
  inventory_trees <- c(
    "R" = "[.]R$",
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
  sort(unique(c(fixed[file.exists(fixed)], active_function_files(root), trees)))
}

#' Resident Analysis Spec.
#'
#' @param root A project-relative or absolute filesystem path.
#' @return An R object described by the function name; its structure is validated before later stages use it.
resident_analysis_spec <- function(root = here::here()) {
  list(
    response = "maximum-duration raw relative-cover Delta AUC (warmed minus origin)",
    exclusions = resident_model_exclusions(),
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
    output = list(
      data = file.path(root, "output", "data"),
      models = file.path(root, "output", "models"),
      tables = file.path(root, "output", "tables"),
      audits = file.path(root, "output", "audits"),
      manifests = file.path(root, "output", "manifests")
    )
  )
}
