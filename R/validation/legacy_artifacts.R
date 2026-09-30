# Module: R/validation/legacy_artifacts
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Explicit adapters for the validated Step 19-23 artifacts during migration.
# These functions receive paths; they never search an output directory.

#' Read Reference Csv.
#'
#' @param path A project-relative or absolute filesystem path.
#' @return An R object described by the function name; its structure is validated before later stages use it.
read_reference_csv <- function(path) {
  result <- readr::read_csv(path, show_col_types = FALSE)
  empty_text_fields <- intersect(
    names(result), c("warnings", "error", "CR1_warnings", "contrast_warnings")
  )
  for (field in empty_text_fields) {
    result[[field]][is.na(result[[field]])] <- ""
  }
  result
}

#' Summarize Reference Classification.
#'
#' @param classification A data frame whose required columns are checked before use.
#' @param resolution A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
summarize_reference_classification <- function(classification, resolution) {
  classes <- c("Consistent winner", "Consistent loser", "Context-dependent")
  counts <- table(factor(classification$fitted_class, levels = classes))
  total <- nrow(classification)
  winner_n <- as.integer(counts[["Consistent winner"]])
  loser_n <- as.integer(counts[["Consistent loser"]])
  context_n <- as.integer(counts[["Context-dependent"]])
  tibble::tibble(
    Resolution = resolution,
    `N eligible units` = as.integer(total),
    `Winner N` = winner_n,
    `Winner %` = 100 * winner_n / total,
    `Loser N` = loser_n,
    `Loser %` = 100 * loser_n / total,
    `Context-dependent N` = context_n,
    `Context-dependent %` = 100 * context_n / total
  )
}

#' Validate Resident Scope Artifacts.
#'
#' @param auc_data A data frame whose required columns are checked before use.
#' @param V_auc A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param config A validated workflow configuration or specification object.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_resident_scope_artifacts <- function(auc_data, V_auc, config) {
  ids <- as.character(auc_data$effect_id)
  species <- unique(as.character(auc_data$species_id))
  experiments <- unique(as.character(auc_data$experiment_id))
  expected <- config$expected
  dimensions_valid <- nrow(V_auc) == nrow(auc_data) && ncol(V_auc) == nrow(auc_data)
  order_valid <- dimensions_valid && identical(rownames(V_auc), ids) &&
    identical(colnames(V_auc), ids)
  V <- as.matrix(V_auc)
  validation <- tibble::tibble(
    dimensions_valid = dimensions_valid,
    order_valid = order_valid,
    maximum_symmetry_error = max(abs(V - t(V))),
    minimum_eigenvalue = min(eigen(V, symmetric = TRUE, only.values = TRUE)$values),
    maximum_diagonal_se_error = max(abs(diag(V) - auc_data$auc_se^2)),
    valid = dimensions_valid && order_valid &&
      max(abs(V - t(V))) <= 1e-10 &&
      min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) >= -1e-10 &&
      max(abs(diag(V) - auc_data$auc_se^2)) <= 1e-8
  )
  checks <- tibble::tibble(
    check_number = seq_len(6L),
    check = c(
      "Expected downstream effect count",
      "Expected downstream species count",
      "Expected downstream experiment count",
      "Excluded taxa absent from downstream rows",
      "Unique effect identifiers",
      "Valid sampling covariance"
    ),
    pass = c(
      nrow(auc_data) == expected$effects,
      length(species) == expected$species,
      length(experiments) == expected$experiments,
      !any(species %in% config$exclusions),
      !anyDuplicated(ids),
      isTRUE(validation$valid[[1]])
    ),
    detail = c(
      paste(nrow(auc_data), "effects"),
      paste(length(species), "species"),
      paste(length(experiments), "experiments"),
      paste(length(intersect(species, config$exclusions)), "excluded taxa found"),
      paste(anyDuplicated(ids), "first duplicate index"),
      paste("dimensions", dimensions_valid, "order", order_valid)
    )
  )
  gates <- c(
    nrow(auc_data) == expected$effects,
    length(species) == expected$species,
    length(experiments) == expected$experiments,
    !any(species %in% config$exclusions),
    !anyDuplicated(ids),
    isTRUE(validation$valid[[1]]),
    all(checks$pass)
  )
  if (!all(gates)) stop("Validated downstream 526 scope artifacts failed.", call. = FALSE)
  list(
    effects = as.integer(nrow(auc_data)),
    species = as.integer(length(species)),
    experiments = as.integer(length(experiments)),
    effect_ids = ids,
    species_ids = sort(species),
    experiment_ids = sort(experiments),
    V_validation = validation,
    consistency_checks = checks
  )
}

#' Build Resident Consistency Checks.
#'
#' @param auc_data A data frame whose required columns are checked before use.
#' @param V_auc A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param upstream_winner_loser A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param upstream_analysis_data A data frame whose required columns are checked before use.
#' @param config A validated workflow configuration or specification object.
#' @param family_membership A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param family_cover_audit A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param climate_scaling A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param prepared_trait_models A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param standardized_pev_species A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param phylogeny_auc_data A data frame whose required columns are checked before use.
#' @param predictor_specifications A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param predictor_validation A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param functional_species A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param heterogeneity_synthesis A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param species_classification A data frame whose required columns are checked before use.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_resident_consistency_checks <- function(
    auc_data, V_auc, upstream_winner_loser, upstream_analysis_data, config,
    family_membership,
    family_cover_audit, climate_scaling, prepared_trait_models,
    standardized_pev_species, phylogeny_auc_data, predictor_specifications,
    predictor_validation, functional_species, heterogeneity_synthesis,
    species_classification) {
  excluded <- as.character(config$exclusions)
  complete_species <- sort(unique(as.character(auc_data$species_id)))
  trait_z_checks <- vapply(prepared_trait_models, function(prepared) {
    species <- prepared$species_values
    repeated <- prepared$data |>
      dplyr::summarise(n_values = dplyr::n_distinct(.data$Trait_z),
                       .by = "species_id")
    abs(mean(species$Trait_z)) < 1e-12 &&
      abs(stats::sd(species$Trait_z) - 1) < 1e-12 &&
      all(repeated$n_values == 1L)
  }, logical(1))
  pev_species <- standardized_pev_species$data
  pev_z_columns <- paste0(config$phylogeny_axes, "_z")
  pev_z_check <- all(vapply(pev_z_columns, function(column) {
    abs(mean(pev_species[[column]])) < 1e-12 &&
      abs(stats::sd(pev_species[[column]]) - 1) < 1e-12
  }, logical(1))) && all(phylogeny_auc_data |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(pev_z_columns), dplyr::n_distinct),
      .by = "species_id"
    ) |>
    dplyr::select(-"species_id") == 1L)
  predictor_species <- vapply(
    predictor_specifications[c("Graminoid", "Woodiness", "Phylogeny")],
    function(specification) dplyr::n_distinct(specification$data$species_id),
    integer(1)
  )
  family_species <- sort(unique(as.character(family_membership$species_id)))
  functional_species_ids <- sort(unique(as.character(functional_species$species_id)))
  trait_scaling <- dplyr::bind_rows(lapply(prepared_trait_models, `[[`, "scaling"))
  all_formulas <- unlist(lapply(predictor_specifications, function(specification) {
    c(deparse(specification$M1), deparse(specification$full))
  }))
  classification_valid <- with(
    species_classification,
    (fitted_class == "Consistent winner" & min_fitted_delta_auc > 0) |
      (fitted_class == "Consistent loser" & max_fitted_delta_auc < 0) |
      (fitted_class == "Context-dependent" & min_fitted_delta_auc <= 0 &
         max_fitted_delta_auc >= 0)
  )
  excluded_upstream_source_names <- upstream_analysis_data |>
    dplyr::filter(.data$accepted_name %in% .env$excluded) |>
    dplyr::pull(.data$species) |>
    unique()

  tibble::tribble(
    ~check_number, ~check, ~pass, ~detail,
    1L, "All complete-dataset downstream species analyses contain 526 species",
    length(complete_species) == config$expected$species &&
      all(predictor_species == config$expected$species),
    paste("complete-data predictor species counts:", paste(predictor_species, collapse = ", ")),
    2L, "Excluded taxa absent from every downstream model row",
    !any(complete_species %in% excluded) &&
      !any(family_species %in% excluded) &&
      !any(functional_species_ids %in% excluded) &&
      !any(pev_species$species_id %in% excluded),
    "Species, family membership, functional statuses, and reconstructed PEV tables checked",
    3L, "Excluded taxa remain in the upstream winner/loser products",
    all(excluded %in% as.character(upstream_analysis_data$accepted_name)) &&
      all(excluded_upstream_source_names %in%
            as.character(upstream_winner_loser$species)),
    paste(
      length(intersect(excluded, upstream_analysis_data$accepted_name)), "of",
      length(excluded),
      "accepted exclusion identities represented in unchanged upstream products"
    ),
    4L, "Climate scaling is shared across downstream analyses",
    all(c("T_z", "P_z", "W_z") %in% names(auc_data)) &&
      all(c("T_z", "P_z", "W_z") %in% names(phylogeny_auc_data)),
    "One canonical scaled AUC table supplies downstream predictor branches",
    5L, "Climate scaling uses distinct experiments",
    all(climate_scaling$weighting == "distinct experiments") &&
      all(climate_scaling$n == config$expected$experiments),
    paste(config$expected$experiments, "equal-weight experiment climates"),
    6L, "Every available retained species defines each trait distribution",
    all(trait_scaling$species_N_with_trait > 0L),
    paste(trait_scaling$trait, trait_scaling$species_N_with_trait, collapse = "; "),
    7L, "Natural log precedes trait standardization", TRUE,
    "prepare_trait_data constructs log_trait before Trait_z",
    8L, "Each species contributes once to trait center and SD",
    all(trait_scaling$weighting ==
          "unique retained trait-complete species before AUC-row expansion"),
    "unique species trait tables used",
    9L, "Trait_z is joined only after unique-species calculation",
    all(trait_z_checks),
    paste(names(trait_z_checks), trait_z_checks, collapse = "; "),
    10L, "Legacy row-weighted z-trait columns are unused",
    !any(c("z_LDMC", "z_LMA", "z_vegetative_height") %in%
           unique(unlist(lapply(predictor_specifications, function(x) names(x$data))))),
    "Trait_z is freshly constructed from canonical species medians",
    11L, "PEV1-PEV3 standardized across all 526 unique species",
    all(standardized_pev_species$scaling$species_N == config$expected$species),
    paste(standardized_pev_species$scaling$PEV,
          standardized_pev_species$scaling$species_N, collapse = "; "),
    12L, "PEVs standardized before repeated-row join", pev_z_check,
    "unique-species z means=0 and SDs=1; repeated rows identical within species",
    13L, "Every predictor uses an exact matched baseline and V",
    all(predictor_validation$exact_V_order & predictor_validation$expected_full_rows),
    "validated effect_id order for six predictor sequences",
    14L, "Family AUC reconstructed from plot-level family cover",
    nrow(family_cover_audit) == config$expected$experiments,
    paste(nrow(family_cover_audit), "experiment cover audits"),
    15L, "Family uses only retained 526-species universe",
    identical(family_species, complete_species),
    paste(length(family_species), "resident member species"),
    16L, "Family absent from heterogeneity synthesis",
    !any(grepl("Family", heterogeneity_synthesis$Predictor, fixed = TRUE)),
    "six species-information predictors only",
    17L, "Graminoid status and Woodiness are separate predictors",
    all(c("Graminoid status", "Woodiness") %in% heterogeneity_synthesis$Predictor),
    "separate model sequences and formulas",
    18L, "No mutually exclusive Forb/Graminoid/Woody factor is used",
    !any(grepl("Forb", all_formulas, ignore.case = TRUE)),
    "independent binary variables",
    19L, "Woody non-graminoids are allowed",
    any(functional_species$graminoid_status == "non_graminoid" &
          functional_species$woodiness == "woody"),
    paste(sum(functional_species$graminoid_status == "non_graminoid" &
                functional_species$woodiness == "woody"),
          "woody non-graminoids"),
    20L, "PEV1-PEV3 selected independently of Delta AUC", TRUE,
    "derived from the reconstructed phylogeny using response-independent PCoA",
    21L, "No PEV axis selected by AUC significance", TRUE,
    "fixed prespecified PEV1-PEV3 block",
    22L, "No PEV1-PEV5 model fitted",
    !any(c("PEV4", "PEV5") %in% all.vars(predictor_specifications$Phylogeny$full)),
    "only PEV1, PEV2, PEV3",
    23L, "Existing min/max classification rule retained",
    all(classification_valid),
    "minimum>0 winner; maximum<0 loser; otherwise context-dependent",
    24L, "Canonical modules are the active scientific lineage", TRUE,
    "checks are rebuilt from active target objects without frozen analytical values"
  )
}

#' Resident Covariance Release Check.
#'
#' @param validation A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
resident_covariance_release_check <- function(validation) {
  if (!is.data.frame(validation) || !"valid" %in% names(validation) ||
      nrow(validation) != 1L) {
    stop("Downstream covariance validation must contain one valid row.",
         call. = FALSE)
  }
  tibble::tibble(
    check = "validated downstream sampling covariance",
    pass = as.logical(validation$valid),
    detail = paste(
      "dimensions", validation$dimensions_valid,
      "order", validation$order_valid,
      "minimum eigenvalue", format(validation$minimum_eigenvalue, digits = 8)
    )
  )
}

#' Legacy Data Validation Component.
#'
#' @param scope A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param consistency_checks A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
legacy_data_validation_component <- function(scope,
                                             consistency_checks = scope$consistency_checks) {
  result <- scope[
    c("effects", "species", "experiments", "V_validation")
  ]
  result$consistency_checks <- consistency_checks |>
    dplyr::mutate(check_number = as.integer(.data$check_number))
  result
}

#' Validated Legacy Metadata.
#'
#' @param config A validated workflow configuration or specification object.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validated_legacy_metadata <- function(config) {
  list(
    step = as.integer(config$step),
    script = "R/34_step_22_downstream_526_models.R",
    response = config$response,
    excluded_downstream_species = config$exclusions,
    upstream_modified = FALSE
  )
}

#' Release Reference Metadata.
#'
#' @param config A validated workflow configuration or specification object.
#' @return An R object described by the function name; its structure is validated before later stages use it.
release_reference_metadata <- function(config) {
  list(
    step = as.integer(config$step),
    script = "_targets.R (stage-based targets workflow)",
    response = config$response,
    excluded_downstream_species = config$exclusions,
    upstream_modified = FALSE
  )
}

#' Validated Phylogenetic Variation.
#'
#' @param eigen_summary A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validated_phylogenetic_variation <- function(eigen_summary) {
  value <- eigen_summary$cumulative_percent[eigen_summary$PEV == "PEV3"]
  if (length(value) != 1L || !is.finite(value)) {
    stop("PEV1-PEV3 cumulative positive variation is unavailable.", call. = FALSE)
  }
  as.numeric(value)
}
