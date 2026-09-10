# Executable audit of target ancestry at the scientific/reference boundary.

target_ancestors <- function(target, edges) {
  visited <- character()
  frontier <- target
  while (length(frontier)) {
    upstream <- unique(edges$from[edges$to %in% frontier])
    upstream <- setdiff(upstream, visited)
    if (!length(upstream)) break
    visited <- c(visited, upstream)
    frontier <- upstream
  }
  unique(visited)
}

audit_canonical_fixture_dependencies <- function(root = NULL) {
  if (is.null(root)) {
    root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    repeat {
      if (file.exists(file.path(root, "_targets.R"))) break
      parent <- dirname(root)
      if (identical(parent, root)) {
        stop("Could not locate _targets.R from the current directory.", call. = FALSE)
      }
      root <- parent
    }
  }
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  previous_directory <- setwd(root)
  on.exit(setwd(previous_directory), add = TRUE)
  network <- targets::tar_network(
    targets_only = TRUE,
    outdated = FALSE,
    callr_function = NULL,
    script = file.path(root, "_targets.R"),
    store = file.path(root, "_targets")
  )
  edges <- network$edges
  scientific <- c(
    "downstream_auc_data",
    "mboth_526_refit",
    "species_classification_526",
    "family_model_bundle_526",
    "family_classification_526",
    "ldmc_model_sequence_526",
    "lma_model_sequence_526",
    "height_model_sequence_526",
    "graminoid_model_sequence_526",
    "woodiness_model_sequence_526",
    "resident_phylogeny_candidate",
    "resident_pcoa_candidate",
    "resident_pev_scores_oriented",
    "standardized_pev_species_526",
    "phylogeny_predictor_specification_526",
    "phylogeny_model_sequence_526",
    "heterogeneity_bundle_526",
    "heterogeneity_components",
    "heterogeneity_synthesis",
    "canonical_reference_results",
    "current_downstream_model_report"
  )
  validation <- c(
    "historical_phylogeny_model_sequence_reference_526",
    "historical_heterogeneity_bundle_reference_526",
    "phylogeny_migration_equivalence",
    "phylogeny_migration_equivalence_file",
    "phylogeny_migration_equivalence_report_file",
    "step22_reference_validation"
  )
  targets_to_check <- c(scientific, validation)
  forbidden_transition <- c(
    "transition_reference_fixture_file", "transition_reference"
  )
  forbidden_step22 <- c(
    "step22_reference_fixture_file", "step22_reference_results"
  )
  historical_step19 <- c(
    "pcoa_eigen_summary_reference", "resident_phylogeny_reference",
    "resident_phylogeny_newick_reference",
    "resident_phylogeny_newick_reference_restored",
    "phylogeny_species_order_reference", "resident_pev_scores_reference",
    "phylogenetic_distance_reference", "zero_distance_audit_reference"
  )
  result <- dplyr::bind_rows(lapply(targets_to_check, function(target) {
    direct <- sort(unique(edges$from[edges$to == target]))
    ancestors <- target_ancestors(target, edges)
    transition_ancestor <- any(ancestors %in% forbidden_transition)
    step22_ancestor <- any(ancestors %in% forbidden_step22)
    historical_step19_ancestor <- any(ancestors %in% historical_step19)
    role <- if (target %in% scientific) "CANONICAL_SCIENTIFIC_COMPUTATION" else
      "REFERENCE_VALIDATION"
    expected_fixture_ancestry <- role == "REFERENCE_VALIDATION"
    observed_fixture_ancestry <- transition_ancestor || step22_ancestor ||
      historical_step19_ancestor
    tibble::tibble(
      target = target,
      role = role,
      direct_upstream_targets = paste(direct, collapse = "; "),
      transition_reference_ancestor = transition_ancestor,
      step22_reference_ancestor = step22_ancestor,
      historical_step19_reference_ancestor = historical_step19_ancestor,
      fixture_ancestor = observed_fixture_ancestry,
      pass = if (expected_fixture_ancestry) {
        observed_fixture_ancestry
      } else {
        !observed_fixture_ancestry
      },
      detail = if (role == "CANONICAL_SCIENTIFIC_COMPUTATION") {
        "Frozen fixtures must not occur anywhere in the upstream ancestry."
      } else {
        "Fixture ancestry is expected and restricted to this validation branch."
      }
    )
  }))
  if (any(!result$pass)) {
    stop(
      "Canonical/reference dependency boundary failed for: ",
      paste(result$target[!result$pass], collapse = ", "),
      call. = FALSE
    )
  }
  result
}
