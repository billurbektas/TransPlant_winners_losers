# Step 20: final phylogeny/PCoA model-readiness audit -----------------------
#
# Audit only. This script reads the completed 526-species phylogeny/PCoA
# artifacts. It does not read AUC responses, rebuild or alter the tree, change
# taxonomy, rerun PCoA, standardize/select PEVs, or fit an ecological model.

source("config/analysis.R")

required_packages <- c("ape", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}

output_dir <- file.path(
  "output", "step_by_step_modeling", "resident_phylogeny_pcoa"
)
input_files <- c(
  tree = file.path(output_dir, "resident_phylogeny.rds"),
  distances = file.path(output_dir, "phylogenetic_distance_matrix.rds"),
  crosswalk = file.path(output_dir, "tnrs_backbone_crosswalk.csv"),
  species_order = file.path(output_dir, "phylogeny_species_order.csv"),
  scores = file.path(output_dir, "species_pev_scores.csv"),
  extremes = file.path(output_dir, "pcoa_axis_extreme_species.csv"),
  report = file.path(output_dir, "phylogeny_pcoa_audit_report.txt")
)
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs)) {
  stop("Missing completed phylogeny audit inputs: ",
       paste(missing_inputs, collapse = ", "), call. = FALSE)
}

output_files <- c(
  zero_pairs = file.path(output_dir, "zero_distance_pair_audit.csv"),
  family_grafts = file.path(output_dir, "family_graft_pev_audit.csv")
)

write_csv_checked <- function(x, path) {
  readr::write_csv(x, path, na = "")
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Failed to write audit table: ", path, call. = FALSE)
  }
  invisible(path)
}

tree <- readRDS(input_files[["tree"]])
phylogenetic_distance <- readRDS(input_files[["distances"]])
crosswalk <- readr::read_csv(input_files[["crosswalk"]], show_col_types = FALSE)
species_order <- readr::read_csv(
  input_files[["species_order"]], show_col_types = FALSE
)
scores <- readr::read_csv(input_files[["scores"]], show_col_types = FALSE)
extremes <- readr::read_csv(input_files[["extremes"]], show_col_types = FALSE)

if (!inherits(tree, "phylo")) {
  stop("resident_phylogeny.rds is not an ape phylo object.", call. = FALSE)
}
if (!is.matrix(phylogenetic_distance)) {
  phylogenetic_distance <- as.matrix(phylogenetic_distance)
}

# -------------------------------------------------------------------------
# 1. Every off-diagonal numerically zero patristic distance
# -------------------------------------------------------------------------

distance_tolerance <- 1e-10 * max(1, max(abs(phylogenetic_distance)))
zero_indices <- which(
  upper.tri(phylogenetic_distance) &
    abs(phylogenetic_distance) <= distance_tolerance,
  arr.ind = TRUE
)

crosswalk_lookup <- crosswalk |>
  dplyr::transmute(
    species_id = as.character(.data$species_id),
    tnrs_name = as.character(.data$`TNRS accepted name`),
    match_type = as.character(.data$match_type),
    parent_binomial = as.character(.data$parent_binomial),
    parent_tip_present = as.logical(.data$parent_in_backbone),
    shared_parent_placement = as.logical(.data$shared_parent_placement),
    infraspecific = as.logical(.data$infraspecific_name),
    backbone_matched_name = as.character(.data$`backbone matched name`)
  )

if (nrow(zero_indices)) {
  zero_pairs <- tibble::tibble(
    species_id_1 = rownames(phylogenetic_distance)[zero_indices[, "row"]],
    species_id_2 = colnames(phylogenetic_distance)[zero_indices[, "col"]],
    phylogenetic_distance = phylogenetic_distance[zero_indices]
  ) |>
    dplyr::left_join(
      crosswalk_lookup |>
        dplyr::rename_with(~ paste0(.x, "_1"), -"species_id") |>
        dplyr::rename(species_id_1 = "species_id"),
      by = "species_id_1"
    ) |>
    dplyr::left_join(
      crosswalk_lookup |>
        dplyr::rename_with(~ paste0(.x, "_2"), -"species_id") |>
        dplyr::rename(species_id_2 = "species_id"),
      by = "species_id_2"
    ) |>
    dplyr::mutate(
      intentional = .data$shared_parent_placement_1 &
        .data$shared_parent_placement_2 &
        !is.na(.data$backbone_matched_name_1) &
        .data$backbone_matched_name_1 == .data$backbone_matched_name_2 &
        (.data$infraspecific_1 | .data$infraspecific_2),
      reason = dplyr::if_else(
        .data$intentional,
        paste0(
          "Ecologically distinct taxa deliberately share parent-species ",
          "placement at ", .data$backbone_matched_name_1,
          "; the backbone does not resolve their infraspecific divergence."
        ),
        "No documented shared parent-species placement explains this zero distance."
      ),
      classification = dplyr::if_else(
        .data$intentional, "intentional shared placement",
        "unexpected zero distance"
      )
    ) |>
    dplyr::transmute(
      species_id_1 = .data$species_id_1,
      `TNRS accepted name 1` = .data$tnrs_name_1,
      species_id_2 = .data$species_id_2,
      `TNRS accepted name 2` = .data$tnrs_name_2,
      phylogenetic_distance = .data$phylogenetic_distance,
      `match_type species 1` = .data$match_type_1,
      `match_type species 2` = .data$match_type_2,
      `parent binomial 1` = .data$parent_binomial_1,
      `parent tip present 1` = .data$parent_tip_present_1,
      `parent binomial 2` = .data$parent_binomial_2,
      `parent tip present 2` = .data$parent_tip_present_2,
      `either taxon infraspecific` =
        .data$infraspecific_1 | .data$infraspecific_2,
      reason = .data$reason, classification = .data$classification
    ) |>
    dplyr::arrange(.data$species_id_1, .data$species_id_2)
} else {
  zero_pairs <- tibble::tibble(
    species_id_1 = character(), `TNRS accepted name 1` = character(),
    species_id_2 = character(), `TNRS accepted name 2` = character(),
    phylogenetic_distance = numeric(),
    `match_type species 1` = character(),
    `match_type species 2` = character(),
    `parent binomial 1` = character(), `parent tip present 1` = logical(),
    `parent binomial 2` = character(), `parent tip present 2` = logical(),
    `either taxon infraspecific` = logical(), reason = character(),
    classification = character()
  )
}
write_csv_checked(zero_pairs, output_files[["zero_pairs"]])

zero_pair_n <- nrow(zero_pairs)
intentional_pair_n <- sum(
  zero_pairs$classification == "intentional shared placement"
)
unexpected_pair_n <- sum(
  zero_pairs$classification == "unexpected zero distance"
)

# -------------------------------------------------------------------------
# 2. Every zero-length tree edge
# -------------------------------------------------------------------------

edge_tolerance <- 1e-12 * max(1, max(abs(tree$edge.length)))
zero_edge_rows <- which(abs(tree$edge.length) <= edge_tolerance)

descendant_tips <- function(node) {
  if (node <= ape::Ntip(tree)) return(tree$tip.label[[node]])
  children <- tree$edge[tree$edge[, 1] == node, 2]
  unlist(lapply(children, descendant_tips), use.names = FALSE)
}

if (length(zero_edge_rows)) {
  zero_edge_audit <- dplyr::bind_rows(lapply(zero_edge_rows, function(i) {
    child <- tree$edge[i, 2]
    descendants <- descendant_tips(child)
    placement <- crosswalk_lookup |>
      dplyr::filter(.data$species_id %in% .env$descendants)
    intentional <- nrow(placement) > 0L &&
      all(placement$shared_parent_placement)
    documented_graft <- !intentional && nrow(placement) > 0L &&
      all(placement$match_type %in%
            c("grafted within genus", "grafted within family"))
    tibble::tibble(
      edge_row = i,
      parent_node = tree$edge[i, 1],
      child_node = child,
      edge_length = tree$edge.length[[i]],
      descendant_species = paste(sort(descendants), collapse = " | "),
      classification = dplyr::case_when(
        intentional ~ "intentional shared placement",
        documented_graft ~ "documented V.PhyloMaker2 placement",
        TRUE ~ "unexplained"
      )
    )
  }))
} else {
  zero_edge_audit <- tibble::tibble(
    edge_row = integer(), parent_node = integer(), child_node = integer(),
    edge_length = numeric(), descendant_species = character(),
    classification = character()
  )
}

zero_edge_n <- nrow(zero_edge_audit)
intentional_edge_n <- sum(
  zero_edge_audit$classification == "intentional shared placement"
)
documented_graft_edge_n <- sum(
  zero_edge_audit$classification == "documented V.PhyloMaker2 placement"
)
unexplained_edge_n <- sum(zero_edge_audit$classification == "unexplained")

# -------------------------------------------------------------------------
# 3. Family-level graft influence on PEV1-PEV3
# -------------------------------------------------------------------------

score_axes <- c("PEV1", "PEV2", "PEV3")
missing_score_axes <- setdiff(score_axes, names(scores))
if (length(missing_score_axes)) {
  stop("Missing leading PEV columns: ", paste(missing_score_axes, collapse = ", "),
       call. = FALSE)
}

score_ranks <- scores |>
  dplyr::transmute(
    species_id = as.character(.data$species_id),
    PEV1 = as.numeric(.data$PEV1), PEV2 = as.numeric(.data$PEV2),
    PEV3 = as.numeric(.data$PEV3),
    abs_PEV1 = abs(.data$PEV1), abs_PEV2 = abs(.data$PEV2),
    abs_PEV3 = abs(.data$PEV3),
    rank_abs_PEV1 = rank(-.data$abs_PEV1, ties.method = "min"),
    rank_abs_PEV2 = rank(-.data$abs_PEV2, ties.method = "min"),
    rank_abs_PEV3 = rank(-.data$abs_PEV3, ties.method = "min")
  )

top10_membership <- extremes |>
  dplyr::filter(.data$PEV %in% .env$score_axes) |>
  dplyr::transmute(
    species_id = as.character(.data$species_id),
    PEV = as.character(.data$PEV), in_top10 = TRUE
  ) |>
  tidyr::pivot_wider(
    names_from = "PEV", values_from = "in_top10",
    names_prefix = "top10_", values_fill = FALSE
  )

family_graft_audit <- crosswalk |>
  dplyr::filter(.data$match_type == "grafted within family") |>
  dplyr::transmute(
    species_id = as.character(.data$species_id),
    `TNRS accepted name` = as.character(.data$`TNRS accepted name`),
    family = as.character(.data$`TNRS family`)
  ) |>
  dplyr::left_join(score_ranks, by = "species_id") |>
  dplyr::left_join(top10_membership, by = "species_id") |>
  dplyr::mutate(
    dplyr::across(
      dplyr::starts_with("top10_"), ~ dplyr::coalesce(.x, FALSE)
    )
  ) |>
  dplyr::arrange(.data$family, .data$species_id)
write_csv_checked(family_graft_audit, output_files[["family_grafts"]])

family_top10_counts <- vapply(score_axes, function(axis) {
  sum(family_graft_audit[[paste0("top10_", axis)]])
}, integer(1))
names(family_top10_counts) <- score_axes
family_graft_disproportionate <- any(family_top10_counts > 0L)

# -------------------------------------------------------------------------
# 4. Final PEV dataset integrity
# -------------------------------------------------------------------------

raw_score_header <- names(utils::read.csv(
  input_files[["scores"]], nrows = 0L, check.names = FALSE
))
raw_pev_header <- raw_score_header[grepl("^PEV[0-9]+$", raw_score_header)]
integrity_checks <- c(
  exactly_526_rows = nrow(scores) == 526L,
  exactly_526_unique_species = dplyr::n_distinct(scores$species_id) == 526L,
  no_duplicated_species_id = !anyDuplicated(scores$species_id),
  score_set_matches_tree = setequal(scores$species_id, tree$tip.label),
  score_set_matches_species_order = setequal(
    scores$species_id, species_order$species_id
  ),
  leading_pevs_complete_finite = all(vapply(
    scores[score_axes],
    function(x) !anyNA(x) && all(is.finite(x)), logical(1)
  )),
  no_duplicate_pev_columns = !anyDuplicated(raw_pev_header),
  no_downstream_excluded_taxa = !any(
    scores$species_id %in% resident_downstream_model_exclusions()
  )
)

ready <- all(integrity_checks) && unexpected_pair_n == 0L &&
  unexplained_edge_n == 0L

# -------------------------------------------------------------------------
# 5. Idempotently append Section 7 to the existing report
# -------------------------------------------------------------------------

report_lines <- readLines(input_files[["report"]], warn = FALSE)
section_start <- grep("^7[.] FINAL MODEL-READINESS AUDIT$", report_lines)
if (length(section_start)) {
  report_lines <- report_lines[seq_len(section_start[[1]] - 1L)]
}

failure_reasons <- c(
  if (unexpected_pair_n > 0L) {
    paste0(unexpected_pair_n, " unexpected zero-distance pair(s)")
  },
  if (unexplained_edge_n > 0L) {
    paste0(unexplained_edge_n, " unexplained zero-length edge(s)")
  },
  names(integrity_checks)[!integrity_checks]
)

final_section <- c(
  "",
  "7. FINAL MODEL-READINESS AUDIT",
  "------------------------------",
  paste0("Zero-distance ecological taxon pairs: ", zero_pair_n),
  paste0("Intentional shared placements: ", intentional_pair_n),
  paste0("Unexpected zero-distance pairs: ", unexpected_pair_n),
  "",
  paste0("Zero-length tree edges: ", zero_edge_n),
  paste0(
    "Explained by documented placements: ",
    intentional_edge_n + documented_graft_edge_n,
    " (intentional shared placement = ", intentional_edge_n,
    "; other documented V.PhyloMaker2 placements = ",
    documented_graft_edge_n, ")"
  ),
  paste0("Unexplained zero-length edges: ", unexplained_edge_n),
  "",
  paste0("Family-level grafts among PEV1 top 10: ",
         family_top10_counts[["PEV1"]]),
  paste0("Family-level grafts among PEV2 top 10: ",
         family_top10_counts[["PEV2"]]),
  paste0("Family-level grafts among PEV3 top 10: ",
         family_top10_counts[["PEV3"]]),
  paste0(
    "Any leading PEV disproportionately driven by family-level grafts: ",
    family_graft_disproportionate
  ),
  "",
  paste0(
    "PEV dataset contains exactly 526 unique resident taxa: ",
    integrity_checks[["exactly_526_rows"]] &&
      integrity_checks[["exactly_526_unique_species"]] &&
      integrity_checks[["no_duplicated_species_id"]]
  ),
  paste0(
    "PEV1-PEV3 complete and finite for all taxa: ",
    integrity_checks[["leading_pevs_complete_finite"]]
  ),
  paste0("PEV/tree/species-order sets agree: ",
         integrity_checks[["score_set_matches_tree"]] &&
           integrity_checks[["score_set_matches_species_order"]]),
  paste0("Accidental duplicate PEV columns absent: ",
         integrity_checks[["no_duplicate_pev_columns"]]),
  paste0("Downstream-excluded taxa absent: ",
         integrity_checks[["no_downstream_excluded_taxa"]]),
  "PEVs remain unstandardized; scaling is deferred to ecological modeling.",
  "",
  paste0(
    "Final phylogeny/PEV dataset ready for downstream ecological modeling: ",
    ready
  ),
  if (!ready) paste0("Reason(s): ", paste(failure_reasons, collapse = "; "))
)
writeLines(c(report_lines, final_section), input_files[["report"]], useBytes = TRUE)

cat("Final phylogeny/PCoA model-readiness audit\n")
cat("Zero-distance pairs:", zero_pair_n,
    "(intentional:", intentional_pair_n,
    "; unexpected:", unexpected_pair_n, ")\n")
cat("Zero-length edges:", zero_edge_n,
    "(explained:", intentional_edge_n + documented_graft_edge_n,
    "; unexplained:", unexplained_edge_n, ")\n")
cat("Family grafts in PEV1/2/3 top 10:",
    paste(family_top10_counts, collapse = "/"), "\n")
cat("Ready for downstream ecological modeling:", ready, "\n")

if (!ready) {
  stop(
    "Final phylogeny/PEV audit failed: ",
    paste(failure_reasons, collapse = "; "),
    ". No tree, taxonomy, or PCoA object was altered.", call. = FALSE
  )
}
