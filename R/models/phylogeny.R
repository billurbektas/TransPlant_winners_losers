# Module: R/models/phylogeny
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

#' Read the positive phylogenetic variation represented by PEV1–PEV3
#'
#' @param eigen_summary PCoA eigenvalue summary with `PEV` and
#'   `cumulative_percent` columns.
#' @return The finite cumulative percentage at PEV3; errors if absent.
validated_phylogenetic_variation <- function(eigen_summary) {
  value <- eigen_summary$cumulative_percent[eigen_summary$PEV == "PEV3"]
  if (length(value) != 1L || !is.finite(value)) {
    stop("PEV1-PEV3 cumulative positive variation is unavailable.", call. = FALSE)
  }
  as.numeric(value)
}

# Response-independent phylogeny/PCoA validation and orientation.

#' Phylogenetic Broken Stick.
#'
#' @param n_positive A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
phylogenetic_broken_stick <- function(n_positive) {
  if (n_positive < 1L) return(numeric())
  vapply(seq_len(n_positive), function(index) {
    sum(1 / seq.int(index, n_positive)) / n_positive
  }, numeric(1))
}

#' Phylogenetic Eigen Summary.
#'
#' @param pcoa_object A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param solution A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param corrected A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
phylogenetic_eigen_summary <- function(pcoa_object, solution,
                                       corrected = FALSE) {
  values <- if (corrected) {
    as.numeric(pcoa_object$values$Corr_eig)
  } else {
    as.numeric(pcoa_object$values$Eigenvalues)
  }
  positive <- values > 0
  positive_values <- values[positive]
  positive_total <- sum(positive_values)
  proportion <- cumulative <- expectation <- rep(NA_real_, length(values))
  exceeds <- rep(NA, length(values))
  if (length(positive_values)) {
    proportion[positive] <- positive_values / positive_total
    cumulative[positive] <- cumsum(positive_values) / positive_total
    expectation[positive] <- phylogenetic_broken_stick(length(positive_values))
    exceeds[positive] <- proportion[positive] > expectation[positive]
  }
  tibble::tibble(
    solution = solution,
    axis_number = seq_along(values),
    PEV = paste0("PEV", seq_along(values)),
    eigenvalue = values,
    positive_eigenvalue = positive,
    percent_positive_variation = 100 * proportion,
    cumulative_percent = 100 * cumulative,
    broken_stick_expectation_percent = 100 * expectation,
    exceeds_broken_stick = exceeds,
    correction = if (corrected) "Lingoes" else "none"
  )
}
#' Snap Zero Distance Coordinates.
#'
#' @param coordinates A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param distance_matrix A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param tolerance A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
snap_zero_distance_coordinates <- function(coordinates, distance_matrix,
                                           tolerance = 1e-12) {
  indices <- which(
    upper.tri(distance_matrix) & abs(distance_matrix) <= tolerance,
    arr.ind = TRUE
  )
  if (!nrow(indices)) return(coordinates)
  for (row in seq_len(nrow(indices))) {
    first <- rownames(distance_matrix)[indices[row, 1]]
    second <- colnames(distance_matrix)[indices[row, 2]]
    coordinates[second, ] <- coordinates[first, ]
  }
  coordinates
}

#' Recompute Phylogenetic Pcoa.
#'
#' @param tree A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param metadata A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param tip_order A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param meaningful_negative_threshold A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
recompute_phylogenetic_pcoa <- function(tree, metadata = NULL, tip_order = NULL,
                                        meaningful_negative_threshold = 1e-6) {
  distance_matrix <- ape::cophenetic.phylo(tree)
  if (!is.null(tip_order)) {
    tip_order <- as.character(tip_order)
    if (!setequal(tip_order, rownames(distance_matrix))) {
      stop("Requested PCoA order does not match the tree tips.", call. = FALSE)
    }
    distance_matrix <- distance_matrix[tip_order, tip_order, drop = FALSE]
  }
  tolerance <- 1e-10 * max(1, max(abs(distance_matrix)))
  checks <- tibble::tibble(
    check = c(
      "square", "unique taxa", "row/column names", "tree order",
      "symmetric", "negligible diagonal", "complete"
    ),
    pass = c(
      nrow(distance_matrix) == ncol(distance_matrix),
      !anyDuplicated(rownames(distance_matrix)) &&
        !anyDuplicated(colnames(distance_matrix)),
      identical(rownames(distance_matrix), colnames(distance_matrix)),
      identical(
        rownames(distance_matrix),
        if (is.null(tip_order)) tree$tip.label else tip_order
      ),
      max(abs(distance_matrix - t(distance_matrix))) <= tolerance,
      max(abs(diag(distance_matrix))) <= tolerance,
      !anyNA(distance_matrix) && all(is.finite(distance_matrix))
    )
  )
  if (any(!checks$pass)) {
    stop("Phylogenetic distance validation failed: ",
         paste(checks$check[!checks$pass], collapse = ", "), call. = FALSE)
  }

  uncorrected <- ape::pcoa(stats::as.dist(distance_matrix), correction = "none")
  eigenvalues <- as.numeric(uncorrected$values$Eigenvalues)
  positive_sum <- sum(eigenvalues[eigenvalues > 0])
  negative <- eigenvalues[eigenvalues < 0]
  negative_ratio <- if (positive_sum > 0) sum(abs(negative)) / positive_sum else Inf
  correction_required <- length(negative) > 0L && is.finite(negative_ratio) &&
    negative_ratio > meaningful_negative_threshold
  corrected <- if (correction_required) {
    ape::pcoa(stats::as.dist(distance_matrix), correction = "lingoes")
  } else {
    NULL
  }
  eigen_summary <- phylogenetic_eigen_summary(
    uncorrected, "uncorrected", corrected = FALSE
  )
  if (correction_required) {
    eigen_summary <- dplyr::bind_rows(
      eigen_summary,
      phylogenetic_eigen_summary(corrected, "corrected", corrected = TRUE)
    )
  }
  selected <- if (correction_required) corrected else uncorrected
  solution <- if (correction_required) "corrected" else "uncorrected"
  coordinates <- if (correction_required) {
    as.matrix(selected$vectors.cor)
  } else {
    as.matrix(selected$vectors)
  }
  selected_summary <- eigen_summary |>
    dplyr::filter(
      .data$solution == .env$solution,
      .data$positive_eigenvalue
    )
  if (ncol(coordinates) != nrow(selected_summary)) {
    stop("PCoA coordinate count differs from the positive-axis count.",
         call. = FALSE)
  }
  colnames(coordinates) <- paste0("PEV", seq_len(ncol(coordinates)))
  if (is.null(rownames(coordinates))) rownames(coordinates) <- rownames(distance_matrix)
  coordinate_order <- if (is.null(tip_order)) tree$tip.label else tip_order
  coordinates <- coordinates[coordinate_order, , drop = FALSE]
  coordinates <- snap_zero_distance_coordinates(coordinates, distance_matrix)
  scores <- tibble::as_tibble(coordinates, rownames = "species_id")
  if (!is.null(metadata)) {
    metadata_columns <- setdiff(names(metadata), grep("^PEV[0-9]+$", names(metadata), value = TRUE))
    scores <- scores |>
      dplyr::left_join(metadata[, metadata_columns, drop = FALSE], by = "species_id") |>
      dplyr::select(dplyr::all_of(metadata_columns), dplyr::starts_with("PEV"))
  }
  scores <- orient_pev_axes(scores)
  list(
    distance = distance_matrix,
    distance_checks = checks,
    eigen_summary = eigen_summary,
    scores = scores,
    correction_required = correction_required,
    negative_ratio = negative_ratio
  )
}

#' Pev Orientation Anchors.
#'
#' @return An R object described by the function name; its structure is validated before later stages use it.
pev_orientation_anchors <- function() {
  tibble::tribble(
    ~axis, ~species_id, ~required_sign,
    "PEV1", "Poa trivialis", 1,
    "PEV2", "Erigeron formosissimus", -1,
    "PEV3", "Euphrasia minima", 1
  )
}
#' Build Phylogeny Species Order.
#'
#' @param tree A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param decisions A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param resident_species Species identifiers in the resident modeling cohort.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_phylogeny_species_order <- function(tree, decisions, resident_species) {
  if (!inherits(tree, "phylo")) {
    stop("Phylogeny species order requires a phylo object.", call. = FALSE)
  }
  resident_species <- sort(unique(as.character(resident_species)))
  tree_species <- as.character(tree$tip.label)
  if (anyDuplicated(tree_species) || !setequal(tree_species, resident_species)) {
    stop(
      "Tree tips do not uniquely match the resident species set.",
      call. = FALSE
    )
  }
  assert_columns(
    decisions,
    c("species_id", "TNRS accepted name", "TNRS genus", "TNRS family"),
    "phylogeny placement decisions"
  )
  taxonomy <- decisions |>
    dplyr::transmute(
      species_id = as.character(.data$species_id),
      tnrs_accepted_name = as.character(.data$`TNRS accepted name`),
      tnrs_genus = as.character(.data$`TNRS genus`),
      tnrs_family = as.character(.data$`TNRS family`)
    )
  result <- tibble::tibble(
    tree_tip_position = seq_along(tree_species),
    species_id = tree_species,
    newick_tip_label = gsub(" ", "_", tree_species, fixed = TRUE)
  ) |>
    dplyr::left_join(taxonomy, by = "species_id")
  if (anyNA(result[c("tnrs_accepted_name", "tnrs_genus", "tnrs_family")])) {
    stop("Active phylogeny species order has missing reviewed taxonomy.", call. = FALSE)
  }
  result
}

#' Build Zero Distance Audit.
#'
#' @param distance_matrix A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param decisions A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param tolerance A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_zero_distance_audit <- function(distance_matrix, decisions,
                                      tolerance = NULL) {
  matrix <- as.matrix(distance_matrix)
  if (is.null(tolerance)) {
    tolerance <- 1e-10 * max(1, max(abs(matrix)))
  }
  indices <- which(
    upper.tri(matrix) & abs(matrix) <= tolerance,
    arr.ind = TRUE
  )
  assert_columns(
    decisions,
    c(
      "species_id", "TNRS accepted name", "match_type", "parent_binomial",
      "parent_in_backbone", "shared_parent_placement", "infraspecific_name",
      "backbone matched name"
    ),
    "phylogeny placement decisions"
  )
  lookup <- decisions |>
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
  if (!nrow(indices)) {
    return(tibble::tibble(
      species_id_1 = character(), `TNRS accepted name 1` = character(),
      species_id_2 = character(), `TNRS accepted name 2` = character(),
      phylogenetic_distance = numeric(), `match_type species 1` = character(),
      `match_type species 2` = character(), `parent binomial 1` = character(),
      `parent tip present 1` = logical(), `parent binomial 2` = character(),
      `parent tip present 2` = logical(),
      `either taxon infraspecific` = logical(), reason = character(),
      classification = character()
    ))
  }
  tibble::tibble(
    species_id_1 = rownames(matrix)[indices[, "row"]],
    species_id_2 = colnames(matrix)[indices[, "col"]],
    phylogenetic_distance = matrix[indices]
  ) |>
    dplyr::left_join(
      lookup |>
        dplyr::rename_with(~ paste0(.x, "_1"), -"species_id") |>
        dplyr::rename(species_id_1 = "species_id"),
      by = "species_id_1"
    ) |>
    dplyr::left_join(
      lookup |>
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
        .data$intentional,
        "intentional shared placement",
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
      reason = .data$reason,
      classification = .data$classification
    )
}

#' Build Phylogeny from Placement Decisions.
#'
#' @param decisions A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param resident_species Species identifiers in the resident modeling cohort.
#' @return An R object described by the function name; its structure is validated before later stages use it.
build_phylogeny_from_placement_decisions <- function(decisions,
                                                     resident_species) {
  required <- c(
    "species_id", "TNRS genus", "TNRS family", "backbone_status",
    "match_type", "phylo_input_label", "shared_parent_placement"
  )
  assert_columns(decisions, required, "phylogeny placement decisions")
  decisions <- decisions |>
    dplyr::mutate(
      species_id = as.character(.data$species_id),
      phylo_input_label = as.character(.data$phylo_input_label),
      shared_parent_placement = as.logical(.data$shared_parent_placement)
    ) |>
    dplyr::arrange(.data$species_id)
  resident_species <- sort(unique(as.character(resident_species)))
  if (anyDuplicated(decisions$species_id) ||
      !identical(decisions$species_id, resident_species)) {
    stop(
      "Phylogeny placement decisions do not exactly match the resident species set.",
      call. = FALSE
    )
  }
  if (any(decisions$match_type == "unresolved") ||
      anyNA(decisions[c("TNRS genus", "TNRS family", "phylo_input_label")])) {
    stop("Phylogeny placement decisions contain unresolved or incomplete rows.",
         call. = FALSE)
  }
  placement_groups <- decisions |>
    dplyr::summarise(
      n_ecological_taxa = dplyr::n(),
      shared_parent_documented = all(.data$shared_parent_placement),
      .by = "phylo_input_label"
    )
  invalid_duplicates <- placement_groups |>
    dplyr::filter(
      .data$n_ecological_taxa > 1L,
      !.data$shared_parent_documented
    )
  if (nrow(invalid_duplicates)) {
    stop("A duplicated tree placement lacks a documented shared parent.",
         call. = FALSE)
  }
  sp_list <- decisions |>
    dplyr::transmute(
      species = .data$phylo_input_label,
      genus = as.character(.data$`TNRS genus`),
      family = as.character(.data$`TNRS family`)
    ) |>
    dplyr::distinct(.data$species, .keep_all = TRUE) |>
    as.data.frame(stringsAsFactors = FALSE)
  phylo_warnings <- character()
  phylo_messages <- character()
  phylo_result <- withCallingHandlers(
    V.PhyloMaker2::phylo.maker(
      sp.list = sp_list,
      tree = V.PhyloMaker2::GBOTB.extended.LCVP,
      nodes = V.PhyloMaker2::nodes.info.1.LCVP,
      output.sp.list = TRUE,
      output.tree = FALSE,
      scenarios = "S3"
    ),
    warning = function(condition) {
      phylo_warnings <<- c(phylo_warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    },
    message = function(condition) {
      phylo_messages <<- c(phylo_messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    }
  )
  if (!all(c("scenario.3", "species.list") %in% names(phylo_result)) ||
      !inherits(phylo_result$scenario.3, "phylo")) {
    stop("V.PhyloMaker2 did not return the deterministic scenario.3 tree.",
         call. = FALSE)
  }
  returned_placement <- tibble::as_tibble(phylo_result$species.list) |>
    dplyr::transmute(
      phylo_input_label = as.character(.data$species),
      returned_backbone_status = as.character(.data$status)
    )
  placement_check <- decisions |>
    dplyr::distinct(
      .data$phylo_input_label, expected_backbone_status = .data$backbone_status
    ) |>
    dplyr::left_join(returned_placement, by = "phylo_input_label")
  if (anyNA(placement_check$returned_backbone_status) ||
      !identical(
        placement_check$expected_backbone_status,
        placement_check$returned_backbone_status
      )) {
    stop("Rebuilt S3 placement statuses differ from the decision ledger.",
         call. = FALSE)
  }

  tree <- phylo_result$scenario.3
  tip_map <- decisions |>
    dplyr::select("phylo_input_label", "species_id")
  tip_groups <- split(tip_map$species_id, tip_map$phylo_input_label)
  duplicate_tip_groups <- tip_groups[lengths(tip_groups) > 1L]
  for (placement_label in names(duplicate_tip_groups)) {
    where <- match(placement_label, tree$tip.label)
    if (is.na(where)) {
      stop("Shared parent placement is absent from rebuilt tree: ",
           placement_label, call. = FALSE)
    }
    daughter <- ape::stree(
      length(duplicate_tip_groups[[placement_label]]), type = "star",
      tip.label = duplicate_tip_groups[[placement_label]]
    )
    daughter$edge.length <- rep(0, nrow(daughter$edge))
    tree <- ape::bind.tree(tree, daughter, where = where, position = 0)
  }
  singleton_map <- tip_map |>
    dplyr::add_count(.data$phylo_input_label, name = "n_ecological_taxa") |>
    dplyr::filter(.data$n_ecological_taxa == 1L)
  map_index <- match(tree$tip.label, singleton_map$phylo_input_label)
  mapped <- !is.na(map_index)
  tree$tip.label[mapped] <- singleton_map$species_id[map_index[mapped]]

  edge_lengths <- as.numeric(tree$edge.length)
  tip_depths <- ape::node.depth.edgelength(tree)[seq_len(ape::Ntip(tree))]
  root_to_tip_discrepancy <- max(tip_depths) - min(tip_depths)
  tree_checks <- tibble::tibble(
    check = c(
      "exact resident species", "unique tree tips", "rooted", "ultrametric",
      "finite branch lengths", "nonnegative branch lengths"
    ),
    pass = c(
      identical(sort(tree$tip.label), resident_species),
      !anyDuplicated(tree$tip.label),
      ape::is.rooted(tree),
      ape::is.ultrametric(tree, option = 1, tol = 1e-8),
      length(edge_lengths) > 0L && all(is.finite(edge_lengths)),
      length(edge_lengths) > 0L && all(edge_lengths >= 0)
    ),
    detail = c(
      paste(length(tree$tip.label), "tips"),
      paste(anyDuplicated(tree$tip.label), "duplicated tip index"),
      "scenario 3 tree",
      paste("root-to-tip discrepancy", format(root_to_tip_discrepancy, digits = 9)),
      paste(length(edge_lengths), "edges"),
      paste(sum(edge_lengths == 0), "zero-length edges")
    )
  )
  if (any(!tree_checks$pass)) {
    stop("Rebuilt phylogeny failed: ",
         paste(tree_checks$check[!tree_checks$pass], collapse = ", "),
         call. = FALSE)
  }
  list(
    tree = tree,
    decisions = decisions,
    sp_list = sp_list,
    returned_placement = returned_placement,
    warnings = unique(phylo_warnings),
    messages = unique(phylo_messages),
    audit = dplyr::bind_rows(
      tree_checks,
      tibble::tibble(
        check = c("placement statuses", "V.PhyloMaker2 warnings"),
        pass = c(TRUE, !length(phylo_warnings)),
        detail = c(
          paste(nrow(placement_check), "unique scenario-3 placement inputs"),
          if (length(phylo_warnings)) {
            paste(unique(phylo_warnings), collapse = " | ")
          } else {
            "No warnings"
          }
        )
      )
    )
  )
}

#' Orient Pev Axes.
#'
#' @param scores A data frame whose required columns are checked before use.
#' @param anchors A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
orient_pev_axes <- function(scores, anchors = pev_orientation_anchors()) {
  result <- scores
  for (index in seq_len(nrow(anchors))) {
    axis <- anchors$axis[[index]]
    species <- anchors$species_id[[index]]
    required_sign <- anchors$required_sign[[index]]
    if (!axis %in% names(result)) stop("Missing PCoA axis: ", axis, call. = FALSE)
    value <- result[[axis]][match(species, result$species_id)]
    if (length(value) != 1L || is.na(value) || value == 0) {
      stop("PCoA orientation anchor unavailable for ", axis, ": ", species,
           call. = FALSE)
    }
    if (sign(value) != required_sign) result[[axis]] <- -result[[axis]]
  }
  result
}




#' Validate Zero Distance Placements.
#'
#' @param distance_matrix A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param audit A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return A structured validation result; fail-closed checks raise an error when an invariant is violated.
validate_zero_distance_placements <- function(distance_matrix, audit) {
  matrix <- as.matrix(distance_matrix)
  indices <- which(upper.tri(matrix) & abs(matrix) <= 1e-12, arr.ind = TRUE)
  pairs <- if (nrow(indices)) {
    apply(indices, 1, function(index) {
      paste(sort(c(rownames(matrix)[index[[1]]], colnames(matrix)[index[[2]]])),
            collapse = "||")
    })
  } else character()
  audited <- apply(audit[, c("species_id_1", "species_id_2")], 1, function(x) {
    paste(sort(as.character(x)), collapse = "||")
  })
  tibble::tibble(
    check = c("zero-distance pair set", "intentional shared placements"),
    pass = c(setequal(pairs, audited), all(audit$classification == "intentional shared placement")),
    detail = c(paste(length(pairs), "matrix pairs and", length(audited), "audit rows"),
               paste(sum(audit$classification == "intentional shared placement"), "documented"))
  )
}
