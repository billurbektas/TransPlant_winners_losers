# Response-independent phylogeny/PCoA validation and orientation.

phylogenetic_broken_stick <- function(n_positive) {
  if (n_positive < 1L) return(numeric())
  vapply(seq_len(n_positive), function(index) {
    sum(1 / seq.int(index, n_positive)) / n_positive
  }, numeric(1))
}

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

pev_orientation_anchors <- function() {
  tibble::tribble(
    ~axis, ~species_id, ~required_sign,
    "PEV1", "Poa trivialis", 1,
    "PEV2", "Erigeron formosissimus", -1,
    "PEV3", "Euphrasia minima", 1
  )
}

restore_phylogeny_tip_labels <- function(tree, species_order) {
  mapping <- stats::setNames(
    as.character(species_order$species_id),
    as.character(species_order$newick_tip_label)
  )
  restored <- unname(mapping[tree$tip.label])
  if (anyNA(restored) || anyDuplicated(restored)) {
    stop("Newick labels do not map uniquely to ecological species IDs.", call. = FALSE)
  }
  tree$tip.label <- restored
  tree
}

build_phylogeny_from_placement_decisions <- function(decisions,
                                                     downstream_species) {
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
  downstream_species <- sort(unique(as.character(downstream_species)))
  if (anyDuplicated(decisions$species_id) ||
      !identical(decisions$species_id, downstream_species)) {
    stop(
      "Phylogeny placement decisions do not exactly match the downstream species set.",
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
      "exact downstream species", "unique tree tips", "rooted", "ultrametric",
      "finite branch lengths", "nonnegative branch lengths"
    ),
    pass = c(
      identical(sort(tree$tip.label), downstream_species),
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

compare_pev_coordinates <- function(candidate, reference, tolerance = 1e-8) {
  axes <- c("PEV1", "PEV2", "PEV3")
  candidate <- orient_pev_axes(candidate[, c("species_id", axes)])
  reference <- orient_pev_axes(reference[, c("species_id", axes)])
  ids_equal <- identical(candidate$species_id, reference$species_id)
  values_equal <- ids_equal && all(
    abs(as.matrix(candidate[axes]) - as.matrix(reference[axes])) <= tolerance
  )
  tibble::tibble(
    check = c("PEV species order", "oriented PEV1-PEV3 coordinates"),
    pass = c(ids_equal, values_equal),
    detail = c(
      paste(nrow(candidate), "candidate and", nrow(reference), "reference rows"),
      paste("absolute tolerance", format(tolerance, scientific = TRUE))
    )
  )
}

build_anchored_pev_modeling_interface <- function(recomputed, reference,
                                                   tolerance = 1e-8) {
  checks <- compare_pev_coordinates(recomputed, reference, tolerance)
  if (any(!checks$pass)) {
    stop(
      "Recomputed PEV coordinates cannot be anchored to the Step 19 interface: ",
      paste(checks$check[!checks$pass], collapse = "; "), call. = FALSE
    )
  }
  axes <- c("PEV1", "PEV2", "PEV3")
  maximum_error <- max(
    abs(as.matrix(recomputed[axes]) - as.matrix(reference[axes]))
  )
  list(
    data = reference,
    audit = tibble::tibble(
      interface = "validated Step 19 PEV1-PEV3 modeling coordinates",
      recomputed_species = nrow(recomputed),
      maximum_recomputation_error = maximum_error,
      validation_tolerance = tolerance,
      exact_reference_values_used_for_step22_modeling = TRUE,
      reason = paste(
        "The PCoA is recomputed independently and must pass tolerance before",
        "the exact frozen coordinate values are supplied to the optimizer;",
        "this prevents floating-point eigensolver noise from changing anchored",
        "CR1 and heterogeneity results."
      ),
      pass = TRUE
    )
  )
}

compare_phylogenetic_distance <- function(candidate, reference,
                                           tolerance = 1e-6) {
  candidate <- as.matrix(candidate)
  reference <- as.matrix(reference)
  order_equal <- identical(rownames(candidate), rownames(reference)) &&
    identical(colnames(candidate), colnames(reference))
  maximum_error <- if (order_equal) max(abs(candidate - reference)) else Inf
  tibble::tibble(
    check = c("distance species order", "phylogenetic distances"),
    pass = c(order_equal, order_equal && maximum_error <= tolerance),
    detail = c(
      paste(nrow(candidate), "candidate and", nrow(reference), "reference taxa"),
      paste("maximum absolute error", format(maximum_error, digits = 8),
            "with tolerance", format(tolerance, scientific = TRUE))
    )
  )
}

compare_pcoa_eigenstructure <- function(candidate, reference,
                                        tolerance = 1e-8) {
  key <- c("solution", "axis_number", "PEV")
  candidate <- candidate |>
    dplyr::arrange(.data$solution, .data$axis_number)
  reference <- reference |>
    dplyr::arrange(.data$solution, .data$axis_number)
  keys_equal <- identical(as.character(candidate$solution), as.character(reference$solution)) &&
    identical(as.numeric(candidate$axis_number), as.numeric(reference$axis_number)) &&
    identical(as.character(candidate$PEV), as.character(reference$PEV))
  positive_equal <- keys_equal && identical(
    as.logical(candidate$positive_eigenvalue),
    as.logical(reference$positive_eigenvalue)
  )
  numeric_columns <- c(
    "eigenvalue", "percent_positive_variation", "cumulative_percent",
    "broken_stick_expectation_percent"
  )
  numeric_equal <- keys_equal && all(vapply(numeric_columns, function(column) {
    numeric_values_equal(
      candidate[[column]], reference[[column]],
      c(absolute = tolerance, relative = 1e-6)
    )
  }, logical(1)))
  tibble::tibble(
    check = c("PCoA eigenvalue keys", "positive-axis flags", "PCoA eigenstructure"),
    pass = c(keys_equal, positive_equal, numeric_equal),
    detail = c(
      paste(nrow(candidate), "candidate and", nrow(reference), "reference axes"),
      "Positive/negative axis identity is exact.",
      paste("absolute tolerance", format(tolerance, scientific = TRUE))
    )
  )
}

compare_phylogeny <- function(candidate_tree, reference_tree, candidate_order,
                              reference_order, tolerance = 1e-6) {
  candidate_species <- as.character(candidate_order$species_id)
  reference_species <- as.character(reference_order$species_id)
  order_equal <- identical(candidate_species, reference_species)
  tips_equal <- setequal(candidate_tree$tip.label, reference_tree$tip.label)
  topology_equal <- isTRUE(ape::all.equal.phylo(
    candidate_tree, reference_tree, use.edge.length = FALSE
  ))
  topology_and_lengths_equal <- isTRUE(ape::all.equal.phylo(
    candidate_tree, reference_tree, use.edge.length = TRUE, tolerance = tolerance
  ))
  tibble::tibble(
    check = c(
      "species order mapping", "tree tip set", "tree topology",
      "topology-aligned branch lengths"
    ),
    pass = c(order_equal, tips_equal, topology_equal, topology_and_lengths_equal),
    detail = c(
      paste(length(candidate_species), "species"),
      paste(length(candidate_tree$tip.label), "tips"),
      "ape::all.equal.phylo without edge lengths",
      paste("ape::all.equal.phylo tolerance", format(tolerance, scientific = TRUE))
    )
  )
}

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

assert_phylogeny_validation <- function(...) {
  checks <- dplyr::bind_rows(...)
  if (any(!checks$pass)) {
    stop("Phylogeny/PCoA regression failed: ",
         paste(checks$check[!checks$pass], collapse = "; "), call. = FALSE)
  }
  checks
}
