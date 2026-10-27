# Phylogeny and its PCoA axes (P3), for every study species. The tree is cut
# from the V.PhyloMaker2 mega-tree GBOTB.extended.LCVP; each species' place
# in it is set by the reviewed ledger prep/config/phylogeny_placements.csv.
# Logic follows legacy/R/models/phylogeny.R (lines named below).

backbone_tips <- function() V.PhyloMaker2::GBOTB.extended.LCVP$tip.label

backbone_genus_family <- function() {
  nodes <- V.PhyloMaker2::nodes.info.1.LCVP
  unique(nodes[nodes$level == "G", c("genus", "family")])
}

# Draft of the placement ledger: legacy reviewed rows where the legacy
# species is the same species, then rules for the rest. Written for review;
# the pipeline itself reads only the reviewed ledger.
placement_draft <- function(study_species, legacy, genus_map, taxonomy, legacy_taxonomy,
                            path = "prep/phylogeny_placements_draft.csv") {
  tips <- backbone_tips()
  gf <- backbone_genus_family()
  label_of <- function(x) gsub(" ", "_", x)
  binomial <- function(x) stringr::word(x, 1, 2)

  # a legacy row is reused only for the same species (binomial of the legacy
  # name equals ours): rows that stood for a homonym-corrected species, such
  # as legacy "Gentianella germanica" for Gentiana campestris, do not match
  # A renamed legacy species (Minuartia verna, now Sabulina verna) keeps its
  # reviewed label: legacy name -> submitted names (legacy taxonomy, copied
  # from step 00 to prep/config/legacy_taxonomy.csv by P5) -> our
  # accepted species, accepted only when the epithet stays the same (first
  # five letters, so gender endings match). Homonym errors change the epithet
  # (legacy Gentianella germanica stood for Gentiana campestris) and are not
  # linked.
  epithet5 <- function(x) substr(tolower(stringr::word(x, 2)), 1, 5)
  chain <- legacy_taxonomy |>
    dplyr::select(submitted_name, legacy_accepted = accepted_name_reviewed) |>
    dplyr::inner_join(dplyr::select(taxonomy, submitted_name = name_original, accepted_species),
                      by = "submitted_name") |>
    dplyr::mutate(legacy_binomial = binomial(legacy_accepted)) |>
    dplyr::filter(!is.na(accepted_species), legacy_binomial != accepted_species,
                  epithet5(legacy_binomial) == epithet5(accepted_species)) |>
    dplyr::distinct(legacy_binomial, renamed_to = accepted_species)
  legacy <- legacy |>
    dplyr::mutate(species = binomial(legacy_species_id),
                  renamed_to = chain$renamed_to[match(species, chain$legacy_binomial)],
                  species = dplyr::if_else(!species %in% study_species$species & !is.na(renamed_to),
                                           renamed_to, species),
                  exact_name = legacy_species_id == species,
                  prune = backbone_status == "prune") |>
    dplyr::arrange(species, dplyr::desc(exact_name), dplyr::desc(prune)) |>
    dplyr::distinct(species, .keep_all = TRUE)

  d <- study_species |>
    dplyr::left_join(dplyr::select(legacy, species, legacy_species_id, legacy_label = phylo_input_label,
                                   legacy_status = backbone_status, legacy_match = match_type),
                     by = "species") |>
    dplyr::mutate(
      mapped_genus = genus_map$backbone_genus[match(genus, genus_map$study_genus)],
      own_tip = label_of(species) %in% tips,
      mapped_tip = !is.na(mapped_genus) & paste(mapped_genus, stringr::word(species, 2), sep = "_") %in% tips,
      legacy_tip = !is.na(legacy_label) & legacy_status == "prune" & legacy_label %in% tips,
      placement = dplyr::case_when(
        own_tip ~ "exact tip",
        legacy_tip ~ "exact tip under an older name (legacy reviewed)",
        mapped_tip ~ "exact tip under an older genus",
        genus %in% gf$genus ~ "grafted within genus",
        !is.na(mapped_genus) & mapped_genus %in% gf$genus ~ "grafted within older genus",
        TRUE ~ "grafted within family"),
      backbone_label = dplyr::case_when(
        own_tip ~ label_of(species),
        legacy_tip ~ legacy_label,
        mapped_tip ~ paste(mapped_genus, stringr::word(species, 2), sep = "_"),
        TRUE ~ NA_character_),
      sp_genus = dplyr::case_when(
        !is.na(backbone_label) ~ stringr::word(backbone_label, 1, sep = "_"),
        placement == "grafted within older genus" ~ mapped_genus,
        TRUE ~ genus),
      sp_family = dplyr::coalesce(gf$family[match(sp_genus, gf$genus)], family),
      source = dplyr::case_when(
        !is.na(legacy_species_id) & legacy_tip & binomial(legacy_species_id) != species ~ "legacy reviewed label (renamed species)",
        !is.na(legacy_species_id) & legacy_tip ~ "legacy reviewed label",
        !is.na(legacy_species_id) ~ "legacy row; placement by rule",
        TRUE ~ "new species; placement by rule"))
  d$check <- placement_checks(d, tips, gf)
  out <- dplyr::select(d, species, genus, family, placement, backbone_label, sp_genus, sp_family,
                       source, legacy_species_id, legacy_label, legacy_match, check)
  utils::write.csv(out, path, row.names = FALSE)
  out
}

# Checks on every ledger row; "ok" or the problems found
placement_checks <- function(d, tips = backbone_tips(), gf = backbone_genus_family()) {
  problems <- vapply(seq_len(nrow(d)), function(i) {
    p <- character()
    if (!is.na(d$backbone_label[i]) && !d$backbone_label[i] %in% tips) p <- c(p, "label is not a backbone tip")
    if (grepl("within genus|within older genus", d$placement[i]) && !d$sp_genus[i] %in% gf$genus)
      p <- c(p, "binding genus not in backbone")
    if (!d$sp_family[i] %in% gf$family) p <- c(p, "family not in backbone")
    fam_of_genus <- gf$family[match(d$sp_genus[i], gf$genus)]
    if (!is.na(fam_of_genus) && fam_of_genus != d$family[i]) p <- c(p, paste0("backbone puts ", d$sp_genus[i], " in ", fam_of_genus))
    if (!is.na(d$backbone_label[i]) && sum(d$backbone_label %in% d$backbone_label[i]) > 1) p <- c(p, "label shared with another species")
    if (length(p)) paste(p, collapse = "; ") else "ok"
  }, character(1))
  problems
}

# The reviewed ledger must cover exactly the study species and pass the checks
read_placements <- function(path, study_species) {
  p <- utils::read.csv(path, colClasses = "character", na.strings = c("", "NA"))
  if (!setequal(p$species, study_species$species) || anyDuplicated(p$species)) {
    stop("phylogeny_placements.csv does not hold exactly the study species: missing ",
         paste(utils::head(setdiff(study_species$species, p$species), 10), collapse = ", "), "; extra ",
         paste(utils::head(setdiff(p$species, study_species$species), 10), collapse = ", "))
  }
  p$check <- placement_checks(p)
  hard <- grepl("not a backbone tip|binding genus not in backbone|family not in backbone|shared with another", p$check)
  if (any(hard)) stop("Placement ledger rows failing checks: ", paste(p$species[hard], collapse = ", "))
  p
}

# legacy/R/models/phylogeny.R:356-504, with one tip per study species
build_phylogeny <- function(placements) {
  sp_list <- data.frame(
    species = dplyr::coalesce(placements$backbone_label, gsub(" ", "_", placements$species)),
    genus = placements$sp_genus, family = placements$sp_family)
  warns <- character()
  res <- withCallingHandlers(
    V.PhyloMaker2::phylo.maker(sp.list = sp_list, tree = V.PhyloMaker2::GBOTB.extended.LCVP,
                               nodes = V.PhyloMaker2::nodes.info.1.LCVP, output.sp.list = TRUE,
                               output.tree = FALSE, scenarios = "S3"),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) invokeRestart("muffleMessage"))
  tree <- res$scenario.3
  # tips back to the accepted species names
  tree$tip.label <- placements$species[match(tree$tip.label, sp_list$species)]
  if (!setequal(tree$tip.label, placements$species) || anyDuplicated(tree$tip.label)) stop("Tree tips are not the study species")
  if (!ape::is.rooted(tree)) stop("Tree is not rooted")
  if (!ape::is.ultrametric(tree, option = 1, tol = 1e-8)) stop("Tree is not ultrametric")
  if (any(!is.finite(tree$edge.length) | tree$edge.length < 0)) stop("Negative or non-finite branch lengths")
  list(tree = tree, species_list = res$species.list, warnings = warns)
}

# legacy/R/models/phylogeny.R:26-196 and :535 (PCoA, eigen summary,
# zero-distance snapping, orientation)
phylogeny_pcoa <- function(phylogeny, negative_threshold = 1e-6) {
  D <- ape::cophenetic.phylo(phylogeny$tree)
  D <- D[sort(rownames(D)), sort(colnames(D))]
  p <- ape::pcoa(stats::as.dist(D), correction = "none")
  eig <- p$values$Eigenvalues
  negative_ratio <- sum(abs(eig[eig < 0])) / sum(eig[eig > 0])
  corrected <- any(eig < 0) && negative_ratio > negative_threshold
  if (corrected) {
    p <- ape::pcoa(stats::as.dist(D), correction = "lingoes")
    eig <- p$values$Corr_eig
    vectors <- p$vectors.cor
  } else {
    vectors <- p$vectors
  }
  positive <- which(eig > 0)
  vectors <- vectors[, positive, drop = FALSE]
  colnames(vectors) <- paste0("PEV", seq_len(ncol(vectors)))
  # zero-distance pairs get identical coordinates (legacy/R/models/phylogeny.R:76)
  zero <- which(abs(D) <= 1e-12 & upper.tri(D), arr.ind = TRUE)
  for (k in seq_len(nrow(zero))) vectors[zero[k, 2], ] <- vectors[zero[k, 1], ]
  n <- length(positive)
  broken_stick <- vapply(seq_len(n), function(i) sum(1 / (i:n)) / n, numeric(1)) * 100
  pct <- eig[positive] / sum(eig[positive]) * 100
  eigenvalues <- tibble::tibble(PEV = colnames(vectors), axis_number = seq_len(n), eigenvalue = eig[positive],
                                percent_positive_variation = pct, cumulative_percent = cumsum(pct),
                                broken_stick_expectation_percent = broken_stick,
                                exceeds_broken_stick = pct > broken_stick,
                                correction = if (corrected) "Lingoes" else "none")
  scores <- tibble::as_tibble(vectors, rownames = "species")
  # orientation anchors as legacy/R/models/phylogeny.R:196
  anchors <- c(PEV1 = "Poa trivialis", PEV2 = "Erigeron formosissimus", PEV3 = "Euphrasia minima")
  signs <- c(PEV1 = 1, PEV2 = -1, PEV3 = 1)
  for (ax in names(anchors)) {
    v <- scores[[ax]][scores$species == anchors[[ax]]]
    if (!length(v) || is.na(v) || v == 0) stop("Orientation anchor ", anchors[[ax]], " unavailable for ", ax)
    if (sign(v) != signs[[ax]]) scores[[ax]] <- -scores[[ax]]
  }
  list(scores = scores, eigenvalues = eigenvalues, negative_ratio = negative_ratio, corrected = corrected)
}
