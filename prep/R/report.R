# Output files and the P1 review packet

write_prep_outputs <- function(community, sites, experiments, taxonomy, species_merges,
                               dir = "data/derived/prep") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  files <- c(community = file.path(dir, "community.rds"), sites = file.path(dir, "sites.rds"),
             experiments = file.path(dir, "experiments.rds"),
             species_taxonomy = file.path(dir, "species_taxonomy.rds"),
             species_merges = file.path(dir, "species_merges.csv"))
  saveRDS(community, files[["community"]])
  saveRDS(sites, files[["sites"]])
  saveRDS(experiments, files[["experiments"]])
  saveRDS(taxonomy, files[["species_taxonomy"]])
  utils::write.csv(species_merges, files[["species_merges"]], row.names = FALSE)
  unname(files)
}

md_table <- function(d) {
  d <- dplyr::mutate(d, dplyr::across(dplyr::everything(), ~ ifelse(is.na(.x), "", as.character(.x))))
  c(paste("|", paste(names(d), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(d)), collapse = " | "), "|"),
    apply(d, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}

write_p1_report <- function(community_raw, community_named, community_dated, community, sites,
                            experiments, taxonomy, name_resolution_study, name_resolution_try, names_all,
                            name_resolution, try_candidate_names, try_names, species_merges,
                            community_edits, accepted_renames, rdata_file, frozen_taxonomy_file,
                            path = "reports/P1_community_taxonomy.md") {
  run <- attr(name_resolution_study, "resolved")

  # Reconciliation with the RData the legacy pipeline read
  env <- new.env()
  load(rdata_file, envir = env)
  old <- dplyr::count(tibble::tibble(region = env$dat$Region, year = as.integer(env$dat$Year)),
                      region, year, name = "rows_rdata")
  new <- dplyr::count(community_raw, region, year, name = "rows_sqlite")
  recon <- dplyr::full_join(new, old, by = c("region", "year")) |>
    dplyr::mutate(dplyr::across(dplyr::starts_with("rows"), ~ dplyr::coalesce(.x, 0L)),
                  difference = rows_sqlite - rows_rdata) |>
    dplyr::arrange(region, year)
  recon_diff <- dplyr::filter(recon, difference != 0)

  # Counts at every step
  dropped <- attr(community_dated, "dropped")
  set_aside <- dplyr::count(dplyr::filter(community, !is.na(drop_reason)), drop_reason, name = "rows")
  analysis <- dplyr::filter(community, is.na(drop_reason))
  steps <- tibble::tibble(
    step = c("SQLite community rows", paste("dropped:", dropped$reason),
             "rows kept (all vegetation records)",
             paste("  of which set aside:", set_aside$drop_reason),
             "analysis records before collapsing", "analysis rows (plot-year x species)"),
    rows = c(nrow(community_raw), dropped$rows, nrow(community_dated), set_aside$rows,
             sum(analysis$n_records), nrow(analysis)))

  # Name resolution
  res <- dplyr::left_join(names_all, dplyr::select(name_resolution, name = submitted_name, accepted_name, needs_review),
                          by = "name")
  n_try_all <- attr(try_candidate_names, "n_try_names")
  res_summary <- tibble::tibble(
    names = c("community", "TRY", "Lynn"),
    submitted = c(sum(res$in_community), sum(res$in_try), sum(res$in_lynn)),
    resolved = c(sum(res$in_community & !is.na(res$accepted_name)), sum(res$in_try & !is.na(res$accepted_name)),
                 sum(res$in_lynn & !is.na(res$accepted_name))),
    needs_review = c(sum(res$in_community & res$needs_review), sum(res$in_try & res$needs_review),
                     sum(res$in_lynn & res$needs_review)))
  levels_tab <- dplyr::count(taxonomy, identification_level = dplyr::coalesce(identification_level, "unidentified"),
                             decided_by, name = "names")

  # Changes against the frozen legacy taxonomy (step 00)
  frozen <- readRDS(frozen_taxonomy_file) |>
    dplyr::select(name_original = submitted_name, legacy_accepted = accepted_name_reviewed)
  changes <- taxonomy |>
    dplyr::filter(identification_level == "species_or_lower") |>
    dplyr::inner_join(frozen, by = "name_original") |>
    dplyr::filter(accepted_species != stringr::word(legacy_accepted, 1, 2)) |>
    dplyr::select(name_original, legacy_accepted, accepted_species, decided_by, n_records, regions) |>
    dplyr::arrange(dplyr::desc(n_records))

  # Community edits
  edits <- community_edits |>
    dplyr::mutate(rows = vapply(edit_id, function(id) {
      if (field[edit_id == id] == "year") return(dropped$rows[dropped$reason == "year edit"])
      sum(community$drop_reason %in% paste0("community edit ", id))
    }, integer(1))) |>
    dplyr::select(edit_id, region, field, value, rows, legacy_source)

  unmapped <- community_named |>
    dplyr::filter(!in_dictionary) |>
    dplyr::count(region, destination_site, local_species_code, name = "rows")

  merge_types <- dplyr::count(species_merges, merge_type, name = "names")

  # TRY names flagged for review: which could belong to a study species?
  study <- unique(stats::na.omit(taxonomy$accepted_species[taxonomy$identification_level %in% "species_or_lower"]))
  try_rows <- dplyr::bind_rows(
    dplyr::transmute(try_names, name = SpeciesName, n = n_rows),
    dplyr::transmute(try_names, name = AccSpeciesName, n = n_rows)) |>
    dplyr::summarise(n = sum(n), .by = name)
  try_review <- name_resolution_try |>
    dplyr::filter(needs_review) |>
    dplyr::mutate(
      reason = sub(":.*$", "", gsub(" \\(score [0-9.]+\\)", "", review_reason)),
      acc_sp = dplyr::coalesce(accepted_species, stringr::word(accepted_name, 1, 2)),
      tied = ifelse(grepl("^tie: ", review_reason), sub("^tie: ", "", review_reason), ""),
      could_be_study = acc_sp %in% study |
        vapply(strsplit(tied, " | ", fixed = TRUE), function(x) any(stringr::word(x, 1, 2) %in% study), logical(1))) |>
    dplyr::left_join(try_rows, by = c(submitted_name = "name"))
  try_review_tab <- try_review |>
    dplyr::summarise(names = dplyr::n(), try_rows = sum(n, na.rm = TRUE), .by = c(reason, could_be_study)) |>
    dplyr::arrange(dplyr::desc(could_be_study), reason)
  exp_tab <- dplyr::select(experiments, experiment_id, first_year, last_year, year_range, included, exclusion_reason)

  report <- c(
    "# P1: Community and taxonomy", "",
    paste0("The prep pipeline read the network SQLite (", nrow(community_raw), " records, ",
           dplyr::n_distinct(community_raw$region), " regions) and the local-code dictionary of the RData, ",
           "resolved every community, TRY and Lynn name once with tidyTRY ", run$tidyTRY, " (TNRS ", run$TNRS,
           ", ", run$date, "; legitimate names only, WCVP before WFO), and wrote the community table with one row ",
           "per plot-year x accepted species for the analysis and every other vegetation record kept with its ",
           "`drop_reason`. No relative cover is computed here: step 01 computes it over all vegetation records of ",
           "a plot-year. ", nrow(experiments), " warmed comparisons, ", sum(experiments$included), " included."),
    "", "## Counts", "", md_table(steps),
    "", "## Reconciliation with the RData", "",
    paste0("Rows per region: SQLite ", nrow(community_raw), ", RData `dat` ", nrow(env$dat),
           ". Region-years where the two differ:"), "",
    if (nrow(recon_diff)) md_table(recon_diff) else "None.",
    "", "## How names were decided", "",
    "Every community name gets one of three origins (`decided_by`) and one identification level.", "",
    paste("- **rule**: chosen automatically by `tidyTRY::resolve_species()` (v", run$tidyTRY, ") from all TNRS",
          "candidates: only legitimate names (taxonomic status Accepted or Synonym, never Illegitimate, Invalid or",
          "No opinion); the highest match score; WCVP before WFO; Accepted before Synonym. A match into another genus",
          "is never chosen, and a name is flagged for review (not chosen automatically) when it comes only from WFO,",
          "matches with a score below 0.9, ties between different accepted names, or has no legitimate candidate."),
    "- **decision**: set by Bills in `prep/config/taxonomy_decisions.csv`, for every community name the rule flagged and for names the rule resolved to a plant from the wrong region (with the reason in the ledger).",
    "- **unidentified**: an unidentified field morphotaxon (`prep/config/unidentified_taxa.csv`, from the legacy overrides); its records stay as vegetation but never enter the analysis.",
    "",
    "Identification level (legacy rule, `legacy/R/data/community.R:252-285`): **species_or_lower** when the recorded and the accepted name both have at least two words; **genus_only** when either has one word or the recorded name ends in sp./spp./spec./indet.; **family_only** when the name ends in -aceae and has one word. Only species_or_lower names enter the analysis; the others keep their records with a `drop_reason`.",
    "",
    paste0("After the rule and the decisions, these accepted names are replaced in every table (`prep/config/accepted_renames.csv`): ",
           paste(sprintf("%s -> %s", accepted_renames$from_accepted_species, accepted_renames$to_accepted_species), collapse = "; "), "."),
    "",
    "**needs_review** in the table below counts the names the rule could not settle on its own (the reasons above). For community names every one of them is decided in the ledger before the pipeline continues; for TRY and Lynn names they are counted here and handled where those names are used (P2, P4).",
    "", "## Name resolution", "",
    paste("TRY names are resolved here, in the same run as the community names, because traits (P2) are attached to a",
          "study species through the TRY record's own name, not through TRY's `AccSpeciesName`, which merges species",
          "(TRY files Festuca pratensis records under Festuca ovina). One run, one rule and one TNRS date serve every join."),
    "",
    paste0("TRY names: ", n_try_all, " distinct names in the export; ", nrow(try_candidate_names),
           " in a genus of a study name or of its accepted name were resolved (the rest cannot attach to a study species ",
           "unless filed under a third genus)."), "", md_table(res_summary), "",
    "Community names by identification level and how they were decided:", "", md_table(levels_tab),
    "", sprintf(paste("TRY names flagged for review: %d. Those whose accepted name, or one of the tied names, is a study",
                      "species could carry trait records for the analysis (`could_be_study`); their handling is decided in P2."),
                nrow(try_review)), "", md_table(try_review_tab),
    "", sprintf("### Accepted species that differ from the frozen legacy taxonomy (%d names)", nrow(changes)), "",
    if (nrow(changes)) md_table(changes) else "None.",
    "", sprintf("### Names merged into one species (%d names in %d species; full table in `data/derived/prep/species_merges.csv`)",
                nrow(species_merges), dplyr::n_distinct(species_merges$accepted_species)), "",
    md_table(merge_types), "", md_table(species_merges),
    "", "## Community edits (legacy hard codes, kept and documented)", "", md_table(edits),
    "", "## Codes not in the dictionary (recorded name kept)", "",
    if (nrow(unmapped)) md_table(unmapped) else "None.",
    "", "## Experiments", "", md_table(exp_tab),
    "", "## Warnings captured", "",
    {
      w <- c(paste("study names:", attr(name_resolution_study, "warnings")),
             paste("TRY names:", attr(name_resolution_try, "warnings")))
      w <- w[!grepl(": $", w)]
      if (length(w)) paste("-", gsub("\n", " ", w)) else "None."
    },
    "", "## Deviations from the plan", "",
    paste("(1) The TNRS calls rerun only when their list of names changes (default targets cue), rather than never;",
          "with `never`, a changed name list would silently keep stale names. (2) TRY names are resolved in a second",
          "call limited to the genera of the study names and their accepted names, because the TRY export lists about",
          "67,000 species; TRY records filed under a third genus are not attached.",
          "(3) This run used tidyTRY", run$tidyTRY, "installed from the branch `fix/tnrs-id-check` (TNRS rows checked",
          "against their submitted name, rows of unrelated names dropped, no match into another genus, retries);",
          "`renv.lock` pins the release and moves to it once it is tagged.",
          "(4) FR_AlpeHuez::HIGH::LOW is included: the SQLite has its years 2015 to 2023, so the legacy reason",
          "\"warmed observations have missing Year\" no longer applies (decision at review).")
  )
  writeLines(report, path)
  path
}

write_trait_outputs <- function(trait_summary, dir = "data/derived/prep") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  files <- file.path(dir, c("traits.rds", "traits_long.rds"))
  saveRDS(trait_summary$traits, files[1])
  saveRDS(dplyr::filter(trait_summary$long, species %in% trait_summary$traits$species), files[2])
  files
}

write_p2_report <- function(try_ids, try_attached, try_clean, try_climate, trait_summary, try_name_map,
                            taxonomy, natural_values, frozen_traits_file, path = "reports/P2_traits.md") {
  traits <- trait_summary$traits
  long <- trait_summary$long
  study <- unique(stats::na.omit(taxonomy$accepted_species[taxonomy$identification_level %in% "species_or_lower"]))
  read <- attr(try_attached, "rows_read")
  counts <- dplyr::bind_rows(
    tibble::tibble(stage = sprintf("read by TRY species ID (%d IDs)", length(try_ids)),
                   rows = read[["rows"]], observations = read[["observations"]]),
    try_climate$counts)

  removed <- try_clean$removed_values
  removed_top <- utils::head(dplyr::arrange(removed, dplyr::desc(n_rows)), 40)

  per_trait <- tibble::tibble(
    trait = c("LDMC", "LMA", "height", "any of the three"),
    species = c(sum(!is.na(traits$LDMC)), sum(!is.na(traits$LMA)), sum(!is.na(traits$height)), nrow(traits)))
  missing <- sort(setdiff(study, traits$species))

  # Trusting TRY's AccSpeciesName against attaching by the record's own name
  ours <- long |>
    dplyr::filter(species %in% study) |>
    dplyr::summarise(n_ours = dplyr::n(), median_ours = stats::median(value), .by = c(species, trait))
  trust <- long |>
    dplyr::filter(species_if_trusting_try %in% study) |>
    dplyr::select(species = species_if_trusting_try, trait, value) |>
    dplyr::summarise(n_trusting_try = dplyr::n(), median_trusting_try = stats::median(value),
                     .by = c(species, trait))
  compare <- dplyr::full_join(ours, trust, by = c("species", "trait")) |>
    dplyr::mutate(relative_difference = round((median_trusting_try - median_ours) / median_ours, 3))
  differ <- compare |>
    dplyr::filter(dplyr::coalesce(n_ours, 0L) != dplyr::coalesce(n_trusting_try, 0L) |
                    dplyr::coalesce(abs(relative_difference) > 0.05, TRUE)) |>
    dplyr::arrange(dplyr::desc(abs(dplyr::coalesce(relative_difference, Inf))))
  focus <- dplyr::filter(compare, species %in% c("Festuca ovina", "Lolium pratense"))
  festuca_pratensis_in_ovina <- sum(long$species %in% "Festuca ovina" & long$SpeciesName %in% "Festuca pratensis")

  # Against the frozen legacy traits (step 00)
  frozen <- readRDS(frozen_traits_file) |>
    dplyr::select(species, LDMC_00 = LDMC, LMA_00 = LMA, height_00 = height)
  vs00 <- dplyr::full_join(dplyr::select(traits, species, LDMC, LMA, height), frozen, by = "species")
  gained <- sum(!traits$species %in% frozen$species)
  lost <- sum(!frozen$species %in% traits$species)
  changed <- dplyr::bind_rows(lapply(c("LDMC", "LMA", "height"), function(t) {
    tibble::tibble(species = vs00$species, trait = t, new = vs00[[t]], legacy = vs00[[paste0(t, "_00")]])
  })) |>
    dplyr::filter(!is.na(new), !is.na(legacy)) |>
    dplyr::mutate(relative_difference = round((new - legacy) / legacy, 3)) |>
    dplyr::filter(abs(relative_difference) > 0.05) |>
    dplyr::arrange(dplyr::desc(abs(relative_difference)))

  flagged <- dplyr::filter(try_name_map, study, needs_review) |>
    dplyr::count(used, reason = sub(":.*$", "", gsub(" \\(score [0-9.]+\\)", "", review_reason)), name = "names")

  report <- c(
    "# P2: TRY traits", "",
    paste0("TRY records were read by TRY species ID (", length(try_ids), " IDs whose contributor or TRY name resolves to a ",
           "study species in P1) and each record was attached to a species through its own resolved name (the contributor's ",
           "`SpeciesName` first, TRY's `AccSpeciesName` only as fallback). Observations under experimental or ambiguous ",
           "treatment values were removed, then duplicates, ErrorRisk >= 4, observations without coordinates and those in ",
           "arid and tropical climate zones, as in the legacy pipeline. ", nrow(traits), " of ", length(study),
           " study species have at least one trait. LDMC in g/g, LMA in mg mm-2 (1 / SLA per observation), height in m."),
    "", "## Counts", "", md_table(counts),
    "", "## Treatment values", "",
    paste0("Counted as natural (`prep/config/try_natural_values.csv`, compared after lower-casing and dropping trailing punctuation): ",
           paste(sprintf("\"%s\"", natural_values$value), collapse = ", "), "."), "",
    sprintf("Every other value marks the observation as experimental: %d distinct values on %d rows. The 40 most frequent:",
            nrow(removed), sum(removed$n_rows)), "",
    md_table(removed_top),
    "", "## Values not used", "",
    sprintf("Non-positive values dropped before the medians (decided 26 Oct): %d observations.",
            sum(trait_summary$non_positive$observations)), "",
    if (nrow(trait_summary$non_positive)) md_table(trait_summary$non_positive) else "None.", "",
    {
      rep <- attr(try_attached, "repaired_utf8")
      if (length(rep)) paste0("Text fields repaired from windows-1252 to UTF-8 after reading: ",
                              paste(sprintf("%s (%d values)", names(rep), rep), collapse = ", "), ".")
      else "No text field needed UTF-8 repair."
    },
    "", "## Species with traits", "", md_table(per_trait), "",
    sprintf("Study species without any of the three traits (%d): %s.", length(missing), paste(missing, collapse = ", ")),
    "", "## Trusting TRY's accepted name against attaching by the record's own name", "",
    sprintf(paste("Records of species x trait whose count or median differs between the two (median differing by more than 5%%):",
                  "%d of %d species x trait combinations. Festuca ovina now holds %d records submitted as Festuca pratensis."),
            nrow(differ), nrow(compare), festuca_pratensis_in_ovina), "",
    md_table(focus), "", "The 30 largest differences:", "", md_table(utils::head(differ, 30)),
    "", "## Against the frozen legacy traits (step 00)", "",
    sprintf("Species with traits now: %d; in step 00: %d; species gained: %d; lost: %d (names that changed count as both); species x trait medians changed by more than 5%%: %d. The 30 largest:",
            nrow(traits), nrow(frozen), gained, lost, nrow(changed)), "",
    md_table(utils::head(changed, 30)),
    "", "## Flagged TRY names of study species", "",
    "A flagged TRY name is used when its accepted species is a study species and the matched name has the submitted genus and epithet (decided 26 Oct).", "",
    md_table(flagged),
    "", "## Warnings captured", "",
    if (length(try_clean$warnings)) paste("-", gsub("\n", " ", unique(try_clean$warnings))) else "None.",
    "", "## Deviations from the plan", "",
    paste("(1) Community name decisions (e.g. Thlaspi montanum -> Noccaea fendleri for US_Montana) are not carried over to TRY names,",
          "because they concern the study records; TRY names are resolved by the rule alone, with the accepted-name renames of P1",
          "(Taraxacum officinale) applied. (2) ErrorRisk < 3 instead of the legacy < 4 (decided 26 Oct): with experiments removed per",
          "observation, datasets that legacy dropped wholesale now contribute, and some of their records 3 to 4 SD from TRY's species",
          "mean gave implausible medians (Viola canina height 1.16 m, Euphrasia officinalis LMA 0.43 mg mm-2). (3) Non-positive values",
          "dropped (decided 26 Oct). (4) TRY text fields repaired to UTF-8 after reading (legacy repaired names only).")
  )
  writeLines(report, path)
  path
}

write_p3_outputs <- function(functional_groups, phylogeny, phylogeny_axes, placements, dir = "data/derived/prep") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  files <- file.path(dir, c("functional_groups.rds", "phylogeny_axes.rds", "phylogeny_tree.rds", "phylogeny_placements.csv"))
  saveRDS(functional_groups, files[1])
  saveRDS(list(scores = phylogeny_axes$scores, eigenvalues = phylogeny_axes$eigenvalues), files[2])
  saveRDS(phylogeny$tree, files[3])
  utils::write.csv(dplyr::select(placements, species, placement, backbone_label, sp_genus, sp_family, source),
                   files[4], row.names = FALSE, na = "")
  files
}

write_p3_report <- function(study_species, functional_groups, placements, placement_draft_tbl, phylogeny,
                            phylogeny_axes, frozen_groups_file, frozen_axes_file,
                            path = "reports/P3_groups_phylogeny.md") {
  fg <- functional_groups
  eig <- phylogeny_axes$eigenvalues
  n_bs <- sum(eig$exceeds_broken_stick)
  woody_tab <- dplyr::select(dplyr::filter(fg, woody), species, family, woody_basis)
  ptab <- dplyr::count(placements, placement, name = "species")
  family_grafts <- placements$species[placements$placement == "grafted within family"]
  edited <- dplyr::inner_join(dplyr::select(placements, species, ledger = backbone_label, ledger_placement = placement),
                              dplyr::select(placement_draft_tbl, species, draft = backbone_label, draft_placement = placement),
                              by = "species") |>
    dplyr::filter(!mapply(identical, ledger, draft) | ledger_placement != draft_placement)

  # against step 00
  frozen_fg <- readRDS(frozen_groups_file)
  fg_change <- dplyr::inner_join(dplyr::select(fg, species, graminoid, woody),
                                 dplyr::select(frozen_fg, species, graminoid_00 = graminoid, woody_00 = woody), by = "species") |>
    dplyr::filter(graminoid != graminoid_00 | woody != woody_00)
  frozen_ax <- readRDS(frozen_axes_file)$scores
  common <- intersect(phylogeny_axes$scores$species, frozen_ax$species)
  cors <- vapply(paste0("PEV", 1:3), function(ax) {
    stats::cor(phylogeny_axes$scores[[ax]][match(common, phylogeny_axes$scores$species)],
               frozen_ax[[ax]][match(common, frozen_ax$species)])
  }, numeric(1))

  report <- c(
    "# P3: Functional groups and phylogeny", "",
    paste0("Graminoid and woody status and a phylogeny with its PCoA axes for all ", nrow(study_species),
           " study species (residents and destination pool), on the P1 accepted species. The tree is cut from the ",
           "V.PhyloMaker2 mega-tree GBOTB.extended.LCVP (scenario S3); each species' place is set by the reviewed ",
           "ledger `prep/config/phylogeny_placements.csv`, which reuses the legacy hand-checked labels for the same species ",
           "(and for renamed species whose epithet is unchanged), and places new species by rule."),
    "", "## Functional groups", "",
    sprintf("Graminoid (Poaceae, Cyperaceae, Juncaceae): %d species. Woody (`prep/config/woodiness.csv`): %d species. Family taken from the backbone genus (no family in P1): %s.",
            sum(fg$graminoid), sum(fg$woody), paste(study_species$species[study_species$family_source == "backbone genus"], collapse = ", ")),
    "", md_table(woody_tab),
    "", sprintf("Species whose graminoid or woody status differs from step 00 (%d):", nrow(fg_change)), "",
    if (nrow(fg_change)) md_table(fg_change) else "None.",
    "", "## Placements", "", md_table(ptab), "",
    sprintf("Grafted at family level, their genus being absent from the mega-tree (%d): %s.", length(family_grafts), paste(family_grafts, collapse = ", ")),
    "", sprintf("Ledger rows edited by hand relative to the rule-based draft: %d.", nrow(edited)),
    if (nrow(edited)) c("", md_table(edited)) else character(),
    "", "## Tree and axes", "",
    sprintf("Tips: %d, rooted and ultrametric. PCoA on cophenetic distances: %d positive axes; %s (negative/positive eigenvalue ratio %.2g). %d axes exceed the broken-stick expectation; PEV1 to PEV3 carry %.1f%% of the positive variation. Orientation anchors: Poa trivialis +PEV1, Erigeron formosissimus -PEV2, Euphrasia minima +PEV3.",
            length(phylogeny$tree$tip.label), nrow(eig), if (phylogeny_axes$corrected) "Lingoes correction" else "no correction needed",
            phylogeny_axes$negative_ratio, n_bs, eig$cumulative_percent[3]),
    "", md_table(dplyr::mutate(utils::head(eig, 6), dplyr::across(dplyr::where(is.numeric), ~ round(.x, 2)))),
    "", sprintf("Against step 00 (%d resident species in both): correlation of scores PEV1 %.3f, PEV2 %.3f, PEV3 %.3f.",
                length(common), cors[1], cors[2], cors[3]),
    "", "## Warnings captured", "",
    if (length(phylogeny$warnings)) paste("-", unique(phylogeny$warnings)) else "None.",
    "", "## Deviations from the plan", "",
    paste("The genus map gained three rows found while checking the draft (Eriocoma -> Stipa instead of Achnatherum, Ptilagrostis -> Stipa,",
          "Cherleria -> Minuartia): V.PhyloMaker2 binds only into genera listed in its node table, and Achnatherum, Ptilagrostis and",
          "Cherleria have tips but no node there.")
  )
  writeLines(report, path)
  path
}

write_p4_outputs <- function(experiment_climate, niche_optima, site_bio10, dir = "data/derived/prep") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  files <- file.path(dir, c("experiment_climate.rds", "niche_optima.rds", "site_bio10_check.csv"))
  saveRDS(experiment_climate, files[1])
  saveRDS(niche_optima, files[2])
  utils::write.csv(site_bio10, files[3], row.names = FALSE)
  files
}

write_p4_report <- function(experiment_climate, experiments, site_bio10, niche_optima, lynn, lynn_disagree,
                            study_species, frozen_climate_file, frozen_niche_file,
                            path = "reports/P4_climate_niche.md") {
  ec <- experiment_climate
  r_P <- stats::cor(ec$P_annual, ec$P_summer)
  clim_summary <- tibble::tibble(
    variable = c("T (destination JJA, degC)", "W (destination - origin JJA, degC)", "P_summer (log mm)", "P_annual (log mm)"),
    min = c(min(ec$T), min(ec$W), min(ec$P_summer), min(ec$P_annual)),
    median = c(stats::median(ec$T), stats::median(ec$W), stats::median(ec$P_summer), stats::median(ec$P_annual)),
    max = c(max(ec$T), max(ec$W), max(ec$P_summer), max(ec$P_annual))) |>
    dplyr::mutate(dplyr::across(c(min, median, max), ~ round(.x, 3)))
  partial <- dplyr::filter(ec, !climate_complete) |>
    dplyr::select(experiment_id, destination_years_available, destination_years_expected, missing_years = destination_missing_years)

  frozen <- readRDS(frozen_climate_file)
  cmp <- dplyr::inner_join(dplyr::select(ec, experiment_id, T, W, P_summer, P_annual, n_years = destination_years_expected),
                           dplyr::select(frozen, experiment_id, T_00 = T, W_00 = W, P_summer_00 = P_summer, P_annual_00 = P_annual,
                                         n_years_00 = destination_years_expected), by = "experiment_id") |>
    dplyr::mutate(max_abs_difference = pmax(abs(T - T_00), abs(W - W_00), abs(P_summer - P_summer_00), abs(P_annual - P_annual_00)))
  unchanged <- dplyr::filter(cmp, n_years == n_years_00)
  changed <- dplyr::filter(cmp, n_years != n_years_00) |>
    dplyr::mutate(dplyr::across(c(T, W, P_summer, P_annual, T_00, W_00, P_summer_00, P_annual_00), ~ round(.x, 3)))
  new_exp <- setdiff(ec$experiment_id, frozen$experiment_id)
  if (any(unchanged$max_abs_difference > 1e-9)) stop("Climate of an unchanged experiment differs from step 00")

  r_bio10 <- stats::cor(site_bio10$T, site_bio10$bio10_site, use = "complete.obs")

  merged <- lynn |>
    dplyr::filter(accepted_species %in% niche_optima$species[niche_optima$n_lynn_names > 1], !is.na(mean_sum_sub)) |>
    dplyr::left_join(dplyr::select(niche_optima, accepted_species = species, lynn_name_used), by = "accepted_species") |>
    dplyr::transmute(species = accepted_species, lynn_name, regions, n_obs_sub, T_opt = round(mean_sum_sub, 2),
                     AP_mm = round(mean_ap_sub), used = lynn_name == lynn_name_used) |>
    dplyr::arrange(species, dplyr::desc(n_obs_sub))
  unused_lynn <- lynn |>
    dplyr::filter(!lynn_name %in% c(niche_optima$lynn_name_used, unlist(strsplit(stats::na.omit(niche_optima$lynn_names_merged), "; ")))) |>
    dplyr::mutate(reason = dplyr::case_when(is.na(mean_sum_sub) ~ "Lynn's niche values are NA",
                                            is.na(accepted_species) ~ "name not resolved in P1",
                                            TRUE ~ "accepted species is not a study species")) |>
    dplyr::count(reason, name = "lynn_names")
  frozen_niche <- readRDS(frozen_niche_file)

  report <- c(
    "# P4: Climate and niche", "",
    paste0("Experiment climate for the ", nrow(ec), " included experiments from the cached CHELSA V2.1 monthly files ",
           "(2008 to 2021; read, not downloaded): T = destination June to August mean temperature (day-weighted), ",
           "W = destination minus origin JJA temperature, P_summer and P_annual = log destination JJA and annual precipitation, ",
           "each averaged over the experiment's years. Lynn's niche optima joined to the P1 accepted species by Lynn's name: ",
           "T_opt = mean_sum_sub, P_opt = log(mean_ap_sub)."),
    "", "## Climate", "", md_table(clim_summary), "",
    sprintf("Correlation of P_annual and P_summer across the %d included experiments: r = %.3f (PLAN.md: about 0.78 on the earlier 41).",
            nrow(ec), r_P),
    "", sprintf("Partial coverage (%d experiments; CHELSA has no 2022 to 2023), averaged over the years available:", nrow(partial)), "",
    md_table(partial),
    "", sprintf(paste("Against step 00: %d experiments with the same years agree exactly (largest difference %.2g).",
                      "Experiments whose years changed: %d. New experiments: %s."),
                nrow(unchanged), max(unchanged$max_abs_difference), nrow(changed),
                if (length(new_exp)) paste(new_exp, collapse = ", ") else "none"), "",
    if (nrow(changed)) md_table(dplyr::select(changed, experiment_id, n_years_00, n_years, T_00, T, W_00, W, P_annual_00, P_annual)) else character(),
    "", "### Summer temperature: CHELSA 1981-2010 climatology against the experiment JJA T", "",
    sprintf(paste("Lynn's site summer temperature (`sum_temp_site`, the bio10 analogue behind T_opt) against our destination JJA T",
                  "averaged over the experiment years: r = %.3f; mean difference (ours minus climatology) %.2f degC, range %.2f to %.2f.",
                  "Table: `data/derived/prep/site_bio10_check.csv`."),
            r_bio10, mean(site_bio10$difference, na.rm = TRUE), min(site_bio10$difference, na.rm = TRUE),
            max(site_bio10$difference, na.rm = TRUE)),
    "", "## Niche", "",
    sprintf(paste("Study species with a niche: %d of %d; with n_obs_sub >= 30: %d. Lynn names: %d, of which %d used",
                  "(%d species from more than one Lynn name). Step 00 had %d Lynn names, keyed by Lynn's own names."),
            nrow(niche_optima), nrow(study_species), sum(niche_optima$n_obs_sub >= 30), nrow(lynn),
            nrow(lynn) - sum(unused_lynn$lynn_names), sum(niche_optima$n_lynn_names > 1), nrow(frozen_niche)), "",
    "Lynn names not used:", "", md_table(unused_lynn),
    "", sprintf("### Species with more than one Lynn name (%d): the name with most occurrences is used (decided 26 Oct)",
                sum(niche_optima$n_lynn_names > 1)), "",
    paste("The second name is a synonym or an infraspecific name recorded in one region; Lynn computed a niche per recorded name,",
          "so the same plant has two samples of its range. The values differ by little except for Solidago virgaurea, whose",
          "alpine subsp. minuta (recorded in CH_Lavey) has a niche about 5 degC colder than the species, which is used."), "",
    md_table(merged),
    "", sprintf("### Rows where Lynn's own code-to-name mapping disagrees with ours (%d)", nrow(lynn_disagree)), "",
    "Lynn's niche for these rows belongs to the name Lynn assigned, not to our species; joining by name keeps it away from our species.", "",
    md_table(dplyr::select(lynn_disagree, region, local_species_code, lynn_name, lynn_accepted, our_species)),
    "", "## Warnings captured", "", "None.",
    "", "## Deviations from the plan", "",
    "Alopecurus gerardii -> Alopecurus gerardi added to `prep/config/accepted_renames.csv` (decided 26 Oct), so that Lynn's and TRY's records of the species join it."
  )
  writeLines(report, path)
  path
}
