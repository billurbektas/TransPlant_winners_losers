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
