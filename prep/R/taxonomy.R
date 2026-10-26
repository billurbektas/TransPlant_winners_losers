# Taxonomy (P1): one resolution run over the community, TRY and Lynn names,
# with tidyTRY >= 0.2.0 (legitimate TNRS names only, WCVP before WFO).

# legacy/R/infrastructure/io.R:39-44
repair_utf8 <- function(x) {
  x <- as.character(x)
  invalid <- !is.na(x) & !stringi::stri_enc_isutf8(x)
  x[invalid] <- stringi::stri_encode(x[invalid], from = "windows-1252", to = "UTF-8")
  x
}

word_count <- function(x) {
  x <- stringr::str_squish(x)
  ifelse(is.na(x) | x == "", NA_integer_, stringr::str_count(x, stringr::fixed(" ")) + 1L)
}

unspecified_species <- function(x) {
  stringr::str_detect(x, stringr::regex("(^|[[:space:]_])(sp|spp|spec|indet)[.]?$", ignore_case = TRUE))
}

# Study names: community names recorded at species level or below (a
# one-word or "sp." name stays genus-level whatever TNRS returns) and Lynn's
# niche names
names_study <- function(community_named, unidentified, lynn_names) {
  comm <- unique(stats::na.omit(community_named$name_original))
  comm <- comm[!comm %in% unidentified$name_original & word_count(comm) >= 2L & !unspecified_species(comm)]
  all <- unique(stringr::str_squish(c(comm, lynn_names)))
  all <- all[!is.na(all) & nzchar(all)]
  tibble::tibble(name = all, in_community = all %in% comm, in_lynn = all %in% lynn_names, in_try = FALSE)
}

# TRY names worth resolving. The export lists every TRY species (about
# 67,000), so only names in a genus of a study name or of its accepted name
# are resolved: this keeps TRY synonyms within a genus and known genus
# changes (TRY "Festuca pratensis" and "Lolium pratense" for a study
# "Festuca pratensis"), and leaves out records filed under a third genus.
names_try <- function(try_names, name_resolution_study) {
  genera <- unique(stringr::word(c(name_resolution_study$submitted_name,
                                   stats::na.omit(name_resolution_study$accepted_name)), 1))
  try_all <- unique(stringr::str_squish(repair_utf8(c(try_names$SpeciesName, try_names$AccSpeciesName))))
  try_all <- try_all[!is.na(try_all) & nzchar(try_all)]
  keep <- tolower(stringr::word(try_all, 1)) %in% tolower(genera)
  cand <- setdiff(try_all[keep], name_resolution_study$submitted_name)
  out <- tibble::tibble(name = cand, in_community = FALSE, in_lynn = FALSE, in_try = TRUE)
  attr(out, "n_try_names") <- length(try_all)
  out
}

# The TNRS calls of the pipeline (study names, then TRY names); each reruns
# only when its list of names changes. The date and package versions travel
# with the result.
# All candidates for 5,000 names in one request fails on the TNRS server
# (HTTP 500), hence smaller batches.
resolve_names <- function(names_all, batch_size = 1000L) {
  # warnings (e.g. TNRS rows of unrelated names dropped) go to the report
  warns <- character()
  res <- withCallingHandlers(
    tidyTRY::resolve_species(names_all$name, method = "tnrs", sources = c("wcvp", "wfo"),
                             batch_size = batch_size),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  attr(res, "warnings") <- warns
  attr(res, "resolved") <- list(date = as.character(Sys.Date()),
                                tidyTRY = as.character(utils::packageVersion("tidyTRY")),
                                TNRS = as.character(utils::packageVersion("TNRS")))
  res
}

# One row per community name: accepted name and species, family, status and
# how the name was decided. Identification level follows
# legacy/R/data/community.R:252-285.
build_taxonomy <- function(community_named, name_resolution, decisions, unidentified,
                           unresolved_path = "prep/unresolved_names.csv") {
  comm <- community_named |>
    dplyr::filter(!is.na(name_original)) |>
    dplyr::summarise(n_records = dplyr::n(), regions = paste(sort(unique(region)), collapse = "; "),
                     .by = name_original)
  tax <- comm |>
    dplyr::left_join(dplyr::select(name_resolution, name_original = submitted_name, accepted_name,
                                   accepted_species, accepted_family, taxonomic_status, source,
                                   needs_review, review_reason, n_legitimate),
                     by = "name_original", relationship = "one-to-one") |>
    dplyr::mutate(
      unidentified = name_original %in% unidentified$name_original,
      # species from the accepted name when TNRS gives none (e.g. WFO matches)
      accepted_species = dplyr::coalesce(
        accepted_species,
        dplyr::if_else(word_count(accepted_name) >= 2L, stringr::word(accepted_name, 1, 2), NA_character_)),
      decided_by = dplyr::if_else(unidentified, "unidentified", "rule"))

  bad <- setdiff(decisions$name_original, tax$name_original)
  if (length(bad)) stop("taxonomy_decisions.csv names not in the community: ", paste(bad, collapse = ", "))
  dec <- match(tax$name_original, decisions$name_original)
  has_dec <- !is.na(dec)
  tax$accepted_species[has_dec] <- ifelse(decisions$decision[dec[has_dec]] == "accept",
                                          decisions$accepted_species[dec[has_dec]], NA_character_)
  tax$accepted_name[has_dec] <- tax$accepted_species[has_dec]
  tax$decided_by[has_dec] <- "decision"
  tax$needs_review[has_dec] <- FALSE

  tax <- tax |>
    dplyr::mutate(
      original_words = word_count(name_original),
      accepted_words = word_count(accepted_name),
      identification_level = dplyr::case_when(
        unidentified ~ NA_character_,
        stringr::str_detect(dplyr::coalesce(accepted_name, name_original),
                            stringr::regex("aceae$", ignore_case = TRUE)) &
          dplyr::coalesce(accepted_words, original_words) == 1L ~ "family_only",
        original_words == 1L | accepted_words %in% 1L | unspecified_species(name_original) ~ "genus_only",
        decided_by == "decision" & is.na(accepted_name) ~ "excluded by decision",
        is.na(accepted_name) ~ "unresolved",
        original_words >= 2L & accepted_words >= 2L ~ "species_or_lower",
        TRUE ~ "unresolved rank"))

  open <- tax |>
    dplyr::filter(!unidentified, decided_by != "decision",
                  original_words >= 2L, !unspecified_species(name_original),
                  is.na(accepted_name) | dplyr::coalesce(needs_review, TRUE)) |>
    dplyr::select(name_original, regions, n_records, review_reason, n_legitimate, accepted_name, source)
  if (nrow(open)) {
    utils::write.csv(open, unresolved_path, row.names = FALSE)
    stop(nrow(open), " community names need a decision; listed in ", unresolved_path,
         ". Add them to prep/config/taxonomy_decisions.csv.")
  }
  if (file.exists(unresolved_path)) file.remove(unresolved_path)

  dplyr::select(tax, name_original, accepted_name, accepted_species, accepted_family, taxonomic_status,
                source, decided_by, identification_level, needs_review_reason = review_reason,
                n_records, regions)
}

# Submitted names that end up as one species, and why
species_merge_table <- function(taxonomy) {
  taxonomy |>
    dplyr::filter(identification_level == "species_or_lower") |>
    dplyr::filter(dplyr::n() > 1L, .by = accepted_species) |>
    dplyr::mutate(merge_type = dplyr::case_when(
      name_original == accepted_species ~ "accepted name",
      accepted_name != accepted_species ~ "infraspecific to species",
      diag(utils::adist(name_original, accepted_species)) <= 2 ~ "spelling",
      TRUE ~ "synonym")) |>
    dplyr::select(accepted_species, name_original, accepted_name, merge_type, regions, n_records) |>
    dplyr::arrange(accepted_species, name_original)
}
