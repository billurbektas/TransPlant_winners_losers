# Study species and functional groups (P3)

# One row per study species (accepted species at species rank). A family
# missing in P1 (names set by decision) is taken from the backbone's genus.
build_study_species <- function(taxonomy) {
  sp <- taxonomy |>
    dplyr::filter(identification_level %in% "species_or_lower") |>
    dplyr::distinct(species = accepted_species, family = accepted_family) |>
    dplyr::filter(!is.na(species)) |>
    dplyr::mutate(genus = stringr::word(species, 1))
  # a species can carry its family on one name and NA on another
  sp <- dplyr::summarise(sp, family = dplyr::first(stats::na.omit(family)), .by = c(species, genus))
  nodes <- V.PhyloMaker2::nodes.info.1.LCVP
  genus_family <- unique(nodes[nodes$level == "G", c("genus", "family")])
  from_backbone <- genus_family$family[match(sp$genus, genus_family$genus)]
  sp$family_source <- ifelse(is.na(sp$family) & !is.na(from_backbone), "backbone genus", "P1 taxonomy")
  sp$family <- dplyr::coalesce(sp$family, from_backbone)
  if (anyNA(sp$family)) stop("Study species without a family: ", paste(sp$species[is.na(sp$family)], collapse = ", "))
  if (anyDuplicated(sp$species)) stop("Duplicate study species")
  dplyr::arrange(sp, species)
}

# graminoid by family (legacy/R/data/analysis_dataset.R:33); woody from
# prep/config/woodiness.csv, a species row overriding its genus row
build_functional_groups <- function(study_species, woodiness) {
  woodiness$woody <- as.logical(woodiness$woody)
  bad <- woodiness |>
    dplyr::filter((level == "genus" & !name %in% study_species$genus) |
                    (level == "species" & !name %in% study_species$species))
  if (nrow(bad)) stop("woodiness.csv names no study taxon: ", paste(bad$name, collapse = ", "))
  by_species <- dplyr::filter(woodiness, level == "species")
  by_genus <- dplyr::filter(woodiness, level == "genus")
  out <- study_species |>
    dplyr::mutate(
      graminoid = family %in% c("Poaceae", "Cyperaceae", "Juncaceae"),
      woody_species = by_species$woody[match(species, by_species$name)],
      woody_genus = by_genus$woody[match(genus, by_genus$name)],
      woody = dplyr::coalesce(woody_species, woody_genus, FALSE),
      woody_basis = dplyr::case_when(
        !is.na(woody_species) ~ paste("species row:", by_species$reason[match(species, by_species$name)]),
        !is.na(woody_genus) ~ paste("genus row:", by_genus$reason[match(genus, by_genus$name)]),
        TRUE ~ "not in woodiness.csv: non-woody")) |>
    dplyr::select(species, family, genus, graminoid, woody, woody_basis)
  if (anyNA(out$graminoid) || anyNA(out$woody)) stop("Missing graminoid or woody status")
  out
}
