# Lynn's niche optima per accepted species (P4). Lynn's names are joined
# through the P1 resolutions by name, never through Lynn's own code-to-name
# mapping, which repeats the homonym errors (Festuca pratensis -> Festuca
# ovina). T_opt = mean_sum_sub (summer temperature, bio10 analogue);
# P_opt = log(mean_ap_sub).

lynn_table <- function(niche_file, name_resolution) {
  raw <- utils::read.csv(niche_file)
  regions <- dplyr::summarise(raw, regions = paste(sort(unique(Region)), collapse = "; "), .by = species)
  raw |>
    dplyr::distinct(lynn_name = species, mean_sum_sub, mean_ap_sub, mean_mat_sub, n_obs_sub) |>
    dplyr::left_join(regions, by = c(lynn_name = "species")) |>
    dplyr::mutate(key = stringr::str_squish(lynn_name),
                  accepted_species = dplyr::coalesce(
                    name_resolution$accepted_species,
                    stringr::word(name_resolution$accepted_name, 1, 2))[match(key, name_resolution$submitted_name)])
}

# One niche per study species; when several Lynn names map to one species
# the name with the most occurrences is used (decided 26 Oct) and all of
# them are recorded
build_niche_optima <- function(lynn, study_species) {
  usable <- lynn |>
    dplyr::filter(accepted_species %in% study_species$species, !is.na(mean_sum_sub), !is.na(mean_ap_sub)) |>
    dplyr::arrange(accepted_species, dplyr::desc(n_obs_sub))
  out <- usable |>
    dplyr::summarise(
      T_opt = dplyr::first(mean_sum_sub), P_opt = log(dplyr::first(mean_ap_sub)),
      mean_mat_sub = dplyr::first(mean_mat_sub), n_obs_sub = dplyr::first(n_obs_sub),
      lynn_name_used = dplyr::first(lynn_name),
      lynn_names_merged = if (dplyr::n() > 1) paste(lynn_name, collapse = "; ") else NA_character_,
      n_lynn_names = dplyr::n(),
      .by = accepted_species) |>
    dplyr::rename(species = accepted_species)
  if (anyDuplicated(out$species)) stop("More than one niche row for a species")
  out
}

# Lynn's own region x code -> species mapping against ours, for the report
lynn_code_disagreements <- function(niche_file, lynn, community_named, taxonomy) {
  ours <- community_named |>
    dplyr::filter(!is.na(name_original)) |>
    dplyr::distinct(region, local_species_code, name_original) |>
    dplyr::left_join(dplyr::select(taxonomy, name_original, our_species = accepted_species), by = "name_original") |>
    dplyr::distinct(region, local_species_code, our_species)
  utils::read.csv(niche_file) |>
    dplyr::distinct(region = Region, local_species_code = SpeciesName, lynn_name = species) |>
    dplyr::left_join(dplyr::select(lynn, lynn_name, lynn_accepted = accepted_species), by = "lynn_name") |>
    dplyr::inner_join(ours, by = c("region", "local_species_code"), relationship = "many-to-many") |>
    dplyr::filter(!is.na(our_species), !mapply(identical, lynn_accepted, our_species))
}
