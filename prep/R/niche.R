# Niche optima per study species (P4) from P6's GBIF niche
# (prep/gbif_niche.R): T_opt = mean CHELSA 1981-2010 bio10 and P_opt = log
# mean bio18 over one value per species x CHELSA cell. Replaces Lynn's niche
# file (decided 28 Oct).

# The GBIF niche is computed outside the pipeline for one species list; if
# the taxonomy has changed since, the script is rerun deliberately
build_niche_optima <- function(gbif_niche_file, study_species) {
  niche <- readRDS(gbif_niche_file)
  if (!setequal(niche$species, study_species$species) || anyDuplicated(niche$species)) {
    stop("data/derived/prep/gbif_niche.rds does not hold exactly the study species; rerun prep/gbif_niche.R")
  }
  niche |>
    dplyr::filter(n_obs_sub > 0) |>
    dplyr::select(species, T_opt, P_opt, n_obs_sub, n_records)
}
