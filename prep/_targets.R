# Data preparation pipeline (prep only; the analysis is in analysis/).
# Run from the repository root:
#   Rscript -e 'targets::tar_make(script = "prep/_targets.R", store = "prep/_targets")'

library(targets)
tar_option_set(packages = c("dplyr", "tibble", "stringr"))
tar_source("prep/R")

read_ledger <- function(path) utils::read.csv(path, colClasses = "character", encoding = "UTF-8")

list(
  # inputs
  tar_target(sqlite_file, "data/transplant_network_clean.sqlite", format = "file"),
  tar_target(rdata_file, "data/TransPlantNetwork_060824.RData", format = "file"),
  tar_target(try_file, "data/try_47163.txt", format = "file"),
  tar_target(niche_file, "data/MacroContextTransplantDat_2.csv", format = "file"),
  tar_target(frozen_taxonomy_file, "data/derived/00_inputs/species_taxonomy.rds", format = "file"),
  tar_target(site_edits_file, "prep/config/site_edits.csv", format = "file"),
  tar_target(community_edits_file, "prep/config/community_edits.csv", format = "file"),
  tar_target(mapping_decisions_file, "prep/config/species_mapping_decisions.csv", format = "file"),
  tar_target(unidentified_file, "prep/config/unidentified_taxa.csv", format = "file"),
  tar_target(taxonomy_decisions_file, "prep/config/taxonomy_decisions.csv", format = "file"),
  tar_target(site_edits, read_ledger(site_edits_file)),
  tar_target(community_edits, read_ledger(community_edits_file)),
  tar_target(mapping_decisions, read_ledger(mapping_decisions_file)),
  tar_target(unidentified, read_ledger(unidentified_file)),
  tar_target(taxonomy_decisions, read_ledger(taxonomy_decisions_file)),

  # P1: community and taxonomy
  tar_target(community_raw, read_community(sqlite_file)),
  tar_target(sites, read_sites(sqlite_file, site_edits)),
  tar_target(dictionary, code_dictionary(rdata_file, mapping_decisions)),
  tar_target(community_named, name_community(community_raw, dictionary)),
  tar_target(try_names, tidyTRY::try_species_names(try_file, progress = FALSE)),
  tar_target(lynn_names, unique(utils::read.csv(niche_file)$species)),
  tar_target(study_names, names_study(community_named, unidentified, lynn_names)),
  tar_target(name_resolution_study, resolve_names(study_names)),
  tar_target(try_candidate_names, names_try(try_names, name_resolution_study)),
  tar_target(name_resolution_try, resolve_names(try_candidate_names)),
  tar_target(names_all, dplyr::bind_rows(study_names, try_candidate_names)),
  tar_target(name_resolution, dplyr::bind_rows(name_resolution_study, name_resolution_try)),
  tar_target(taxonomy, build_taxonomy(community_named, name_resolution_study, taxonomy_decisions, unidentified)),
  tar_target(community_dated, date_community(community_named, sites, community_edits)),
  tar_target(community, finish_community(community_dated, taxonomy, community_edits)),
  tar_target(experiments, build_experiments(community, sites)),
  tar_target(species_merges, species_merge_table(taxonomy)),
  tar_target(prep_outputs,
             write_prep_outputs(community, sites, experiments, taxonomy, species_merges),
             format = "file"),
  tar_target(p1_report,
             write_p1_report(community_raw, community_named, community_dated, community, sites,
                             experiments, taxonomy, name_resolution_study, name_resolution_try, names_all,
                             name_resolution, try_candidate_names, species_merges,
                             community_edits, rdata_file, frozen_taxonomy_file),
             format = "file")
)
