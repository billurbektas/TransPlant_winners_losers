# 00_freeze_inputs.R
#
# PLAN.md section: Step A, "Data rules"; Appendix, Fig S1.
# Exports the validated data preparation of the legacy targets pipeline once,
# as plain files, so that no later step depends on `targets`.
#
# Reads
#   _targets/ (legacy store, read only): experiments, sites, experiment_climate,
#     analysis_community, community_taxonomy_species, try_trait_outputs,
#     functional_groups, resident_pcoa_candidate
#   data/MacroContextTransplantDat_2.csv (J. Lynn, 11.2024; niche summaries)
#
# Writes, in data/derived/00_inputs/
#   experiments.rds, sites.rds, experiment_climate.rds, community.rds,
#   species_taxonomy.rds, traits.rds, functional_groups.rds,
#   phylogeny_axes.rds, niche_optima.rds, manifest.csv
#   and reports/00_freeze_inputs.md
#
# Expected counts checked
#   experiments: 44 rows, 42 included, all with origin and destination coordinates
#   experiment_climate: 42 rows, one per included experiment, T W P finite
#   phylogeny: one score column per eigenvalue, all eigenvalues positive
#   every output has a row in manifest.csv

library(targets)
library(dplyr)
library(tidyr)

store <- "_targets"
out_dir <- "data/derived/00_inputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
niche_csv <- "data/MacroContextTransplantDat_2.csv"

# warnings are collected for the report, never suppressed
warns <- character()
collect <- function(expr) withCallingHandlers(expr, warning = function(w) {
  warns <<- c(warns, conditionMessage(w))
  invokeRestart("muffleWarning")
})
read_target <- function(name) collect(tar_read_raw(name, store = store))

counts <- tibble(object = character(), step = character(), rows = integer())
count_row <- function(object, step, rows) {
  counts <<- add_row(counts, object = object, step = step, rows = as.integer(rows))
}

# ---- experiments and sites -------------------------------------------------
sites <- read_target("sites")
stopifnot(!anyDuplicated(sites[, c("region", "site_id")]))  # site key is region x site_id
count_row("sites", "read and written", nrow(sites))

experiments <- read_target("experiments")
count_row("experiments", "read", nrow(experiments))
# the exclusion reason is carried for Fig S1; every excluded experiment must have one
stopifnot("exclusion_reason" %in% names(experiments),
          !anyNA(experiments$exclusion_reason[!experiments$included]))
experiments <- experiments |>
  left_join(select(sites, region, origin_site = site_id,
                   origin_longitude = longitude, origin_latitude = latitude),
            by = c("region", "origin_site"), relationship = "many-to-one") |>
  left_join(select(sites, region, destination_site = site_id,
                   destination_longitude = longitude, destination_latitude = latitude),
            by = c("region", "destination_site"), relationship = "many-to-one")
# all 44 legacy rows are kept; step 01 filters to the included ones for Fig S1
stopifnot(nrow(experiments) == 44, sum(experiments$included) == 42,
          !anyDuplicated(experiments$experiment_id))
# every experiment has origin and destination site coordinates
stopifnot(all(is.finite(c(experiments$origin_longitude, experiments$origin_latitude,
                          experiments$destination_longitude, experiments$destination_latitude))))
count_row("experiments", "written (all rows)", nrow(experiments))
count_row("experiments", "of which included", sum(experiments$included))

# ---- experiment climate ------------------------------------------------------
# destination_P_annual is the legacy sum of monthly CHELSA precipitation per
# year, averaged over the experiment's years (legacy/R/data/climate.R:346)
experiment_climate <- read_target("experiment_climate") |>
  mutate(T = destination_T_JJA,
         W = diff_T_JJA,
         P_summer = log(destination_P_JJA),
         P_annual = log(destination_P_annual)) |>
  select(experiment_id, T, W, P_summer, P_annual,
         destination_T_JJA, origin_T_JJA, diff_T_JJA,
         destination_P_JJA, destination_P_annual,
         origin_years_expected, origin_years_available, origin_missing_years,
         destination_years_expected, destination_years_available,
         destination_missing_years, climate_complete, climate_coverage_note)
# one row per included experiment
stopifnot(nrow(experiment_climate) == 42, !anyDuplicated(experiment_climate$experiment_id),
          setequal(experiment_climate$experiment_id,
                   experiments$experiment_id[experiments$included]))
stopifnot(all(is.finite(c(experiment_climate$T, experiment_climate$W,
                          experiment_climate$P_summer, experiment_climate$P_annual))))
count_row("experiment_climate", "read and written", nrow(experiment_climate))
r_P <- cor(experiment_climate$P_annual, experiment_climate$P_summer)
partial_climate <- filter(experiment_climate, !climate_complete) |>
  select(experiment_id, origin_years_available, origin_years_expected,
         destination_years_available, destination_years_expected,
         missing_years = destination_missing_years)

# ---- community ---------------------------------------------------------------
community <- read_target("analysis_community")
count_row("community", "read", nrow(community))
# legacy relative_cover was computed before the deep-lineage removal; step 01
# recomputes relative_cover from cover, so the legacy column must not travel
community <- select(community, -relative_cover)
# columns step 01 builds cover_long from
stopifnot(all(c("region", "origin_site", "destination_site", "treatment", "plot_year_id",
                "stable_plot_id", "year", "relative_year", "species", "cover",
                "supplied_relative_cover") %in% names(community)))
stopifnot(!anyNA(community$species))
# US_Arizona supplies relative cover only, so cover is NA there and
# supplied_relative_cover carries the record; nowhere else may cover be missing
na_cover <- filter(community, is.na(cover))
stopifnot(all(na_cover$region == "US_Arizona"), !anyNA(na_cover$supplied_relative_cover))
# year is missing only in FR_AlpeHuez, the reason the legacy pipeline excluded
# FR_AlpeHuez::HIGH::LOW; those rows carry no relative_year either
na_year <- filter(community, is.na(year))
stopifnot(all(na_year$region == "FR_AlpeHuez"), all(is.na(na_year$relative_year)),
          !anyNA(community$relative_year[!is.na(community$year)]))
count_row("community", "written", nrow(community))
count_row("community", "of which cover NA (US_Arizona, supplied relative cover only)", nrow(na_cover))
count_row("community", "of which year and relative_year NA (FR_AlpeHuez)", nrow(na_year))
community_summary <- tibble(
  regions = n_distinct(community$region),
  plots = n_distinct(community$stable_plot_id),
  plot_years = n_distinct(community$plot_year_id),
  species = n_distinct(community$species),
  relative_year = paste(range(community$relative_year, na.rm = TRUE), collapse = " to "),
  cover_NA_rows = nrow(na_cover),
  cover_NA_plot_years = n_distinct(na_cover$plot_year_id),
  year_NA_rows = nrow(na_year)
)

species_taxonomy <- read_target("community_taxonomy_species")
stopifnot(!anyDuplicated(species_taxonomy$submitted_name))  # one row per submitted name
count_row("species_taxonomy", "read and written", nrow(species_taxonomy))

# ---- species properties --------------------------------------------------------
trait_summary <- read_target("try_trait_outputs")$summary
count_row("traits", "read (species x trait rows)", nrow(trait_summary))
n_blank_trait <- sum(trait_summary$accepted_name == "")
# a TRY record without an accepted name cannot be joined to any species
trait_summary <- filter(trait_summary, accepted_name != "")
stopifnot(n_blank_trait == 1)
count_row("traits", "after dropping blank-name rows", nrow(trait_summary))
trait_summary <- mutate(trait_summary, trait = recode(trait, vegetative_height = "height"))
trait_units <- distinct(trait_summary, trait, unit)
traits <- trait_summary |>
  select(species = accepted_name, trait, median, n_observations) |>
  pivot_wider(names_from = trait, values_from = c(median, n_observations)) |>
  rename(LDMC = median_LDMC, LMA = median_LMA, height = median_height,
         n_obs_LDMC = n_observations_LDMC, n_obs_LMA = n_observations_LMA,
         n_obs_height = n_observations_height)
stopifnot(!anyDuplicated(traits$species))  # one row per species
count_row("traits", "written (species)", nrow(traits))

functional_groups <- read_target("functional_groups") |>
  transmute(species = accepted_name, family,
            graminoid = graminoid_status == "graminoid",
            woody = woody_status == "woody",
            functional_group)
stopifnot(!anyDuplicated(functional_groups$species),
          !anyNA(functional_groups$graminoid), !anyNA(functional_groups$woody))
count_row("functional_groups", "read and written", nrow(functional_groups))

pcoa <- read_target("resident_pcoa_candidate")
eigenvalues <- pcoa$eigen_summary
scores <- rename(pcoa$scores, species = species_id)
# all axes are exported; k is chosen where the axes are first used (step 05)
stopifnot(!pcoa$correction_required, all(eigenvalues$positive_eigenvalue),
          !anyDuplicated(scores$species),
          identical(names(scores)[-1], eigenvalues$PEV))
phylogeny_axes <- list(scores = scores, eigenvalues = eigenvalues)
count_row("phylogeny_axes", "read and written (species)", nrow(scores))
n_broken_stick <- sum(eigenvalues$exceeds_broken_stick)
cum_PEV3 <- eigenvalues$cumulative_percent[eigenvalues$PEV == "PEV3"]

# Lynn's table has one row per turf x species; the niche columns are species-level
niche_raw <- collect(read.csv(niche_csv))
count_row("niche_optima", "read (turf x species rows)", nrow(niche_raw))
niche_optima <- niche_raw |>
  distinct(species, mean_sum_sub, mean_ap_sub, mean_mat_sub, n_obs_sub) |>
  transmute(species, T_opt = mean_sum_sub, P_opt = log(mean_ap_sub),
            mean_mat_sub, n_obs_sub)
stopifnot(!anyDuplicated(niche_optima$species))  # one niche row per Lynn species name
count_row("niche_optima", "written (species)", nrow(niche_optima))
# !is.finite also catches log(0) = -Inf in P_opt
n_niche_na <- sum(!is.finite(niche_optima$T_opt) | !is.finite(niche_optima$P_opt))
n_niche_30 <- sum(niche_optima$n_obs_sub >= 30, na.rm = TRUE)

# ---- write and manifest ---------------------------------------------------------
outputs <- list(
  experiments = list(x = experiments, from = c("experiments", "sites")),
  sites = list(x = sites, from = "sites"),
  experiment_climate = list(x = experiment_climate, from = "experiment_climate"),
  community = list(x = community, from = "analysis_community"),
  species_taxonomy = list(x = species_taxonomy, from = "community_taxonomy_species"),
  traits = list(x = traits, from = "try_trait_outputs"),
  functional_groups = list(x = functional_groups, from = "functional_groups"),
  phylogeny_axes = list(x = phylogeny_axes, from = "resident_pcoa_candidate"),
  niche_optima = list(x = niche_optima, from = niche_csv)
)
meta <- tar_meta(store = store, fields = c("name", "data", "time"))
sha256_file <- function(path) digest::digest(file = path, algo = "sha256")
units_note <- c(traits = paste(sprintf("%s: %s", trait_units$trait,
                                       sub("1/(mm2 mg-1)", "mg mm-2", trait_units$unit, fixed = TRUE)),
                               collapse = "; "))

manifest <- bind_rows(lapply(names(outputs), function(nm) {
  o <- outputs[[nm]]
  path <- file.path(out_dir, paste0(nm, ".rds"))
  saveRDS(o$x, path)
  tab <- if (is.data.frame(o$x)) o$x else o$x$scores
  from_store <- o$from %in% meta$name
  m <- meta[match(o$from[from_store], meta$name), ]
  tibble(
    file = basename(path), rows = nrow(tab), columns = ncol(tab),
    sha256 = sha256_file(path),
    source = paste(o$from, collapse = "; "),
    target_data_hash = if (any(from_store)) paste(m$data, collapse = "; ") else NA_character_,
    target_time_built = if (any(from_store)) paste(format(m$time, "%Y-%m-%d %H:%M:%S"), collapse = "; ") else NA_character_,
    source_sha256 = if (any(from_store)) NA_character_ else sha256_file(o$from),
    source_modified = if (any(from_store)) NA_character_ else format(file.mtime(o$from), "%Y-%m-%d %H:%M:%S"),
    units = if (nm %in% names(units_note)) units_note[[nm]] else NA_character_
  )
}))
write.csv(manifest, file.path(out_dir, "manifest.csv"), row.names = FALSE)
# every output has a row in manifest.csv
stopifnot(setequal(manifest$file, setdiff(list.files(out_dir), "manifest.csv")))

# ---- report -----------------------------------------------------------------------
md_table <- function(d) {
  d <- mutate(d, across(everything(), ~ ifelse(is.na(.x), "", as.character(.x))))
  c(paste("|", paste(names(d), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(d)), collapse = " | "), "|"),
    apply(d, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
climate_summary <- experiment_climate |>
  summarise(across(c(T, W, P_summer, P_annual),
                   list(min = min, median = median, max = max))) |>
  pivot_longer(everything(), names_to = c("variable", ".value"),
               names_pattern = "^(.*)_(min|median|max)$") |>
  mutate(across(c(min, median, max), ~ round(.x, 3)))

report <- c(
  "# Step 00: Freeze inputs", "",
  sprintf(paste(
    "The script read %d objects from the legacy targets store and Lynn's niche file and",
    "wrote them unchanged in content to `%s`, one `.rds` per object, with `manifest.csv`.",
    "The only changes are: the four site coordinates joined onto `experiments`; T, W and",
    "log precipitation derived in `experiment_climate`; the legacy `relative_cover` dropped",
    "from `community`, because PLAN.md binds that name to the cover step 01 recomputes after",
    "the deep-lineage removal; traits pivoted to one row per species; and niche columns",
    "renamed to `T_opt` and `P_opt`. No filter changes the analysis set in this step."),
    sum(unique(unlist(lapply(outputs, `[[`, "from"))) %in% meta$name), out_dir),
  "", sprintf("Run with %s, targets %s.", R.version.string, packageVersion("targets")),
  "", "## Counts", "", md_table(counts),
  "", sprintf("Experiments excluded by the legacy pipeline (kept in `experiments.rds` with `included == FALSE`): %s.",
              paste(sprintf("%s (%s)", experiments$experiment_id[!experiments$included],
                            experiments$exclusion_reason[!experiments$included]), collapse = "; ")),
  "", "## Key results", "", "### Manifest", "", md_table(manifest),
  "", "### Community (`analysis_community`, legacy `relative_cover` dropped)", "",
  md_table(community_summary),
  "", sprintf(paste(
    "The %d rows without cover are US_Arizona, which recorded relative cover only; it is carried",
    "in `supplied_relative_cover` (%d of these rows filled), not in `cover`. US_Arizona is removed",
    "in step 01 by the pool-size rule."),
    nrow(na_cover), sum(!is.na(na_cover$supplied_relative_cover))),
  "", "### Climate of the 42 included experiments",
  "", paste("T = destination JJA temperature (degC), W = destination minus origin JJA temperature,",
            "P_summer = log destination JJA precipitation, P_annual = log destination annual precipitation",
            "(legacy sum of monthly CHELSA precipitation per year, averaged over the experiment's years;",
            "`legacy/R/data/climate.R:346`, not recomputed here). Not standardised."),
  "", md_table(climate_summary),
  "", sprintf("Correlation of P_annual and P_summer across the %d included experiments: r = %.3f (PLAN.md expects about 0.78).",
              nrow(experiment_climate), r_P),
  "", sprintf("Experiments with partial climate coverage (%d), averaged over the years available in the CHELSA cache:", nrow(partial_climate)),
  "", md_table(partial_climate),
  "", "### Species properties", "",
  sprintf("- Traits: %d species with at least one of LDMC, LMA, height; %d with LDMC, %d with LMA, %d with height. Legacy per-species medians; %d blank-name LMA row dropped.",
          nrow(traits), sum(!is.na(traits$LDMC)), sum(!is.na(traits$LMA)), sum(!is.na(traits$height)), n_blank_trait),
  sprintf("- Functional groups: %d species, %d graminoid, %d woody.",
          nrow(functional_groups), sum(functional_groups$graminoid), sum(functional_groups$woody)),
  sprintf("- Phylogeny: %d resident species, %d PCoA axes with eigenvalues, all exported. %d axes exceed the broken-stick expectation; PEV1 to PEV3 carry %.1f%% of the positive variation.",
          nrow(scores), nrow(eigenvalues), n_broken_stick, cum_PEV3),
  sprintf("- Niche optima: %d species in Lynn's file, keyed by Lynn's name (no join to community names here); %d without a finite T_opt or P_opt; %d with `n_obs_sub >= 30`.",
          nrow(niche_optima), n_niche_na, n_niche_30),
  "", paste("Niche provenance: `MacroContextTransplantDat_2.csv` was built by J. Lynn (11.2024) with",
            "`DataProcessing_MacroContext.R` in TransPlant-Updates_oct2021, following the method of",
            "Lynn et al. 2021 (Global Change Biology). Species niche summaries are CHELSA 1981 to 2010",
            "climatologies over GBIF occurrences; the `_sub` columns use occurrences thinned to one per",
            "CHELSA (30 arc-second) cell, and `n_obs_sub` is the number of thinned occurrences behind them.",
            "None of Lynn's processing is rerun; only columns are taken. Source of this description:",
            "`legacy/R/additional_checks/thermal_niche/random_forest_thermal_niche.R`, lines 31 to 37."),
  "", paste("`T_opt` = `mean_sum_sub`, taken as summer temperature (the bio10 analogue); this reading of",
            "the `sum` columns still needs confirming. `P_opt` = log(`mean_ap_sub`), annual precipitation",
            "(bio12). `mean_mat_sub` (bio1) is carried only for the bio10 against JJA check of Fig S5."),
  "", sprintf("All %d columns of the CSV: %s.", ncol(niche_raw),
              paste(sprintf("`%s`", names(niche_raw)), collapse = ", ")),
  "", "## Diagnostic figures", "", "None in this step.",
  "", "## Warnings captured", "",
  if (length(warns)) paste("-", unique(warns)) else "None.",
  "", "## Deviations from PLAN.md", "",
  paste("None in the analysis. Open question for step 07: the PCoA covers resident species only,",
        "so whether colonizers get a phylogeny property, and from which PCoA, is not settled.")
)
writeLines(report, "reports/00_freeze_inputs.md")
