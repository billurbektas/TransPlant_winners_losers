# P5: compare the prep inputs (data/derived/prep/) with the frozen legacy
# inputs of step 00 (data/derived/00_inputs/) and name a cause for every
# difference.
#
# A one-off script, not a target: it needs the step 00 files, which are
# archived after this run (~/Desktop/TransPlant_archive/2026-10-legacy/). It
# reruns only with data/derived/00_inputs/ restored.
#
# Reads   data/derived/prep/*.rds; data/derived/00_inputs/*.rds;
#         data/TransPlantNetwork_060824.RData (`dat`, for record-level
#         reconciliation); the prep store (community_raw, name_resolution,
#         lynn, try_attached); prep/config/ ledgers
# Writes  reports/P5_compare_to_frozen.md; reports/P5_differences/<file>.csv;
#         prep/config/legacy_taxonomy.csv
# Checks  every difference has a named cause (stop otherwise); frozen and
#         prep community rows add up over the causes

suppressMessages(library(dplyr))
source("prep/R/report.R") # md_table()

f00 <- function(x) readRDS(file.path("data/derived/00_inputs", paste0(x, ".rds")))
fp <- function(x) readRDS(file.path("data/derived/prep", paste0(x, ".rds")))
store <- "prep/_targets"
binomial <- function(x) stringr::word(x, 1, 2)
dir.create("reports/P5_differences", showWarnings = FALSE)

# every cause is named here with its kind; a difference whose cause is not
# in this list stops the script
causes <- tribble(
  ~cause, ~kind,
  "species renamed or merged by the P1 taxonomy", "correction",
  "US_Arizona: absolute cover from the SQLite (step 00 had supplied relative cover only)", "new data",
  "SQLite value differs from the RData (data cleaning repository)", "new data",
  "SQLite recodes the record (FR_Lautaret Pilosella officinarum -> Pilosella)", "new data",
  "cover summed over names merged into one species", "correction",
  "FR_AlpeHuez: survey years in the SQLite only (the RData had no year)", "new data",
  "CH_Calanda2 2020 kept (legacy drop C6 removed, decided 26 Oct)", "correction",
  "species in one plot-year of its experiment x treatment: legacy filter not applied (legacy/R/data/community.R:242-245)", "correction",
  "vegetation record kept with its drop_reason (relative-cover denominator, step 01)", "correction",
  "excluded by decision (prep/config/taxonomy_decisions.csv)", "correction",
  "taxonomy decision (prep/config/taxonomy_decisions.csv)", "correction",
  "accepted-name rename (prep/config/accepted_renames.csv)", "correction",
  "legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match", "correction",
  "name below species level, now kept with its identification level", "correction",
  "name new in the SQLite (FR_AlpeHuez)", "new data",
  "unidentified morphotaxon, now listed (prep/config/unidentified_taxa.csv)", "correction",
  "name recorded only in establishment years (dropped by both pipelines; P1 resolves every recorded name)", "correction",
  "records attached by the contributor's own name, not TRY's AccSpeciesName", "correction",
  "P2 cleaning rules (experiments removed per observation, ErrorRisk < 3, non-positive values dropped)", "correction",
  "not a P1 study species", "correction",
  "no TRY record left after the P2 rules", "correction",
  "destination-pool species, now included (P3 covers both pools)", "correction",
  "study species new to the community (P1 taxonomy or FR_AlpeHuez)", "correction",
  "woodiness ledger (prep/config/woodiness.csv)", "correction",
  "family from the P1 taxonomy", "correction",
  "axes recomputed on the tree of all 904 study species", "correction",
  "experiment years changed (CH_Calanda2 2020 kept)", "correction",
  "experiment included: FR_AlpeHuez years in the SQLite", "new data",
  "Lynn name resolves to another accepted species", "correction",
  "one of several Lynn names of one species; the one with most occurrences is used", "correction",
  "Lynn name is not a P1 study species", "correction",
  "Lynn's niche values are NA", "correction",
  "Lynn name not resolved by the P1 rule (needs review; community decisions are not applied to Lynn's names)", "correction")

diffs <- list()
add_diffs <- function(file, d) {
  d <- mutate(d, file = file, .before = 1)
  # a combined cause ("a + b") is named when each of its parts is
  named <- vapply(strsplit(coalesce(d$cause, ""), " + ", fixed = TRUE), function(x) length(x) && all(x %in% causes$cause), TRUE)
  bad <- filter(d, !named)
  if (nrow(bad)) {
    print(as.data.frame(utils::head(bad, 30)))
    stop(nrow(bad), " differences in ", file, " without a named cause")
  }
  d <- mutate(d, kind = vapply(strsplit(cause, " + ", fixed = TRUE),
                               function(x) paste(sort(unique(causes$kind[match(x, causes$cause)])), collapse = " + "), ""))
  utils::write.csv(d, file.path("reports/P5_differences", paste0(file, ".csv")), row.names = FALSE, na = "")
  diffs[[file]] <<- d
  invisible(d)
}
cause_table <- function(d, n = "rows") {
  out <- count(d, kind, cause, name = n)
  arrange(out, kind, desc(.data[[n]]))
}
columns_row <- function(file, old, new) {
  cols <- function(x) {
    pev <- grepl("^PEV[0-9]+$", x)
    if (sum(pev) > 3) x <- c(x[!pev], sprintf("PEV%s to PEV%s", min(as.integer(sub("PEV", "", x[pev]))), max(as.integer(sub("PEV", "", x[pev])))))
    paste(x, collapse = ", ")
  }
  tibble(file = file, rows_step00 = nrow(old), rows_prep = nrow(new),
         columns_added = cols(setdiff(names(new), names(old))), columns_dropped = cols(setdiff(names(old), names(new))))
}

# ---- legacy names -> P1 species ---------------------------------------------
ot <- f00("species_taxonomy")
nt <- fp("species_taxonomy")
renames <- utils::read.csv("prep/config/accepted_renames.csv")
name_resolution <- targets::tar_read(name_resolution, store = store)

# the placement draft (P3) reuses legacy labels through these two columns
utils::write.csv(select(ot, submitted_name, accepted_name_reviewed), "prep/config/legacy_taxonomy.csv",
                 row.names = FALSE)

# a legacy accepted name maps to the P1 species of the names it stood for; a
# legacy name absent from the community (TRY only) maps through the P1
# resolution of that name
legacy_map <- ot |>
  select(name_original, legacy = accepted_name_reviewed) |>
  inner_join(select(nt, name_original, species = accepted_species), by = "name_original") |>
  filter(!is.na(species))
# keyed on the full legacy name (traits, groups) and on its binomial (community)
legacy_map <- distinct(bind_rows(legacy_map, mutate(legacy_map, legacy = binomial(legacy))), legacy, species)
map_legacy <- function(x) {
  m <- filter(legacy_map, legacy %in% x)
  rest <- setdiff(x, m$legacy)
  r <- name_resolution |>
    filter(submitted_name %in% rest) |>
    transmute(legacy = submitted_name, species = coalesce(accepted_species, binomial(accepted_name))) |>
    filter(!is.na(species))
  distinct(bind_rows(m, r))
}

# species status of a legacy species set against a prep species set
species_status <- function(old, new) {
  m <- map_legacy(old) |> filter(species %in% new)
  legacy <- tibble(species_00 = old) |>
    mutate(status = case_when(species_00 %in% new ~ "same name",
                              species_00 %in% m$legacy ~ "renamed or merged",
                              TRUE ~ "lost"),
           species_prep = ifelse(status == "same name", species_00,
                                 vapply(species_00, function(s) paste(m$species[m$legacy == s], collapse = "; "), "")))
  gained <- tibble(species_prep = setdiff(new, c(old, m$species)), status = "gained")
  bind_rows(legacy, gained)
}
species_summary <- function(file, st) {
  tibble(file = file, same_name = sum(st$status == "same name"), renamed_or_merged = sum(st$status == "renamed or merged"),
         lost = sum(st$status == "lost"), gained = sum(st$status == "gained"))
}

# ---- community ------------------------------------------------------------------
o <- f00("community")
n <- fp("community")
keys <- c("region", "year", "destination_plot", "origin_site", "destination_site", "treatment", "species")

# SQLite against the RData, record by record (FR_AlpeHuez aside: the RData
# lacks its years)
raw <- targets::tar_read(community_raw, store = store)
env <- new.env()
load("data/TransPlantNetwork_060824.RData", envir = env)
rdata <- transmute(env$dat, region = Region, year = as.integer(Year), destination_plot = destPlotID,
                   treatment = Treatment, local_species_code = SpeciesName, cover_rdata = Cover)
recon <- full_join(transmute(raw, region, year, destination_plot, treatment, local_species_code, cover_sqlite = cover),
                   rdata, by = c("region", "year", "destination_plot", "treatment", "local_species_code")) |>
  filter(region != "FR_AlpeHuez", !mapply(identical, cover_sqlite, cover_rdata))
value_changed <- filter(recon, !is.na(cover_sqlite), !is.na(cover_rdata), region != "US_Arizona")
recoded <- filter(recon, is.na(cover_sqlite))

# frozen rows carry the recorded name; they are re-keyed to the P1 species
# through region x local code
codes <- n |>
  filter(!is.na(species)) |>
  tidyr::separate_rows(local_species_code, sep = "; ") |>
  distinct(region, local_species_code, p1_species = species)
stopifnot(!anyDuplicated(codes[c("region", "local_species_code")])) # one P1 species per local code
o <- left_join(o, codes, by = c("region", "local_species_code"))
o_unmapped <- filter(o, is.na(p1_species)) |>
  left_join(tidyr::separate_rows(filter(n, is.na(species)), local_species_code, sep = "; ") |>
              distinct(region, year, destination_plot, local_species_code, p1_drop = drop_reason),
            by = c("region", "year", "destination_plot", "local_species_code")) |>
  mutate(cause = case_when(
    paste(region, year, destination_plot, local_species_code) %in%
      paste(recoded$region, recoded$year, recoded$destination_plot, recoded$local_species_code) ~
      "SQLite recodes the record (FR_Lautaret Pilosella officinarum -> Pilosella)",
    p1_drop == "excluded by decision" ~ "excluded by decision (prep/config/taxonomy_decisions.csv)"))
oc <- o |>
  filter(!is.na(p1_species)) |>
  mutate(legacy_name = species, species = p1_species) |>
  summarise(cover_00 = sum(cover), rows_00 = n(), cover_00_na = all(is.na(cover)),
            codes_00 = paste(sort(unique(local_species_code)), collapse = "; "), .by = all_of(keys))
nk <- filter(n, !is.na(species))
region_years_00 <- distinct(o, region, year) |> mutate(region_year_in_00 = TRUE)
plot_years_00 <- distinct(o, region, year, destination_plot) |> mutate(plot_year_in_00 = TRUE)

# the legacy filter kept a species x experiment x treatment only with more
# than one plot-year, counted per recorded name over the years legacy had
legacy_rows <- nk |>
  semi_join(region_years_00, by = c("region", "year")) |>
  filter(!(region == "CH_Calanda2" & year == 2020L)) |>
  tidyr::separate_rows(name_original, sep = "; ")
legacy_single <- legacy_rows |>
  mutate(n_plot_years = n_distinct(paste(destination_plot, year)),
         .by = c(region, origin_site, destination_site, treatment, name_original)) |>
  summarise(legacy_single = all(n_plot_years == 1L), .by = all_of(keys))

changed_keys <- paste(value_changed$region, value_changed$year, value_changed$destination_plot)
comm <- full_join(oc, select(nk, all_of(keys), cover_prep = cover, n_records, local_species_code), by = keys) |>
  left_join(region_years_00, by = c("region", "year")) |>
  left_join(plot_years_00, by = c("region", "year", "destination_plot")) |>
  left_join(legacy_single, by = keys) |>
  mutate(
    in_00 = !is.na(rows_00), in_prep = !is.na(cover_prep),
    same = in_00 & in_prep & !cover_00_na & abs(cover_00 - cover_prep) < 1e-9,
    cause = case_when(
      same ~ NA_character_,
      in_00 & in_prep & cover_00_na & region == "US_Arizona" ~
        "US_Arizona: absolute cover from the SQLite (step 00 had supplied relative cover only)",
      in_00 & in_prep & n_records > rows_00 ~ "cover summed over names merged into one species",
      in_00 & in_prep & paste(region, year, destination_plot) %in% changed_keys &
        mapply(function(r, y, p, c) any(value_changed$region == r & value_changed$year == y &
                                          value_changed$destination_plot == p &
                                          value_changed$local_species_code %in% strsplit(c, "; ")[[1]]),
               region, year, destination_plot, coalesce(local_species_code, "")) ~
        "SQLite value differs from the RData (data cleaning repository)",
      in_00 & is.na(year) & region == "FR_AlpeHuez" ~ "FR_AlpeHuez: survey years in the SQLite only (the RData had no year)",
      in_prep & region == "FR_AlpeHuez" & (is.na(region_year_in_00) | is.na(plot_year_in_00)) ~
        "FR_AlpeHuez: survey years in the SQLite only (the RData had no year)",
      in_prep & region == "CH_Calanda2" & year == 2020L ~ "CH_Calanda2 2020 kept (legacy drop C6 removed, decided 26 Oct)",
      in_prep & !in_00 & legacy_single %in% TRUE ~
        "species in one plot-year of its experiment x treatment: legacy filter not applied (legacy/R/data/community.R:242-245)"))
# rows set aside in prep (no species): every vegetation record is kept
comm_aside <- filter(n, is.na(species)) |>
  transmute(region, year, destination_plot, treatment, local_species_code, cover_prep = cover, drop_reason,
            cause = case_when(region == "FR_AlpeHuez" ~ "FR_AlpeHuez: survey years in the SQLite only (the RData had no year)",
                              region == "CH_Calanda2" & year == 2020L ~ "CH_Calanda2 2020 kept (legacy drop C6 removed, decided 26 Oct)",
                              TRUE ~ "vegetation record kept with its drop_reason (relative-cover denominator, step 01)"))

comm_diff <- bind_rows(
  filter(comm, !same) |>
    transmute(key = paste(region, year, destination_plot, treatment, species, sep = " | "),
              difference = case_when(in_00 & in_prep ~ sprintf("cover %s -> %s", ifelse(cover_00_na, "NA", format(cover_00)), format(cover_prep)),
                                     in_00 ~ sprintf("step 00 only (%d rows)", rows_00), TRUE ~ "prep only"),
              cause, rows_00 = coalesce(rows_00, 0L), rows_prep = as.integer(in_prep)),
  transmute(o_unmapped, key = paste(region, year, destination_plot, treatment, local_species_code, sep = " | "),
            difference = "step 00 only (1 row)", cause, rows_00 = 1L, rows_prep = 0L),
  transmute(comm_aside, key = paste(region, year, destination_plot, treatment, local_species_code, sep = " | "),
            difference = paste("prep only, drop_reason:", drop_reason), cause, rows_00 = 0L, rows_prep = 1L))
add_diffs("community", comm_diff)
comm_same <- filter(comm, same)
# frozen and prep rows add up over the causes
stopifnot(sum(comm_same$rows_00) + sum(diffs$community$rows_00) == nrow(o),
          nrow(comm_same) + sum(diffs$community$rows_prep) == nrow(n))
comm_causes <- diffs$community |>
  summarise(rows_step00 = sum(rows_00), rows_prep = sum(rows_prep), .by = c(kind, cause)) |>
  arrange(kind, desc(rows_prep + rows_step00))
comm_causes <- bind_rows(tibble(kind = "", cause = "same row, same cover", rows_step00 = sum(comm_same$rows_00),
                                rows_prep = nrow(comm_same)), comm_causes,
                         tibble(kind = "", cause = "total", rows_step00 = nrow(o), rows_prep = nrow(n)))

# species of the analysis rows: the legacy accepted name of each recorded name
legacy_species <- unique(binomial(ot$accepted_name_reviewed[ot$name_original %in% o$species]))
st_comm <- species_status(legacy_species, unique(nk$species))
excluded <- binomial(ot$accepted_name_reviewed[ot$name_original %in% nt$name_original[nt$identification_level %in% "excluded by decision"]])
comm_species <- st_comm |>
  filter(status %in% c("lost", "gained")) |>
  mutate(cause = ifelse(status == "lost" & species_00 %in% excluded,
                        "excluded by decision (prep/config/taxonomy_decisions.csv)", NA_character_))
# a gained species takes the cause(s) of its prep rows
gained_rows <- filter(diffs$community, rows_prep == 1L) |>
  mutate(species_prep = sub("^.*\\| ", "", key)) |>
  filter(species_prep %in% comm_species$species_prep) |>
  summarise(cause = paste(unique(cause), collapse = " + "), .by = species_prep)
comm_species <- rows_update(comm_species, gained_rows, by = "species_prep", unmatched = "ignore")

# ---- sites and experiments -------------------------------------------------------
value_diffs <- function(old, new, by) {
  common <- setdiff(intersect(names(old), names(new)), by)
  j <- inner_join(old, new, by = by, suffix = c("_00", "_prep"))
  bind_rows(lapply(common, function(col) {
    a <- j[[paste0(col, "_00")]]
    b <- j[[paste0(col, "_prep")]]
    differ <- !mapply(identical, a, b)
    if (!any(differ)) return(NULL)
    tibble(key = do.call(paste, c(j[differ, by], sep = "::")), column = col,
           step00 = as.character(a[differ]), prep = as.character(b[differ]))
  }))
}
os <- f00("sites")
ns <- fp("sites")
site_diff <- value_diffs(os, ns, c("region", "site_id")) |>
  mutate(cause = case_when(startsWith(key, "FR_AlpeHuez::") ~ "FR_AlpeHuez: survey years in the SQLite only (the RData had no year)"),
         difference = sprintf("%s: %s -> %s", column, step00, prep))
add_diffs("sites", select(site_diff, key, difference, cause))

oe <- f00("experiments")
ne <- fp("experiments")
exp_diff <- value_diffs(oe, ne, "experiment_id") |>
  mutate(cause = case_when(key == "FR_AlpeHuez::HIGH::LOW" ~ "experiment included: FR_AlpeHuez years in the SQLite",
                           key == "CH_Calanda2::Cal::Nes" ~ "experiment years changed (CH_Calanda2 2020 kept)"),
         difference = sprintf("%s: %s -> %s", column, step00, prep))
add_diffs("experiments", select(exp_diff, key, difference, cause))

# ---- species_taxonomy (one row per recorded name) -----------------------------------
alpe_names <- unique(n$name_original[n$region == "FR_AlpeHuez"])
other_names <- unique(n$name_original[n$region != "FR_AlpeHuez"])
community_names <- unique(unlist(strsplit(n$name_original, "; ")))
analysis_names <- unique(unlist(strsplit(nk$name_original, "; ")))
# legacy built its taxonomy from the filtered community, so a name whose
# records all fell to a legacy filter never reached it: such a name takes
# the cause of its rows
row_cause <- comm |>
  filter(in_prep) |>
  left_join(select(nk, all_of(keys), name_original), by = keys) |>
  tidyr::separate_rows(name_original, sep = "; ") |>
  summarise(row_cause = if (any(is.na(cause))) NA_character_ else paste(sort(unique(cause)), collapse = " + "),
            .by = name_original)
tax <- full_join(select(ot, name_original, legacy = accepted_name_reviewed),
                 select(nt, name_original, accepted_species, decided_by, identification_level), by = "name_original") |>
  mutate(
    in_00 = !is.na(legacy), in_prep = !is.na(decided_by),
    same = in_00 & in_prep & coalesce(binomial(legacy) == accepted_species, FALSE),
    cause = case_when(
      same ~ NA_character_,
      in_00 & !in_prep ~ "SQLite recodes the record (FR_Lautaret Pilosella officinarum -> Pilosella)",
      in_00 & identification_level %in% "excluded by decision" ~ "excluded by decision (prep/config/taxonomy_decisions.csv)",
      in_00 & decided_by == "decision" ~ "taxonomy decision (prep/config/taxonomy_decisions.csv)",
      in_00 & accepted_species %in% renames$to_accepted_species ~ "accepted-name rename (prep/config/accepted_renames.csv)",
      in_00 ~ "legitimate-name rule, WCVP first (tidyTRY 0.2.1), replaces the legacy TNRS best match",
      decided_by == "unidentified" ~ "unidentified morphotaxon, now listed (prep/config/unidentified_taxa.csv)",
      identification_level %in% c("genus_only", "family_only") ~ "name below species level, now kept with its identification level",
      name_original %in% alpe_names & !name_original %in% other_names ~ "name new in the SQLite (FR_AlpeHuez)",
      !name_original %in% community_names ~
        "name recorded only in establishment years (dropped by both pipelines; P1 resolves every recorded name)",
      !name_original %in% analysis_names ~ "vegetation record kept with its drop_reason (relative-cover denominator, step 01)",
      TRUE ~ row_cause$row_cause[match(name_original, row_cause$name_original)]),
    difference = case_when(in_00 & in_prep ~ sprintf("%s -> %s", legacy, coalesce(accepted_species, identification_level)),
                           in_00 ~ "step 00 only", TRUE ~ sprintf("prep only (%s)", coalesce(identification_level, decided_by))))
add_diffs("species_taxonomy", transmute(filter(tax, !same), key = name_original, difference, cause))

# ---- traits -------------------------------------------------------------------------
of <- f00("traits")
nf <- fp("traits")
tl <- fp("traits_long")
attached <- targets::tar_read(try_attached, store = store)
study <- unique(na.omit(nt$accepted_species[nt$identification_level %in% "species_or_lower"]))
legacy_study <- unique(map_legacy(legacy_species)$species)
st_tr <- species_status(of$species, nf$species)
own_name_moved <- tl |>
  filter(species != coalesce(species_if_trusting_try, "")) |>
  tidyr::pivot_longer(c(species, species_if_trusting_try), values_to = "sp") |>
  distinct(sp, trait)
pairs <- st_tr |>
  filter(status %in% c("same name", "renamed or merged")) |>
  tidyr::separate_rows(species_prep, sep = "; ")
long_vals <- function(d, suffix) {
  tidyr::pivot_longer(select(d, species, LDMC, LMA, height), -species, names_to = "trait", values_to = paste0("value", suffix)) |>
    left_join(tidyr::pivot_longer(select(d, species, n_obs_LDMC, n_obs_LMA, n_obs_height), -species,
                                  names_to = "trait", names_prefix = "n_obs_", values_to = paste0("n_obs", suffix)),
              by = c("species", "trait"))
}
trait_cmp <- pairs |>
  left_join(long_vals(of, "_00"), by = c(species_00 = "species"), relationship = "many-to-many") |>
  left_join(long_vals(nf, "_prep"), by = c(species_prep = "species", "trait")) |>
  filter(!is.na(value_00) | !is.na(value_prep)) |>
  mutate(relative_difference = (value_prep - value_00) / value_00,
         same = coalesce(abs(relative_difference) < 1e-9 & n_obs_00 == n_obs_prep, FALSE),
         cause = case_when(
           same ~ NA_character_,
           status == "renamed or merged" ~ "species renamed or merged by the P1 taxonomy",
           paste(species_prep, trait) %in% paste(own_name_moved$sp, own_name_moved$trait) ~
             "records attached by the contributor's own name, not TRY's AccSpeciesName",
           TRUE ~ "P2 cleaning rules (experiments removed per observation, ErrorRisk < 3, non-positive values dropped)"),
         difference = sprintf("%s: %s (n %s) -> %s (n %s)", trait, signif(value_00, 3), n_obs_00, signif(value_prep, 3), n_obs_prep))
trait_species <- filter(st_tr, status %in% c("lost", "gained")) |>
  mutate(
    mapped = vapply(species_00, function(s) paste(map_legacy(s)$species, collapse = "; "), ""),
    cause = case_when(
      status == "lost" & (mapped == "" | !vapply(strsplit(mapped, "; "), function(x) any(x %in% study), TRUE)) ~ "not a P1 study species",
      status == "lost" & !vapply(strsplit(mapped, "; "), function(x) any(x %in% attached$species), TRUE) ~
        "records attached by the contributor's own name, not TRY's AccSpeciesName",
      status == "lost" ~ "no TRY record left after the P2 rules",
      !species_prep %in% legacy_study ~ "study species new to the community (P1 taxonomy or FR_AlpeHuez)",
      species_prep %in% own_name_moved$sp ~ "records attached by the contributor's own name, not TRY's AccSpeciesName",
      TRUE ~ "P2 cleaning rules (experiments removed per observation, ErrorRisk < 3, non-positive values dropped)"),
    difference = ifelse(status == "lost", paste("no traits now; P1 species:", ifelse(mapped == "", "none", mapped)), "traits now, none in step 00"))
add_diffs("traits", bind_rows(
  transmute(filter(trait_cmp, !same), key = paste(species_00, "->", species_prep), difference, cause,
            relative_difference = round(relative_difference, 4)),
  transmute(trait_species, key = coalesce(species_00, species_prep), difference, cause)))

# ---- functional groups ---------------------------------------------------------------
og <- f00("functional_groups")
ng <- fp("functional_groups")
wood <- utils::read.csv("prep/config/woodiness.csv")
st_fg <- species_status(og$species, ng$species)
fg_cmp <- st_fg |>
  filter(status %in% c("same name", "renamed or merged")) |>
  tidyr::separate_rows(species_prep, sep = "; ") |>
  left_join(select(og, species_00 = species, family_00 = family, graminoid_00 = graminoid, woody_00 = woody), by = "species_00") |>
  left_join(select(ng, species_prep = species, genus, family, graminoid, woody), by = "species_prep") |>
  mutate(cause = case_when(
    woody != woody_00 & (species_prep %in% wood$name[wood$level == "species"] | genus %in% wood$name[wood$level == "genus"] |
                           woody_00) ~ "woodiness ledger (prep/config/woodiness.csv)",
    graminoid != graminoid_00 | family != family_00 ~ "family from the P1 taxonomy",
    status == "renamed or merged" ~ "species renamed or merged by the P1 taxonomy"),
    same = status == "same name" & woody == woody_00 & graminoid == graminoid_00 & family == family_00,
    difference = sprintf("family %s -> %s; graminoid %s -> %s; woody %s -> %s", family_00, family, graminoid_00, graminoid, woody_00, woody))
fg_species <- filter(st_fg, status %in% c("lost", "gained")) |>
  mutate(mapped = vapply(species_00, function(s) paste(map_legacy(s)$species, collapse = "; "), ""),
         cause = case_when(status == "lost" & species_00 %in% excluded ~ "excluded by decision (prep/config/taxonomy_decisions.csv)",
                           status == "lost" ~ NA_character_,
                           !species_prep %in% legacy_study ~ "study species new to the community (P1 taxonomy or FR_AlpeHuez)"),
         difference = ifelse(status == "lost", "not in prep", "prep only"))
add_diffs("functional_groups", bind_rows(
  transmute(filter(fg_cmp, !same), key = paste(species_00, "->", species_prep), difference, cause),
  transmute(fg_species, key = coalesce(species_00, species_prep), difference, cause)))

# ---- phylogeny axes ------------------------------------------------------------------
oa <- f00("phylogeny_axes")
na_ <- fp("phylogeny_axes")
st_ax <- species_status(oa$scores$species, na_$scores$species)
ax_pairs <- st_ax |>
  filter(status %in% c("same name", "renamed or merged")) |>
  tidyr::separate_rows(species_prep, sep = "; ") |>
  left_join(select(oa$scores, species_00 = species, PEV1_00 = PEV1, PEV2_00 = PEV2, PEV3_00 = PEV3), by = "species_00") |>
  left_join(select(na_$scores, species_prep = species, PEV1, PEV2, PEV3), by = "species_prep")
pev_cor <- with(filter(ax_pairs, status == "same name"),
                c(PEV1 = cor(PEV1, PEV1_00), PEV2 = cor(PEV2, PEV2_00), PEV3 = cor(PEV3, PEV3_00)))
ax_species <- filter(st_ax, status %in% c("lost", "gained")) |>
  mutate(cause = case_when(status == "lost" & species_00 %in% excluded ~ "excluded by decision (prep/config/taxonomy_decisions.csv)",
                           status == "gained" & species_prep %in% legacy_study ~ "destination-pool species, now included (P3 covers both pools)",
                           status == "gained" ~ "study species new to the community (P1 taxonomy or FR_AlpeHuez)"),
         difference = ifelse(status == "lost", "not in prep", "prep only"))
add_diffs("phylogeny_axes", bind_rows(
  transmute(ax_pairs, key = paste(species_00, "->", species_prep),
            difference = sprintf("PEV1 %.3f -> %.3f; PEV2 %.3f -> %.3f; PEV3 %.3f -> %.3f", PEV1_00, PEV1, PEV2_00, PEV2, PEV3_00, PEV3),
            cause = ifelse(status == "renamed or merged", "species renamed or merged by the P1 taxonomy",
                           "axes recomputed on the tree of all 904 study species"),
            abs_difference_PEV1 = abs(PEV1 - PEV1_00)),
  transmute(ax_species, key = coalesce(species_00, species_prep), difference, cause)))

# ---- experiment climate ---------------------------------------------------------------
oc_ <- f00("experiment_climate")
nc_ <- fp("experiment_climate")
clim <- full_join(select(oc_, experiment_id, T_00 = T, W_00 = W, P_summer_00 = P_summer, P_annual_00 = P_annual,
                         years_00 = destination_years_expected),
                  select(nc_, experiment_id, T, W, P_summer, P_annual, years = destination_years_expected), by = "experiment_id") |>
  mutate(max_abs_difference = pmax(abs(T - T_00), abs(W - W_00), abs(P_summer - P_summer_00), abs(P_annual - P_annual_00)),
         same = !is.na(max_abs_difference) & max_abs_difference < 1e-9,
         cause = case_when(same ~ NA_character_,
                           is.na(T_00) & experiment_id == "FR_AlpeHuez::HIGH::LOW" ~ "experiment included: FR_AlpeHuez years in the SQLite",
                           years != years_00 & experiment_id == "CH_Calanda2::Cal::Nes" ~ "experiment years changed (CH_Calanda2 2020 kept)"),
         difference = sprintf("years %s -> %s; T %s -> %s; W %s -> %s; P_annual %s -> %s", years_00, years,
                              round(T_00, 3), round(T, 3), round(W_00, 3), round(W, 3), round(P_annual_00, 3), round(P_annual, 3)))
add_diffs("experiment_climate", transmute(filter(clim, !same), key = experiment_id, difference, cause))

# ---- niche optima ----------------------------------------------------------------------
on <- f00("niche_optima")
nn <- fp("niche_optima")
lynn <- targets::tar_read(lynn, store = store)
merged_names <- unlist(strsplit(na.omit(nn$lynn_names_merged), "; "))
niche <- on |>
  select(lynn_name = species, T_opt_00 = T_opt, P_opt_00 = P_opt) |>
  left_join(select(lynn, lynn_name, accepted_species), by = "lynn_name") |>
  left_join(select(nn, accepted_species = species, T_opt, P_opt, lynn_name_used), by = "accepted_species") |>
  mutate(used = lynn_name %in% nn$lynn_name_used,
         values_equal = used & abs(T_opt - T_opt_00) < 1e-9 & abs(P_opt - P_opt_00) < 1e-9,
         same = values_equal & lynn_name == accepted_species,
         cause = case_when(
           same ~ NA_character_,
           !is.finite(T_opt_00) | !is.finite(P_opt_00) ~ "Lynn's niche values are NA",
           is.na(accepted_species) ~ "Lynn name not resolved by the P1 rule (needs review; community decisions are not applied to Lynn's names)",
           values_equal ~ "Lynn name resolves to another accepted species",
           !accepted_species %in% study ~ "Lynn name is not a P1 study species",
           lynn_name %in% merged_names & !used ~ "one of several Lynn names of one species; the one with most occurrences is used"),
         difference = case_when(values_equal ~ sprintf("%s -> %s (values unchanged)", lynn_name, accepted_species),
                                is.na(accepted_species) ~ "not resolved",
                                TRUE ~ sprintf("not used (P1 species %s uses %s)", accepted_species, coalesce(lynn_name_used, "none"))))
add_diffs("niche_optima", transmute(filter(niche, !same), key = lynn_name, difference, cause))
# Lynn's names come from the same community records: a Lynn name that is a
# decided community name but resolves otherwise leaves the decided species
# without its niche
decisions <- utils::read.csv("prep/config/taxonomy_decisions.csv")
lynn_decided <- lynn |>
  filter(lynn_name %in% decisions$name_original) |>
  left_join(select(decisions, lynn_name = name_original, decision_species = accepted_species), by = "lynn_name") |>
  filter(decision_species != "", !mapply(identical, accepted_species, decision_species)) |>
  transmute(lynn_name, regions, n_obs_sub, lynn_resolves_to = coalesce(accepted_species, "unresolved"), decision_species,
            decision_species_has_niche = decision_species %in% nn$species)

# ---- what later steps name, and what still depends on step 00 or legacy ------------------
# columns dropped from step 00 that steps 01-10 or PLAN.md name in backticks
spec <- readLines("MILESTONES.md")
spec <- spec[grep("^# Part III", spec):(grep("^## 15", spec) - 1)]
spec <- c(spec, readLines("PLAN.md", warn = FALSE))
olds <- list(community = o, sites = os, experiments = oe, species_taxonomy = ot, traits = of, functional_groups = og,
             phylogeny_axes = oa$scores, experiment_climate = oc_, niche_optima = on)
news <- list(community = n, sites = ns, experiments = ne, species_taxonomy = nt, traits = nf, functional_groups = ng,
             phylogeny_axes = na_$scores, experiment_climate = nc_, niche_optima = nn)
dropped <- bind_rows(lapply(names(olds), function(f) {
  cols <- setdiff(names(olds[[f]]), c(names(news[[f]]), "p1_species"))
  tibble(file = f, column = cols)
})) |>
  mutate(named_later = vapply(column, function(cl) any(grepl(paste0("`", cl, "`"), spec, fixed = TRUE)), TRUE))
# where the content of a dropped column that a later step names now lives
moved_to <- c(supplied_relative_cover = "`cover` (US_Arizona absolute cover from the SQLite); step 01 computes relative cover for every region")
named_later <- filter(dropped, named_later)
if (!all(named_later$column %in% names(moved_to))) stop("A dropped column named by a later step has no new place: ",
                                                       paste(setdiff(named_later$column, names(moved_to)), collapse = ", "))
named_later$now_in <- moved_to[named_later$column]

grep_code <- function(pattern, paths) {
  hits <- suppressWarnings(system2("grep", c("-rnE", shQuote(pattern), paths), stdout = TRUE, stderr = FALSE))
  hits[!grepl("^[^:]+:[0-9]+:\\s*#", hits)] # comments do not count
}
analysis_hits <- if (dir.exists("analysis") && length(list.files("analysis"))) grep_code("_targets|legacy/", "analysis") else character()
prep_code <- c("prep/_targets.R", "prep/R")
prep_hits <- grep_code("00_inputs|store = \"_targets\"|\"legacy/", prep_code)

# counts marked "(to confirm)" in steps 01-10 that the prep data settle now
inc <- filter(ne, included) |> left_join(select(ns, region, destination_site = site_id, gradient), by = c("region", "destination_site"))
confirm <- tibble(
  count = c("included experiments", "included experiments without US_Arizona", "gradients of the included experiments",
            "gradients without US_Arizona"),
  expected = c(NA, 42L, NA, 20L),
  prep = c(nrow(inc), sum(inc$region != "US_Arizona"), n_distinct(paste(inc$region, inc$gradient)),
           n_distinct(paste(inc$region, inc$gradient)[inc$region != "US_Arizona"])))

# ---- report -------------------------------------------------------------------------------
all_diffs <- bind_rows(lapply(diffs, function(d) select(d, file, key, difference, cause, kind)))
summary_tab <- bind_rows(lapply(names(olds), function(f) columns_row(f, olds[[f]], news[[f]]))) |>
  left_join(bind_rows(species_summary("community", st_comm), species_summary("traits", st_tr),
                      species_summary("functional_groups", st_fg), species_summary("phylogeny_axes", st_ax)), by = "file") |>
  left_join(count(all_diffs, file, name = "differences"), by = "file")
summary_tab$rows_prep[summary_tab$file == "community"] <- sprintf("%d (%d analysis rows)", nrow(n), nrow(nk))
kind_tab <- count(all_diffs, file, kind, name = "differences")

top <- function(d, col, k = 10) utils::head(arrange(d, desc(abs(.data[[col]]))), k)
tax_changes <- filter(tax, in_00, in_prep, !same) |>
  left_join(select(nt, name_original, n_records, regions), by = "name_original") |>
  transmute(name_original, legacy_accepted = legacy, accepted_species = coalesce(accepted_species, identification_level),
            cause, n_records) |>
  arrange(desc(n_records))
# the counts P2 gave against step 00, by name
traits_by_name <- full_join(select(nf, species, LDMC, LMA, height),
                            select(of, species, LDMC_00 = LDMC, LMA_00 = LMA, height_00 = height), by = "species")
medians_moved <- bind_rows(lapply(c("LDMC", "LMA", "height"), function(t) {
  tibble(new = traits_by_name[[t]], legacy = traits_by_name[[paste0(t, "_00")]])
})) |> filter(!is.na(new), !is.na(legacy), abs(round((new - legacy) / legacy, 3)) > 0.05) # rounded as P2 did
trait_top <- trait_cmp |>
  filter(!same, !is.na(relative_difference)) |>
  top("relative_difference") |>
  transmute(species_00, species_prep, trait, step00 = signif(value_00, 3), n_00 = n_obs_00, prep = signif(value_prep, 3),
            n_prep = n_obs_prep, relative_difference = round(relative_difference, 3), cause)
fg_changes <- filter(fg_cmp, !same, woody != woody_00 | graminoid != graminoid_00) |>
  transmute(species_00, species_prep, graminoid_00, graminoid, woody_00, woody, cause)
lines_cause <- function(d) md_table(cause_table(d, "differences"))

report <- c(
  "# P5: Comparison with the frozen legacy inputs", "",
  paste0("Every file the analysis reads from `data/derived/prep/` was compared with its counterpart in the frozen legacy ",
         "inputs of step 00, row by row, after mapping legacy names to the P1 species. Each of the ", nrow(all_diffs),
         " differences has a named cause, classed as a correction (a rule, ledger or decision of P1 to P4) or new data ",
         "(the network SQLite against the RData); the script stops on any difference without one. The differences are ",
         "listed in full in `reports/P5_differences/`. Step 00 is retired with this step: its script and report remain in ",
         "the git history, its data are archived outside the repository, and steps 01 to 10 read `data/derived/prep/` only. ",
         "This report carries every comparison with step 00 that the P1 to P4 reports made before."),
  "", "## Summary", "", md_table(summary_tab), "",
  "Differences by kind:", "", md_table(tidyr::pivot_wider(kind_tab, names_from = kind, values_from = differences, values_fill = 0L)),
  "", "## Community", "",
  paste0("Frozen rows (recorded names) were re-keyed to the P1 species through region x local code and summed per plot-year; ",
         "prep rows are one per plot-year x species, plus every other vegetation record with its `drop_reason`. ",
         "Rows of each cause (they add up to both tables):"), "",
  md_table(comm_causes), "",
  sprintf(paste("Record by record, the SQLite equals the RData outside FR_AlpeHuez except %d US_Arizona covers (absent from the",
                "RData), %d CN_Heibei covers and %d FR_Lautaret records recoded."),
          sum(recon$region == "US_Arizona"), nrow(value_changed), nrow(recoded)), "",
  md_table(transmute(value_changed, region, year, destination_plot, local_species_code, cover_rdata, cover_sqlite)), "",
  sprintf("Species of the analysis rows: %d legacy species, %d P1 species; same name %d, renamed or merged %d, lost %d, gained %d. Lost and gained species by cause:",
          length(legacy_species), n_distinct(nk$species), sum(st_comm$status == "same name"),
          sum(st_comm$status == "renamed or merged"), sum(st_comm$status == "lost"), sum(st_comm$status == "gained")), "",
  md_table(count(comm_species, status, cause, name = "species")),
  "", "## Sites and experiments", "",
  md_table(select(diffs$sites, key, difference, cause)), "", md_table(select(diffs$experiments, key, difference, cause)),
  "", "## Species taxonomy", "",
  sprintf("Recorded names: %d in step 00 (species level only), %d in prep (every recorded name). By cause:", nrow(ot), nrow(nt)), "",
  lines_cause(diffs$species_taxonomy), "",
  sprintf("### Accepted species that differ from the frozen legacy taxonomy (%d names)", nrow(tax_changes)), "",
  md_table(tax_changes),
  "", "## Traits", "",
  sprintf(paste("Species with traits now: %d; in step 00: %d. By name, as P2 counted: gained %d, lost %d (names that changed count",
                "as both); species x trait medians changed by more than 5%%: %d. After mapping legacy names: same name %d,",
                "renamed or merged %d, lost %d, gained %d."),
          nrow(nf), nrow(of), sum(!nf$species %in% of$species), sum(!of$species %in% nf$species), nrow(medians_moved),
          sum(st_tr$status == "same name"), sum(st_tr$status == "renamed or merged"), sum(st_tr$status == "lost"),
          sum(st_tr$status == "gained")), "",
  lines_cause(diffs$traits), "", "The 10 largest relative differences of a median:", "", md_table(trait_top), "",
  "Legacy species without traits now:", "",
  md_table(transmute(filter(trait_species, status == "lost"), species_00, difference, cause)),
  "", "## Functional groups", "",
  sprintf("Species: %d in step 00, %d in prep. By cause:", nrow(og), nrow(ng)), "", lines_cause(diffs$functional_groups), "",
  sprintf("Species whose graminoid or woody status differs from step 00 (%d):", nrow(fg_changes)), "", md_table(fg_changes),
  "", "## Phylogeny axes", "",
  sprintf(paste("Species: %d in step 00 (residents), %d in prep (both pools). Against step 00 (%d species under the same name):",
                "correlation of scores PEV1 %.3f, PEV2 %.3f, PEV3 %.3f; the destination-pool species rotate the axes after the",
                "first, so PEV2 and PEV3 are read again in step 10 (Table S6). Eigenvalues: %d axes in step 00, %d in prep."),
          nrow(oa$scores), nrow(na_$scores), sum(st_ax$status == "same name"), pev_cor[1], pev_cor[2], pev_cor[3],
          nrow(oa$eigenvalues), nrow(na_$eigenvalues)), "",
  lines_cause(diffs$phylogeny_axes),
  "", "## Experiment climate", "",
  sprintf("%d experiments with the same years agree exactly (largest difference %.2g). Differences:",
          sum(clim$same), max(clim$max_abs_difference[clim$same])), "",
  md_table(select(diffs$experiment_climate, key, difference, cause)),
  "", "## Niche optima", "",
  sprintf(paste("Step 00 had %d Lynn names, keyed by Lynn's own names; prep has %d study species, keyed by the P1 species",
                "through Lynn's name. Unchanged (same name, same values): %d. Differences:"),
          nrow(on), nrow(nn), sum(niche$same)), "",
  lines_cause(diffs$niche_optima), "", md_table(select(diffs$niche_optima, key, difference, cause)),
  "", "## Found while comparing", "",
  sprintf(paste("Lynn names that are decided community names (`prep/config/taxonomy_decisions.csv`) but resolve otherwise in P4,",
                "because community decisions are not applied to Lynn's names (%d):"), nrow(lynn_decided)), "",
  if (nrow(lynn_decided)) md_table(lynn_decided) else "None.", "",
  "Lynn's names come from the same community records, so the decision applies to them as well; changing this is a P4 change (open question for Bills).",
  "", "## Columns", "",
  "Columns dropped from step 00 that a later step (MILESTONES steps 01 to 10, PLAN.md) names:", "",
  if (nrow(named_later)) md_table(select(named_later, file, column, now_in)) else "None.", "",
  sprintf("Other dropped columns (named by no later step): %s.",
          paste(sprintf("%s `%s`", dropped$file[!dropped$named_later], dropped$column[!dropped$named_later]), collapse = "; ")),
  "", "## Counts to confirm for steps 01 to 10", "",
  md_table(confirm), "",
  "The gradient count uses the SQLite's `gradient` per destination site; step 01 defines `gradient_id`. The number of effects (about 1,596) is confirmed in step 02.",
  "", "## Switch-over", "",
  sprintf("`grep -rnE \"_targets|legacy/\" analysis/` outside comments: %s. Steps 01 to 10 run the same check on their own script.",
          if (length(analysis_hits)) paste(analysis_hits, collapse = "; ") else "no hits (`analysis/` holds no script yet)"),
  sprintf("Code in `prep/` (outside this script) that reads step 00 files, the legacy store or `legacy/`: %s.",
          if (length(prep_hits)) paste(prep_hits, collapse = "; ") else "none"),
  paste("The P3 placement draft reuses legacy labels through `prep/config/legacy_taxonomy.csv` (step 00's submitted name and",
        "legacy accepted name, copied by this script). `legacy/` is read by nothing; prep comments cite its lines, and it stays",
        "until step 15 is closed (also on the tag `targets-pipeline-2026-10`)."),
  "", "## Warnings captured", "", "None.",
  "", "## Deviations from the plan", "",
  paste("(1) The comparison is a one-off script (`prep/compare_to_frozen.R`), not a target: it needs the step 00 files, which",
        "are archived after this run. (2) The approved P1 to P4 reports no longer carry their step 00 sections (decided 27 Oct);",
        "their content is in this report.")
)
writeLines(report, "reports/P5_compare_to_frozen.md")
