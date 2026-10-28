# P6: species climatic niche from our own GBIF download. One script, outside
# the targets pipeline: the download is a frozen input with a DOI, not a
# target that reruns when a ledger changes. The pipeline (P4) reads its result.
#
# Method: Lynn's SpeciesClimate&RangeAreaDet_TransPlant.R (MacroContext_Transplant,
# 2024; lines cited below), with the changes decided 28 Oct (Bills):
#   - GBIF keys chosen from the P1 species with their WCVP author and the
#     reviewed ledger prep/config/gbif_keys.csv (Lynn: names from the Global
#     Names Resolver), and records grouped by our species (Lynn: GBIF's
#     `species` column);
#   - P_opt from CHELSA bio18, precipitation of the warmest quarter (Lynn: bio12);
#   - LIVING_SPECIMEN records dropped as well as fossils;
#   - thinning on the CHELSA grid itself, one value per species x climate cell,
#     which is what Lynn et al. 2021 (GCB, Methods 2.2) describe; Lynn's 0.0083
#     degree grid laid from the point extent drifts against the 30" grid.
#
# Run from the repository root, after tar_make has built P1-P3:
#   Rscript prep/gbif_niche.R keys   # stage 1 only (no GBIF credentials needed)
#   Rscript prep/gbif_niche.R        # all stages; every stage resumes from its cache
#
# Reads:  the prep store (study_species, taxonomy, name_resolution);
#         prep/config/gbif_keys.csv; data/MacroContextTransplantDat_2.csv (comparison only)
# Writes: data/cache/gbif/, data/cache/chelsa_bioclim/, data/cache/naturalearth/;
#         data/derived/prep/gbif_keys_used.csv, gbif_download.csv, gbif_niche.rds and .csv;
#         reports/P6_gbif_niche.md, reports/P6_gbif_download_request.json
# Checks: 904 study species; every non-automatic match decided in the ledger; every
#         key claimed by one species; every species keyed or excluded; GBIF records
#         in a study country or confirmed; rows read = GBIF's total; every record
#         assigned to a species; bio10 within -30 to 40 degC; 903 species out

suppressPackageStartupMessages(library(dplyr))
source("prep/R/report.R") # md_table()

keys_only <- identical(commandArgs(trailingOnly = TRUE), "keys")
backbone <- "d7dddbf4-2cf0-4f39-9b2a-bb099caae36c" # GBIF Backbone, as Lynn's download; rgbif 3.9 defaults to COL XR
cache <- "data/cache/gbif"
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
warns <- character()
capture <- function(expr) withCallingHandlers(expr, warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
})
store <- function(x) targets::tar_read_raw(x, store = "prep/_targets")
# GBIF's API slows down or times out when busy; retry with growing pauses
gbif_retry <- function(call) {
  for (attempt in 1:5) {
    out <- tryCatch(capture(call()), error = function(e) e)
    if (!inherits(out, "error")) return(out)
    message("GBIF call failed (", conditionMessage(out), "); retry ", attempt)
    Sys.sleep(30 * attempt)
  }
  stop("GBIF call failed five times; rerun to resume")
}

# ---- Stage 1: names -> GBIF keys ---------------------------------------------
study <- store("study_species")
taxonomy <- store("taxonomy")
res <- store("name_resolution")
stopifnot(nrow(study) == 903) # P3: residents and destination pool (904 before the 28 Oct homonym fixes)

ledger_file <- "prep/config/gbif_keys.csv"
ledger <- if (file.exists(ledger_file)) utils::read.csv(ledger_file, colClasses = "character", encoding = "UTF-8") else
  tibble::tibble(species = character(), name = character(), author = character(), gbif_key = character(),
                 gbif_name = character(), decision = character(), reason = character())
stopifnot(all(ledger$decision %in% c("accepted", "synonym_only", "gbif_subtaxon", "include", "exclude", "author", "geography_ok")))

# the WCVP author TNRS returned with the accepted name; the most frequent one
# where TNRS gave several, unless the ledger names it
authors <- res |>
  filter(accepted_name %in% study$species, !is.na(accepted_author)) |>
  count(species = accepted_name, author = accepted_author) |>
  arrange(species, desc(n), author) |>
  distinct(species, .keep_all = TRUE) |>
  select(species, author)
authors <- left_join(select(study, species), authors, by = "species")
fix <- filter(ledger, decision == "author")
authors$author[match(fix$species, authors$species)] <- fix$author
n_legit <- res |>
  filter(accepted_name %in% study$species, !is.na(n_legitimate)) |>
  summarise(n_legitimate = max(n_legitimate), .by = accepted_name) |>
  rename(species = accepted_name)

community <- taxonomy |>
  filter(accepted_species %in% study$species) |>
  mutate(name = stringr::str_squish(name_original))
others <- community |> filter(name != accepted_species) |> distinct(species = accepted_species, name)
countries <- community |>
  tidyr::separate_rows(regions, sep = "; ") |>
  distinct(species = accepted_species, country = substr(regions, 1, 2))

# matched once and cached; rematched only when the names or authors change.
# Each accepted name is matched twice: with its author, and bare, because GBIF
# can hold one name under two usages that split its records
query <- bind_rows(transmute(authors, species, name = species, author, role = "accepted"),
                   transmute(authors, species, name = species, author = NA_character_, role = "bare"),
                   mutate(others, author = NA_character_, role = "other"))
match_file <- file.path(cache, "name_matches.rds")
if (file.exists(match_file) && identical(readRDS(match_file)$query, query)) {
  matches <- readRDS(match_file)$matches
} else {
  m <- capture(rgbif::name_backbone_checklist(
    data.frame(scientificName = query$name, scientificNameAuthorship = query$author),
    kingdom = "Plantae", strict = TRUE, checklistKey = backbone))
  stopifnot(nrow(m) == nrow(query))
  matches <- bind_cols(query, transmute(m, matchType, status, rank, usageKey = as.character(usageKey),
                                        gbif_name = scientificName, gbif_species = species,
                                        speciesKey = as.character(speciesKey)))
  saveRDS(list(query = query, matches = matches), match_file)
}
ledger_rows <- filter(ledger, !decision %in% c("author", "geography_ok")) |>
  select(species, name, ledger_key = gbif_key, decision, reason)

# Our accepted name is automatic when GBIF files it under a species of the
# same name; every other result is decided in the ledger
bare <- filter(matches, role == "bare") |>
  transmute(species, bare_key = speciesKey, bare_name = gbif_name,
            bare_ok = matchType == "EXACT" & !is.na(speciesKey) & gbif_species == species)
acc <- filter(matches, role == "accepted") |>
  left_join(bare, by = "species") |>
  mutate(class = case_when(
    matchType == "EXACT" & !is.na(speciesKey) & gbif_species == species ~ paste("automatic:", tolower(status)),
    matchType == "EXACT" & !is.na(speciesKey) ~ "ledger: GBIF files it under another species",
    TRUE ~ "ledger: no exact match")) |>
  left_join(ledger_rows, by = c("species", "name"))
auto <- filter(acc, startsWith(class, "automatic"))
# whole-species keys take every record GBIF files under that species; where
# GBIF holds our author's usage only provisionally or as a synonym, the
# same-named usage it accepts holds records of the same name and is added
# (Bistorta officinalis Delarbre 99,494 records, "Raf." 139,030)
species_keys <- bind_rows(
  transmute(auto, species, name, gbif_key = speciesKey, how = class),
  filter(auto, status != "ACCEPTED", bare_ok, bare_key != speciesKey) |>
    transmute(species, name, gbif_key = bare_key, how = "automatic: same-named accepted usage"),
  filter(acc, decision == "accepted") |>
    transmute(species, name, gbif_key = coalesce(na_if(ledger_key, ""), speciesKey), how = "ledger: accepted"))

oth <- filter(matches, role == "other") |>
  mutate(class = case_when(
    !matchType %in% c("EXACT", "VARIANT") | is.na(speciesKey) ~ "ignored: no exact match",
    paste(species, speciesKey) %in% paste(species_keys$species, species_keys$gbif_key) ~ "nothing to do: same GBIF species",
    TRUE ~ "ledger: GBIF keeps it as another species")) |>
  left_join(ledger_rows, by = c("species", "name"))

decide <- bind_rows(filter(acc, startsWith(class, "ledger")), filter(oth, startsWith(class, "ledger")))
undecided <- filter(decide, is.na(decision))
if (nrow(undecided)) {
  utils::write.csv(transmute(undecided, species, name, author, gbif_key = usageKey, gbif_name, gbif_species,
                             gbif_status = status, match = matchType, decision = "", reason = ""),
                   file.path(cache, "ledger_needed.csv"), row.names = FALSE)
  stop(nrow(undecided), " GBIF matches need a decision in ", ledger_file, "; listed in ", cache, "/ledger_needed.csv")
}

# name-only keys take only the records labelled with that exact name
keys <- bind_rows(
  mutate(species_keys, mode = "species"),
  filter(oth, decision == "include") |>
    transmute(species, name, gbif_key = coalesce(na_if(ledger_key, ""), speciesKey), mode = "species", how = "ledger: include"),
  filter(decide, decision == "synonym_only") |>
    transmute(species, name, gbif_key = coalesce(na_if(ledger_key, ""), usageKey), mode = "name", how = "ledger: synonym_only"))

# gbif_subtaxon (decided 28 Oct for Ranunculus tuberosus). GBIF files our
# species as a taxon below another study species: R. tuberosus Lapeyr. sits in
# R. polyanthemos subsp. nemorosus, together with R. nemorosus DC. and the other
# names GBIF lists as its synonyms, while WCVP gives every one of those names
# to R. tuberosus. Taking only the records labelled "R. tuberosus" would leave
# the rest of the plant's records in R. polyanthemos, which takes its whole
# GBIF species. So the GBIF subtaxon is requested, and every record filed under
# it, by any of its names, is assigned to our species by the record's own name
# key; that carves it out of the parent species.
sub <- filter(decide, decision == "gbif_subtaxon")
sub_names <- bind_rows(lapply(seq_len(nrow(sub)), function(i) {
  k <- as.numeric(sub$ledger_key[i])
  syn <- capture(rgbif::name_usage(key = k, data = "synonyms", limit = 1000))$data
  kids <- capture(rgbif::name_usage(key = k, data = "children", limit = 1000))$data
  stopifnot(NROW(kids) == 0) # names ranked below the subtaxon would need their own synonyms
  # every record filed under the subtaxon must carry one of these name keys
  usage <- capture(rgbif::name_usage(key = k))$data
  higher <- as.character(unlist(usage[, intersect(c("kingdomKey", "phylumKey", "classKey", "orderKey", "familyKey",
                                                    "genusKey", "speciesKey"), names(usage))]))
  by_key <- gbif_retry(function() rgbif::occ_count(taxonKey = k, hasCoordinate = TRUE, hasGeospatialIssue = FALSE,
                                                   year = "1950,2026", facet = "taxonKey", facetLimit = 200,
                                                   curlopts = list(http_version = 2, timeout = 120)))
  names_k <- c(as.character(k), as.character(syn$key))
  stray <- setdiff(by_key$taxonKey, c(higher, names_k))
  if (length(stray)) stop("Records under GBIF key ", k, " carry name keys not listed as its synonyms: ", paste(stray, collapse = ", "))
  tibble::tibble(species = sub$species[i], name = c(sub$name[i], syn$scientificName), gbif_key = names_k, mode = "name",
                 how = c("ledger: gbif_subtaxon", rep("ledger: gbif_subtaxon, name filed under it", nrow(syn))),
                 records = by_key$count[match(names_k, by_key$taxonKey)])
}))
keys <- bind_rows(keys, select(sub_names, -records))
stopifnot(!anyNA(keys$gbif_key)) # a match without a GBIF key needs the key in the ledger
excluded <- filter(decide, decision == "exclude", name == species)$species
stopifnot(setequal(c(keys$species, excluded), study$species))   # every species keyed or excluded
claimed <- distinct(keys, gbif_key, mode, species) |> count(gbif_key, mode)
if (any(claimed$n > 1)) stop("GBIF keys claimed by two species: ", paste(claimed$gbif_key[claimed$n > 1], collapse = ", "))

# Geographic check: the chosen keys must have at least 10 records in the
# countries the species was recorded in, so a wrong homonym cannot pass on a
# stray record (Primula auricula passed with 1 Chinese record of 26,539)
count_file <- file.path(cache, "country_counts.rds")
key_set <- sort(unique(keys$gbif_key[keys$how != "ledger: gbif_subtaxon, name filed under it"]))
# cached per key, so an interrupted run (GBIF throttles repeated requests) resumes
by_country <- if (file.exists(count_file)) readRDS(count_file) else
  tibble::tibble(gbif_key = character(), country = character(), count = integer())
if (!is.data.frame(by_country)) by_country <- by_country$counts # cache written before 28 Oct
for (k in setdiff(key_set, by_country$gbif_key)) {
  x <- gbif_retry(function() rgbif::occ_count(taxonKey = as.numeric(k), hasCoordinate = TRUE, hasGeospatialIssue = FALSE,
                                              year = "1950,2026", facet = "country", facetLimit = 300,
                                              curlopts = list(http_version = 2, timeout = 120)))
  by_country <- bind_rows(by_country, if (nrow(x)) tibble::tibble(gbif_key = k, country = x$country, count = as.integer(x$count))
                          else tibble::tibble(gbif_key = k, country = NA_character_, count = 0L))
  saveRDS(by_country, count_file)
}
by_country <- filter(by_country, gbif_key %in% key_set)
geo <- keys |>
  distinct(species, gbif_key) |>
  left_join(by_country, by = "gbif_key", relationship = "many-to-many") |>
  summarise(records = sum(count, na.rm = TRUE),
            in_study_countries = sum(count[country %in% countries$country[countries$species == first(species)]], na.rm = TRUE),
            .by = species)
far <- filter(geo, in_study_countries < 10, !species %in% ledger$species[ledger$decision == "geography_ok"])
if (nrow(far)) {
  utils::write.csv(transmute(far, species, name = species, author = authors$author[match(species, authors$species)],
                             gbif_key = "", gbif_name = "", decision = "geography_ok", reason = "", records),
                   file.path(cache, "geography_needed.csv"), row.names = FALSE)
  stop(nrow(far), " species have fewer than 10 GBIF records in their study countries; confirm in ", ledger_file,
       " (listed in ", cache, "/geography_needed.csv)")
}

keys_used <- keys |>
  mutate(gbif_name = if_else(mode == "name", matches$gbif_name[match(gbif_key, matches$usageKey)],
                             matches$gbif_species[match(gbif_key, matches$speciesKey)]),
         gbif_name = coalesce(gbif_name, ledger$gbif_name[match(gbif_key, ledger$gbif_key)])) |>
  left_join(summarise(by_country, records = sum(count), .by = gbif_key), by = "gbif_key") |>
  arrange(species, name)
utils::write.csv(keys_used, "data/derived/prep/gbif_keys_used.csv", row.names = FALSE)
message("Stage 1: ", n_distinct(keys$species), " species keyed with ", length(key_set), " GBIF keys; ",
        length(excluded), " excluded")
if (keys_only) quit(save = "no")

# ---- Stage 2: the GBIF download -----------------------------------------------
key_file <- file.path(cache, "download_key.txt")
record_file <- "data/derived/prep/gbif_download.csv"
request_file <- "reports/P6_gbif_download_request.json"
if (!file.exists(key_file)) {
  stopifnot(all(nzchar(Sys.getenv(c("GBIF_USER", "GBIF_PWD", "GBIF_EMAIL"))))) # from ~/.Renviron
  # Lynn's download predicates (download 0068537-240626123714530): coordinates,
  # no geospatial issue, 1950 onwards; simple CSV
  dl <- rgbif::occ_download(rgbif::pred_in("taxonKey", as.numeric(key_set)), rgbif::pred("hasCoordinate", TRUE),
                            rgbif::pred("hasGeospatialIssue", FALSE), rgbif::pred_gte("year", 1950),
                            format = "SIMPLE_CSV", checklistKey = backbone)
  writeLines(as.character(dl), key_file) # saved before waiting: a rerun never requests twice
}
download_key <- readLines(key_file)
meta <- rgbif::occ_download_meta(download_key)
jsonlite::write_json(meta$request, request_file, auto_unbox = TRUE, pretty = TRUE)
requested <- sort(as.character(unlist(Filter(function(p) identical(p$key, "TAXON_KEY"), meta$request$predicate$predicates)[[1]]$values)))
if (!identical(requested, key_set)) stop("The keys now chosen differ from the download ", download_key,
                                         "; remove ", key_file, " only if a new download is intended")
if (meta$status != "SUCCEEDED") meta <- rgbif::occ_download_wait(download_key, status_ping = 60, quiet = FALSE)
meta <- rgbif::occ_download_meta(download_key)
zip <- file.path(cache, paste0(download_key, ".zip"))
if (!file.exists(zip)) rgbif::occ_download_get(download_key, path = cache)
sha <- digest::digest(file = zip, algo = "sha256")
if (file.exists(record_file) && utils::read.csv(record_file)$sha256 != sha) stop(zip, " differs from the recorded download")
record <- tibble::tibble(
  download_key, doi = meta$doi, url = paste0("https://doi.org/", meta$doi),
  citation = sprintf("GBIF.org (%s) GBIF Occurrence Download https://doi.org/%s",
                     format(as.Date(substr(meta$created, 1, 10)), "%d %B %Y"), meta$doi),
  created = meta$created, status = meta$status, total_records = meta$totalRecords, size_bytes = file.size(zip),
  sha256 = sha, erase_after = meta$eraseAfter, checklist = paste("GBIF Backbone", backbone),
  rgbif = as.character(utils::packageVersion("rgbif")), n_keys = length(key_set))
utils::write.csv(record, record_file, row.names = FALSE)

# ---- Stage 3: CHELSA climatologies and the land polygons ----------------------
clim_dir <- "data/cache/chelsa_bioclim"
dir.create(clim_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = 7200)
for (v in c("bio10", "bio18")) {
  url <- sprintf("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_%s_1981-2010_V.2.1.tif", v)
  dest <- file.path(clim_dir, basename(url))
  size <- as.numeric(sub(".*: ", "", grep("^content-length", curlGetHeaders(url), ignore.case = TRUE, value = TRUE)))
  if (!file.exists(dest) || file.size(dest) != size) utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  stopifnot(file.size(dest) == size)
}
clim <- terra::rast(file.path(clim_dir, c("CHELSA_bio10_1981-2010_V.2.1.tif", "CHELSA_bio18_1981-2010_V.2.1.tif")))
names(clim) <- c("bio10", "bio18")
# clean_coordinates() defaults Lynn ran (lines 40-46): seas against Natural
# Earth 1:50 m land, zeros within 0.5 degree, equal coordinates absolute
ne_dir <- "data/cache/naturalearth"
dir.create(ne_dir, showWarnings = FALSE)
if (!length(list.files(ne_dir, pattern = "land[.]shp$")))
  rnaturalearth::ne_download(scale = 50, type = "land", category = "physical", destdir = ne_dir, load = FALSE)
land <- terra::vect(rnaturalearth::ne_load(scale = 50, type = "land", category = "physical", destdir = ne_dir,
                                           returnclass = "sf"))

# ---- Stage 4: stream the archive, clean, thin --------------------------------
cells_file <- file.path(cache, paste0("species_cells_", download_key, ".rds"))
if (!file.exists(cells_file)) {
  csv <- utils::unzip(zip, list = TRUE)$Name[1]
  header <- strsplit(readLines(pipe(sprintf("unzip -p %s %s | head -1", shQuote(zip), shQuote(csv)))), "\t")[[1]]
  cols <- c("taxonKey", "speciesKey", "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters",
            "basisOfRecord", "occurrenceStatus", "individualCount")
  idx <- sort(match(cols, header))
  stopifnot(!anyNA(idx))
  # read from the zip without unpacking it (about 55 GB unpacked)
  con <- pipe(sprintf("unzip -p %s %s | cut -f %s", shQuote(zip), shQuote(csv), paste(idx, collapse = ",")), "r")
  invisible(readLines(con, n = 1))
  name_keys <- filter(keys, mode == "name")
  species_keys <- filter(keys, mode == "species")
  steps <- c("read", "unassigned", "fossil or living specimen", "absent", "uncertainty > 10 km",
             "individual count 0", "equal coordinates", "zero coordinates", "sea", "kept")
  flow <- setNames(numeric(length(steps)), steps)
  parts <- list()
  repeat {
    lines <- readLines(con, n = 5e6)
    if (!length(lines)) break
    d <- data.table::fread(text = lines, sep = "\t", header = FALSE, quote = "", col.names = header[idx],
                           colClasses = list(character = c("taxonKey", "speciesKey", "basisOfRecord", "occurrenceStatus")),
                           na.strings = "")
    flow["read"] <- flow["read"] + nrow(d)
    # records under a name-only key go to that species, all others by GBIF species
    d$species <- coalesce(name_keys$species[match(d$taxonKey, name_keys$gbif_key)],
                          species_keys$species[match(d$speciesKey, species_keys$gbif_key)])
    flow["unassigned"] <- flow["unassigned"] + sum(is.na(d$species))
    # Lynn's filters (lines 29-37), plus living specimens (decided 28 Oct)
    keep <- function(d, ok, step) { flow[step] <<- flow[step] + sum(!ok); d[ok, ] }
    d <- d[!is.na(d$species) & !is.na(d$decimalLatitude) & !is.na(d$decimalLongitude), ]
    d <- keep(d, !d$basisOfRecord %in% c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN"), "fossil or living specimen")
    d <- keep(d, d$occurrenceStatus == "PRESENT", "absent")
    d <- keep(d, is.na(d$coordinateUncertaintyInMeters) | d$coordinateUncertaintyInMeters / 1000 <= 10, "uncertainty > 10 km")
    d <- keep(d, is.na(d$individualCount) | d$individualCount > 0, "individual count 0")
    d <- as.data.frame(d)
    d <- keep(d, CoordinateCleaner::cc_equ(d, lon = "decimalLongitude", lat = "decimalLatitude", test = "absolute",
                                           value = "flagged", verbose = FALSE), "equal coordinates")
    d <- keep(d, CoordinateCleaner::cc_zero(d, lon = "decimalLongitude", lat = "decimalLatitude", buffer = 0.5,
                                            value = "flagged", verbose = FALSE), "zero coordinates")
    d <- keep(d, capture(CoordinateCleaner::cc_sea(d, lon = "decimalLongitude", lat = "decimalLatitude", ref = land,
                                                   value = "flagged", verbose = FALSE)), "sea")
    flow["kept"] <- flow["kept"] + nrow(d)
    # one value per species x CHELSA cell: Lynn's thinning (lines 91-109) on the climate grid itself
    d$cell <- terra::cellFromXY(clim, cbind(d$decimalLongitude, d$decimalLatitude))
    parts[[length(parts) + 1]] <- count(d, species, cell, name = "n_records")
    message("read ", format(flow[["read"]], big.mark = ","))
  }
  close(con)
  cells <- bind_rows(parts) |> summarise(n_records = sum(n_records), .by = c(species, cell))
  saveRDS(list(cells = cells, flow = flow), cells_file)
}
thinned <- readRDS(cells_file)
flow <- thinned$flow
stopifnot(flow[["read"]] == record$total_records) # every row of the archive read
stopifnot(flow[["unassigned"]] == 0)               # every record came in under one of our keys

# ---- Stage 5: climate per cell, niche per species ---------------------------
cell_ids <- unique(stats::na.omit(thinned$cells$cell)) # NA: outside the CHELSA extent
vals <- cbind(cell = cell_ids, terra::extract(clim, cell_ids))
cells <- left_join(thinned$cells, vals, by = "cell")
no_value <- filter(cells, is.na(bio10) | is.na(bio18))
cells <- filter(cells, !is.na(bio10), !is.na(bio18))
stopifnot(all(cells$bio10 > -30 & cells$bio10 < 40)) # degC, scale and offset applied
niche <- cells |>
  summarise(T_opt = mean(bio10), P_opt = log(mean(bio18)), mean_bio18 = mean(bio18),
            n_obs_sub = n(), n_records = sum(n_records), .by = species)
niche <- left_join(select(study, species), niche, by = "species") |>
  mutate(n_obs_sub = coalesce(n_obs_sub, 0L), n_records = coalesce(n_records, 0L))
stopifnot(nrow(niche) == nrow(study), !anyDuplicated(niche$species))
saveRDS(niche, "data/derived/prep/gbif_niche.rds")
utils::write.csv(niche, "data/derived/prep/gbif_niche.csv", row.names = FALSE)

# ---- Comparison with Lynn's niche (joined through the community taxonomy) ---
lynn <- utils::read.csv("data/MacroContextTransplantDat_2.csv") |>
  distinct(lynn_name = species, lynn_T = mean_sum_sub, lynn_cells = n_obs_sub) |>
  filter(!is.na(lynn_T)) |>
  mutate(species = community$accepted_species[match(stringr::str_squish(lynn_name), community$name)])
n_lynn_only <- sum(is.na(lynn$species))
lynn <- lynn |> filter(!is.na(species)) |> arrange(species, desc(lynn_cells)) |> distinct(species, .keep_all = TRUE)
cmp <- inner_join(lynn, filter(niche, n_obs_sub > 0), by = "species") |>
  mutate(difference = round(T_opt - lynn_T, 2),
         cause = if_else(lynn_name == species, "same name: new records, keys and cleaning",
                         "Lynn's niche from another GBIF name"))
gained <- filter(niche, n_obs_sub > 0, !species %in% lynn$species)

# ---- Report --------------------------------------------------------------
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
class_table <- bind_rows(mutate(count(acc, class), names = "accepted name"),
                         mutate(count(oth, class), names = "other community names")) |>
  select(names, class, n)
ledger_table <- filter(ledger, decision != "geography_ok") |> select(species, name, author, gbif_key, gbif_name, decision, reason)
multi_author <- n_legit |> filter(n_legitimate > 1) |> left_join(authors, by = "species") |>
  left_join(geo, by = "species") |> arrange(species)
flow_table <- tibble::tibble(step = names(flow), records = fmt(unname(flow)))
report <- c(
  "# P6: Niche from GBIF", "",
  paste0("Species climatic optima from our own GBIF download (", record$citation, "), keyed by the ",
         nrow(study), " P1 study species with their WCVP author and the reviewed ledger `prep/config/gbif_keys.csv`, ",
         "cleaned with Lynn's rules plus living specimens, thinned to one value per species x CHELSA 30\" cell, and ",
         "summarised as T_opt = mean CHELSA 1981-2010 bio10 and P_opt = log mean bio18. Replaces Lynn's niche file ",
         "(decided 28 Oct). The download record is `data/derived/prep/gbif_download.csv` and the exact request ",
         "`reports/P6_gbif_download_request.json`."),
  "", "## The download", "",
  md_table(transmute(record, download_key, doi, total_records = fmt(total_records), size_GB = round(size_bytes / 1e9, 1),
                     created = substr(created, 1, 10), erase_after = substr(erase_after, 1, 10), sha256 = substr(sha256, 1, 16))),
  "", paste("GBIF deletes downloads nobody cites after `erase_after`; the DOI goes into the methods and data",
            "availability text, so the download is kept."),
  "", "## Keys", "",
  md_table(class_table), "",
  sprintf("%d species keyed with %d GBIF keys; %d excluded.", n_distinct(keys$species), length(key_set), length(excluded)),
  "", sprintf(paste("Same-named usages added (%d species): GBIF holds our author's usage only provisionally or as a",
                    "synonym and accepts another usage of the same name, so the records of the name are split between",
                    "the two keys; both are requested:"), sum(species_keys$how == "automatic: same-named accepted usage")),
  "", md_table(filter(keys_used, how == "automatic: same-named accepted usage") |> select(species, gbif_key, gbif_name, records)),
  "", sprintf("### Ledger decisions (%d)", nrow(ledger_table)), "", md_table(ledger_table),
  "", "### Ranunculus tuberosus and Ranunculus polyanthemos (`gbif_subtaxon`, decided 28 Oct)", "",
  paste("GBIF and WCVP cut these two plants differently. GBIF files R. tuberosus Lapeyr. as a synonym of",
        "R. polyanthemos subsp. nemorosus, under which it also files R. nemorosus DC. and the other names below;",
        "WCVP (our P1 taxonomy) gives every one of these names to R. tuberosus and keeps R. polyanthemos for the rest.",
        "Taking only the records labelled \"R. tuberosus\" would have given R. tuberosus a minority of its own records",
        "and left the rest in R. polyanthemos, which takes its whole GBIF species. So GBIF's subspecies nemorosus is",
        "requested, and every record filed under it, by any of its names, is assigned to R. tuberosus by the record's",
        "own name key and so carved out of R. polyanthemos. The script stops if a record under the subspecies carries a",
        "name key not listed here. Records with coordinates, no geospatial issue, 1950 onwards, before cleaning:"),
  "", md_table(filter(sub_names, !is.na(records)) |> arrange(desc(records)) |>
                 transmute(species, gbif_name = name, gbif_key, how, records = fmt(records))),
  "", sprintf(paste("The first row counts all names together: %s records go to R. tuberosus, of which %s are labelled",
                    "\"Ranunculus tuberosus Lapeyr.\" (the only ones a name-only key would have taken)."),
              fmt(sum(sub_names$records[sub_names$how == "ledger: gbif_subtaxon"], na.rm = TRUE)),
              fmt(sum(sub_names$records[grepl("^Ranunculus tuberosus", sub_names$name)], na.rm = TRUE))),
  "", sprintf("### Species whose names had more than one legitimate author in TNRS (%d)", nrow(multi_author)), "",
  md_table(transmute(multi_author, species, author, records = fmt(records), in_study_countries = fmt(in_study_countries))),
  "", sprintf("Geographic check: %d of %d keyed species have at least 10 GBIF records in their study countries; %d confirmed in the ledger.",
              sum(geo$in_study_countries >= 10), nrow(geo), sum(ledger$decision == "geography_ok")),
  "", "## Cleaning and thinning", "", md_table(flow_table), "",
  sprintf("%s species x cell values kept; %s cells without a CHELSA value dropped.",
          fmt(nrow(cells)), fmt(nrow(no_value))),
  "", "## Niche", "",
  sprintf("Study species with a niche: %d of %d; with n_obs_sub >= 30: %d. Median n_obs_sub %s.",
          sum(niche$n_obs_sub > 0), nrow(niche), sum(niche$n_obs_sub >= 30), fmt(stats::median(niche$n_obs_sub))),
  "", "## Comparison with Lynn's niche", "",
  sprintf(paste("Lynn's names joined through the community taxonomy (%d names; %d hand-mapped Lynn-only names left out).",
                "T_opt against Lynn's mean_sum_sub on %d shared species: r = %.3f. Species with a niche now and none",
                "in Lynn's file: %d."), nrow(lynn) + n_lynn_only, n_lynn_only, nrow(cmp), stats::cor(cmp$T_opt, cmp$lynn_T), nrow(gained)),
  "", "Largest differences:", "",
  md_table(cmp |> arrange(desc(abs(difference))) |> head(15) |>
             transmute(species, lynn_name, lynn_cells, our_cells = n_obs_sub, lynn_T = round(lynn_T, 2),
                       T_opt = round(T_opt, 2), difference, cause)),
  "", sprintf("Koeleria macrantha: n_obs_sub = %s, T_opt = %.2f.", fmt(niche$n_obs_sub[niche$species == "Koeleria macrantha"]),
              niche$T_opt[niche$species == "Koeleria macrantha"]),
  "", "## Warnings captured", "",
  if (length(warns)) paste("-", unique(warns)) else "None.",
  "", "## Deviations from PLAN.md", "",
  paste("Against Lynn's script, all decided 28 Oct: keys from our species with author and ledger; records grouped by our",
        "species; bio18 for P_opt; living specimens dropped; thinning on the CHELSA grid (Lynn et al. 2021: \"we created a",
        "raster grid with cell length equal to 0.0083 deg and sampled a single representative point from the occurrence",
        "data per grid cell. The grid resolution approximates the 30 s resolution of the climate data\");",
        "n_obs_sub counts cells with a climate value.")
)
writeLines(report, "reports/P6_gbif_niche.md")
message("P6 written: ", sum(niche$n_obs_sub > 0), " species with a niche")
