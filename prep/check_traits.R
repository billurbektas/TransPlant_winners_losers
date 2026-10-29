# P7: checks of the TRY trait records behind P2's medians, within TRY only.
# A one-off script outside the targets pipeline: it reports and changes
# nothing (CLAUDE.md: no audit-only targets in prep/). The rules Bills accepts
# after reading the report become P2 cleaning targets (decided 29 Oct).
#
# Run from the repository root:  Rscript prep/check_traits.R
# Reads:  data/derived/prep/traits_long.rds, traits.rds, functional_groups.rds
# Writes: reports/P7_trait_checks.md; reports/P7_trait_checks/flagged_records.csv,
#         flagged_datasets.csv, flagged_species.csv;
#         figures/diagnostics/P7_dataset_offsets.png
# Checks, all on the log scale (traits are positive and right-skewed):
#   1. records: robust z within species x trait (>= 5 records), |z| > 3.5; and
#      unit signatures, a value 10, 100 or 1000 times (or 1/...) the species
#      median within 15%
#   2. datasets: median over shared species of log(dataset median / median of
#      the other datasets); off by >= 1.5x on >= 3 shared species
#   3. species: median against the other species of its genus (family if the
#      genus has < 3 others), |robust z| > 3.5; plausibility rules for woody
#      height, graminoid LDMC and LMA-LDMC discordance
#   4. fragile medians: < 3 records or a single dataset
# Thresholds decided 29 Oct (Bills).

suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
source("prep/R/report.R") # md_table()

warns <- character()
globalCallingHandlers(warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
})
z_cut <- 3.5
dataset_factor <- 1.5
dataset_min_species <- 3L
unit_tolerance <- 0.15
out_dir <- "reports/P7_trait_checks"
dir.create(out_dir, showWarnings = FALSE)
dir.create("figures/diagnostics", recursive = TRUE, showWarnings = FALSE)

records <- readRDS("data/derived/prep/traits_long.rds") |>
  mutate(log_value = log(value))
medians <- readRDS("data/derived/prep/traits.rds")
groups <- readRDS("data/derived/prep/functional_groups.rds")
stopifnot(all(records$value > 0)) # P2 drops non-positive values
robust_z <- function(x, centre, spread) ifelse(spread > 0, (x - centre) / spread, NA_real_)

# ---- 1. Records within species x trait --------------------------------------
records <- records |>
  mutate(n_records = n(), species_median = median(log_value), mad = stats::mad(log_value),
         .by = c(species, trait)) |>
  mutate(z = if_else(n_records >= 5, robust_z(log_value, species_median, mad), NA_real_),
         ratio = exp(log_value - species_median),
         unit_signature = n_records >= 3 &
           Reduce(`|`, lapply(c(10, 100, 1000, 0.1, 0.01, 0.001),
                              function(f) abs(ratio / f - 1) <= unit_tolerance)),
         outlier = !is.na(z) & abs(z) > z_cut,
         flagged = outlier | unit_signature)
flagged_records <- records |>
  filter(flagged) |>
  transmute(species, trait, dataset_id, observation_id, SpeciesName, value = signif(value, 4),
            species_median = signif(exp(species_median), 4), ratio = signif(ratio, 3), z = round(z, 2),
            reason = case_when(outlier & unit_signature ~ "outlier and unit signature",
                               outlier ~ "outlier", TRUE ~ "unit signature")) |>
  arrange(desc(abs(log(ratio))))
# how far each species median moves without its flagged records
median_shift <- records |>
  summarise(median_all = median(value), median_without = median(value[!flagged]), n_flagged = sum(flagged),
            .by = c(species, trait)) |>
  filter(n_flagged > 0) |>
  mutate(shift_pct = round(100 * (median_without / median_all - 1), 1),
         across(c(median_all, median_without), ~ signif(.x, 3))) |>
  arrange(desc(abs(shift_pct)))
mad_zero <- records |> distinct(species, trait, n_records, mad) |> filter(n_records >= 5, mad == 0) |> nrow()

# ---- 2. Datasets against the other datasets of the same species -----------
by_dataset <- records |>
  summarise(dataset_median = median(log_value), n = n(), .by = c(trait, species, dataset_id))
others <- by_dataset |>
  select(trait, species, other_dataset = dataset_id) |>
  distinct()
offsets <- by_dataset |>
  inner_join(records |> select(trait, species, other = dataset_id, log_value), by = c("trait", "species"),
             relationship = "many-to-many") |>
  filter(other != dataset_id) |>
  summarise(other_median = median(log_value), .by = c(trait, species, dataset_id, dataset_median, n)) |>
  mutate(log_ratio = dataset_median - other_median)
dataset_summary <- offsets |>
  summarise(shared_species = n(), records = sum(n), offset = exp(median(log_ratio)),
            share_species_above_1.5x = round(mean(abs(log_ratio) >= log(dataset_factor)), 2),
            .by = c(trait, dataset_id)) |>
  mutate(flagged = shared_species >= dataset_min_species & abs(log(offset)) >= log(dataset_factor)) |>
  arrange(desc(abs(log(offset))))
feeds <- by_dataset |>
  left_join(summarise(records, total = n(), .by = c(trait, species)), by = c("trait", "species")) |>
  mutate(share_of_records = round(n / total, 2))
flagged_datasets <- dataset_summary |>
  filter(flagged) |>
  mutate(offset = signif(offset, 3)) |>
  left_join(feeds |> summarise(species_fed = n(), species_mostly_this_dataset = sum(share_of_records >= 0.5),
                               .by = c(trait, dataset_id)), by = c("trait", "dataset_id"))

# ---- 3. Species against their genus (family) and plausibility rules ---------
long_medians <- medians |>
  select(species, LMA, LDMC, height) |>
  tidyr::pivot_longer(-species, names_to = "trait", values_to = "median") |>
  filter(!is.na(median)) |>
  left_join(groups |> select(species, genus, family, graminoid, woody), by = "species") |>
  mutate(log_median = log(median))
relatives <- function(level) {
  long_medians |>
    inner_join(long_medians |> select(trait, other = species, other_log = log_median, !!level := all_of(level)),
               by = c("trait", level), relationship = "many-to-many") |>
    filter(other != species) |>
    summarise(n_relatives = n(), relatives_median = median(other_log), relatives_mad = stats::mad(other_log),
              .by = c(species, trait))
}
# a genus of 3 to 5 species can have almost no spread, which inflates z; the
# spread used is the larger of the group's own and the typical within-genus
# spread of that trait across all genera with at least 3 species
typical_spread <- long_medians |>
  filter(n() >= 3, .by = c(trait, genus)) |>
  summarise(spread = stats::mad(log_median), .by = c(trait, genus)) |>
  summarise(typical_spread = stats::median(spread), .by = trait)
species_check <- long_medians |>
  left_join(relatives("genus"), by = c("species", "trait")) |>
  left_join(rename_with(relatives("family"), ~ paste0("fam_", .x), -c(species, trait)), by = c("species", "trait")) |>
  left_join(typical_spread, by = "trait") |>
  mutate(level = if_else(coalesce(n_relatives, 0L) >= 3, "genus", "family"),
         z = if_else(level == "genus", robust_z(log_median, relatives_median, pmax(relatives_mad, typical_spread)),
                     if_else(coalesce(fam_n_relatives, 0L) >= 3,
                             robust_z(log_median, fam_relatives_median, pmax(fam_relatives_mad, typical_spread)), NA_real_)),
         relatives = if_else(level == "genus", n_relatives, fam_n_relatives))
herb_height_p05 <- stats::quantile(long_medians$median[long_medians$trait == "height" & !long_medians$woody], 0.05)
herb_ldmc_median <- stats::median(long_medians$median[long_medians$trait == "LDMC" & !long_medians$woody &
                                                        !long_medians$graminoid])
wide <- long_medians |> select(species, trait, log_median) |> tidyr::pivot_wider(names_from = trait, values_from = log_median)
q <- function(x, p) stats::quantile(x, p, na.rm = TRUE)
discordant <- wide |>
  filter((LMA > q(LMA, 0.9) & LDMC < q(LDMC, 0.1)) | (LMA < q(LMA, 0.1) & LDMC > q(LDMC, 0.9))) |>
  pull(species)
flagged_species <- bind_rows(
  species_check |> filter(!is.na(z), abs(z) > z_cut) |>
    transmute(species, trait, median = signif(median, 3), reason = paste0("far from its ", level, " (robust z ", round(z, 1),
                                                                         ", ", relatives, " relatives)")),
  long_medians |> filter(trait == "height", woody, median < herb_height_p05) |>
    transmute(species, trait, median = signif(median, 3), reason = "woody, shorter than 95% of herbs"),
  long_medians |> filter(trait == "LDMC", graminoid, median < herb_ldmc_median) |>
    transmute(species, trait, median = signif(median, 3), reason = "graminoid LDMC below the herb median"),
  long_medians |> filter(species %in% discordant, trait %in% c("LMA", "LDMC")) |>
    transmute(species, trait, median = signif(median, 3), reason = "LMA and LDMC at opposite extremes")) |>
  arrange(species, trait)

# ---- 4. Fragile medians ------------------------------------------------------
fragile <- records |>
  summarise(records = n(), datasets = n_distinct(dataset_id), .by = c(species, trait)) |>
  filter(records < 3 | datasets == 1)

# ---- Write -------------------------------------------------------------------
utils::write.csv(flagged_records, file.path(out_dir, "flagged_records.csv"), row.names = FALSE)
utils::write.csv(flagged_datasets, file.path(out_dir, "flagged_datasets.csv"), row.names = FALSE)
utils::write.csv(flagged_species, file.path(out_dir, "flagged_species.csv"), row.names = FALSE)
fig <- ggplot(dataset_summary, aes(x = log2(offset), y = trait, size = shared_species, colour = flagged)) +
  geom_vline(xintercept = c(-1, 1) * log2(dataset_factor), linetype = "dashed") +
  geom_jitter(height = 0.15, width = 0, alpha = 0.6) +
  scale_colour_manual(values = c(`FALSE` = "grey50", `TRUE` = "firebrick")) +
  labs(x = "log2(dataset median / other datasets' median), median over shared species", y = NULL,
       size = "Shared species", colour = "Flagged") +
  theme_classic()
ggsave("figures/diagnostics/P7_dataset_offsets.png", fig, width = 7, height = 3.5, dpi = 150)

# ---- Report --------------------------------------------------------------
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
euphrasia <- list(
  record = any(flagged_records$species == "Euphrasia officinalis" & flagged_records$trait == "LMA"),
  dataset = any(flagged_datasets$trait == "LMA" & flagged_datasets$dataset_id == 524),
  species = any(flagged_species$species == "Euphrasia officinalis" & flagged_species$trait == "LMA"))
counts <- tibble::tibble(
  check = c("records: outlier (|z| > 3.5)", "records: unit signature", "species x trait medians moving > 10% without flagged records",
            "datasets off by >= 1.5x", "species flagged against relatives or plausibility", "fragile medians (< 3 records or one dataset)"),
  n = c(sum(records$outlier), sum(records$unit_signature), sum(abs(median_shift$shift_pct) > 10),
        nrow(flagged_datasets), n_distinct(paste(flagged_species$species, flagged_species$trait)), nrow(fragile)))
report <- c(
  "# P7: Trait checks", "",
  paste0("Checks of the ", fmt(nrow(records)), " TRY records behind P2's medians (", n_distinct(records$species),
         " species; LDMC, LMA, height), within TRY only and on the log scale. Nothing is changed: the rules Bills accepts become ",
         "P2 cleaning targets (decided 29 Oct). Thresholds: robust z ", z_cut, "; datasets off by ", dataset_factor,
         "x on at least ", dataset_min_species, " shared species; unit signatures at 10, 100 and 1000x within ",
         100 * unit_tolerance, "%. Full lists in `reports/P7_trait_checks/`."),
  "", "## Counts", "", md_table(counts),
  "", sprintf("Species x trait groups with at least 5 records and zero spread (z undefined): %d.", mad_zero),
  sprintf("Species check: spread = the larger of the relatives' own spread and the typical within-genus spread (log scale; %s).",
          paste(sprintf("%s %.2f", typical_spread$trait, typical_spread$typical_spread), collapse = ", ")),
  "", sprintf("Known case (P2 notes): Euphrasia officinalis LMA from dataset 524 flagged at record level: %s; dataset level: %s; species level: %s.",
              euphrasia$record, euphrasia$dataset, euphrasia$species),
  "", "## Datasets off by 1.5x or more against the other datasets", "",
  md_table(flagged_datasets |> select(trait, dataset_id, offset, shared_species, records, species_fed, species_mostly_this_dataset)),
  "", "![Dataset offsets](../figures/diagnostics/P7_dataset_offsets.png)",
  "", "## Records: the 20 largest departures from their species median", "",
  md_table(head(flagged_records, 20)),
  "", "## Species medians that move most without their flagged records", "",
  md_table(head(median_shift, 20)),
  "", "## Species flagged against their relatives or by plausibility rules", "",
  md_table(head(flagged_species, 40)),
  if (nrow(flagged_species) > 40) sprintf("... %d more in `flagged_species.csv`.", nrow(flagged_species) - 40) else character(),
  "", "## What to decide", "",
  paste("1. Datasets: exclude a flagged dataset for that trait if its offset reflects a unit or protocol difference",
        "(the table above lists how many species medians it feeds and for how many it is the main source)."),
  "2. Records: drop outliers and unit-signature records as a rule, or only those that move a species median by more than 10%.",
  "3. Species: species flagged against their relatives are inspected by name; the plausibility flags are soft and not exclusions by themselves.",
  "", "## Warnings captured", "", if (length(warns)) paste("-", unique(warns)) else "None.",
  "", "## Deviations from PLAN.md", "", "None; a check of P2's inputs, requested 29 Oct."
)
writeLines(report, "reports/P7_trait_checks.md")
message("P7 written: ", sum(records$flagged), " records, ", nrow(flagged_datasets), " datasets, ",
        n_distinct(paste(flagged_species$species, flagged_species$trait)), " species x trait flagged")
