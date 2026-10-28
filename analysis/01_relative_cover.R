# Step 01: relative cover table and species pools.
#
# PLAN.md: Step A "Data rules"; The colonizer story (destination pool);
#   Appendix Fig S1 (flow counts).
# Reads:  data/derived/prep/community.rds, experiments.rds, species_taxonomy.rds
# Writes: data/derived/01_relative_cover/cover_long.rds, a0.rds,
#         experiments_kept.csv, pool_sizes.csv, flow_counts.csv;
#         figures/diagnostics/01_pool_sizes.png; reports/01_relative_cover.md
# Checks: relative cover sums to 1 per plot-year over all vascular records;
#         the six deep-lineage taxa absent afterwards; 43 included experiments;
#         no species in both pools; 43 experiments and 21 gradients after the
#         pool-size rule (k = 5; US_Arizona kept, decided 28 Oct); no duplicate experiment x plot x year x
#         species x treatment; every pool species has both A0 values.
#
# Decided 28 Oct (Bills): the observation filter (>= 2 plot-year records)
# decides pool membership only and all observed cover is kept; pools over all
# survey years; k = 5 resident species; A0 from the first available survey of
# each control group, log of the plot mean in % + 1; gradient_id = region;
# relative cover over all vascular records, deep-lineage taxa included, and
# those taxa removed afterwards; destination pool as PLAN.md (present in
# destination controls, absent from origin controls, whether or not ever
# found in warmed plots).

suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
source("R/md_table.R")

warns <- character()
globalCallingHandlers(warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
})
out_dir <- "data/derived/01_relative_cover"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("figures/diagnostics", recursive = TRUE, showWarnings = FALSE)

deep_lineage <- c("Gymnocarpium dryopteris", "Athyrium spinulosum", "Botrychium lunaria",
                  "Equisetum arvense", "Selaginella selaginoides", "Juniperus communis")
k <- 5L       # minimum resident pool size (decided 28 Oct)
min_obs <- 2L # plot-year records for a species to count as present in a group

community <- readRDS("data/derived/prep/community.rds")
experiments <- readRDS("data/derived/prep/experiments.rds")
taxonomy <- readRDS("data/derived/prep/species_taxonomy.rds")
stopifnot(all(community$is_vegetation), all(community$cover > 0))
stopifnot(all(deep_lineage %in% taxonomy$accepted_species)) # names as in the P1 taxonomy

# ---- Relative cover per plot-year, over every vascular record ---------------
# The denominator is all vascular cover of the plot-year: identified, genus-level
# and unidentified records, the deep-lineage taxa included (decided 28 Oct).
community <- community |>
  mutate(relative_cover = cover / sum(cover), .by = plot_year_id)
sums <- summarise(community, total = sum(relative_cover), .by = plot_year_id)
stopifnot(all(abs(sums$total - 1) < 1e-9)) # relative cover sums to 1 per plot-year

# analysis rows: identified species, no drop_reason, deep-lineage taxa removed
is_analysis <- !is.na(community$species) & is.na(community$drop_reason)
community$not_analysed <- !is_analysis | community$species %in% deep_lineage
dropped_share <- summarise(community, share = sum(relative_cover[not_analysed]), .by = c(region, plot_year_id))
dropped_by_region <- dropped_share |>
  summarise(plot_years = n(), median_share = round(stats::median(share), 3), max_share = round(max(share), 3),
            plot_years_over_25pct = sum(share > 0.25), .by = region) |>
  arrange(desc(median_share))
analysis <- community[is_analysis & !community$species %in% deep_lineage, ]
stopifnot(!any(analysis$species %in% deep_lineage)) # deep-lineage taxa absent
stopifnot(nrow(analysis) == sum(is_analysis) - sum(community$species[is_analysis] %in% deep_lineage))

# every plot-year with any record counts as surveyed, whatever the record
surveyed <- distinct(community, region, origin_site, destination_site, treatment, plot_id, plot_year_id,
                     year = relative_year)

# ---- Experiments and their three groups of plots ----------------------------
included <- filter(experiments, included)
stopifnot(nrow(included) == 43) # 44 experiments, IN_Kashmir excluded for insufficient duration (P1)
groups <- bind_rows(
  transmute(included, experiment_id, region, origin_site, destination_site, treatment = "Warm", group = "warmed"),
  transmute(included, experiment_id, region, origin_site, destination_site = origin_site,
            treatment = "LocalControl", group = "origin_control"),
  transmute(included, experiment_id, region, origin_site = destination_site, destination_site,
            treatment = "LocalControl", group = "destination_control"))
plots <- inner_join(groups, surveyed, by = c("region", "origin_site", "destination_site", "treatment"),
                    relationship = "many-to-many") |>
  select(experiment_id, region, group, plot_id, plot_year_id, year)
records <- inner_join(plots, select(analysis, plot_year_id, species, relative_cover), by = "plot_year_id",
                      relationship = "many-to-many")

# ---- Both-groups rule: a year enters a pool's table only if both of its groups were surveyed
years <- distinct(plots, experiment_id, group, year)
pool_groups <- tibble::tibble(species_pool = c("resident", "resident", "destination_pool", "destination_pool"),
                              group = c("origin_control", "warmed", "destination_control", "warmed"))
pool_years <- inner_join(pool_groups, years, by = "group", relationship = "many-to-many") |>
  summarise(n_groups = n_distinct(group), groups = paste(sort(group), collapse = " + "),
            .by = c(experiment_id, species_pool, year))
year_fails <- filter(pool_years, n_groups < 2) |> arrange(experiment_id, species_pool, year)
pool_years <- filter(pool_years, n_groups == 2) |> select(experiment_id, species_pool, year)

# ---- Pools, from plot-year records over all survey years --------------------
n_rec <- records |>
  summarise(n = n_distinct(plot_year_id), .by = c(experiment_id, species, group)) |>
  tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0)
for (g in c("origin_control", "destination_control", "warmed")) if (!g %in% names(n_rec)) n_rec[[g]] <- 0L
membership <- n_rec |>
  mutate(species_pool = case_when(
    origin_control >= min_obs ~ "resident",
    destination_control >= min_obs ~ "destination_pool",
    origin_control == 0 & destination_control == 0 ~ "warmed only",
    TRUE ~ "below the filter in both control groups"))
pools <- filter(membership, species_pool %in% c("resident", "destination_pool"))
stopifnot(!anyDuplicated(pools[c("experiment_id", "species")])) # no species in both pools

# ---- Pool-size rule ------------------------------------------------------
pool_sizes <- pools |>
  count(experiment_id, species_pool) |>
  tidyr::pivot_wider(names_from = species_pool, values_from = n, values_fill = 0) |>
  right_join(select(included, experiment_id, region), by = "experiment_id") |>
  mutate(across(c(resident, destination_pool), ~ coalesce(.x, 0L)), k = k, rule = "resident pool >= k",
         kept = resident >= k) |>
  select(experiment_id, region, resident, destination_pool, rule, k, kept) |>
  arrange(resident)
kept <- filter(pool_sizes, kept)
# US_Arizona has 8 residents on the SQLite cover and stays (decided 28 Oct;
# PLAN.md expected about 4 from the legacy relative cover, and 42 experiments)
stopifnot(nrow(kept) == 43)                    # every included experiment passes k = 5
stopifnot(n_distinct(kept$region) == 21)       # 21 gradients, gradient_id = region

# ---- Zero-filled long table ------------------------------------------------
pool_plots <- pools |>
  filter(experiment_id %in% kept$experiment_id) |>
  select(experiment_id, species, species_pool) |>
  inner_join(pool_groups, by = "species_pool", relationship = "many-to-many") |>
  inner_join(plots, by = c("experiment_id", "group"), relationship = "many-to-many") |>
  semi_join(pool_years, by = c("experiment_id", "species_pool", "year"))
cover_long <- pool_plots |>
  left_join(select(analysis, plot_year_id, species, observed = relative_cover), by = c("plot_year_id", "species")) |>
  transmute(experiment_id, gradient_id = region, plot_id, year = factor(year, levels = sort(unique(year))),
            treatment = group, species, species_pool, relative_cover = coalesce(observed, 0))
stopifnot(!anyDuplicated(cover_long[c("experiment_id", "plot_id", "year", "species", "treatment")]))
stopifnot(n_distinct(cover_long$gradient_id) == 21, n_distinct(cover_long$experiment_id) == 43)

# spot check, counted independently: rows = pool species x surveyed plot-years
# of the pool's groups in its years; zeros = rows - observed records
spot <- kept$experiment_id[1]
spot_expected <- pool_groups |>
  inner_join(filter(plots, experiment_id == spot), by = "group", relationship = "many-to-many") |>
  semi_join(filter(pool_years, experiment_id == spot), by = c("experiment_id", "species_pool", "year")) |>
  count(species_pool, name = "plot_years") |>
  left_join(count(filter(pools, experiment_id == spot), species_pool, name = "species"), by = "species_pool") |>
  summarise(rows = sum(plot_years * species))
spot_check <- cover_long |>
  filter(experiment_id == spot) |>
  summarise(rows = n(), zero_rows = sum(relative_cover == 0), observed = sum(relative_cover > 0))
stopifnot(spot_check$rows == spot_expected$rows) # zero-filling covers every pool species x plot-year

# ---- A0: first available survey of each control group ----------------------
first_year <- plots |>
  filter(group %in% c("origin_control", "destination_control"), experiment_id %in% kept$experiment_id) |>
  summarise(first_year = min(year), .by = c(experiment_id, group))
a0_part <- function(g) {
  base <- plots |>
    filter(group == g) |>
    inner_join(filter(first_year, group == g), by = c("experiment_id", "group")) |>
    filter(year == first_year)
  pools |>
    filter(experiment_id %in% kept$experiment_id) |>
    select(experiment_id, species) |>
    inner_join(base, by = "experiment_id", relationship = "many-to-many") |>
    left_join(select(analysis, plot_year_id, species, rc = relative_cover), by = c("plot_year_id", "species")) |>
    summarise(A0 = log(100 * mean(coalesce(rc, 0)) + 1), .by = c(experiment_id, species))
}
a0 <- full_join(rename(a0_part("origin_control"), A0_origin = A0), rename(a0_part("destination_control"), A0_dest = A0),
                by = c("experiment_id", "species"))
stopifnot(nrow(a0) == nrow(distinct(cover_long, experiment_id, species)), !anyNA(a0$A0_origin), !anyNA(a0$A0_dest))
a0_not_year1 <- first_year |>
  filter(first_year != 1) |>
  tidyr::pivot_wider(names_from = group, values_from = first_year) |>
  arrange(experiment_id)

# ---- Flow counts (Fig S1) --------------------------------------------------
combos <- function(d) tibble::tibble(n_experiments = n_distinct(d$experiment_id), n_species = n_distinct(d$species),
                                     n_combinations = nrow(distinct(d, experiment_id, species)))
all_groups <- bind_rows(
  transmute(experiments, experiment_id, region, origin_site, destination_site, treatment = "Warm", group = "warmed"),
  transmute(experiments, experiment_id, region, origin_site, destination_site = origin_site, treatment = "LocalControl", group = "origin_control"),
  transmute(experiments, experiment_id, region, origin_site = destination_site, destination_site, treatment = "LocalControl", group = "destination_control"))
input_records <- inner_join(all_groups, filter(community, !is.na(species), is.na(drop_reason)),
                            by = c("region", "origin_site", "destination_site", "treatment"), relationship = "many-to-many")
flow_counts <- bind_rows(
  mutate(combos(input_records), rule = "input (P1 community, identified species)", branch = "all"),
  mutate(combos(filter(input_records, experiment_id %in% included$experiment_id)), rule = "experiment exclusions (P1)", branch = "all"),
  mutate(combos(records), rule = "deep-lineage taxa removed", branch = "all"),
  mutate(combos(pools), rule = "observation filter (>= 2 plot-year records in a control group)", branch = "all"),
  mutate(combos(filter(pools, species_pool == "resident")), rule = "pool assignment", branch = "resident"),
  mutate(combos(filter(pools, species_pool == "destination_pool")), rule = "pool assignment", branch = "destination_pool"),
  mutate(combos(distinct(cover_long, experiment_id, species)), rule = "pool-size rule (resident pool >= 5)", branch = "all"),
  mutate(combos(distinct(filter(cover_long, species_pool == "resident"), experiment_id, species)),
         rule = "pool-size rule (resident pool >= 5)", branch = "resident"),
  mutate(combos(distinct(filter(cover_long, species_pool == "destination_pool"), experiment_id, species)),
         rule = "pool-size rule (resident pool >= 5)", branch = "destination_pool")) |>
  mutate(order = row_number()) |>
  select(order, rule, branch, n_experiments, n_species, n_combinations)
not_pooled <- count(filter(membership, !species_pool %in% c("resident", "destination_pool")), reason = species_pool,
                    name = "species_x_experiment")

# ---- Write -------------------------------------------------------------------
saveRDS(cover_long, file.path(out_dir, "cover_long.rds"))
saveRDS(a0, file.path(out_dir, "a0.rds"))
utils::write.csv(select(filter(included, experiment_id %in% kept$experiment_id), experiment_id, region, origin_site,
                        destination_site, first_year, last_year), file.path(out_dir, "experiments_kept.csv"), row.names = FALSE)
utils::write.csv(pool_sizes, file.path(out_dir, "pool_sizes.csv"), row.names = FALSE)
utils::write.csv(flow_counts, file.path(out_dir, "flow_counts.csv"), row.names = FALSE)
fig <- ggplot(pool_sizes, aes(resident)) +
  geom_histogram(binwidth = 2, boundary = 0, fill = "grey60", colour = "white") +
  geom_vline(xintercept = k - 0.5, linetype = "dashed") +
  annotate("text", x = k - 0.5, y = Inf, label = paste0(" k = ", k), hjust = 0, vjust = 1.5) +
  labs(x = "Resident pool size (species)", y = "Experiments") +
  theme_classic()
ggsave("figures/diagnostics/01_pool_sizes.png", fig, width = 5, height = 3.5, dpi = 150)

# ---- Report --------------------------------------------------------------
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
share_q <- stats::quantile(dropped_share$share, c(0, 0.5, 0.9, 1))
report <- c(
  "# Step 01: Relative cover table and pools", "",
  paste0("Relative cover was computed per plot-year over all ", fmt(nrow(community)), " vascular records (identified, ",
         "genus-level and unidentified, the six deep-lineage taxa included), so it sums to 1 in every one of ",
         fmt(nrow(sums)), " plot-years. The analysis keeps identified species without a `drop_reason` and removes ",
         "the deep-lineage taxa afterwards. Each of the ", nrow(included), " included experiments gets its warmed ",
         "turfs, origin controls and destination controls; a species is resident with at least ", min_obs,
         " plot-year records in origin controls, in the destination pool with at least ", min_obs,
         " in destination controls and fewer than ", min_obs, " in origin controls, whether or not it was ever found in ",
         "warmed turfs. Experiments with fewer than k = ", k, " resident species leave, all observed cover is kept, and ",
         "every pool species gets a row, zero where unrecorded, in every surveyed plot-year of its pool's two groups."),
  "", "## Counts at every rule (`flow_counts.csv`, Fig S1)", "", md_table(select(flow_counts, -order)),
  "", "Species x experiment combinations in neither pool:", "", md_table(not_pooled),
  "", "## Relative cover", "",
  sprintf("Sum of relative cover per plot-year: min %.12f, max %.12f. Share of plot-year cover in rows not analysed (dropped or deep-lineage): min %.3f, median %.3f, 90%% %.3f, max %.3f.",
          min(sums$total), max(sums$total), share_q[1], share_q[2], share_q[3], share_q[4]),
  "", "Share of plot-year cover in rows not analysed, by region:", "", md_table(dropped_by_region),
  "", "## Pool sizes and the rule on k (`pool_sizes.csv`)", "",
  sprintf("Rule: an experiment enters only if its resident pool has at least k = %d species. %d of %d experiments kept%s. Smallest resident pools: %s.",
          k, nrow(kept), nrow(pool_sizes),
          if (any(!pool_sizes$kept)) paste0("; left: ", paste(pool_sizes$experiment_id[!pool_sizes$kept], collapse = ", ")) else "",
          paste(sprintf("%s %d", head(pool_sizes$experiment_id, 3), head(pool_sizes$resident, 3)), collapse = "; ")),
  paste("US_Arizona stays: on the SQLite's absolute cover it has", pool_sizes$resident[pool_sizes$region == "US_Arizona"],
        "resident species, where PLAN.md expected about four from the legacy relative cover (decided 28 Oct: kept)."),
  sprintf("Gradients (`gradient_id` = region): %d.", n_distinct(kept$region)),
  "", "![Resident pool size per experiment](../figures/diagnostics/01_pool_sizes.png)", "",
  md_table(head(pool_sizes, 10)),
  "", sprintf("## Survey years failing the both-groups rule (%d experiment-years)", nrow(year_fails)), "",
  "A year enters a pool's table only if both of that pool's groups were surveyed in it.", "",
  md_table(year_fails),
  "", "## A0 from the first available control survey", "",
  "A0 = log(100 x mean relative cover over the control plots of the group's first surveyed year + 1). Experiments whose first control survey is not year 1:", "",
  md_table(a0_not_year1),
  "", "## Table and spot check", "",
  sprintf("`cover_long`: %s rows, %d experiments, %d gradients, %s species x experiment combinations (%s resident, %s destination pool); %s rows observed, %s zeros.",
          fmt(nrow(cover_long)), n_distinct(cover_long$experiment_id), n_distinct(cover_long$gradient_id),
          fmt(nrow(a0)), fmt(sum(pools$species_pool == "resident" & pools$experiment_id %in% kept$experiment_id)),
          fmt(sum(pools$species_pool == "destination_pool" & pools$experiment_id %in% kept$experiment_id)),
          fmt(sum(cover_long$relative_cover > 0)), fmt(sum(cover_long$relative_cover == 0))),
  sprintf("Spot check, %s: pool species x surveyed plot-years of their groups = %s, rows written = %s (%s zeros + %s observed records).",
          spot, fmt(spot_expected$rows), fmt(spot_check$rows), fmt(spot_check$zero_rows), fmt(spot_check$observed)),
  "", "## Warnings captured", "", if (length(warns)) paste("-", unique(warns)) else "None.",
  "", "## Deviations from PLAN.md", "",
  paste("Relative cover includes the six deep-lineage taxa in its denominator and removes them afterwards (decided 28 Oct;",
        "PLAN.md Data rules updated). The observation filter decides pool membership only; observed cover of pool species",
        "is kept whatever its count (decided 28 Oct; legacy removed it per treatment).")
)
writeLines(report, "reports/01_relative_cover.md")
message("Step 01 written: ", nrow(kept), " experiments, ", fmt(nrow(cover_long)), " rows")
