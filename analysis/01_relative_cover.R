# Step 01: relative cover table and species pools.
#
# PLAN.md: Step A "Data rules"; The colonizer story (destination pool);
#   Appendix Fig S1 (flow counts).
# Reads:  data/derived/prep/community.rds, experiments.rds, species_taxonomy.rds
# Writes: data/derived/01_relative_cover/cover_long.rds, a0.rds,
#         experiments_kept.csv, pool_sizes.csv, flow_counts.csv;
#         figures/diagnostics/01_pool_rule.png, 01_pool_sizes.png;
#         reports/01_relative_cover.md
# Checks: relative cover sums to 1 per plot-year over all vascular records;
#         the six deep-lineage taxa absent afterwards; 43 included experiments;
#         no species in both pools; pool rule 1,008 resident effects (349
#         species) and 407 destination-pool effects; the pool-size rule removes
#         US_Arizona::MC::PP only; then 42 experiments, 20 gradients, 1,004
#         resident effects (346 species), 402 destination-pool effects (236
#         species); no duplicate experiment x plot x year x species x treatment;
#         every pool species has both A0 values.
#
# Pool rule (decided 2 Nov, Bills; PLAN.md Data rules): a species is resident
# if, in the origin controls, it was recorded in >= 2 plots in at least half of
# the survey years (A) or in at least half of the plots in >= 2 survey years
# (B), halves rounded up, plots counted among those surveyed that year, years
# among those the group was surveyed; destination pool = the same rule in the
# destination controls and no record in the origin controls. It replaces the
# 28 Oct filter (>= 2 plot-year records). Pool-size rule (reinstated 2 Nov):
# at least 5 residents per experiment. Decided 28 Oct and kept: the pool rule
# decides membership only and all observed cover is kept; A0 from the first
# available survey of each control group, log of the plot mean in % + 1;
# gradient_id = region; relative cover over all vascular records, deep-lineage
# taxa included, and those taxa removed afterwards; the destination pool
# whether or not ever found in warmed plots.

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
k_min <- 5L # residents for an experiment to enter (pool-size rule, reinstated 2 Nov)

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

# ---- Pools: the 2 Nov pool rule, per control group over its survey years ----
ctrl <- c("origin_control", "destination_control")
plots_year <- plots |> filter(group %in% ctrl) |>
  summarise(n_plots = n_distinct(plot_id), .by = c(experiment_id, group, year)) # plots surveyed that year
years_group <- summarise(plots_year, n_years = n_distinct(year), .by = c(experiment_id, group))
rule <- records |> filter(group %in% ctrl) |>
  summarise(plots_with = n_distinct(plot_id), .by = c(experiment_id, group, species, year)) |>
  left_join(plots_year, by = c("experiment_id", "group", "year")) |>
  summarise(years_2plots = sum(plots_with >= 2), years_half = sum(plots_with >= ceiling(n_plots / 2)),
            .by = c(experiment_id, group, species)) |>
  left_join(years_group, by = c("experiment_id", "group")) |>
  mutate(A = years_2plots >= ceiling(n_years / 2), B = years_half >= 2, pass = A | B) # halves rounded up
n_rec <- records |>
  summarise(n = n_distinct(plot_year_id), .by = c(experiment_id, species, group)) |>
  tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0)
for (g in c(ctrl, "warmed")) if (!g %in% names(n_rec)) n_rec[[g]] <- 0L
passes <- function(g) select(filter(rule, group == g), experiment_id, species, A, B, pass) |>
  rename_with(~ paste0(.x, "_", g), c(A, B, pass))
membership <- n_rec |>
  left_join(passes("origin_control"), by = c("experiment_id", "species")) |>
  left_join(passes("destination_control"), by = c("experiment_id", "species")) |>
  mutate(across(starts_with(c("A_", "B_", "pass_")), ~ coalesce(.x, FALSE)),
         species_pool = case_when(
           pass_origin_control ~ "resident",
           pass_destination_control & origin_control == 0 ~ "destination_pool", # no record in origin controls
           origin_control == 0 & destination_control == 0 ~ "warmed only",
           pass_destination_control ~ "destination rule met, but recorded in origin controls",
           TRUE ~ "below the pool rule in both control groups"))
pools <- filter(membership, species_pool %in% c("resident", "destination_pool"))
stopifnot(!anyDuplicated(pools[c("experiment_id", "species")])) # no species in both pools
stopifnot(sum(pools$species_pool == "resident") == 1008,        # pool rule (2 Nov), before the pool-size rule
          n_distinct(pools$species[pools$species_pool == "resident"]) == 349,
          sum(pools$species_pool == "destination_pool") == 407)

# ---- Pool-size rule (reinstated 2 Nov) ---------------------------------------
pool_sizes <- pools |>
  count(experiment_id, species_pool) |>
  tidyr::pivot_wider(names_from = species_pool, values_from = n, values_fill = 0) |>
  right_join(select(included, experiment_id, region), by = "experiment_id") |>
  mutate(across(c(resident, destination_pool), ~ coalesce(.x, 0L)), kept = resident >= k_min) |>
  select(experiment_id, region, resident, destination_pool, kept) |>
  arrange(resident)
stopifnot(nrow(pool_sizes) == 43) # every included experiment, kept or not
dropped_k <- filter(pool_sizes, !kept)
stopifnot(identical(dropped_k$experiment_id, "US_Arizona::MC::PP")) # fewer than 5 residents: US_Arizona only (4)
pools_before_k <- pools
pools <- filter(pools, !experiment_id %in% dropped_k$experiment_id)
kept_experiments <- filter(included, !experiment_id %in% dropped_k$experiment_id)
stopifnot(n_distinct(pools$experiment_id) == 42, n_distinct(sub("::.*", "", pools$experiment_id)) == 20, # gradient_id = region
          sum(pools$species_pool == "resident") == 1004, n_distinct(pools$species[pools$species_pool == "resident"]) == 346,
          sum(pools$species_pool == "destination_pool") == 402,
          n_distinct(pools$species[pools$species_pool == "destination_pool"]) == 236)

# ---- Zero-filled long table ------------------------------------------------
pool_plots <- pools |>
  select(experiment_id, species, species_pool) |>
  inner_join(pool_groups, by = "species_pool", relationship = "many-to-many") |>
  inner_join(plots, by = c("experiment_id", "group"), relationship = "many-to-many") |>
  semi_join(pool_years, by = c("experiment_id", "species_pool", "year"))
cover_long <- pool_plots |>
  left_join(select(analysis, plot_year_id, species, observed = relative_cover), by = c("plot_year_id", "species")) |>
  transmute(experiment_id, gradient_id = region, plot_id, year = factor(year, levels = sort(unique(year))),
            treatment = group, species, species_pool, relative_cover = coalesce(observed, 0))
stopifnot(!anyDuplicated(cover_long[c("experiment_id", "plot_id", "year", "species", "treatment")]))
stopifnot(n_distinct(cover_long$gradient_id) == 20, n_distinct(cover_long$experiment_id) == 42,
          !any(grepl("^US_Arizona", cover_long$experiment_id))) # pool-size rule

# spot check, counted independently: rows = pool species x surveyed plot-years
# of the pool's groups in its years; zeros = rows - observed records
spot <- filter(pool_sizes, kept)$experiment_id[1] # smallest kept resident pool
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
  filter(group %in% c("origin_control", "destination_control")) |>
  summarise(first_year = min(year), .by = c(experiment_id, group))
a0_part <- function(g) {
  base <- plots |>
    filter(group == g) |>
    inner_join(filter(first_year, group == g), by = c("experiment_id", "group")) |>
    filter(year == first_year)
  pools |>
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
  mutate(combos(pools_before_k), rule = "pool rule (2 Nov: >= 2 plots in half the years, or half the plots in >= 2 years)", branch = "all"),
  mutate(combos(pools), rule = "pool-size rule (at least 5 residents per experiment)", branch = "all"),
  mutate(combos(filter(pools, species_pool == "resident")), rule = "pool assignment", branch = "resident"),
  mutate(combos(filter(pools, species_pool == "destination_pool")), rule = "pool assignment", branch = "destination_pool"),
) |>
  mutate(order = row_number()) |>
  select(order, rule, branch, n_experiments, n_species, n_combinations)
not_pooled <- count(filter(membership, !species_pool %in% c("resident", "destination_pool")), reason = species_pool,
                    name = "species_x_experiment")

# ---- Write -------------------------------------------------------------------
saveRDS(cover_long, file.path(out_dir, "cover_long.rds"))
saveRDS(a0, file.path(out_dir, "a0.rds"))
utils::write.csv(select(kept_experiments, experiment_id, region, origin_site, destination_site, first_year, last_year),
                 file.path(out_dir, "experiments_kept.csv"), row.names = FALSE)
utils::write.csv(pool_sizes, file.path(out_dir, "pool_sizes.csv"), row.names = FALSE)
utils::write.csv(flow_counts, file.path(out_dir, "flow_counts.csv"), row.names = FALSE)
fig <- ggplot(pool_sizes, aes(resident)) +
  geom_histogram(binwidth = 2, boundary = 0, fill = "grey60", colour = "white") +
  geom_vline(xintercept = k_min, colour = "#B2182B", linetype = "dashed") +
  labs(x = "Resident pool size (species; dashed: pool-size rule)", y = "Experiments") +
  theme_classic(base_size = 14)
ggsave("figures/diagnostics/01_pool_sizes.png", fig, width = 6, height = 4, dpi = 150)

# ---- Pool-rule diagnostic: what the rule keeps and leaves behind (2 Nov) -----
# every species x experiment recorded in a control group, by its plot-year
# records there (the old cut-off axis), split by what the pool rule did with it
removed_k <- dropped_k$experiment_id
status <- membership |>
  mutate(resident_status = case_when(origin_control == 0 ~ NA_character_,
                                     species_pool == "resident" & experiment_id %in% removed_k ~ "left behind: experiment under 5 residents",
                                     species_pool == "resident" ~ "kept",
                                     TRUE ~ "left behind: fails the pool rule"),
         colonizer_status = case_when(destination_control == 0 | species_pool == "resident" ~ NA_character_,
                                      species_pool == "destination_pool" & experiment_id %in% removed_k ~ "left behind: experiment under 5 residents",
                                      species_pool == "destination_pool" ~ "kept",
                                      pass_destination_control ~ "left behind: recorded in origin controls",
                                      TRUE ~ "left behind: fails the pool rule"))
status_long <- bind_rows(
  transmute(filter(status, !is.na(resident_status)), pool = "Residents (origin controls)", plot_years = origin_control, status = resident_status),
  transmute(filter(status, !is.na(colonizer_status)), pool = "Colonizers (destination controls)", plot_years = destination_control,
            status = colonizer_status)) |>
  mutate(pool = factor(pool, levels = c("Residents (origin controls)", "Colonizers (destination controls)")),
         status = factor(status, levels = c("kept", "left behind: fails the pool rule", "left behind: recorded in origin controls",
                                            "left behind: experiment under 5 residents")))
stopifnot(sum(status_long$status == "kept") == nrow(pools)) # every kept combination appears once
status_cols <- c("kept" = "#1B9E77", "left behind: fails the pool rule" = "grey70",
                 "left behind: recorded in origin controls" = "#E6A15C", "left behind: experiment under 5 residents" = "#B2182B")
big <- theme_classic(base_size = 21) + theme(plot.title = element_text(face = "bold", size = 29), axis.title = element_text(face = "bold"),
  strip.text = element_text(face = "bold", size = 19), strip.background = element_blank(), legend.position = "top",
  legend.text = element_text(size = 17), plot.background = element_rect(fill = "white", colour = NA))
pa <- ggplot(status_long, aes(pmin(plot_years, 40), fill = status)) +
  geom_histogram(binwidth = 1, boundary = 0.5, colour = "white", linewidth = 0.2) +
  facet_wrap(~ pool, ncol = 1, scales = "free_y") + scale_fill_manual(values = status_cols, name = NULL, drop = FALSE) +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 20, 30, 40)) + guides(fill = guide_legend(ncol = 2)) +
  labs(title = "a)", subtitle = "Each species x experiment recorded in the control group, by its plot-year records there",
       x = "Plot-year records in the control group (40 or more in the last bar)", y = "Species x experiment") + big
pb_dat <- count(status_long, pool, status)
pb <- ggplot(pb_dat, aes(status, n, fill = status)) + geom_col(width = 0.7) +
  geom_text(aes(label = format(n, big.mark = ",")), vjust = -0.3, size = 6.6, fontface = "bold") +
  facet_wrap(~ pool, ncol = 1, scales = "free") + scale_fill_manual(values = status_cols, guide = "none") +
  scale_x_discrete(labels = function(x) gsub("(.{1,16})(\\s|$)", "\\1\n", sub("left behind: ", "left behind: ", x))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "b)", subtitle = paste0("Kept: ", n_distinct(pools$species[pools$species_pool == "resident"]), " resident and ",
                                       n_distinct(pools$species[pools$species_pool == "destination_pool"]), " colonizer species"),
       x = NULL, y = "Species x experiment") + big + theme(axis.text.x = element_text(size = 15))
pc_dat <- status |> filter(!is.na(resident_status)) |> count(experiment_id, status = resident_status) |>
  mutate(status = factor(status, levels = levels(status_long$status)),
         experiment = factor(sub("::", " ", experiment_id), levels = rev(sort(unique(sub("::", " ", experiment_id))))))
pc <- ggplot(pc_dat, aes(n, experiment, fill = status)) + geom_col(width = 0.75) +
  geom_vline(xintercept = k_min, colour = "#B2182B", linetype = "dashed") +
  scale_fill_manual(values = status_cols, guide = "none") +
  labs(title = "c)", subtitle = paste0("Species recorded in each experiment's origin controls\n(dashed: ", k_min, " residents)"),
       x = "Species", y = NULL) + big + theme(axis.text.y = element_text(size = 14))
kept_grad <- pools |> mutate(gradient_id = sub("::.*", "", experiment_id),
                             pool = ifelse(species_pool == "resident", "Residents", "Colonizers")) |>
  summarise(gradients = n_distinct(gradient_id), .by = c(pool, species)) |> mutate(pool = factor(pool, levels = c("Residents", "Colonizers")))
pd <- ggplot(kept_grad, aes(gradients, fill = pool, alpha = gradients >= 2)) +
  geom_bar(position = position_dodge(preserve = "single"), width = 0.8) +
  geom_text(stat = "count", aes(label = after_stat(count)), position = position_dodge(width = 0.8), vjust = -0.3, size = 5.6,
            fontface = "bold", alpha = 1) +
  geom_vline(xintercept = 1.5, colour = "#D95F2B", linetype = "dashed", linewidth = 0.9) +
  scale_fill_manual(values = c(Residents = "#1B9E77", Colonizers = "#7442A0"), name = NULL) +
  scale_alpha_manual(values = c(`FALSE` = 0.4, `TRUE` = 1), guide = "none") + scale_x_continuous(breaks = 1:21) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "d)", subtitle = "Kept species by the number of gradients they occur in\n(dashed: at least 2, as in steps 04 and 05)",
       x = "Gradients", y = "Species") + big
ggsave("figures/diagnostics/01_pool_rule.png",
       # each panel drawn on its own, so its axis titles sit close instead of aligning to c)'s long labels (Bills, 2 Nov)
       patchwork::wrap_plots(lapply(list(pa, pb, pc, pd), function(p) patchwork::wrap_elements(full = p)),
                             design = "AB\nCD", widths = c(1.3, 1), heights = c(1, 1.3)),
       width = 24, height = 24, dpi = 110, bg = "white")
criteria <- pools |> transmute(pool = species_pool,
  A = ifelse(pool == "resident", A_origin_control, A_destination_control), B = ifelse(pool == "resident", B_origin_control, B_destination_control)) |>
  summarise(`A only` = sum(A & !B), `B only` = sum(!A & B), `A and B` = sum(A & B), .by = pool)

# ---- Report --------------------------------------------------------------
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
share_q <- stats::quantile(dropped_share$share, c(0, 0.5, 0.9, 1))
report <- c(
  "# Step 01: Relative cover table and pools", "",
  paste0("Relative cover was computed per plot-year over all ", fmt(nrow(community)), " vascular records (identified, ",
         "genus-level and unidentified, the six deep-lineage taxa included), so it sums to 1 in every one of ",
         fmt(nrow(sums)), " plot-years. The analysis keeps identified species without a `drop_reason` and removes ",
         "the deep-lineage taxa afterwards. Each of the ", nrow(included), " included experiments gets its warmed ",
         "turfs, origin controls and destination controls. Pool rule (2 Nov): a species is resident if, in the origin ",
         "controls, it was recorded in at least two plots in at least half of the survey years, or in at least half of the ",
         "plots in at least two survey years (halves rounded up; plots surveyed that year; years the group was surveyed); ",
         "it is in the destination pool if it meets the same rule in the destination controls and has no record in the ",
         "origin controls, whether or not it was ever found in warmed turfs. Pool-size rule (reinstated 2 Nov): an ",
         "experiment needs at least ", k_min, " residents. All observed cover of pool species is kept, and every pool ",
         "species gets a row, zero where unrecorded, in every surveyed plot-year of its pool's two groups."),
  "", "## Counts at every rule (`flow_counts.csv`, Fig S1)", "", md_table(select(flow_counts, -order)),
  "", "Species x experiment combinations in neither pool (before the pool-size rule):", "", md_table(not_pooled),
  "", "## What the pool rule keeps and leaves behind", "",
  paste0("Every species x experiment recorded in a control group, by its plot-year records there and by what the rule ",
         "did with it (a); totals (b); per experiment, the species recorded in the origin controls (c); kept species by the ",
         "number of gradients they occur in, with the at-least-2 cut-off that steps 04 and 05 apply (d). The rule is not a ",
         "cut-off on records: combinations with many records but in few plots or few years are left behind, and some with ",
         "few records are kept where the plots are few."), "",
  "![Kept and left behind](../figures/diagnostics/01_pool_rule.png)", "",
  md_table(tidyr::pivot_wider(count(status_long, pool, status), names_from = status, values_from = n, values_fill = 0)),
  "", "Kept effects by criterion (A: at least 2 plots in at least half of the survey years; B: at least half of the plots in at least 2 survey years):", "",
  md_table(mutate(criteria, pool = ifelse(pool == "resident", "residents", "destination pool"))),
  "", "## Relative cover", "",
  sprintf("Sum of relative cover per plot-year: min %.12f, max %.12f. Share of plot-year cover in rows not analysed (dropped or deep-lineage): min %.3f, median %.3f, 90%% %.3f, max %.3f.",
          min(sums$total), max(sums$total), share_q[1], share_q[2], share_q[3], share_q[4]),
  "", "Share of plot-year cover in rows not analysed, by region:", "", md_table(dropped_by_region),
  "", "## Pool sizes (`pool_sizes.csv`)", "",
  sprintf(paste("Pool-size rule: %d of %d experiments have at least %d residents. Left out: %s, with %d residents: %d species",
                "were ever recorded in its origin controls, a median of %.1f per plot-year. Smallest kept pools: %s."),
          sum(pool_sizes$kept), nrow(pool_sizes), k_min, paste(dropped_k$experiment_id, collapse = ", "), dropped_k$resident[1],
          n_distinct(filter(records, experiment_id %in% removed_k, group == "origin_control")$species),
          stats::median(summarise(filter(records, experiment_id %in% removed_k, group == "origin_control"),
                                  n = n_distinct(species), .by = plot_year_id)$n),
          paste(sprintf("%s %d", head(filter(pool_sizes, kept)$experiment_id, 3), head(filter(pool_sizes, kept)$resident, 3)),
                collapse = "; ")),
  sprintf("Gradients (`gradient_id` = region) after the rule: %d.", n_distinct(filter(pool_sizes, kept)$region)),
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
          fmt(nrow(a0)), fmt(sum(pools$species_pool == "resident")), fmt(sum(pools$species_pool == "destination_pool")),
          fmt(sum(cover_long$relative_cover > 0)), fmt(sum(cover_long$relative_cover == 0))),
  sprintf("Spot check, %s: pool species x surveyed plot-years of their groups = %s, rows written = %s (%s zeros + %s observed records).",
          spot, fmt(spot_expected$rows), fmt(spot_check$rows), fmt(spot_check$zero_rows), fmt(spot_check$observed)),
  "", "## Warnings captured", "", if (length(warns)) paste("-", unique(warns)) else "None.",
  "", "## Deviations from PLAN.md", "",
  paste("Relative cover includes the six deep-lineage taxa in its denominator and removes them afterwards (decided 28 Oct;",
        "PLAN.md Data rules updated). The pool rule decides pool membership only; observed cover of pool species",
        "is kept whatever its count (decided 28 Oct; legacy removed it per treatment).")
)
writeLines(report, "reports/01_relative_cover.md")
message("Step 01 written: ", n_distinct(cover_long$experiment_id), " experiments, ", fmt(nrow(cover_long)), " rows")
