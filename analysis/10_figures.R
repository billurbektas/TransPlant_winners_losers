# Step 10: figures, appendix and report.
#
# PLAN.md: Figure map, Captions, Appendix.
# Reads:  data/derived/prep/ (experiments, sites, experiment_climate, site_bioclim_check,
#         traits, niche_optima, functional_groups, phylogeny_axes);
#         data/derived/01_relative_cover/ (flow_counts.csv, cover_long.rds, a0.rds,
#         pool_sizes.csv); data/derived/03b_strength_y3/strength_df.rds;
#         data/derived/07a_colonizers/strength_df_y3.rds; the species summaries of steps
#         04 and 07b; the tables of steps 05 and 07c; the reviewed PNGs of steps 03 to 07f.
# Writes: figures/main/fig1.png ... fig7.png; figures/appendix/figS1.png ... figS8b.png;
#         tables/tableS1.csv ... tableS5.csv; data/derived/10_figures/fig4_pie_cover.csv,
#         fig4_pie_counts.csv;
#         reports/10_figures.md
# Checks: 42 experiments; 1004 resident effects (703 losers, 301 winners); 262
#         colonizer effects with warmed cover above zero in 28 experiments; every pie's
#         slices sum to its species; every figure number exists exactly once.
#
# Decided 2 Nov (Bills): one script for the whole step, run without stopping. Seven main
# figures: 1 Sankey, 2 residents' classification, 3 colonizers' classification, 4
# climate space, 5 traits (pathways), 6 the residents' ladder; the results-dependent
# Fig 6 and the old Fig S8 are dropped. The climate-space pies (Fig 4 as cover, Fig S2 as counts) use
# Fig 5's classes per experiment: residents by the sign of their net strength from
# year 3 on, colonizers with net warmed cover above zero (07a, so the 14 experiments
# below the colonizer pool-size rule have no purple); one radius, pies sharing a
# destination climate fanned out around it. Reviewed figures are copied byte for byte.

suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
source("R/md_table.R")
warns <- character()
globalCallingHandlers(warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
})
out_dir <- "data/derived/10_figures"
for (d in c(out_dir, "figures/appendix", "tables")) dir.create(d, recursive = TRUE, showWarnings = FALSE)
cols <- c("Resident losers" = "#D95F2B", "Resident winners" = "#1B9E77", "Colonizers" = "#5E3C99") # Fig 5
big <- theme_classic(base_size = 30) +
  theme(plot.title = element_text(face = "bold", size = 40), axis.title = element_text(face = "bold"),
        plot.background = element_rect(fill = "white", colour = NA))

# ---- Shared inputs ------------------------------------------------------------------
res <- readRDS("data/derived/03b_strength_y3/strength_df.rds")
stopifnot(nrow(res) == 1004, n_distinct(res$experiment_id) == 42) # step 03b
col <- readRDS("data/derived/07a_colonizers/strength_df_y3.rds")
stopifnot(nrow(col) == 379, n_distinct(col$experiment_id) == 28) # step 07a, pool-size rule
sites <- readRDS("data/derived/prep/sites.rds") |> select(region, site_id, elevation, country)
experiments <- readRDS("data/derived/prep/experiments.rds") |> filter(experiment_id %in% res$experiment_id) |>
  left_join(rename(sites, origin_site = site_id, origin_elevation = elevation), by = c("region", "origin_site")) |>
  left_join(select(rename(sites, destination_site = site_id, destination_elevation = elevation), -country), by = c("region", "destination_site")) |>
  left_join(readRDS("data/derived/prep/experiment_climate.rds"), by = "experiment_id")
stopifnot(nrow(experiments) == 42, !anyNA(experiments$origin_elevation), !anyNA(experiments$destination_elevation),
          !anyNA(experiments$destination_P_JJA)) # every experiment placed

# ---- Fig 4: the experiments in climate space ---------------------------------------
pies <- bind_rows(count(res, experiment_id, group = ifelse(strength < 0, "Resident losers", "Resident winners")),
                  count(filter(col, strength > 0), experiment_id, group = "Colonizers")) # net warmed cover above zero
stopifnot(sum(res$strength < 0) == 703, sum(res$strength > 0) == 301, sum(col$strength > 0) == 262) # Fig 5's classes per effect
pies <- tidyr::pivot_wider(pies, names_from = group, values_from = n, values_fill = 0)
pies <- left_join(select(experiments, experiment_id, origin_elevation, destination_elevation, destination_T_JJA, destination_P_JJA, diff_T_JJA),
                  pies, by = "experiment_id") |> mutate(across(all_of(names(cols)), ~ coalesce(.x, 0L)))
stopifnot(sum(pies$`Resident losers` + pies$`Resident winners`) == 1004, sum(pies$Colonizers) == 262,
          sum(pies$Colonizers > 0) == n_distinct(col$experiment_id[col$strength > 0])) # one denominator: counts, not cover
write.csv(pies, file.path(out_dir, "fig4_pie_counts.csv"), row.names = FALSE)
# main Fig 4 shows the same classes as shares of the warmed turf's relative cover, years 3 and
# later (Bills, 2 Nov); the counts move to Fig S2. Relative cover is over all vascular cover of
# the plot-year (step 01), so what no class holds is a slice of its own: species outside the
# pools, unidentified and genus-level records
cv <- readRDS("data/derived/01_relative_cover/cover_long.rds") |> mutate(year = as.integer(as.character(year))) |>
  filter(experiment_id %in% res$experiment_id, treatment == "warmed", year >= 3)
cls <- bind_rows(transmute(res, experiment_id, species, group = ifelse(strength < 0, "Resident losers", "Resident winners")),
                 transmute(filter(col, strength > 0), experiment_id, species, group = "Colonizers"))
cover_shares <- cv |> inner_join(cls, by = c("experiment_id", "species")) |>
  summarise(cover = sum(relative_cover), .by = c(experiment_id, group)) |>
  left_join(summarise(cv, n_py = n_distinct(paste(plot_id, year)), .by = experiment_id), by = "experiment_id") |>
  mutate(cover = cover / n_py) |> select(-n_py) |> tidyr::pivot_wider(names_from = group, values_from = cover, values_fill = 0)
cover_shares <- left_join(select(experiments, experiment_id), cover_shares, by = "experiment_id") |>
  mutate(across(all_of(names(cols)), ~ coalesce(.x, 0)), `Not analysed` = 1 - `Resident losers` - `Resident winners` - Colonizers)
stopifnot(nrow(cover_shares) == 42, all(cover_shares$`Not analysed` > -1e-9)) # shares of one turf: the classes never exceed it
write.csv(cover_shares, file.path(out_dir, "fig4_pie_cover.csv"), row.names = FALSE)
cols_cover <- c(cols, "Not analysed" = "grey80")
# pies need equal units on both axes: log precipitation rescaled to the temperature range
k <- diff(range(pies$destination_T_JJA)) / diff(range(log(pies$destination_P_JJA))) * 1.1
r <- 0.21
pies <- mutate(pies, x0 = destination_T_JJA, y0 = log(destination_P_JJA) * k)
# experiments sharing a destination climate fan out around the shared point on short lines
# (Bills, 2 Nov); every other pie sits at its own climate
pies <- pies |> mutate(n_shared = n(), i = row_number(), .by = c(x0, y0)) |>
  mutate(d = ifelse(n_shared > 1, pmax(1.7 * r, 1.1 * r / sin(pi / n_shared)), 0),
         angle = pi / 2 + 2 * pi * (i - 1) / n_shared, x = x0 + d * cos(angle), y = y0 + d * sin(angle))
stopifnot(sum(pies$n_shared > 1) == 15, n_distinct(paste(pies$x0, pies$y0)) == 33) # 5 pairs and one five sharing a destination
# pies of different climates that touch are pushed apart; each moved pie keeps a dot at its
# true climate and a line to it, the same mark as the shared destinations (Bills, 2 Nov)
for (it in 1:500) for (i in seq_len(nrow(pies) - 1)) for (j in (i + 1):nrow(pies)) {
  dx <- pies$x[j] - pies$x[i]; dy <- pies$y[j] - pies$y[i]; d <- sqrt(dx^2 + dy^2)
  if (d < 2.1 * r) { m <- (2.1 * r - d) / 2; pies$x[i] <- pies$x[i] - m * dx / d; pies$y[i] <- pies$y[i] - m * dy / d
                     pies$x[j] <- pies$x[j] + m * dx / d; pies$y[j] <- pies$y[j] + m * dy / d }
}
pies <- mutate(pies, moved = sqrt((x - x0)^2 + (y - y0)^2) > 0.01)
centre_d <- as.matrix(dist(cbind(pies$x, pies$y)))
stopifnot(all(centre_d[upper.tri(centre_d)] >= 2 * r * 0.99)) # no two pies overlap
# each pie is named by its gradient, with a letter where the gradient has several experiments
# (highest origin first); the key beside the figure gives each one's elevations (Bills, 2 Nov)
pies <- pies |> left_join(distinct(res, experiment_id, gradient_id), by = "experiment_id") |>
  arrange(gradient_id, desc(origin_elevation), desc(destination_elevation)) |>
  mutate(letter = if (n() > 1) letters[row_number()] else "", .by = gradient_id) |>
  mutate(label = trimws(paste(sub("_", " ", gradient_id), letter)),
         elevations = paste0(round(origin_elevation), " → ", round(destination_elevation), " m"))
stopifnot(n_distinct(pies$label) == 42, n_distinct(pies$gradient_id) == 20) # one name per pie, 20 gradients
warm_scale <- scale_colour_gradient(low = "#C99A1C", high = "#7B0D1E", limits = range(pies$diff_T_JJA), name = "Warming\n(°C, JJA)")
p_breaks <- c(150, 200, 300, 400, 600, 800, 1000)
climate_map <- function(slices, fill_values) { # the same layout, slices as counts or as cover shares
  d <- left_join(select(pies, -any_of(names(cols_cover))), slices, by = "experiment_id")
  ggplot(d) +
  geom_segment(data = filter(d, moved), aes(x0, y0, xend = x, yend = y), colour = "grey25", linewidth = 0.9) +
  geom_point(data = distinct(filter(d, moved), x0, y0), aes(x0, y0), size = 4) + # the true climate of a moved pie
  scatterpie::geom_scatterpie(aes(x = x, y = y, r = r), data = mutate(d, r = r), cols = names(fill_values), colour = "white", linewidth = 0.3) +
  ggrepel::geom_text_repel(aes(x, y, label = label, colour = diff_T_JJA), size = 7.5, fontface = "bold", bg.colour = "white", bg.r = 0.15,
                           point.size = 24, point.padding = 0, box.padding = 0.25,
                           nudge_x = ifelse(d$n_shared > 1, 2.1 * r * cos(d$angle), 0), # fanned pies: label outward, off the shared point
                           nudge_y = ifelse(d$n_shared > 1, 2.1 * r * sin(d$angle), -1.9 * r), force = 5, force_pull = 3, # point.size = the pie's size
                           max.iter = 1e6, min.segment.length = 0.15, segment.colour = "grey45", segment.size = 0.6, max.overlaps = Inf, seed = 1) +
  scale_fill_manual(values = fill_values, breaks = names(fill_values), name = NULL) + warm_scale +
  scale_y_continuous(breaks = log(p_breaks) * k, labels = p_breaks) + coord_equal() +
  labs(x = "Destination June–August mean temperature (°C)", y = "Destination June–August\nprecipitation (mm, log scale)") +
  big + theme(legend.position = "bottom", legend.box = "vertical", legend.key.size = unit(1.4, "cm"), legend.key.width = unit(3, "cm"),
              legend.title = element_text(face = "bold"), panel.grid.major = element_line(colour = "grey90"))
}
key <- arrange(pies, desc(diff_T_JJA)) |> mutate(row = row_number(), name = sub("_", " ", gradient_id)) # strongest warming first (Bills, 2 Nov)
fig4key <- ggplot(key, aes(y = -row)) +
  geom_text(aes(x = 0, label = name, colour = diff_T_JJA), hjust = 0, size = 10, fontface = "bold") +
  geom_text(aes(x = 3.7, label = letter, colour = diff_T_JJA), hjust = 0.5, size = 10, fontface = "bold") +
  geom_text(aes(x = 4.0, label = elevations, colour = diff_T_JJA), hjust = 0, size = 10) +
  warm_scale + guides(colour = "none") + scale_x_continuous(limits = c(0, 6.6), expand = c(0, 0)) +
  labs(title = "Origin → destination elevation,\nstrongest warming first") + theme_void() +
  theme(plot.title = element_text(face = "bold", size = 34, margin = margin(b = 10)), plot.background = element_rect(fill = "white", colour = NA))
fig4 <- patchwork::wrap_plots(climate_map(cover_shares, cols_cover), fig4key, ncol = 2, widths = c(2.2, 0.7))
figS2 <- patchwork::wrap_plots(climate_map(select(pies, experiment_id, all_of(names(cols))), cols), fig4key, ncol = 2, widths = c(2.2, 0.7))

# ---- Fig 7: what the models find, residents beside colonizers --------------------------
# Climate + abundance, Height and LDMC models (Bills, 2 Nov). Both responses are relative cover
# and W, T, P and the traits share one scale (steps 02, 07c), so those terms sit on one axis;
# the abundance terms are standardised within each pool and are drawn apart
coefs <- bind_rows(mutate(read.csv("data/derived/05_ladder_fit/tableS5.csv"), pool = "Residents"),
                   mutate(read.csv("data/derived/07c_ladder_fit/coefficients.csv"), pool = "Colonizers")) |>
  filter(model %in% c("Climate + start", "Height", "LDMC"), term != "intrcpt") |>
  mutate(pool = factor(pool, levels = c("Residents", "Colonizers")), excludes_zero = ci_lower > 0 | ci_upper < 0)
climate_terms <- c(W = "Warming", T = "Temperature", P = "Precipitation", `W:T` = "Warming × temperature", `W:P` = "Warming × precipitation",
                   `T:P` = "Temperature × precipitation", `W:T:P` = "Warming × temperature × precipitation")
abund_terms <- c(A0_origin = "Abundance at origin", A0_dest = "Abundance at destination", `W:A0_origin` = "Warming × abundance at origin",
                 `W:A0_dest` = "Warming × abundance at destination")
trait_terms <- function(z, name) setNames(paste(name, c("", "× warming", "× temperature", "× precipitation", "× temperature × precipitation")),
                                         paste0(c("", "W:", "T:", "P:", ""), z, c("", "", "", "", ":TP")))
prop_terms <- c(trait_terms("height_z", "Height"), trait_terms("ldmc_z", "LDMC"))
f7_clim <- filter(coefs, model == "Climate + start", term %in% names(climate_terms)) |> mutate(name = factor(climate_terms[term], levels = rev(climate_terms)))
f7_abund <- filter(coefs, model == "Climate + start", term %in% names(abund_terms)) |> mutate(name = factor(abund_terms[term], levels = rev(abund_terms)))
f7_prop <- filter(coefs, model %in% c("Height", "LDMC"), term %in% names(prop_terms)) |> mutate(name = factor(trimws(prop_terms[term]), levels = rev(trimws(prop_terms))))
stopifnot(nrow(f7_clim) == 14, nrow(f7_abund) == 6, nrow(f7_prop) == 20) # 7 climate terms and 5 terms per trait in each pool; A0_origin residents only
pool_cols <- c(Residents = "#3D3D3D", Colonizers = "#5E3C99")
x_shared <- range(c(f7_clim$ci_lower, f7_clim$ci_upper, f7_prop$ci_lower, f7_prop$ci_upper)) # a, c, d on one axis
coef_panel <- function(d, tag, xlim = NULL) {
  ggplot(d, aes(estimate, name, colour = pool, group = pool)) + geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_linerange(aes(xmin = ci_lower, xmax = ci_upper), linewidth = 1.6, position = position_dodge(width = 0.6)) +
    geom_point(aes(fill = ifelse(excludes_zero, as.character(pool), "none")), shape = 21, size = 6, stroke = 1.6, position = position_dodge(width = 0.6)) +
    scale_colour_manual(values = pool_cols, name = NULL, breaks = names(pool_cols)) +
    scale_fill_manual(values = c(pool_cols, none = "white"), guide = "none") +
    (if (!is.null(xlim)) scale_x_continuous(limits = xlim)) +
    labs(title = tag, x = "Effect on relative cover per SD of the predictor", y = NULL) + big +
    theme(legend.position = "bottom", legend.text = element_text(size = 30), legend.key.size = unit(1.6, "cm"))
}
fig7 <- patchwork::wrap_plots(
  coef_panel(f7_clim, "a) Climate (climate + abundance model)", x_shared),
  coef_panel(f7_abund, "b) Abundance (same model; scaled within each pool)") + facet_wrap(~ pool, ncol = 1, scales = "free") +
    theme(strip.text = element_text(face = "bold", size = 28), strip.background = element_blank()),
  coef_panel(filter(f7_prop, model == "Height"), "c) Height model", x_shared),
  coef_panel(filter(f7_prop, model == "LDMC"), "d) LDMC model", x_shared),
  ncol = 2, guides = "collect") & theme(legend.position = "bottom", plot.title = element_text(size = 34))

# ---- Fig S1: data flow ------------------------------------------------------------------
flow <- read.csv("data/derived/01_relative_cover/flow_counts.csv")
stopifnot(nrow(flow) == 7) # step 01's filters
s04 <- read.csv("data/derived/04_classification/species_summary.csv")
s07b <- read.csv("data/derived/07b_classification/species_summary.csv")
h05 <- read.csv("data/derived/05_ladder_fit/heterogeneity.csv") |> filter(rowset == "all", model == "Null")
lad_col <- readRDS("data/derived/07c_ladder_fit/ladder_rows.rds")
box <- function(e, s, n) sprintf("%d experiments · %d species · %d effects", e, s, n)
excl <- filter(readRDS("data/derived/prep/experiments.rds"), !included)
stopifnot(nrow(excl) == flow$n_experiments[1] - flow$n_experiments[2]) # the experiments removed between the first two boxes
excl_rule <- paste(sprintf("%s (%s → %s) removed: %s, %d survey years (%d–%d)", sub("_", " ", excl$region), excl$origin_site, excl$destination_site,
                           excl$exclusion_reason, excl$year_range, excl$first_year, excl$last_year), collapse = "\n")
nodes <- tibble::tribble(
  ~id, ~x, ~y, ~title, ~counts,
  1, 0, 9, "Community data (P1)", box(flow$n_experiments[1], flow$n_species[1], flow$n_combinations[1]),
  2, 0, 8, "Experiments kept", box(flow$n_experiments[2], flow$n_species[2], flow$n_combinations[2]),
  3, 0, 7, "Deep lineages removed", box(flow$n_experiments[3], flow$n_species[3], flow$n_combinations[3]),
  4, 0, 6, "Pool rule", box(flow$n_experiments[4], flow$n_species[4], flow$n_combinations[4]),
  5, 0, 5, "Pool-size rule", box(flow$n_experiments[5], flow$n_species[5], flow$n_combinations[5]),
  6, -1, 4, "Residents: strength from year 3", box(n_distinct(res$experiment_id), n_distinct(res$species), nrow(res)),
  7, -1, 3, "Residents: Step C (Fig 2)", box(NA_integer_, nrow(s04), sum(s04$n_experiments)),
  8, -1, 2, "Residents: the ladder (Fig 6)", box(h05$n_experiments, h05$n_species, h05$n_effects),
  9, 1, 4, "Colonizers: warmed cover from year 3", box(n_distinct(col$experiment_id), n_distinct(col$species), nrow(col)),
  10, 1, 3, "Colonizers: Step C (Fig 3)", box(NA_integer_, nrow(s07b), sum(s07b$n_experiments)),
  11, 1, 2, "Colonizers: the ladder (Fig S7)", box(n_distinct(lad_col$experiment_id), n_distinct(lad_col$species), nrow(lad_col))) |>
  mutate(counts = sub("NA experiments · ", "", counts))
edges <- tibble::tribble(
  ~from, ~to, ~rule,
  1, 2, excl_rule,
  2, 3, "six deep-lineage taxa",
  3, 4, ">= 2 plots in half the years, or half the plots in >= 2 years",
  4, 5, ">= 5 residents per experiment",
  5, 6, "resident pool",
  6, 7, "species in >= 3 gradients",
  7, 8, "species with every property",
  5, 9, "destination pool; >= 5 colonizers per experiment",
  9, 10, "species in >= 3 gradients",
  10, 11, "species with every property") |>
  left_join(select(nodes, from = id, x, y), by = "from") |> left_join(select(nodes, to = id, xend = x, yend = y), by = "to")
figS1 <- ggplot() +
  geom_segment(data = edges, aes(x, y - 0.28, xend = xend, yend = yend + 0.28), arrow = arrow(length = unit(0.4, "cm")), linewidth = 0.9, colour = "grey40") +
  geom_text(data = edges, aes((x + xend) / 2 + ifelse(x == xend, 0.04, 0), (y + yend) / 2 + ifelse(x == xend, 0, 0.22), label = rule,
                              hjust = ifelse(x == xend, 0, 0.5)), size = 6.5, colour = "grey30", fontface = "italic") + # branch labels above their arrows
  geom_label(data = nodes, aes(x, y, label = paste0(title, "\n", counts),
                               fill = ifelse(x < 0, "res", ifelse(x > 0, "col", "all"))), size = 7, label.padding = unit(0.5, "lines"), lineheight = 1) +
  scale_fill_manual(values = c(all = "grey95", res = "#FBE3D6", col = "#E6DDF1"), guide = "none") +
  scale_x_continuous(limits = c(-1.9, 2.3)) + scale_y_continuous(limits = c(1.5, 9.5)) +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))

# ---- Fig S3 and Table S3: phylogenetic axes -------------------------------------------
phy <- readRDS("data/derived/prep/phylogeny_axes.rds")
fg <- readRDS("data/derived/prep/functional_groups.rds")
analysed <- union(res$species, col$species)
pev <- phy$scores |> select(species, PEV1, PEV2, PEV3) |> filter(species %in% analysed) |> left_join(select(fg, species, family, graminoid), by = "species")
no_pev <- setdiff(analysed, pev$species) # species not placed on the tree (P3), listed in the report
eig <- phy$eigenvalues |> slice(1:10) |> mutate(PEV = factor(PEV, levels = PEV))
ext <- bind_rows(lapply(c("PEV2", "PEV3"), function(a) bind_rows(slice_min(pev, .data[[a]], n = 2), slice_max(pev, .data[[a]], n = 2), slice_min(pev, PEV1, n = 2), slice_max(pev, PEV1, n = 2)) |> mutate(axis = a))) |> distinct()
pev_panel <- function(a, tag) {
  ggplot(pev, aes(PEV1, .data[[a]])) + geom_point(aes(colour = graminoid), size = 3, alpha = 0.6) +
    ggrepel::geom_text_repel(data = filter(ext, axis == a), aes(label = species), size = 6, fontface = "italic", max.overlaps = Inf, seed = 1) +
    scale_colour_manual(values = c(`TRUE` = "#B8860B", `FALSE` = "grey45"), labels = c(`TRUE` = "graminoid", `FALSE` = "other"), name = NULL) +
    labs(title = tag) + big + theme(legend.position = "bottom")
}
figS3 <- patchwork::wrap_plots(pev_panel("PEV2", "a)"), pev_panel("PEV3", "b)"),
                               ggplot(eig, aes(PEV, percent_positive_variation, fill = axis_number <= 3)) + geom_col() +
                                 scale_fill_manual(values = c(`TRUE` = "grey25", `FALSE` = "grey75"), guide = "none") +
                                 labs(title = "c)", x = NULL, y = "Variation (%)") + big, nrow = 1, widths = c(1, 1, 0.8))
ends <- bind_rows(lapply(c("PEV1", "PEV2", "PEV3"), function(a) {
  bind_rows(mutate(slice_min(pev, .data[[a]], n = 15), end = "negative"), mutate(slice_max(pev, .data[[a]], n = 15), end = "positive")) |>
    count(end, family, sort = TRUE) |> slice_head(n = 3, by = end) |>
    summarise(families = paste0(family, " (", n, ")", collapse = ", "), .by = end) |> mutate(axis = a)
}))
tableS3 <- tidyr::pivot_wider(ends, names_from = end, values_from = families, names_prefix = "most common families of the 15 most extreme species, ") |>
  left_join(transmute(phy$eigenvalues, axis = PEV, `variation (%)` = round(percent_positive_variation, 1)), by = "axis") |> relocate(axis, `variation (%)`)

# ---- Fig S5: the predictors before standardisation -------------------------------------
traits <- readRDS("data/derived/prep/traits.rds")
niche <- readRDS("data/derived/prep/niche_optima.rds") |> filter(n_obs_sub >= 30)
a0 <- readRDS("data/derived/01_relative_cover/a0.rds") |> semi_join(res, by = c("experiment_id", "species"))
stopifnot(nrow(a0) == 1004) # A0 for every resident effect
chk <- read.csv("data/derived/prep/site_bioclim_check.csv") |> filter(experiment_id %in% res$experiment_id)
res_sp <- unique(res$species)
hist_panel <- function(v, xlab, tag, log_x = FALSE) {
  d <- tibble::tibble(v = v[!is.na(v)])
  ggplot(d, aes(v)) + geom_histogram(bins = 25, fill = "grey55", colour = "white") +
    (if (log_x) scale_x_log10(labels = scales::label_number(big.mark = "")) else scale_x_continuous()) +
    labs(title = tag, x = paste0(xlab, " (n = ", nrow(d), ")"), y = NULL) + big + theme(axis.title.x = element_text(size = 24))
}
figS5 <- patchwork::wrap_plots(
  hist_panel(experiments$destination_T_JJA, "Destination JJA temperature (°C)", "a)"),
  hist_panel(experiments$diff_T_JJA, "Warming W (°C)", "b)"),
  hist_panel(experiments$destination_P_JJA, "Destination JJA precipitation (mm)", "c)", TRUE),
  hist_panel(a0$A0_origin, "A0 origin (relative cover)", "d)"),
  hist_panel(a0$A0_dest, "A0 destination (relative cover)", "e)"),
  hist_panel(1000 * traits$LDMC[traits$species %in% res_sp], "LDMC (mg g⁻¹)", "f)", TRUE),
  hist_panel(1000 * traits$LMA[traits$species %in% res_sp], "LMA (g m⁻²)", "g)", TRUE),
  hist_panel(100 * traits$height[traits$species %in% res_sp], "Vegetative height (cm)", "h)", TRUE),
  hist_panel(niche$T_opt[niche$species %in% res_sp], "T optimum (°C)", "i)"),
  hist_panel(exp(niche$P_opt[niche$species %in% res_sp]), "P optimum (mm)", "j)", TRUE),
  ggplot(chk, aes(bio10, T)) + geom_abline(linetype = "dashed", colour = "grey50") + geom_point(size = 4, alpha = 0.7) +
    labs(title = "k)", x = "CHELSA bio10 (°C)", y = "Ladder JJA T (°C)") + big,
  ggplot(chk, aes(bio18, destination_P_JJA)) + geom_abline(linetype = "dashed", colour = "grey50") + geom_point(size = 4, alpha = 0.7) +
    scale_x_log10() + scale_y_log10() + labs(title = "l)", x = "CHELSA bio18 (mm)", y = "Ladder JJA P (mm)") + big,
  ncol = 4)

# ---- Tables S1, S2, S4, S5 -----------------------------------------------------------
cover <- readRDS("data/derived/01_relative_cover/cover_long.rds") |> filter(experiment_id %in% res$experiment_id)
plots <- distinct(cover, experiment_id, treatment, plot_id) |> count(experiment_id, treatment) |>
  tidyr::pivot_wider(names_from = treatment, values_from = n, names_prefix = "plots ")
pools <- read.csv("data/derived/01_relative_cover/pool_sizes.csv")
tableS1 <- experiments |>
  transmute(experiment_id, country, origin_site, destination_site, origin_elevation, destination_elevation, first_year, last_year,
            `destination JJA T (°C)` = round(destination_T_JJA, 2), `warming W (°C)` = round(diff_T_JJA, 2), `destination JJA P (mm)` = round(destination_P_JJA)) |>
  left_join(plots, by = "experiment_id") |> left_join(select(pools, experiment_id, `resident pool` = resident, `destination pool` = destination_pool), by = "experiment_id")
stopifnot(nrow(tableS1) == 42, !anyNA(tableS1$`resident pool`))
tableS2 <- tibble::tibble(species = sort(analysed)) |>
  left_join(select(fg, species, family, graminoid, woody), by = "species") |>
  left_join(summarise(res, `gradients as resident` = n_distinct(gradient_id), .by = species), by = "species") |>
  left_join(summarise(col, `gradients as colonizer` = n_distinct(gradient_id), .by = species), by = "species") |>
  left_join(select(s04, species, `Step C class (resident)` = class), by = "species") |>
  left_join(select(s07b, species, `Step C class (colonizer)` = class), by = "species") |>
  mutate(LDMC = species %in% traits$species[!is.na(traits$LDMC)], LMA = species %in% traits$species[!is.na(traits$LMA)],
         height = species %in% traits$species[!is.na(traits$height)], PEV = species %in% phy$scores$species, niche = species %in% niche$species)
stopifnot(nrow(tableS2) == length(analysed), !anyNA(tableS2$graminoid))
tableS4 <- bind_rows(mutate(read.csv("data/derived/05_ladder_fit/tableS4.csv"), pool = "residents"),
                     mutate(read.csv("data/derived/07c_ladder_fit/heterogeneity.csv"), pool = "colonizers")) |> relocate(pool)
tableS5 <- bind_rows(mutate(read.csv("data/derived/05_ladder_fit/tableS5.csv"), pool = "residents"),
                     mutate(read.csv("data/derived/07c_ladder_fit/coefficients.csv"), pool = "colonizers")) |> relocate(pool)
for (t in c("S1", "S2", "S3", "S4", "S5")) write.csv(get(paste0("table", t)), sprintf("tables/table%s.csv", t), row.names = FALSE)

# ---- Write the new figures; copy the reviewed ones under their final numbers ------------
ggsave("figures/main/fig4.png", fig4, width = 38, height = 28, dpi = 100, bg = "white")
ggsave("figures/main/fig7.png", fig7, width = 40, height = 26, dpi = 100, bg = "white")
ggsave("figures/appendix/figS1.png", figS1, width = 34, height = 22, dpi = 100, bg = "white")
ggsave("figures/appendix/figS2.png", figS2, width = 38, height = 28, dpi = 100, bg = "white")
ggsave("figures/appendix/figS3.png", figS3, width = 40, height = 14, dpi = 100, bg = "white")
ggsave("figures/appendix/figS5.png", figS5, width = 40, height = 26, dpi = 100, bg = "white")
copies <- c("figures/main/fig3.png" = "figures/diagnostics/07b_figS2.png",
            "figures/main/fig6.png" = "figures/diagnostics/06_fig4.png",
            "figures/appendix/figS4a.png" = "figures/diagnostics/05_figS4.png",
            "figures/appendix/figS4b.png" = "figures/diagnostics/07c_figS4b.png",
            "figures/appendix/figS6a.png" = "figures/diagnostics/06_figS7.png",
            "figures/appendix/figS6b.png" = "figures/diagnostics/07d_figS7.png",
            "figures/appendix/figS7.png" = "figures/diagnostics/07d_figS9.png",
            "figures/appendix/figS8a.png" = "figures/diagnostics/06b_figS11.png",
            "figures/appendix/figS8b.png" = "figures/diagnostics/07e_figS11.png")
stopifnot(all(file.copy(copies, names(copies), overwrite = TRUE)))
fig_index <- tibble::tribble(
  ~figure, ~file, ~content, ~from,
  "Fig 1", "figures/main/fig1.png", "Sankey of year-specific classes, residents and colonizers", "steps 03, 07a",
  "Fig 2", "figures/main/fig2.png", "Residents' classification (Step C)", "step 04",
  "Fig 3", "figures/main/fig3.png", "Colonizers' classification (Step C)", "step 07b",
  "Fig 4", "figures/main/fig4.png", "The experiments in climate space, pies of warmed-turf cover", "this step",
  "Fig 5", "figures/main/fig5.png", "Traits of losers, winners and colonizers (pathways)", "step 07f",
  "Fig 6", "figures/main/fig6.png", "The residents' ladder", "step 06",
  "Fig 7", "figures/main/fig7.png", "Coefficients of the climate + abundance, Height and LDMC models, residents beside colonizers", "this step",
  "Fig S1", "figures/appendix/figS1.png", "Data flow", "this step",
  "Fig S2", "figures/appendix/figS2.png", "Fig 4 with species counts in the pies", "this step",
  "Fig S3", "figures/appendix/figS3.png", "Phylogenetic axes (with Table S3)", "this step",
  "Fig S4a, b", "figures/appendix/figS4a.png", "Coefficients, residents (a) and colonizers (b)", "steps 05, 07c",
  "Fig S5", "figures/appendix/figS5.png", "Predictors before standardisation, niche checks", "this step",
  "Fig S6a, b", "figures/appendix/figS6a.png", "Held-out metrics per gradient, residents (a) and colonizers (b)", "steps 06, 07d",
  "Fig S7", "figures/appendix/figS7.png", "The colonizers' ladder", "step 07d",
  "Fig S8a, b", "figures/appendix/figS8a.png", "Random forest beside metafor, residents (a) and colonizers (b)", "steps 06b, 07e")
main_files <- sort(list.files("figures/main", pattern = "^fig[0-9]+\\.png$"))
stopifnot(identical(main_files, paste0("fig", 1:7, ".png"))) # seven main figures, each exactly once
stopifnot(all(file.exists(fig_index$file)), all(file.exists(sub("a\\.png$", "b.png", fig_index$file[grepl("a, b", fig_index$figure)]))))

# ---- Report --------------------------------------------------------------------------
count_tab <- tibble::tibble(
  item = c("experiments", "resident effects (03b)", "net strength below zero", "net strength above zero", "colonizer effects (07a)",
           "net warmed cover above zero", "experiments with colonizers", "experiments without purple", "species analysed (residents or colonizers)"),
  n = c(nrow(experiments), nrow(res), sum(res$strength < 0), sum(res$strength > 0), nrow(col), sum(col$strength > 0),
        sum(pies$Colonizers > 0), sum(pies$Colonizers == 0), length(analysed)))
report <- c(
  "# Step 10: Figures, appendix and report", "",
  paste0("Draws what no earlier step produced: Fig 4, the ", nrow(experiments), " experiments in climate space with one pie per experiment ",
         "splitting the warmed turf's relative cover, years 3 and later, among resident losers and winners (by the sign of their net strength from year 3 on), ",
         "colonizers (net warmed cover above zero) and what was not analysed (species outside the pools, unidentified and genus-level records) (Bills, 2 Nov); ",
         "Fig S2 the same map with the pies counting species; Fig S1 the data flow; Fig S3 and Table S3 the phylogenetic axes; Fig S5 the predictors before ",
         "standardisation with the niche checks; Tables S1, S2, S4 and S5. Copies the reviewed figures of steps 03 to 07f under their final numbers ",
         "(seven main figures, Bills 2 Nov; Fig 7 the coefficients of the climate + abundance, Height and LDMC models, residents beside colonizers). Pies are drawn to one radius at their destination climate; the ", sum(pies$n_shared > 1), " experiments that share a destination climate ",
         "(", n_distinct(paste(pies$x0, pies$y0)[pies$n_shared > 1]), " points) fan out around a black dot at that climate on short lines; no two pies overlap (", sum(pies$moved & pies$n_shared == 1), " pies of different climates had to be pushed apart). ",
         "Each pie is named by its gradient (with a letter where a gradient has several experiments), coloured by warming; the key beside the figure gives each one's origin and destination elevation."),
  "", "## Counts", "", md_table(count_tab),
  "", "## Fig 4 and Fig S2 pies", "", "Cover shares (Fig 4) in per cent of the warmed turf; species counts (Fig S2) in brackets.", "",
  md_table(transmute(left_join(pies, rename_with(cover_shares, ~ paste0("cover_", .x), -experiment_id), by = "experiment_id"),
                     experiment_id, name = label, elevations, `T (°C)` = round(destination_T_JJA, 1), `P (mm)` = round(destination_P_JJA),
                     `W (°C)` = round(diff_T_JJA, 1),
                     `Resident losers` = sprintf("%.0f%% (%d)", 100 * `cover_Resident losers`, `Resident losers`),
                     `Resident winners` = sprintf("%.0f%% (%d)", 100 * `cover_Resident winners`, `Resident winners`),
                     Colonizers = sprintf("%.0f%% (%d)", 100 * cover_Colonizers, Colonizers),
                     `Not analysed` = sprintf("%.0f%%", 100 * `cover_Not analysed`))),
  "", paste0("Mean over experiments, per cent of the warmed turf: resident losers ", round(100 * mean(cover_shares$`Resident losers`)),
             ", resident winners ", round(100 * mean(cover_shares$`Resident winners`)), ", colonizers ", round(100 * mean(cover_shares$Colonizers)),
             ", not analysed ", round(100 * mean(cover_shares$`Not analysed`)), "."),
  "", "## Figures", "", md_table(fig_index), "",
  "![Fig 4](../figures/main/fig4.png)", "", "![Fig 7](../figures/main/fig7.png)", "", "![Fig S1](../figures/appendix/figS1.png)", "", "![Fig S2](../figures/appendix/figS2.png)", "",
  "![Fig S3](../figures/appendix/figS3.png)", "", "![Fig S5](../figures/appendix/figS5.png)",
  "", "## Tables", "",
  md_table(tibble::tibble(table = paste("Table", c("S1", "S2", "S3", "S4", "S5")), file = sprintf("tables/table%s.csv", c("S1", "S2", "S3", "S4", "S5")),
                          content = c("Experiments", "Species list", "Phylogenetic axes", "Variance components per model", "Fixed-effect coefficients per model"),
                          rows = c(nrow(tableS1), nrow(tableS2), nrow(tableS3), nrow(tableS4), nrow(tableS5)))),
  "", "## Table S3, phylogenetic axes", "", md_table(tableS3),
  "", paste0(length(no_pev), " analysed species have no phylogenetic scores", if (length(no_pev)) paste0(": ", paste(no_pev, collapse = ", ")) else "", "."),
  "", "## Warnings captured", "", if (length(warns)) paste("-", unique(warns)) else "None.",
  "", "## Deviations from PLAN.md", "",
  "- One script for the whole step (Bills, 2 Nov).",
  "- Seven main figures (Bills, 2 Nov): the colonizers' classification becomes Fig 3, the climate space Fig 4, the pathways Fig 5, the ladder Fig 6, the coefficients of three models Fig 7 (abundance terms drawn per pool, not comparable); the results-dependent Fig 6 and the old Fig S8 dropped; appendix renumbered with the count version of Fig 4 as Fig S2.",
  "- Fig 4's pies split the warmed turf's relative cover (years 3 and later) among Fig 5's classes per experiment, not counts of last-year significance classes; the counts are Fig S2; no purple in the experiments below the colonizer pool-size rule (Bills, 2 Nov).",
  "- Coefficients, per-gradient metrics and the random forest stay as residents (a) and colonizers (b) files.",
  "- `REPORT.html` not rendered: Quarto is not installed.")
writeLines(report, "reports/10_figures.md")
message("Step 10 written: ", nrow(experiments), " pies, ", nrow(fig_index), " figure entries, 5 tables")
