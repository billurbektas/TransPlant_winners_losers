# Step 10: figures, appendix and report.
#
# PLAN.md: Figure map, Captions, Appendix.
# Reads:  data/derived/prep/ (experiments, community, sites, experiment_climate, site_bioclim_check,
#         traits, niche_optima, functional_groups, phylogeny_axes);
#         data/derived/01_relative_cover/ (flow_counts.csv, cover_long.rds, a0.rds,
#         pool_sizes.csv); data/derived/03b_strength_y3/strength_df.rds;
#         data/derived/07a_colonizers/strength_df_y3.rds; the species summaries of steps
#         04 and 07b; the tables of steps 05 and 07c (weights.csv and ladder_rows.rds name the
#         ladder rows); the appendix figures written by the steps.
# Writes: figures/main/fig4.png, fig7.png (the other main figures are written by their steps);
#         figures/appendix/figS1.png, figS3.png, figS4.png, figS10.png (the other appendix
#         figures are written by their steps);
#         tables/tableS1.csv ... tableS6.csv; data/derived/10_figures/fig4_balance.csv,
#         fig4_pie_counts.csv;
#         reports/10_figures.md
# Checks: 42 experiments; 1004 resident effects (703 losers, 301 winners); 264
#         colonizer effects with warmed cover above zero in 28 experiments; every pie's
#         slices sum to its species; Fig 4 on the 512 resident and 110 colonizer ladder
#         effects (39 experiments with modelled residents, 29 with a net loss, 22 with a
#         colonizer gain of the 24 with modelled colonizers), no two discs overlapping;
#         every figure number exists exactly once.
#
# Decided 2 Nov (Bills): one script for the whole step, run without stopping. Seven main
# figures: 1 Sankey, 2 residents' classification, 3 colonizers' classification, 4
# climate space, 5 traits (pathways), 6 the residents' ladder; the results-dependent
# Fig 6 and the old Fig S8 are dropped. Fig S10's pies count Fig 5's classes per experiment:
# residents by the sign of their net strength from year 3 on, colonizers with net warmed
# cover above zero (07a, so the 14 experiments below the colonizer pool-size rule have no
# purple); one radius, pies sharing a destination climate fanned out around it. Fig 4 (3 Nov)
# draws the same places as discs for the modelled effects: the residents' net change and
# the colonizers' gain. Each appendix figure is written by its own step.

suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
source("R/md_table.R"); source("R/term_labels.R")
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
stopifnot(sum(res$strength < 0) == 703, sum(res$strength > 0) == 301, sum(col$strength > 0) == 264) # Fig 5's classes per effect (262 until CH_Calanda2's year 4 was regained, 3 Nov)
pies <- tidyr::pivot_wider(pies, names_from = group, values_from = n, values_fill = 0)
pies <- left_join(select(experiments, experiment_id, origin_elevation, destination_elevation, destination_T_JJA, destination_P_JJA, diff_T_JJA),
                  pies, by = "experiment_id") |> mutate(across(all_of(names(cols)), ~ coalesce(.x, 0L)))
stopifnot(sum(pies$`Resident losers` + pies$`Resident winners`) == 1004, sum(pies$Colonizers) == 264,
          sum(pies$Colonizers > 0) == n_distinct(col$experiment_id[col$strength > 0])) # one denominator: counts, not cover
write.csv(pies, file.path(out_dir, "fig4_pie_counts.csv"), row.names = FALSE)
# main Fig 4 (Bills, 3 Nov) shows, for the effects the ladders model (steps 05 and 07c), one disc
# per experiment: its fill the residents' net change, the sum of their strengths with the sign
# kept, so that noise cancels instead of piling up as it does in sums of |strength|; its purple
# ring the colonizers' total warmed cover, measured cover rather than a magnitude. Both are
# percentage points of the warmed turf from year 3 on, and the colonizers' strength is also their
# change from the origin controls, where they were never recorded. The pies of Fig S10 count species
lad_res <- read.csv("data/derived/05_ladder_fit/weights.csv") |> filter(rowset == "all", model == "Null")
lad_col <- readRDS("data/derived/07c_ladder_fit/ladder_rows.rds")
res_lad <- filter(res, effect_id %in% lad_res$effect_id)
col_lad <- filter(col, effect_id %in% lad_col$effect_id)
stopifnot(nrow(res_lad) == 512, nrow(col_lad) == 110) # the ladder rows of steps 05 and 07c
balance <- left_join(select(experiments, experiment_id),
                     summarise(res_lad, net = 100 * sum(strength), .by = experiment_id), by = "experiment_id") |>
  left_join(summarise(col_lad, gain = 100 * sum(strength), .by = experiment_id), by = "experiment_id") |>
  mutate(gain = coalesce(gain, 0))
write.csv(balance, file.path(out_dir, "fig4_balance.csv"), row.names = FALSE)
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
warm_scale <- scale_colour_gradient(low = "#C99A1C", high = "#7B0D1E", limits = range(pies$diff_T_JJA), name = "Warming\n(°C, JJA)", guide = guide_colourbar(order = 1)) # order as reviewed
p_breaks <- c(150, 200, 300, 400, 600, 800, 1000)
map_under <- function(d) list( # the same layout in Fig 4 and Fig S10: lines and dots under the marks
  geom_segment(data = filter(d, moved), aes(x0, y0, xend = x, yend = y), colour = "grey25", linewidth = 0.9),
  geom_point(data = distinct(filter(d, moved), x0, y0), aes(x0, y0), size = 4)) # the true climate of a moved pie
map_over <- function(d) list( # labels, scales and theme over the marks
  ggrepel::geom_text_repel(aes(x, y, label = label, colour = diff_T_JJA), size = 7.5, fontface = "bold", bg.colour = "white", bg.r = 0.15,
                           point.size = 24, point.padding = 0, box.padding = 0.25,
                           nudge_x = ifelse(d$n_shared > 1, 2.1 * r * cos(d$angle), 0), # fanned pies: label outward, off the shared point
                           nudge_y = ifelse(d$n_shared > 1, 2.1 * r * sin(d$angle), -1.9 * r), force = 5, force_pull = 3, # point.size = the pie's size
                           max.iter = 2e5, max.time = 120, # the iteration count, not the clock, ends the search: same labels on every machine
                           min.segment.length = 0.15, segment.colour = "grey45", segment.size = 0.6, max.overlaps = Inf, seed = 1),
  warm_scale, scale_y_continuous(breaks = log(p_breaks) * k, labels = p_breaks), coord_equal(),
  labs(x = "Destination June–August mean temperature (°C)", y = "Destination June–August\nprecipitation (mm, log scale)"),
  big, theme(legend.position = "bottom", legend.key.size = unit(1.4, "cm"), legend.title = element_text(face = "bold"),
             panel.grid.major = element_line(colour = "grey90")))
# Fig S10: pies counting species, all classified effects
figS10_map <- ggplot(pies) + map_under(pies) +
  scatterpie::geom_scatterpie(aes(x = x, y = y, r = r), data = mutate(pies, r = r), cols = names(cols), colour = "white", linewidth = 0.3) +
  map_over(pies) + scale_fill_manual(values = cols, breaks = names(cols), name = NULL, guide = guide_legend(order = 2)) + # legend order fixed: ggplot2 otherwise orders them differently between runs
  theme(legend.box = "vertical", legend.key.width = unit(3, "cm"))
# Fig 4: disc fill capped at ±30 points (beyond drawn at full colour), ring width up to w_max for the
# largest gain; disc and ring stay inside the pie radius, so the pies' spacing holds
r_in <- 0.16; w_max <- 0.07; net_lim <- 30
d4 <- left_join(pies, balance, by = "experiment_id") |> mutate(w = w_max * gain / max(gain))
ro <- r_in + d4$w; ring_d <- as.matrix(dist(cbind(d4$x, d4$y)))
stopifnot(sum(!is.na(d4$net)) == 39, sum(d4$net < 0, na.rm = TRUE) == 29, sum(d4$gain > 0) == 22, # experiments with modelled residents, net losses, colonizer gains
          n_distinct(col_lad$experiment_id) == 24, # 24 with modelled colonizers; in 2 every one stays at zero warmed cover
          !any(ring_d[upper.tri(ring_d)] < outer(ro, ro, "+")[upper.tri(ring_d)])) # no two discs overlap
ring_key <- tibble::tibble(gain = c(5, 15, 30), w = w_max * gain / max(d4$gain))
fig4_map <- function(mm_per_degree) { # the ring key is drawn in mm, the rings in °C: matched through the panel's scale
  ggplot(d4) + map_under(d4) +
    ggforce::geom_circle(aes(x0 = x, y0 = y, r = r_in + w), data = filter(d4, gain > 0), fill = "#5E3C99", colour = NA) +
    ggforce::geom_circle(aes(x0 = x, y0 = y, r = r_in, fill = net), colour = "grey30", linewidth = 0.4) +
    # invisible points carry the ring key: a pch 21 circle is 0.753 mm across per unit of size and its
    # line 0.5 mm wide per unit of stroke, drawn through the middle of the ring
    geom_point(aes(x = min(d4$x), y = min(d4$y), size = factor(gain)), data = ring_key, inherit.aes = FALSE, alpha = 0) +
    scale_size_manual(values = setNames(2 * (r_in + ring_key$w / 2) * mm_per_degree / 0.753, ring_key$gain),
                      name = "Colonizers'\ngain (% points)",
                      guide = guide_legend(order = 3, override.aes = list(alpha = 1, shape = 21, fill = "white", colour = "#5E3C99",
                                                                          stroke = ring_key$w * mm_per_degree / 0.5))) +
    map_over(d4) +
    scale_fill_gradient2(low = cols[["Resident losers"]], mid = "white", high = cols[["Resident winners"]], midpoint = 0,
                         limits = c(-net_lim, net_lim), oob = scales::squish, breaks = c(-30, -15, 0, 15, 30),
                         labels = c("≤ −30", "−15", "0", "15", "≥ 30"), na.value = "grey85", name = "Residents' net\nchange (% points)",
                         guide = guide_colourbar(order = 2, barwidth = unit(10, "cm"))) +
    guides(colour = guide_colourbar(order = 1, barwidth = unit(10, "cm"))) +
    theme(legend.box = "horizontal", legend.spacing.x = unit(2, "cm"), legend.title.position = "top")
}
panel_mm_per_degree <- function(map, fig) { # width of the drawn map panel over its temperature range
  f <- tempfile(fileext = ".png"); ragg::agg_png(f, width = 38, height = 28, units = "in", res = 100)
  grid::grid.newpage(); grid::grid.draw(patchwork::patchworkGrob(fig)); grid::grid.force(); grid::seekViewport("panel.1-1-1-1")
  w <- grid::convertWidth(grid::unit(1, "npc"), "mm", valueOnly = TRUE)
  x_range <- ggplot_build(map)$layout$panel_params[[1]]$x.range # built while the device is open: no Rplots.pdf
  invisible(dev.off()); unlink(f)
  w / diff(x_range)
}
key <- arrange(pies, desc(diff_T_JJA)) |> mutate(row = row_number(), name = sub("_", " ", gradient_id)) # strongest warming first (Bills, 2 Nov)
fig4key <- ggplot(key, aes(y = -row)) +
  geom_text(aes(x = 0, label = name, colour = diff_T_JJA), hjust = 0, size = 10, fontface = "bold") +
  geom_text(aes(x = 3.7, label = letter, colour = diff_T_JJA), hjust = 0.5, size = 10, fontface = "bold") +
  geom_text(aes(x = 4.0, label = elevations, colour = diff_T_JJA), hjust = 0, size = 10) +
  warm_scale + guides(colour = "none") + scale_x_continuous(limits = c(0, 6.6), expand = c(0, 0)) +
  labs(title = "Origin → destination elevation,\nstrongest warming first") + theme_void() +
  theme(plot.title = element_text(face = "bold", size = 34, margin = margin(b = 10)), plot.background = element_rect(fill = "white", colour = NA))
figS10 <- patchwork::wrap_plots(figS10_map, fig4key, ncol = 2, widths = c(2.2, 0.7))
# the legend row sets the panel's size and the ring key depends on it: redraw until the scale holds
mm_per_degree <- 40
for (pass in 1:5) {
  fig4 <- patchwork::wrap_plots(fig4_map(mm_per_degree), fig4key, ncol = 2, widths = c(2.2, 0.7))
  measured <- panel_mm_per_degree(fig4_map(mm_per_degree), fig4)
  if (abs(measured - mm_per_degree) < 0.01) break
  mm_per_degree <- measured
}
stopifnot(abs(measured - mm_per_degree) < 0.01) # the ring key matches the rings on the page

# ---- Fig 7: what the models find, residents beside colonizers --------------------------
# Climate + abundance, Height and LDMC models (Bills, 2 Nov). Both responses are relative cover
# and W, T, P and the traits share one scale (steps 02, 07c), so those terms sit on one axis;
# the abundance terms are standardised within each pool and are drawn apart
coefs <- bind_rows(mutate(read.csv("data/derived/05_ladder_fit/tableS6.csv"), pool = "Residents"),
                   mutate(read.csv("data/derived/07c_ladder_fit/coefficients.csv"), pool = "Colonizers")) |>
  filter(model %in% c("Climate + start", "Height", "LDMC"), term != "intrcpt") |>
  mutate(pool = factor(pool, levels = c("Residents", "Colonizers")), excludes_zero = ci_lower > 0 | ci_upper < 0)
climate_terms <- c(W = "Warming", T = "Temperature", P = "Precipitation", `W:T` = "Warming × temperature", `W:P` = "Warming × precipitation",
                   `T:P` = "Temperature × precipitation", `W:T:P` = "Warming × temperature × precipitation")
# residents' abundances are initial (first survey of the origin or destination controls, step 01);
# the colonizers' is the mean destination-control cover over every plot and year (07a)
abund_names <- c("Residents A0_origin" = "Initial abundance at origin", "Residents A0_dest" = "Initial abundance at destination",
                 "Colonizers A0_dest" = "Mean abundance at destination", "Residents W:A0_origin" = "Warming × initial abundance at origin",
                 "Residents W:A0_dest" = "Warming × initial abundance at destination", "Colonizers W:A0_dest" = "Warming × mean abundance at destination")
trait_terms <- function(z, name) setNames(paste(name, c("", "× warming", "× temperature", "× precipitation", "× temperature × precipitation")),
                                         paste0(c("", "W:", "T:", "P:", ""), z, c("", "", "", "", ":TP")))
prop_terms <- c(trait_terms("height_z", "Height"), trait_terms("ldmc_z", "LDMC"))
f7_clim <- filter(coefs, model == "Climate + start", term %in% names(climate_terms)) |> mutate(name = factor(climate_terms[term], levels = rev(climate_terms)))
f7_abund <- filter(coefs, model == "Climate + start", grepl("A0", term)) |> mutate(name = factor(abund_names[paste(pool, term)], levels = rev(abund_names)))
f7_prop <- filter(coefs, model %in% c("Height", "LDMC"), term %in% names(prop_terms)) |> mutate(name = factor(trimws(prop_terms[term]), levels = rev(trimws(prop_terms))))
stopifnot(nrow(f7_clim) == 14, nrow(f7_abund) == 6, nrow(f7_prop) == 20) # 7 climate terms and 5 terms per trait in each pool; A0_origin residents only
pool_cols <- c(Residents = "#3D3D3D", Colonizers = "#5E3C99")
x_shared <- range(c(f7_clim$ci_lower, f7_clim$ci_upper, f7_abund$ci_lower, f7_abund$ci_upper, f7_prop$ci_lower, f7_prop$ci_upper)) # every panel on one axis
# a) and b): the climate + abundance model, one graph per pool, climate terms above the abundance terms (Bills)
f7_model <- bind_rows(f7_clim, f7_abund) |> mutate(name = factor(as.character(name), levels = rev(c(climate_terms, abund_names)))) # the same rows in a and b
coef_panel <- function(d, tag, xlim = NULL) {
  ggplot(d, aes(estimate, name, colour = pool, group = pool)) + geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_linerange(aes(xmin = ci_lower, xmax = ci_upper), linewidth = 1.6, position = position_dodge(width = 0.6)) +
    geom_point(aes(fill = ifelse(excludes_zero, as.character(pool), "none")), shape = 21, size = 6, stroke = 1.6, position = position_dodge(width = 0.6)) +
    scale_colour_manual(values = pool_cols, name = NULL, breaks = names(pool_cols)) +
    scale_fill_manual(values = c(pool_cols, none = "white"), guide = "none") +
    (if (!is.null(xlim)) scale_x_continuous(limits = xlim)) +
    labs(title = tag, x = "Effect on relative cover per SD of the predictor", y = NULL) + big +
    theme(legend.position = "bottom", legend.text = element_text(size = 30), legend.key.size = unit(1.6, "cm"),
          panel.grid.major.y = element_line(colour = "grey45", linewidth = 1.4, linetype = "dotted")) # row guides from name to estimate (Bills)
}
fig7 <- patchwork::wrap_plots(
  coef_panel(filter(f7_model, pool == "Residents"), "a) Residents: climate + abundance model", x_shared) + guides(colour = "none") + # the title names the pool
    scale_y_discrete(drop = FALSE), # a row stays empty where the model has no such term
  coef_panel(filter(f7_model, pool == "Colonizers"), "b) Colonizers: climate + abundance model", x_shared) + guides(colour = "none") +
    scale_y_discrete(drop = FALSE),
  coef_panel(filter(f7_prop, model == "Height"), "c) Height model", x_shared),
  coef_panel(filter(f7_prop, model == "LDMC"), "d) LDMC model", x_shared),
  ncol = 2, guides = "collect") & theme(legend.position = "bottom", plot.title = element_text(size = 34))

# ---- Fig S1: data flow ------------------------------------------------------------------
flow <- read.csv("data/derived/01_relative_cover/flow_counts.csv")
stopifnot(nrow(flow) == 7) # step 01's filters
s04 <- read.csv("data/derived/04_classification/species_summary.csv")
s07b <- read.csv("data/derived/07b_classification/species_summary.csv")
h05 <- read.csv("data/derived/05_ladder_fit/heterogeneity.csv") |> filter(rowset == "all", model == "Null")
box <- function(e, s, n) sprintf("%d experiments · %d species · %d combinations", e, s, n) # species × experiment, as the captions
excl <- filter(readRDS("data/derived/prep/experiments.rds"), !included)
stopifnot(nrow(excl) == flow$n_experiments[1] - flow$n_experiments[2]) # the experiments removed between the first two boxes
excl_rule <- paste(sprintf("%s (%s → %s) removed: %s, %d survey years (%d–%d)", sub("_", " ", excl$region), excl$origin_site, excl$destination_site,
                           excl$exclusion_reason, excl$year_range, excl$first_year, excl$last_year), collapse = "\n")
nodes <- tibble::tribble(
  ~id, ~x, ~y, ~title, ~counts,
  1, 0, 9, "Community data", box(flow$n_experiments[1], flow$n_species[1], flow$n_combinations[1]),
  2, 0, 8, "Experiments kept", box(flow$n_experiments[2], flow$n_species[2], flow$n_combinations[2]),
  3, 0, 7, "Deep lineages removed", box(flow$n_experiments[3], flow$n_species[3], flow$n_combinations[3]),
  4, 0, 6, "Pool rule", box(flow$n_experiments[4], flow$n_species[4], flow$n_combinations[4]),
  5, 0, 5, "Pool-size rule", box(flow$n_experiments[5], flow$n_species[5], flow$n_combinations[5]),
  6, -1, 4, "Residents: strength from year 3", box(n_distinct(res$experiment_id), n_distinct(res$species), nrow(res)),
  7, -1, 3, "Residents: classification (Fig 2)", box(NA_integer_, nrow(s04), sum(s04$n_experiments)),
  8, -1, 2, "Residents: model comparison (Fig 6)", box(h05$n_experiments, h05$n_species, h05$n_effects),
  9, 1, 4, "Colonizers: warmed cover from year 3", box(n_distinct(col$experiment_id), n_distinct(col$species), nrow(col)),
  10, 1, 3, "Colonizers: classification (Fig 3)", box(NA_integer_, nrow(s07b), sum(s07b$n_experiments)),
  11, 1, 2, "Colonizers: model comparison (Fig S17)", box(n_distinct(lad_col$experiment_id), n_distinct(lad_col$species), nrow(lad_col))) |>
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

# ---- Fig S3 and Table S4: phylogenetic axes -------------------------------------------
phy <- readRDS("data/derived/prep/phylogeny_axes.rds")
fg <- readRDS("data/derived/prep/functional_groups.rds")
analysed <- union(res$species, col$species)
pev <- phy$scores |> select(species, PEV1, PEV2, PEV3) |> filter(species %in% analysed) |> left_join(select(fg, species, family, graminoid), by = "species")
no_pev <- setdiff(analysed, pev$species) # species not placed on the tree (P3), listed in the report
eig <- phy$eigenvalues |> slice(1:10) |> mutate(PEV = factor(PEV, levels = PEV))
ext <- bind_rows(lapply(c("PEV2", "PEV3"), function(a) bind_rows(slice_min(pev, .data[[a]], n = 2), slice_max(pev, .data[[a]], n = 2), slice_min(pev, PEV1, n = 2), slice_max(pev, PEV1, n = 2)) |> mutate(axis = a))) |> distinct()
pev_panel <- function(a, tag) {
  ggplot(pev, aes(PEV1, .data[[a]])) + geom_point(aes(colour = graminoid), size = 3, alpha = 0.6) +
    ggrepel::geom_text_repel(data = filter(ext, axis == a), aes(label = species), size = 6, fontface = "italic", max.overlaps = Inf, seed = 1, max.time = 120) +
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
tableS4 <- tidyr::pivot_wider(ends, names_from = end, values_from = families, names_prefix = "most common families of the 15 most extreme species, ") |>
  left_join(transmute(phy$eigenvalues, axis = PEV, `variation (%)` = round(percent_positive_variation, 1)), by = "axis") |> relocate(axis, `variation (%)`)

# ---- Fig S4: the predictors before standardisation -------------------------------------
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
figS4 <- patchwork::wrap_plots(
  hist_panel(experiments$destination_T_JJA, "Destination June–August temperature (°C)", "a)"),
  hist_panel(experiments$diff_T_JJA, "Warming (°C)", "b)"),
  hist_panel(experiments$destination_P_JJA, "Destination June–August precipitation (mm)", "c)", TRUE),
  hist_panel(a0$A0_origin, "Initial abundance at origin, log(% cover + 1)", "d)"),
  hist_panel(a0$A0_dest, "Initial abundance at destination, log(% cover + 1)", "e)"),
  hist_panel(1000 * traits$LDMC[traits$species %in% res_sp], "LDMC (mg g⁻¹)", "f)", TRUE),
  hist_panel(1000 * traits$LMA[traits$species %in% res_sp], "LMA (g m⁻²)", "g)", TRUE),
  hist_panel(100 * traits$height[traits$species %in% res_sp], "Vegetative height (cm)", "h)", TRUE),
  hist_panel(niche$T_opt[niche$species %in% res_sp], "Thermal optimum (°C)", "i)"),
  hist_panel(exp(niche$P_opt[niche$species %in% res_sp]), "Precipitation optimum (mm)", "j)", TRUE),
  ggplot(chk, aes(bio10, T)) + geom_abline(linetype = "dashed", colour = "grey50") + geom_point(size = 4, alpha = 0.7) +
    labs(title = "k)", x = "CHELSA bio10, warmest quarter (°C)", y = "June–August temperature (°C)") + big,
  ggplot(chk, aes(bio18, destination_P_JJA)) + geom_abline(linetype = "dashed", colour = "grey50") + geom_point(size = 4, alpha = 0.7) +
    scale_x_log10() + scale_y_log10() + labs(title = "l)", x = "CHELSA bio18, warmest quarter (mm)", y = "June–August precipitation (mm)") + big,
  ncol = 4)

# ---- Tables S1 to S3, S5 and S6 (S4 is drawn with Fig S3) ---------------------------
# Table S1 (Bills, 4 Nov): what each selection removes, in the order the analysis applies it.
# Residents and colonizers split after the pool rule; their later rows count from the row above
# on their own branch, and nothing is removed at the split itself
st_res <- filter(res, species %in% s04$species); st_col <- filter(col, species %in% s07b$species)
left3 <- function(d) c(n_distinct(d$experiment_id), n_distinct(d$species), nrow(d))
sel <- tibble::tibble(
  step = c("P1", "P1", "01", "01", "01", "01", "01", "07a", "04", "07b", "05", "07c"),
  code = c("prep/R/community.R", "prep/R/community.R", rep("analysis/01_relative_cover.R", 5), "analysis/07a_colonizers.R",
           "analysis/04_classification.R", "analysis/07b_classification.R", "analysis/05_ladder_fit.R", "analysis/07c_ladder_fit.R"),
  rule = c("Network data, species identified to species level",
           "Experiments with too short a record are left out (IN_Kashmir, two survey years)",
           "Ferns, horsetails, clubmosses and juniper are left out: lineages far from the flowering plants the questions concern (their cover still counts in the relative-cover denominator)",
           "A species must be established in the controls: recorded in at least two plots in half of the years, or in half of the plots in at least two years",
           "An experiment needs at least 5 residents for a species to be compared within it (US_Arizona has 4)",
           "Residents: established in the origin controls",
           "Destination pool: established in the destination controls, never recorded in the origin controls",
           "Colonizers: an experiment needs at least 5 destination-pool species, as for residents",
           "Residents in at least 3 gradients, so that their variation among experiments can be estimated",
           "Colonizers in at least 3 gradients, as for residents",
           "Residents with every property (traits, growth form, phylogeny, niche), so all twelve models use the same rows",
           "Colonizers with every property, as for residents"),
  parent = c(NA, 1, 2, 3, 4, 5, 5, 7, 6, 8, 9, 10),
  split = c(rep(FALSE, 5), TRUE, TRUE, rep(FALSE, 5)))
sel_n <- rbind(as.matrix(flow[1:7, c("n_experiments", "n_species", "n_combinations")]), left3(col), left3(st_res), left3(st_col),
               c(h05$n_experiments, h05$n_species, h05$n_effects), left3(lad_col))
removed <- function(k) ifelse(is.na(sel$parent) | sel$split, NA, sel_n[sel$parent, k] - sel_n[, k])
tableS1 <- transmute(sel, step, code, rule) |>
  mutate(experiments_removed = removed(1), experiments_left = sel_n[, 1], species_removed = removed(2), species_left = sel_n[, 2],
         combinations_removed = removed(3), combinations_left = sel_n[, 3])
stopifnot(nrow(tableS1) == 12, all(c(tableS1$experiments_removed, tableS1$species_removed, tableS1$combinations_removed) >= 0, na.rm = TRUE),
          tableS1$combinations_left[11] == 512, tableS1$combinations_left[12] == 110) # a rule only removes; the ladder rows
# plots per group counted in the community data with step 01's three groups, so experiments whose
# destination pool is empty keep their destination-control plots (they have no rows in cover_long)
community <- readRDS("data/derived/prep/community.rds")
plot_groups <- bind_rows(
  transmute(experiments, experiment_id, region, origin_site, destination_site, treatment = "Warm", group = "warmed plots"),
  transmute(experiments, experiment_id, region, origin_site, destination_site = origin_site, treatment = "LocalControl", group = "origin-control plots"),
  transmute(experiments, experiment_id, region, origin_site = destination_site, destination_site, treatment = "LocalControl", group = "destination-control plots"))
plots <- inner_join(plot_groups, distinct(community, region, origin_site, destination_site, treatment, plot_id),
                    by = c("region", "origin_site", "destination_site", "treatment"), relationship = "many-to-many") |>
  count(experiment_id, group) |> tidyr::pivot_wider(names_from = group, values_from = n)
pools <- read.csv("data/derived/01_relative_cover/pool_sizes.csv")
tableS2 <- experiments |>
  transmute(experiment_id, country, origin_site, destination_site, origin_elevation, destination_elevation, first_year, last_year,
            `destination June–August temperature (°C)` = round(destination_T_JJA, 2), `warming (°C)` = round(diff_T_JJA, 2),
            `destination June–August precipitation (mm)` = round(destination_P_JJA)) |>
  left_join(select(plots, experiment_id, `destination-control plots`, `origin-control plots`, `warmed plots`), by = "experiment_id") |>
  left_join(select(pools, experiment_id, `resident pool` = resident, `destination pool` = destination_pool), by = "experiment_id")
stopifnot(nrow(tableS2) == 42, !anyNA(tableS2)) # every experiment has its three groups of plots and its pools
tableS3 <- tibble::tibble(species = sort(analysed)) |>
  left_join(select(fg, species, family, graminoid, woody), by = "species") |>
  left_join(summarise(res, `gradients as resident` = n_distinct(gradient_id), .by = species), by = "species") |>
  left_join(summarise(col, `gradients as colonizer` = n_distinct(gradient_id), .by = species), by = "species") |>
  left_join(select(s04, species, `classification (resident)` = class), by = "species") |>
  left_join(select(s07b, species, `classification (colonizer)` = class), by = "species") |>
  mutate(LDMC = species %in% traits$species[!is.na(traits$LDMC)], LMA = species %in% traits$species[!is.na(traits$LMA)],
         height = species %in% traits$species[!is.na(traits$height)], PEV = species %in% phy$scores$species, niche = species %in% niche$species)
stopifnot(nrow(tableS3) == length(analysed), !anyNA(tableS3$graminoid))
# step 07c's heterogeneity table lacks the counts, parameters and W slope variance of Table S5;
# they are read from its saved fits so both pools carry the same columns
col_fits <- lapply(list.files("data/derived/07c_ladder_fit/models", full.names = TRUE), readRDS)
col_extra <- tibble::tibble(model = sapply(col_fits, `[[`, "model"), n_parameters = sapply(col_fits, function(f) f$m$parms),
                            tau2_W = sapply(col_fits, function(f) if (isTRUE(f$m$withG)) f$m$tau2 else NA_real_),
                            struct = sapply(col_fits, function(f) if (isTRUE(f$m$withG)) paste(f$m$struct, collapse = ", ") else "none (no slopes)"))
col_het <- read.csv("data/derived/07c_ladder_fit/heterogeneity.csv") |> left_join(col_extra, by = "model") |>
  mutate(rowset = "all", fitted_on = "all", n_species = n_distinct(lad_col$species), n_experiments = n_distinct(lad_col$experiment_id))
stopifnot(nrow(col_het) == 12, !anyNA(col_het$n_parameters)) # every colonizer model found
# step 05 lists the five floor and ceiling fits again under each property's row set; since 2 Nov
# every row set is the ladder rows, so those repeats are dropped here
res_t4 <- read.csv("data/derived/05_ladder_fit/tableS5.csv")
res_t4 <- res_t4[res_t4$rowset == "all" | !res_t4$model %in% c("Null", "Start", "Climate", "Climate + start", "Species-dependent"), ]
stopifnot(nrow(res_t4) == 12, !anyDuplicated(res_t4$model), all(res_t4$n_effects == res_t4$n_effects[1])) # one row per resident model, one row set
tableS5 <- bind_rows(mutate(res_t4, pool = "residents"),
                     mutate(col_het, pool = "colonizers")) |> relocate(pool) |>
  mutate(struct = ifelse(struct == "none (no slopes)", struct, "warming: one variance; temperature, precipitation, their product: diagonal"), # one label for both pools
         model = model_label(model, pool)) # the names of Figs 6 and S17
tableS6 <- bind_rows(mutate(read.csv("data/derived/05_ladder_fit/tableS6.csv"), pool = "residents"),
                     mutate(read.csv("data/derived/07c_ladder_fit/coefficients.csv"), pool = "colonizers")) |> relocate(pool) |>
  mutate(model = model_label(model, pool), term = term_label(term, pool)) # words, as Fig 7
for (t in paste0("S", 1:6)) write.csv(get(paste0("table", t)), sprintf("tables/table%s.csv", t), row.names = FALSE)

# ---- Write the new figures; copy the reviewed ones under their final numbers ------------
ggsave("figures/main/fig4.png", fig4, width = 38, height = 28, dpi = 100, bg = "white")
ggsave("figures/main/fig7.png", fig7, width = 40, height = 26, dpi = 100, bg = "white")
ggsave("figures/appendix/figS1.png", figS1, width = 34, height = 22, dpi = 100, bg = "white")
ggsave("figures/appendix/figS3.png", figS3, width = 40, height = 14, dpi = 100, bg = "white")
ggsave("figures/appendix/figS4.png", figS4, width = 40, height = 26, dpi = 100, bg = "white")
ggsave("figures/appendix/figS10.png", figS10, width = 38, height = 28, dpi = 100, bg = "white")
# the appendix follows the order of the analysis (Bills, 2 Nov); every other appendix figure is
# written under its number by the step that makes it (S2 by prep P2, S5 to S19 by steps 01 to 07e);
# S11, the variance box, is drawn by hand and committed as it is (Bills, 3 Nov)
fig_index <- tibble::tribble(
  ~figure, ~file, ~content, ~from,
  "Fig 1", "figures/main/fig1.png", "Sankey of year-specific classes, residents and colonizers", "steps 03, 07a",
  "Fig 2", "figures/main/fig2.png", "Residents' classification (Step C)", "step 04",
  "Fig 3", "figures/main/fig3.png", "Colonizers' classification (Step C)", "step 07b",
  "Fig 4", "figures/main/fig4.png", "The experiments in climate space: residents' net change and colonizers' gain, modelled effects", "this step",
  "Fig 5", "figures/main/fig5.png", "Traits of losers, winners and colonizers (pathways)", "step 07f",
  "Fig 6", "figures/main/fig6.png", "The residents' ladder", "step 06",
  "Fig 7", "figures/main/fig7.png", "Coefficients of the climate + abundance, Height and LDMC models, residents beside colonizers", "this step",
  "Fig S1", "figures/appendix/figS1.png", "Data flow", "this step",
  "Fig S2", "figures/appendix/figS2.png", "TRY trait screening", "prep P2",
  "Fig S3", "figures/appendix/figS3.png", "Phylogenetic axes (with Table S4)", "this step (prep P3)",
  "Fig S4", "figures/appendix/figS4.png", "Predictors before standardisation, niche checks", "this step (prep P4, P6)",
  "Fig S5", "figures/appendix/figS5.png", "The pool rule", "step 01",
  "Fig S6", "figures/appendix/figS6.png", "Resident pool sizes", "step 01",
  "Fig S7", "figures/appendix/figS7.png", "Losing-winning strength, all years", "step 02",
  "Fig S8", "figures/appendix/figS8.png", "Year-specific classes, residents in at least 3 gradients", "step 03",
  "Fig S9", "figures/appendix/figS9.png", "Strength from year 3 against all years", "step 03b",
  "Fig S10", "figures/appendix/figS10.png", "The experiments in climate space: pies counting losers, winners and colonizers, all effects", "this step",
  "Fig S11", "figures/appendix/figS11.png", "Where each variance component of the ladder comes from (toy table)", "drawn by hand",
  "Fig S12", "figures/appendix/figS12.png", "Coefficients, residents", "step 05",
  "Fig S13", "figures/appendix/figS13.png", "Metrics per gradient, residents", "step 06",
  "Fig S14", "figures/appendix/figS14.png", "Random forest beside metafor, residents", "step 06b",
  "Fig S15", "figures/appendix/figS15.png", "Colonizer strength from year 3 against all years", "step 07a",
  "Fig S16", "figures/appendix/figS16.png", "Coefficients, colonizers", "step 07c",
  "Fig S17", "figures/appendix/figS17.png", "The colonizers' ladder", "step 07d",
  "Fig S18", "figures/appendix/figS18.png", "Metrics per gradient, colonizers", "step 07d",
  "Fig S19", "figures/appendix/figS19.png", "Random forest beside metafor, colonizers", "step 07e")
main_files <- sort(list.files("figures/main", pattern = "^fig[0-9]+\\.png$"))
stopifnot(identical(main_files, paste0("fig", 1:7, ".png"))) # seven main figures, each exactly once
app_files <- list.files("figures/appendix", pattern = "^figS[0-9]+\\.png$")
stopifnot(all(file.exists(fig_index$file)), setequal(app_files, paste0("figS", 1:19, ".png"))) # S1 to S19, each exactly once

# ---- Report --------------------------------------------------------------------------
count_tab <- tibble::tibble(
  item = c("experiments", "resident effects (03b)", "net strength below zero", "net strength above zero", "colonizer effects (07a)",
           "net warmed cover above zero", "experiments with colonizers", "experiments without purple", "species analysed (residents or colonizers)",
           "Fig 4: resident ladder effects", "Fig 4: colonizer ladder effects", "Fig 4: experiments with modelled residents",
           "Fig 4: experiments with a net resident loss", "Fig 4: experiments with a colonizer gain"),
  n = c(nrow(experiments), nrow(res), sum(res$strength < 0), sum(res$strength > 0), nrow(col), sum(col$strength > 0),
        sum(pies$Colonizers > 0), sum(pies$Colonizers == 0), length(analysed), nrow(res_lad), nrow(col_lad), sum(!is.na(balance$net)),
        sum(balance$net < 0, na.rm = TRUE), sum(balance$gain > 0)))
report <- c(
  "# Step 10: Figures, appendix and report", "",
  paste0("Draws what no earlier step produced: Fig 4, the ", nrow(experiments), " experiments in climate space with one disc per experiment for the effects the ladders model ",
         "(", nrow(res_lad), " resident and ", nrow(col_lad), " colonizer effects, steps 05 and 07c): its fill the residents' net change, the sum of their strengths with the sign kept, ",
         "its purple ring the colonizers' total warmed cover, both in percentage points of the warmed turf from year 3 on (Bills, 3 Nov); ",
         "Fig S10 the same map with pies counting resident losers and winners (by the sign of their net strength from year 3 on) and colonizers (net warmed cover above zero), all classified effects; Fig S1 the data flow; Fig S3 and Table S4 the phylogenetic axes; Fig S4 the predictors before ",
         "standardisation with the niche checks; Tables S1, S2, S3, S5 and S6. Copies the reviewed figures of steps 03 to 07f under their final numbers ",
         "(seven main figures, Bills 2 Nov; Fig 7 the coefficients of the climate + abundance, Height and LDMC models, residents beside colonizers). Pies and discs sit at their destination climate; the ", sum(pies$n_shared > 1), " experiments that share a destination climate ",
         "(", n_distinct(paste(pies$x0, pies$y0)[pies$n_shared > 1]), " points) fan out around a black dot at that climate on short lines; no two pies or discs overlap (", sum(pies$moved & pies$n_shared == 1), " pies of different climates had to be pushed apart). ",
         "Each pie is named by its gradient (with a letter where a gradient has several experiments), coloured by warming; the key beside the figure gives each one's origin and destination elevation. ",
         "The ring key is drawn in millimetres and matched to the rings through the drawn panel's scale (", sprintf("%.2f", mm_per_degree), " mm per °C)."),
  "", "## Counts", "", md_table(count_tab),
  "", "## Fig 4 discs and Fig S10 pies", "", "Fig 4: the modelled residents' net change and the colonizers' gain, percentage points of the warmed turf (empty: no modelled residents). Fig S10: species counts.", "",
  md_table(transmute(left_join(pies, balance, by = "experiment_id"),
                     experiment_id, name = label, elevations, `T (°C)` = round(destination_T_JJA, 1), `P (mm)` = round(destination_P_JJA),
                     `W (°C)` = round(diff_T_JJA, 1),
                     `Net change` = ifelse(is.na(net), "", sprintf("%.1f", net)), `Colonizers' gain` = sprintf("%.1f", gain),
                     `Resident losers`, `Resident winners`, Colonizers)),
  "", paste0("Median over the experiments with modelled residents: net change ", sprintf("%.1f", median(balance$net, na.rm = TRUE)),
             " points; median colonizers' gain where colonizers are modelled ", sprintf("%.1f", median(balance$gain[balance$experiment_id %in% col_lad$experiment_id])), " points."),
  "", "## Figures", "", md_table(fig_index), "",
  "![Fig 4](../figures/main/fig4.png)", "", "![Fig 7](../figures/main/fig7.png)", "", "![Fig S1](../figures/appendix/figS1.png)", "", "![Fig S3](../figures/appendix/figS3.png)", "",
  "![Fig S4](../figures/appendix/figS4.png)", "", "![Fig S10](../figures/appendix/figS10.png)",
  "", "## Tables", "",
  md_table(tibble::tibble(table = paste("Table", paste0("S", 1:6)), file = sprintf("tables/table%s.csv", paste0("S", 1:6)),
                          content = c("What each selection removes", "Experiments", "Species list", "Phylogenetic axes", "Variance components per model", "Fixed-effect coefficients per model"),
                          rows = c(nrow(tableS1), nrow(tableS2), nrow(tableS3), nrow(tableS4), nrow(tableS5), nrow(tableS6)))),
  "", "## Table S4, phylogenetic axes", "", md_table(tableS4),
  "", paste0(length(no_pev), " analysed species have no phylogenetic scores", if (length(no_pev)) paste0(": ", paste(no_pev, collapse = ", ")) else "", "."),
  "", "## Warnings captured", "", if (length(warns)) paste("-", unique(warns)) else "None.",
  "", "## Deviations from PLAN.md", "",
  "- One script for the whole step (Bills, 2 Nov).",
  "- Seven main figures (Bills, 2 Nov): the colonizers' classification becomes Fig 3, the climate space Fig 4, the pathways Fig 5, the ladder Fig 6, the coefficients of three models Fig 7 (abundance terms drawn per pool, not comparable); the results-dependent Fig 6 and the old Fig S8 dropped; appendix S1 to S19 in the order of the analysis (Bills, 2 Nov); the main figures are written by their own steps.",
  "- Fig 4 shows the modelled effects as discs, the residents' net change and the colonizers' gain (Bills, 3 Nov; the cover pies of 2 Nov replaced); Fig S10's pies count Fig 5's classes per experiment, not last-year significance classes; no purple in the experiments below the colonizer pool-size rule (Bills, 2 Nov).",
  "- Coefficients, per-gradient metrics and the random forest stay as residents (a) and colonizers (b) files.",
  "- `REPORT.html` not rendered: Quarto is not installed.")
writeLines(report, "reports/10_figures.md")
# the repository README is the realized plan; its numbers are read from data/derived/ (Bills, 2 Nov)
invisible(knitr::knit("README.Rmd", "README.md", quiet = TRUE))
message("Step 10 written: ", nrow(experiments), " pies, ", nrow(fig_index), " figure entries, 6 tables")
