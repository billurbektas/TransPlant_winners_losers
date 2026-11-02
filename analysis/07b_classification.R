# Step 07b: consistent colonizers (Step C, colonizers).
#
# PLAN.md: The colonizer story, "Consistent colonizers and context-dependent
#   arrival"; Fig S2 (mirrors the anchored Fig 2, Bills 2 Nov).
# Reads:  data/derived/07a_colonizers/strength_df_y3.rds (warmed cover, years >= 3);
#         data/derived/01_relative_cover/cover_long.rds, experiments_kept.csv and
#         data/derived/prep/sites.rds (non-colonizers; Fig S2 b, c)
# Writes: data/derived/07b_classification/m_C.rds, species_summary.csv;
#         figures/diagnostics/07b_figS2.png; reports/07b_classification.md
# Checks: 379 colonizer effects; 29 species in at least 3 gradients (117 effects,
#         24 experiments); 1 consistent non-colonizer set aside, 28 modelled; one
#         tau2_s per modelled species; the prediction interval equals
#         mu +/- 1.96 sqrt(tau2 + se(mu)^2); classes sum to the 29 eligible.
#
# Decided 2 Nov (Bills): strength from years 3 and later, as for residents; the
# residents' Step C model (one fixed mean and one among-experiment variance per
# species, fitted per species, bobyqa where nlminb fails); a species with no warmed
# cover in years >= 3 in any of its experiments is a consistent non-colonizer,
# read from the raw data and set aside (PLAN.md); Fig S2 mirrors the anchored
# resident Fig 2 (code copied from analysis/04_classification.R), in purple.

suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(metafor)})
source("R/md_table.R")

warns <- character()
globalCallingHandlers(warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
})
out_dir <- "data/derived/07b_classification"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
class_levels <- c("Consistent non-colonizer", "Context-dependent", "Consistent colonizer")
cols <- c("Consistent non-colonizer" = "#CFC6DD", "Context-dependent" = "#A68BC8", "Consistent colonizer" = "#5E3C99")

strength_df <- readRDS("data/derived/07a_colonizers/strength_df_y3.rds")
stopifnot(nrow(strength_df) == 379) # step 07a, years >= 3

# ---- Eligible species: at least 3 gradients (2 Nov) -------------------------------
eligible <- summarise(strength_df, n_gradients = n_distinct(gradient_id), .by = species) |>
  filter(n_gradients >= 3) |> pull(species)
d_c <- filter(strength_df, species %in% eligible)
stopifnot(length(eligible) == 29, nrow(d_c) == 117, n_distinct(d_c$experiment_id) == 24) # at least 3 gradients

# ---- Consistent non-colonizers: no warmed cover in years >= 3 anywhere (raw data) ---
cover_long <- readRDS("data/derived/01_relative_cover/cover_long.rds")
warmed_y3 <- cover_long |>
  filter(species_pool == "destination_pool", treatment == "warmed", as.integer(as.character(year)) >= 3,
         paste(experiment_id, species) %in% paste(d_c$experiment_id, d_c$species)) |>
  summarise(arrived = any(relative_cover > 0), .by = species)
non_colonizers <- warmed_y3$species[!warmed_y3$arrived]
stopifnot(nrow(warmed_y3) == 29, identical(non_colonizers, "Ranunculus acris"), # Bills, 2 Nov
          all(d_c$strength[d_c$species %in% non_colonizers] < 1e-12))
d_m <- filter(d_c, !species %in% non_colonizers) |> mutate(species = factor(species))
stopifnot(nlevels(d_m$species) == 28)

# ---- Model and prediction interval per species (as step 04) ------------------------
# one fixed mean and one among-experiment variance per species, fitted per species
# (the same REML likelihood as the joint ~ 0 + species, DIAG model); nlminb can fail
# when tau2_s sits at 0, then bobyqa
fallback <- character()
fit_species <- function(d) {
  m <- try(rma.mv(yi = strength, V = strength_se^2, random = ~ 1 | effect_id, data = d, method = "REML"), silent = TRUE)
  if (inherits(m, "try-error")) {
    fallback <<- c(fallback, as.character(d$species[1]))
    m <- rma.mv(yi = strength, V = strength_se^2, random = ~ 1 | effect_id, data = d, method = "REML",
                control = list(optimizer = "bobyqa"))
  }
  m
}
m_C <- lapply(split(d_m, d_m$species, drop = TRUE), fit_species)
sp <- levels(d_m$species)
tau2 <- sapply(m_C, function(m) m$sigma2)
stopifnot(length(tau2) == 28, identical(names(m_C), sp)) # one tau2_s per modelled species, in order
pred <- lapply(m_C, predict)
modelled <- tibble::tibble(species = sp, n_experiments = as.integer(table(d_m$species)[sp]),
                           n_gradients = as.integer(tapply(d_m$gradient_id, d_m$species, n_distinct)[sp]),
                           mu = sapply(pred, `[[`, "pred"), se_mu = sapply(pred, `[[`, "se"), tau = sqrt(tau2),
                           lower = sapply(pred, `[[`, "pi.lb"), upper = sapply(pred, `[[`, "pi.ub"),
                           n_measured_absences = as.integer(tapply(d_m$se_source == "measured absence", d_m$species, sum)[sp]))
stopifnot(max(abs(modelled$upper - (modelled$mu + qnorm(0.975) * sqrt(modelled$tau^2 + modelled$se_mu^2)))) < 1e-10,
          max(abs(modelled$lower - (modelled$mu - qnorm(0.975) * sqrt(modelled$tau^2 + modelled$se_mu^2)))) < 1e-10)
species_summary <- bind_rows(
  mutate(modelled, class = ifelse(lower > 0, "Consistent colonizer", "Context-dependent")),
  filter(d_c, species %in% non_colonizers) |>
    summarise(n_experiments = n(), n_gradients = n_distinct(gradient_id), n_measured_absences = sum(se_source == "measured absence"), .by = species) |>
    mutate(class = "Consistent non-colonizer")) |>
  mutate(class = factor(class, levels = class_levels))
stopifnot(nrow(species_summary) == 29, sum(table(species_summary$class)) == 29) # classes sum to the eligible species

saveRDS(m_C, file.path(out_dir, "m_C.rds"))
write.csv(mutate(species_summary, class = as.character(class)), file.path(out_dir, "species_summary.csv"), row.names = FALSE)

# ---- Fig S2: mirrors Fig 2 (analysis/04_classification.R:96-183), colonizers --------
# b) c): highest and lowest pooled strength among species in the most gradients
# reached by at least two species
top_g <- sort(unique(modelled$n_gradients), decreasing = TRUE)
top_g <- top_g[sapply(top_g, function(g) sum(modelled$n_gradients >= g) >= 2)][1]
top <- filter(modelled, n_gradients >= top_g)
panel_species <- c(top$species[which.max(top$mu)], top$species[which.min(top$mu)])
# every modelled species is shown: Fig 2's +/- 0.005 rule would hide almost all colonizers
ranges <- filter(species_summary, !is.na(mu)) |> arrange(mu) |> mutate(species = factor(species, levels = species))
x_lim <- range(c(ranges$lower, ranges$upper)) + c(-1, 1) * 0.03 * diff(range(c(ranges$lower, ranges$upper)))
ranges$x_lab <- x_lim[1] - 0.01 * diff(x_lim)
pa <- ggplot(ranges, aes(y = species)) +
  geom_segment(aes(x = x_lim[1], xend = lower, yend = species), colour = cols[as.character(ranges$class)], linewidth = 0.9, alpha = 0.45) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_segment(aes(x = lower, xend = upper, yend = species, colour = n_gradients), linewidth = 4) +
  geom_point(aes(x = mu, fill = class), shape = 21, colour = "black", size = 10, stroke = 0.8) +
  geom_label(aes(x = pmax(upper, 0) + 0.006 * diff(x_lim), label = paste0(n_gradients, " (", n_experiments, ")")), hjust = 0,
             size = 10.4, colour = "grey20", fill = "white", label.size = 0, label.padding = unit(0.05, "lines")) +
  geom_text(aes(x = x_lab, label = species), colour = cols[as.character(ranges$class)], hjust = 1, size = 12.6,
            fontface = ifelse(ranges$species %in% panel_species, "bold.italic", "italic")) +
  coord_cartesian(xlim = x_lim, ylim = c(0.2, nrow(ranges) + 0.8), expand = FALSE, clip = "off") +
  # legend labels carry the count of every class, the non-colonizers (no interval) included
  scale_fill_manual(values = cols, name = NULL, drop = FALSE,
                    labels = paste0(names(cols), " (", as.integer(table(species_summary$class)[names(cols)]), ")")) +
  scale_colour_gradient(name = "Number of gradients", low = "grey85", high = "black") +
  guides(fill = guide_legend(override.aes = list(size = 18), nrow = 1)) +
  labs(title = "a)   Colonizing species", x = "Warmed relative cover", y = NULL) +
  theme_classic(base_size = 32) +
  theme(plot.title = element_text(face = "bold", size = 42), legend.position = "bottom", legend.box = "vertical",
        legend.text = element_text(size = 38, face = "bold"), legend.title = element_text(size = 38, face = "bold"),
        legend.key.width = unit(3.5, "cm"), legend.spacing.x = unit(0.1, "cm"), legend.key.spacing.x = unit(0.6, "cm"),
        legend.justification = "center", legend.location = "plot", axis.title.x = element_text(face = "bold", size = 46),
        axis.text.x = element_text(size = 36), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        plot.margin = margin(15, 40, 15, 565))
sites <- read.csv("data/derived/01_relative_cover/experiments_kept.csv")
site_elev <- readRDS("data/derived/prep/sites.rds") |> select(region, site_id, elevation)
panel_cover <- cover_long |>
  filter(species %in% panel_species, species_pool == "destination_pool", treatment %in% c("destination_control", "warmed"),
         experiment_id %in% strength_df$experiment_id, as.integer(as.character(year)) >= 3) |>
  summarise(cover = 100 * mean(relative_cover), .by = c(species, experiment_id, gradient_id, treatment)) |>
  tidyr::pivot_wider(names_from = treatment, values_from = cover) |>
  left_join(select(sites, experiment_id, origin_site, destination_site), by = "experiment_id") |>
  left_join(transmute(site_elev, gradient_id = region, origin_site = site_id, elev_origin = elevation), by = c("gradient_id", "origin_site")) |>
  left_join(transmute(site_elev, gradient_id = region, destination_site = site_id, elev_dest = elevation), by = c("gradient_id", "destination_site")) |>
  mutate(row = paste0(substr(gradient_id, 1, 2), "  ", elev_origin, " → ", elev_dest, " m"))
corner <- function(sp, tag) {
  d <- filter(panel_cover, species == sp) |> arrange(gradient_id, row) |>
    mutate(row_id = factor(seq_len(n()), levels = rev(seq_len(n()))))
  info <- filter(species_summary, species == sp)
  ggplot(d, aes(y = row_id)) +
    geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.8) +
    geom_segment(aes(x = destination_control, xend = warmed, yend = row_id), colour = "grey40", linewidth = 1.2) + # no arrow: not a change (Bills, 2 Nov)
    geom_point(aes(x = destination_control, shape = "Destination control"), fill = "white", colour = "grey20", size = 13, stroke = 2) +
    geom_point(aes(x = warmed, shape = "Warmed"), colour = "grey50", size = 13, alpha = 0.55) +
    scale_y_discrete(labels = setNames(d$row, d$row_id)) +
    scale_shape_manual(values = c("Destination control" = 21, "Warmed" = 16), name = NULL) +
    guides(shape = guide_legend(override.aes = list(size = 18, stroke = 1.6, colour = c("grey20", "grey50"), fill = c("white", "grey50"), alpha = c(1, 0.55)))) +
    labs(tag = tag, title = sp, x = "Mean relative cover (%)", y = NULL) +
    theme_classic(base_size = 32) +
    theme(plot.tag = element_text(face = "bold", size = 42), plot.tag.position = c(0, 0.985),
          plot.title = element_text(face = "bold.italic", size = 44, colour = cols[[as.character(info$class)]], margin = margin(l = 70, b = 15)),
          plot.title.position = "plot", axis.title.x = element_text(face = "bold", size = 46),
          axis.text.y = element_text(size = 31), axis.text.x = element_text(size = 36),
          legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 42, face = "bold"))
}
bc <- patchwork::wrap_plots(corner(panel_species[1], "b)"), corner(panel_species[2], "c)"), ncol = 1)
fig <- patchwork::wrap_plots(pa, bc, ncol = 2, widths = c(2.2, 1)) +
  patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "vertical",
    legend.text = element_text(size = 42, face = "bold"), legend.title = element_text(size = 42, face = "bold"))
ggsave("figures/diagnostics/07b_figS2.png", fig, width = 44, height = 30, dpi = 100, bg = "white")

# ---- Report ----------------------------------------------------------------------
f4 <- function(x) ifelse(is.na(x), "", sprintf("%.4f", x))
tau0 <- filter(modelled, tau < 1e-4)
report <- c(
  "# Step 07b: Consistent colonizers (Step C, colonizers)", "",
  paste0("For the ", length(eligible), " colonizer species found in at least three gradients (", nrow(d_c), " effects in ",
         n_distinct(d_c$experiment_id), " experiments), on warmed cover from years 3 and later (step 07a): a species with no warmed cover ",
         "in years 3 and later in any of its experiments is a consistent non-colonizer, read from the data and set aside; for the others ",
         "one rma.mv per species estimated a pooled warmed cover and a species-specific among-experiment variance, as the residents' Step C ",
         "(step 04). Each species' prediction interval for a new experiment places it as a consistent colonizer (interval above zero) or ",
         "context-dependent. Measured absences (strength 0, SE 0) enter with their zero SE; their weight is set by the among-experiment variance."),
  "", "## Counts", "",
  md_table(tibble::tibble(rule = c("colonizer effects (step 07a, years 3 and later)", "species in >= 3 gradients", "their effects",
                                   "experiments they occur in", "of their effects, measured absences (SE 0)",
                                   "consistent non-colonizers set aside", "species modelled", "species with tau2_s at about 0 (< 1e-8)"),
                          n = c(nrow(strength_df), length(eligible), nrow(d_c), n_distinct(d_c$experiment_id),
                                sum(d_c$se_source == "measured absence"), length(non_colonizers), nlevels(d_m$species), nrow(tau0)))),
  "", "## Classes", "", md_table(as.data.frame(table(class = species_summary$class), responseName = "species")),
  "", paste0("Consistent non-colonizers: ", paste(non_colonizers, collapse = ", "), "."),
  "", "Every species, by pooled warmed cover:", "",
  md_table(arrange(species_summary, desc(mu)) |> transmute(species, class = as.character(class), n_gradients, n_experiments, n_measured_absences,
                                                           mu = f4(mu), se_mu = f4(se_mu), tau = f4(tau), lower = f4(lower), upper = f4(upper))),
  "", "## Fig S2", "", "![Fig S2](../figures/diagnostics/07b_figS2.png)", "",
  paste0("Fig S2 mirrors Fig 2 (Bills, 2 Nov). a): point = pooled warmed cover, filled by class; line = 95% prediction interval for a new ",
         "experiment, shaded by the number of gradients and labelled gradients (experiments); every modelled species is shown. b) and c): ",
         panel_species[1], " and ", panel_species[2], ", the highest and lowest pooled warmed cover among species in at least ", top_g,
         " gradients; one row per experiment, mean relative cover over years 3 and later in destination controls (hollow) and warmed plots (filled)."),
  "", "## Warnings captured", "",
  if (length(fallback)) paste0("- Refitted with the bobyqa optimizer after nlminb did not declare convergence: ", paste(fallback, collapse = ", "), ".") else character(),
  if (length(warns)) paste("-", unique(warns)) else if (!length(fallback)) "None." else character(),
  "", "## Deviations from PLAN.md", "",
  "- Fig S2 mirrors the anchored Fig 2 (two corner panels by Fig 2's rule, destination control to warmed), not PLAN.md's mock with three panels (Bills, 2 Nov; PLAN.md updated).",
  "- Consistent non-colonizers are read in the years 3 and later window used for strength (Bills, 2 Nov).",
  "- Fig 2's rule hiding context-dependent species with |pooled strength| below 0.005 is not applied: almost every colonizer would be hidden.")
writeLines(report, "reports/07b_classification.md")
message("Step 07b written: ", paste(sprintf("%s %d", names(table(species_summary$class)), table(species_summary$class)), collapse = "; "))
