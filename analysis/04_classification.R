# Step 04: consistency across experiments (Step C, residents).
#
# PLAN.md: Step C (model as revised 30 Oct); Fig 2 (redesigned 1 Nov).
# Reads:  data/derived/03b_strength_y3/strength_df.rds (strength from year 3);
#         data/derived/01_relative_cover/cover_long.rds, experiments_kept.csv and
#         data/derived/prep/sites.rds (Fig 2 b, c)
# Writes: data/derived/04_classification/m_C.rds, species_summary.csv;
#         figures/main/fig2.png;
#         reports/04_classification.md
# Checks: 1,004 effects; 74 species in at least 3 gradients (523 effects);
#         one tau2_s per eligible species, in species order; the prediction
#         interval equals mu +/- 1.96 sqrt(tau2 + se(mu)^2); class counts sum to
#         the eligible count; no interval or whisker on tau anywhere.
#
# Decided 30 Oct (Bills): one fixed pooled mean per species (mods = ~ 0 +
# species) and one among-experiment variance per species (DIAG); the class
# comes from the prediction interval of predict(), which carries tau2_s and the
# uncertainty of mu_s (the species random intercept of the first version had
# its variance estimated at about zero, pulling every mean onto the grand
# mean). Approved with the plan: se_mu in species_summary. The all-years
# sensitivity run was dropped on 30 Oct (Bills: not needed for now). On 1 Nov
# (Bills) the species overview with two species panels became Fig 2, replacing
# the lower-vs-upper scatter with corner panels of the original plan.

suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(metafor)})
source("R/md_table.R")

warns <- character()
globalCallingHandlers(warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
})
out_dir <- "data/derived/04_classification"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
class_levels <- c("Consistent loser", "Context-dependent", "Consistent winner")
cols <- c("Consistent loser" = "#D95F2B", "Context-dependent" = "#3A86D6", "Consistent winner" = "#1B9E77")

strength_df <- readRDS("data/derived/03b_strength_y3/strength_df.rds")
stopifnot(nrow(strength_df) == 1004) # step 03b

# ---- Eligible species (PLAN.md Step C) ------------------------------------
eligible <- strength_df |>
  dplyr::summarise(n_gradients = dplyr::n_distinct(gradient_id), .by = species) |>
  dplyr::filter(n_gradients >= 3) |>   # 2 Nov (Bills): at least 3 gradients
  dplyr::pull(species)
d_c <- dplyr::filter(strength_df, species %in% eligible) |> mutate(species = factor(species))
stopifnot(length(eligible) == 74, nrow(d_c) == 523) # at least 3 gradients

# ---- Model and prediction interval per species ------------------------------
# The model (mods = ~ 0 + species, random = ~ species | effect_id, DIAG) shares
# no parameter between species: its REML likelihood is the sum of one REML
# likelihood per species. It is therefore fitted as one rma.mv per species
# (intercept = mu_s, random ~ 1 | effect_id = tau2_s), with the same optimizer;
# checked 30 Oct against the joint fit (see MILESTONES notes), seconds instead
# of minutes.
# nlminb can fail to declare convergence when tau2_s sits at its boundary of 0;
# those species are refitted with bobyqa and listed in the report
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
m_C <- lapply(split(d_c, d_c$species, drop = TRUE), fit_species)
sp <- levels(d_c$species)
tau2 <- sapply(m_C, function(m) m$sigma2)
stopifnot(length(tau2) == length(eligible), identical(names(m_C), sp)) # one tau2_s per eligible species, in order
pred <- lapply(m_C, predict)
species_summary <- tibble::tibble(species = sp, n_experiments = as.integer(table(d_c$species)[sp]),
                                  n_gradients = as.integer(tapply(d_c$gradient_id, d_c$species, dplyr::n_distinct)[sp]),
                                  mu = sapply(pred, `[[`, "pred"), se_mu = sapply(pred, `[[`, "se"), tau = sqrt(tau2),
                                  lower = sapply(pred, `[[`, "pi.lb"), upper = sapply(pred, `[[`, "pi.ub"))
stopifnot(max(abs(species_summary$upper - (species_summary$mu + qnorm(0.975) * sqrt(species_summary$tau^2 + species_summary$se_mu^2)))) < 1e-10,
          max(abs(species_summary$lower - (species_summary$mu - qnorm(0.975) * sqrt(species_summary$tau^2 + species_summary$se_mu^2)))) < 1e-10)
species_summary <- mutate(species_summary, class = factor(case_when(lower > 0 ~ "Consistent winner", upper < 0 ~ "Consistent loser",
                                                                    TRUE ~ "Context-dependent"), levels = class_levels))
stopifnot(sum(table(species_summary$class)) == 74) # classes sum to the eligible count


# ---- Write -------------------------------------------------------------------
saveRDS(m_C, file.path(out_dir, "m_C.rds"))
utils::write.csv(mutate(species_summary, class = as.character(class)), file.path(out_dir, "species_summary.csv"),
                 row.names = FALSE)

named <- filter(species_summary, class != "Context-dependent")

# ---- Fig 2 (main figure, Bills 1 Nov): a) species overview, b) c) two species --
# one row per species: pooled mean (point) and the model's prediction interval
# for a new experiment (line), the interval that sets the class (Bills, 30 Oct)
# species shown in panels b) and c), chosen here so that a) can mark them:
# highest and lowest pooled strength among species in the most gradients
# reached by at least two species
top_g <- sort(unique(species_summary$n_gradients), decreasing = TRUE)
top_g <- top_g[sapply(top_g, function(g) sum(species_summary$n_gradients >= g) >= 2)][1]
top <- filter(species_summary, n_gradients >= top_g)
panel_species <- c(top$species[which.max(top$mu)], top$species[which.min(top$mu)])
# readable version (Bills, 30 Oct): species with a pooled strength of at least
# half a percentage point of relative cover, plus every consistent species
min_shown <- 0.005
ranges <- species_summary |>
  filter(abs(mu) >= min_shown | class != "Context-dependent") |>
  arrange(mu) |> mutate(species = factor(species, levels = species))
n_hidden <- nrow(species_summary) - nrow(ranges)
x_lim <- range(c(ranges$lower, ranges$upper)) + c(-1, 1) * 0.03 * diff(range(c(ranges$lower, ranges$upper)))
ranges$x_lab <- x_lim[1] - 0.01 * diff(x_lim)
figS6 <- ggplot(ranges, aes(y = species)) +
  # guide from each species name to its interval, in the class colour of the name (Bills, 2 Nov)
  geom_segment(aes(x = x_lim[1], xend = lower, yend = species), colour = cols[as.character(ranges$class)],
               linewidth = 0.9, alpha = 0.45) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_segment(aes(x = lower, xend = upper, yend = species, colour = n_gradients), linewidth = 4) + # thicker (Bills, 2 Nov)
  geom_point(aes(x = mu, fill = class), shape = 21, colour = "black", size = 10, stroke = 0.8) +
  # gradients (experiments), right of each interval's upper end, never left of the zero line
  geom_label(aes(x = pmax(upper, 0) + 0.006 * diff(x_lim), label = paste0(n_gradients, " (", n_experiments, ")")), hjust = 0,
             size = 10.4, colour = "grey20", fill = "white", label.size = 0, label.padding = unit(0.05, "lines")) +
  # species names in their class colour, set directly (one colour scale per plot)
  geom_text(aes(x = x_lab, label = species), colour = cols[as.character(ranges$class)], hjust = 1, size = 12.6,
            fontface = ifelse(ranges$species %in% panel_species, "bold.italic", "italic")) +
  coord_cartesian(xlim = x_lim, ylim = c(0.2, nrow(ranges) + 0.8), expand = FALSE, clip = "off") +
  # legend labels carry the total count of each class (all eligible species, shown or not)
  scale_fill_manual(values = cols, name = NULL,
                    labels = paste0(names(cols), " (", as.integer(table(species_summary$class)[names(cols)]), ")")) +
  scale_colour_gradient(name = "Number of gradients", low = "grey85", high = "black", breaks = c(1, 4, 7, 10, 13)) +
  guides(colour = guide_colourbar(order = 1), fill = guide_legend(override.aes = list(size = 18), nrow = 1, order = 2)) + # fixed legend order: ggplot2 otherwise swaps them between runs
  labs(title = "a)   Resident species", x = "Losing-winning strength", y = NULL) +
  theme_classic(base_size = 32) +
  theme(plot.title = element_text(face = "bold", size = 42), legend.position = "bottom", legend.box = "vertical",
        legend.text = element_text(size = 38, face = "bold"), legend.title = element_text(size = 38, face = "bold"),
        legend.key.width = unit(3.5, "cm"), legend.spacing.x = unit(0.1, "cm"), legend.key.spacing.x = unit(0.6, "cm"), legend.justification = "center", legend.location = "plot", axis.title.x = element_text(face = "bold", size = 46), axis.text.x = element_text(size = 36),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        plot.margin = margin(15, 40, 15, 565))
# ---- Fig 2 panels b) and c): best winner and worst loser among the most widely sampled species
# (Bills, 1 Nov): among species in the largest number of gradients reached by at least two
# species (here 12), the highest and the lowest pooled strength; one row per experiment,
# grouped by gradient, grey arrow from origin-control cover (hollow) to warmed cover (filled; no class colours, Bills),
# mean relative cover over years >= 3 as the strength
sites <- utils::read.csv("data/derived/01_relative_cover/experiments_kept.csv")
site_elev <- readRDS("data/derived/prep/sites.rds") |> select(region, site_id, elevation)
panel_cover <- readRDS("data/derived/01_relative_cover/cover_long.rds") |>
  filter(species %in% panel_species, species_pool == "resident", treatment %in% c("origin_control", "warmed"),
         as.integer(as.character(year)) >= 3) |>
  summarise(cover = 100 * mean(relative_cover), .by = c(species, experiment_id, gradient_id, treatment)) |>
  tidyr::pivot_wider(names_from = treatment, values_from = cover) |>
  left_join(select(sites, experiment_id, origin_site, destination_site), by = "experiment_id") |>
  left_join(transmute(site_elev, gradient_id = region, origin_site = site_id, elev_origin = elevation),
            by = c("gradient_id", "origin_site")) |>
  left_join(transmute(site_elev, gradient_id = region, destination_site = site_id, elev_dest = elevation),
            by = c("gradient_id", "destination_site")) |>
  mutate(row = paste0(substr(gradient_id, 1, 2), "  ", elev_origin, " \u2192 ", elev_dest, " m"),
         direction = ifelse(warmed >= origin_control, "Cover gained under warming", "Cover lost under warming"))
corner <- function(sp, tag) {
  d <- filter(panel_cover, species == sp) |> arrange(gradient_id, row) |>
    mutate(row_id = factor(seq_len(n()), levels = rev(seq_len(n()))))  # one row per experiment, labels may repeat
  info <- filter(species_summary, species == sp)
  ggplot(d, aes(y = row_id)) +
    geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.8) +   # zero cover
    geom_segment(aes(x = origin_control, xend = warmed, yend = row_id), colour = "black", linewidth = 1.6,
                 arrow = grid::arrow(length = grid::unit(0.35, "cm"), type = "closed")) +
    geom_point(aes(x = origin_control, shape = "Origin control"), fill = "white", colour = "grey20", size = 13, stroke = 2) +
    geom_point(aes(x = warmed, shape = "Warmed"), colour = "grey50", size = 13, alpha = 0.55) +
    scale_y_discrete(labels = setNames(d$row, d$row_id)) +
    scale_shape_manual(values = c("Origin control" = 21, "Warmed" = 16), name = NULL) +
    guides(shape = guide_legend(override.aes = list(size = 18, stroke = 1.6, colour = c("grey20", "grey50"), fill = c("white", "grey50"), alpha = c(1, 0.55)))) +
    labs(tag = tag, title = sp, x = "Mean relative cover (%)", y = NULL) +
    theme_classic(base_size = 32) +
    theme(plot.tag = element_text(face = "bold", size = 42), plot.tag.position = c(0, 0.985),
          plot.title = element_text(face = "bold.italic", size = 44, colour = cols[[as.character(info$class)]],
                                    margin = margin(l = 70, b = 15)), plot.title.position = "plot",
          axis.title.x = element_text(face = "bold", size = 46),
          axis.text.y = element_text(size = 31), axis.text.x = element_text(size = 36),
          legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 42, face = "bold"))
}
bc <- patchwork::wrap_plots(corner(panel_species[1], "b)"), corner(panel_species[2], "c)"), ncol = 1)
fig_s6 <- patchwork::wrap_plots(figS6, bc,
                                ncol = 2, widths = c(2.2, 1)) +
  # all legends collected in one row under the whole figure
  patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "vertical",
    legend.text = element_text(size = 42, face = "bold"), legend.title = element_text(size = 42, face = "bold")) # one legend size (Bills, 2 Nov)
dir.create("figures/main", recursive = TRUE, showWarnings = FALSE)
ggsave("figures/main/fig2.png", fig_s6, width = 44, height = 40, dpi = 100, bg = "white")

# ---- Report --------------------------------------------------------------
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
f4 <- function(x) sprintf("%.4f", x)
consistent <- named |> arrange(class, mu) |>
  transmute(species, class = as.character(class), n_experiments, mu = f4(mu), se_mu = f4(se_mu), tau = f4(tau),
            lower = f4(lower), upper = f4(upper))
tau0 <- filter(species_summary, tau < 1e-4)
report <- c(
  "# Step 04: Consistency across experiments (Step C, residents)", "",
  paste0("For the ", length(eligible), " resident species found in at least three gradients (", fmt(nrow(d_c)),
         " effects in ", n_distinct(d_c$experiment_id), " experiments), one `rma.mv` estimated a pooled losing-winning ",
         "strength per species and a species-specific among-experiment variance, on the strength from year 3 (step 03b). ",
         "Each species' prediction interval for a new experiment, which carries both its among-experiment variance and the ",
         "uncertainty of its pooled strength, places it as a consistent winner (interval above zero), a consistent loser ",
         "(below zero) or context-dependent (spanning zero)."),
  "", "## Model as fitted (for reviewers)", "", "```r",
  "# the model: one fixed mean and one among-experiment variance per species",
  "# rma.mv(strength, strength_se^2, mods = ~ 0 + species, random = ~ species | effect_id, struct = \"DIAG\")",
  "# fitted as one rma.mv per species, which maximises the same REML likelihood (no shared parameter):",
  paste0("m_s <- ", paste(trimws(deparse(m_C[[1]]$call, width.cutoff = 500)), collapse = " "), "   # d = one species' effects"),
  "p <- predict(m_s)",
  "# lower, upper = p$pi.lb, p$pi.ub = mu_s -/+ 1.96 * sqrt(tau2_s + se(mu_s)^2)",
  "```", "",
  "## Counts", "",
  md_table(tibble::tibble(rule = c("effects (step 03b)", "species in >= 3 gradients", "their effects",
                                   "experiments they occur in", "species with tau2_s estimated at about 0 (< 1e-8)"),
                          n = fmt(c(nrow(strength_df), length(eligible), nrow(d_c), n_distinct(d_c$experiment_id),
                                    nrow(tau0))))),
  "", "## Classes", "", md_table(as.data.frame(table(class = species_summary$class), responseName = "species")),
  "", "Consistent winners and losers:", "", md_table(consistent),
  "", paste0("Species with tau2_s at about zero (their interval is the uncertainty of the mean alone): ",
             paste(tau0$species, collapse = ", "), "."),
  "", "## Fig 2", "",
  "![Fig 2](../figures/main/fig2.png)", "",
  paste0("Fig 2 b) and c): ", panel_species[1], " and ", panel_species[2], ", the highest and lowest pooled strength among ",
         "species in at least ", top_g, " gradients; one row per experiment (country, origin to destination elevation), mean relative cover over years 3 and later in ",
         "origin controls (hollow) and warmed plots (filled)."),
  paste0("Fig 2 a): point = pooled strength, filled by class; line = 95% interval including the variation among experiments ",
         "(a line that does not cross zero marks a consistent loser or winner), shaded by the number of gradients and labelled gradients (experiments). ",
         n_hidden, " context-dependent species with a pooled strength between -", min_shown, " and +", min_shown,
         " are not shown."),
  "", "## Warnings captured", "",
  if (length(fallback)) paste0("- Refitted with the bobyqa optimizer after nlminb did not declare convergence (tau2_s at 0): ",
                               paste(fallback, collapse = ", "), ".") else character(),
  if (length(warns)) paste("-", unique(warns)) else if (!length(fallback)) "None." else character(),
  "", "## Deviations from PLAN.md", "",
  paste("None: PLAN.md Step C was revised on 30 Oct (Bills) to fixed species means with predict() intervals, and Fig 2",
        "was redesigned on 1 Nov (Bills: species overview with two species panels; PLAN.md updated). Decision approved",
        "with the plan: se_mu added to species_summary.")
)
writeLines(report, "reports/04_classification.md")
message("Step 04 written: ", paste(sprintf("%s %d", names(table(species_summary$class)), table(species_summary$class)),
                                   collapse = "; "))
