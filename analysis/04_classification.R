# Step 04: consistency across experiments (Step C, residents).
#
# PLAN.md: Step C (model as revised 30 Oct); Fig 2 caption; Appendix Fig S6.
# Reads:  data/derived/03b_strength_y3/strength_df.rds (strength from year 3)
# Writes: data/derived/04_classification/m_C.rds, species_summary.csv;
#         figures/diagnostics/04_fig2_draft.png, 04_figS6_draft.png;
#         reports/04_classification.md
# Checks: 1,655 effects; 137 species with at least 4 experiments (998 effects);
#         one tau2_s per eligible species, in species order; the prediction
#         interval equals mu +/- 1.96 sqrt(tau2 + se(mu)^2); class counts sum to
#         the eligible count; no interval or whisker on tau anywhere.
#
# Decided 30 Oct (Bills): one fixed pooled mean per species (mods = ~ 0 +
# species) and one among-experiment variance per species (DIAG); the class
# comes from the prediction interval of predict(), which carries tau2_s and the
# uncertainty of mu_s (the species random intercept of the first version had
# its variance estimated at about zero, pulling every mean onto the grand
# mean). Approved with the plan: se_mu in species_summary; draft figures in
# figures/diagnostics/; Fig 2 corner panels once four species are named. The
# all-years sensitivity run was dropped on 30 Oct (Bills: not needed for now).

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
stopifnot(nrow(strength_df) == 1655) # step 03b

# ---- Eligible species (PLAN.md Step C) ------------------------------------
eligible <- strength_df |>
  dplyr::count(species) |>
  dplyr::filter(n >= 4) |>
  dplyr::pull(species)
d_c <- dplyr::filter(strength_df, species %in% eligible) |> mutate(species = factor(species))
stopifnot(length(eligible) == 137, nrow(d_c) == 998) # at least 4 experiments

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
                                  mu = sapply(pred, `[[`, "pred"), se_mu = sapply(pred, `[[`, "se"), tau = sqrt(tau2),
                                  lower = sapply(pred, `[[`, "pi.lb"), upper = sapply(pred, `[[`, "pi.ub"))
stopifnot(max(abs(species_summary$upper - (species_summary$mu + qnorm(0.975) * sqrt(species_summary$tau^2 + species_summary$se_mu^2)))) < 1e-10,
          max(abs(species_summary$lower - (species_summary$mu - qnorm(0.975) * sqrt(species_summary$tau^2 + species_summary$se_mu^2)))) < 1e-10)
species_summary <- mutate(species_summary, class = factor(case_when(lower > 0 ~ "Consistent winner", upper < 0 ~ "Consistent loser",
                                                                    TRUE ~ "Context-dependent"), levels = class_levels))
stopifnot(sum(table(species_summary$class)) == 137) # classes sum to the eligible count


# ---- Write -------------------------------------------------------------------
saveRDS(m_C, file.path(out_dir, "m_C.rds"))
utils::write.csv(mutate(species_summary, class = as.character(class)), file.path(out_dir, "species_summary.csv"),
                 row.names = FALSE)

# ---- Fig 2 draft (corner panels once four species are named) --------------
lim_x <- range(species_summary$lower); lim_y <- range(species_summary$upper)
named <- filter(species_summary, class != "Context-dependent")
fig2 <- ggplot(species_summary, aes(lower, upper)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = cols[["Context-dependent"]], alpha = 0.07) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = cols[["Consistent winner"]], alpha = 0.09) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = cols[["Consistent loser"]], alpha = 0.09) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
  geom_abline(linetype = "dashed", colour = "grey70") +
  geom_point(aes(colour = class), size = 3, alpha = 0.85) +
  ggrepel::geom_text_repel(data = named, aes(label = species, colour = class), size = 4.6, fontface = "bold.italic",
                           min.segment.length = 0, seed = 1, max.overlaps = Inf, show.legend = FALSE) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(title = "a)", x = "Lower bound of the prediction interval",
       y = "Upper bound of the prediction interval") +
  theme_classic(base_size = 17) +
  theme(plot.title = element_text(face = "bold", size = 26), legend.position = "top",
        legend.text = element_text(size = 15, face = "bold"))
ggsave("figures/diagnostics/04_fig2_draft.png", fig2, width = 10, height = 9, dpi = 150, bg = "white")

# ---- Fig S6 draft: species overview -------------------------------------
# one row per species: pooled mean (point) and the model's prediction interval
# for a new experiment (line), the interval that sets the class (Bills, 30 Oct)
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
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_segment(aes(x = lower, xend = upper, yend = species, colour = class), linewidth = 0.9) +
  geom_point(aes(x = mu, fill = n_experiments), shape = 21, colour = "black", size = 3.4, stroke = 0.5) +
  # species names coloured by class, so classes stay visible where the line is shorter than the point
  geom_text(aes(x = x_lab, label = species, colour = class), hjust = 1, size = 3.8, fontface = "italic",
            show.legend = FALSE) +
  coord_cartesian(xlim = x_lim, ylim = c(0.2, nrow(ranges) + 0.8), expand = FALSE, clip = "off") +
  scale_colour_manual(values = cols, name = NULL) +
  scale_fill_gradient(name = "Experiments", low = "white", high = "black", breaks = c(4, 8, 12, 16, 20)) +
  labs(title = "a)",
       subtitle = paste0("Point: pooled strength across experiments.\nLine: 95% interval including the variation among ",
                         "experiments;\na line that does not cross zero marks a consistent loser or winner.\n",
                         n_hidden, " context-dependent species with a pooled strength between -", min_shown, " and +",
                         min_shown, " are not shown."),
       x = "Losing-winning strength (warmed minus origin-control relative cover)", y = NULL) +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 24), legend.position = "top", legend.box = "vertical",
        legend.text = element_text(size = 14, face = "bold"), legend.title = element_text(size = 14, face = "bold"),
        legend.key.width = unit(1.6, "cm"), plot.subtitle = element_text(size = 13, colour = "grey25"), axis.title.x = element_text(face = "bold", size = 16),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        plot.margin = margin(10, 10, 10, 150))
ggsave("figures/diagnostics/04_figS6_draft.png", figS6, width = 12, height = 14, dpi = 150, bg = "white")

# ---- Report --------------------------------------------------------------
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
f4 <- function(x) sprintf("%.4f", x)
consistent <- named |> arrange(class, mu) |>
  transmute(species, class = as.character(class), n_experiments, mu = f4(mu), se_mu = f4(se_mu), tau = f4(tau),
            lower = f4(lower), upper = f4(upper))
tau0 <- filter(species_summary, tau < 1e-4)
report <- c(
  "# Step 04: Consistency across experiments (Step C, residents)", "",
  paste0("For the ", length(eligible), " resident species found in at least four experiments (", fmt(nrow(d_c)),
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
  md_table(tibble::tibble(rule = c("effects (step 03b)", "species with >= 4 experiments", "their effects",
                                   "experiments they occur in", "species with tau2_s estimated at about 0 (< 1e-8)"),
                          n = fmt(c(nrow(strength_df), length(eligible), nrow(d_c), n_distinct(d_c$experiment_id),
                                    nrow(tau0))))),
  "", "## Classes", "", md_table(as.data.frame(table(class = species_summary$class), responseName = "species")),
  "", "Consistent winners and losers:", "", md_table(consistent),
  "", paste0("Species with tau2_s at about zero (their interval is the uncertainty of the mean alone): ",
             paste(tau0$species, collapse = ", "), "."),
  "", "## Figures (drafts)", "",
  "![Fig 2 draft](../figures/diagnostics/04_fig2_draft.png)", "",
  "Fig 2 corner panels follow once four species are named (MILESTONES).", "",
  "![Fig S6 draft](../figures/diagnostics/04_figS6_draft.png)",
  "", "## Warnings captured", "",
  if (length(fallback)) paste0("- Refitted with the bobyqa optimizer after nlminb did not declare convergence (tau2_s at 0): ",
                               paste(fallback, collapse = ", "), ".") else character(),
  if (length(warns)) paste("-", unique(warns)) else if (!length(fallback)) "None." else character(),
  "", "## Deviations from PLAN.md", "",
  paste("None: PLAN.md Step C was revised on 30 Oct (Bills) to fixed species means with predict() intervals. Decisions",
        "approved with the plan: se_mu added to species_summary; draft figures in figures/diagnostics/ rather than",
        "data/derived/04_classification/.")
)
writeLines(report, "reports/04_classification.md")
message("Step 04 written: ", paste(sprintf("%s %d", names(table(species_summary$class)), table(species_summary$class)),
                                   collapse = "; "))
