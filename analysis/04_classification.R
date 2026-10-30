# Step 04: consistency across experiments (Step C, residents).
#
# PLAN.md: Step C (model as revised 30 Oct); Fig 2 caption; Appendix Fig S6.
# Reads:  data/derived/03b_strength_y3/strength_df.rds (main, years >= 3);
#         data/derived/02_strength/strength_df.rds (all-years sensitivity)
# Writes: data/derived/04_classification/m_C.rds, species_summary.csv;
#         figures/diagnostics/04_fig2_draft.png, 04_figS6_draft.png;
#         reports/04_classification.md
# Checks: 995 effects; 77 species with at least 4 experiments (552 effects);
#         length(m_C$tau2) equals the eligible count and its species order
#         matches m_C$g.levels.f[[1]] and the mods columns; the prediction
#         interval equals mu +/- 1.96 sqrt(tau2 + se(mu)^2); class counts sum to
#         the eligible count; no interval or whisker on tau anywhere.
#
# Decided 30 Oct (Bills): one fixed pooled mean per species (mods = ~ 0 +
# species) and one among-experiment variance per species (DIAG); the class
# comes from the prediction interval of predict(), which carries tau2_s and the
# uncertainty of mu_s (the species random intercept of the first version had
# its variance estimated at about zero, pulling every mean onto the grand
# mean). Approved with the plan: se_mu in species_summary; the all-years
# strength of step 02 as a sensitivity run; draft figures in
# figures/diagnostics/; Fig 2 corner panels once four species are named.

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
stopifnot(nrow(strength_df) == 995) # step 03b

# ---- Eligible species (PLAN.md Step C) ------------------------------------
eligible <- strength_df |>
  dplyr::count(species) |>
  dplyr::filter(n >= 4) |>
  dplyr::pull(species)
d_c <- dplyr::filter(strength_df, species %in% eligible) |> mutate(species = factor(species))
stopifnot(length(eligible) == 77, nrow(d_c) == 552) # at least 4 experiments

# ---- Model and prediction interval per species ------------------------------
fit_c <- function(d) {
  rma.mv(yi = strength, V = strength_se^2,
         mods   = ~ 0 + species,          # one mu_s per species, fixed
         random = ~ species | effect_id,  # one tau2_s per species
         struct = "DIAG", data = d, method = "REML", sparse = TRUE)
}
summarise_c <- function(m, d) {
  sp <- levels(d$species)
  stopifnot(length(m$tau2) == length(sp), identical(m$g.levels.f[[1]], sp),
            identical(sub("^species", "", rownames(m$beta)), sp)) # tau2 and mu in the same species order
  p <- predict(m, newmods = diag(length(sp)), tau2.levels = sp)
  out <- tibble::tibble(species = sp, n_experiments = as.integer(table(d$species)[sp]),
                        mu = p$pred, se_mu = p$se, tau = sqrt(m$tau2), lower = p$pi.lb, upper = p$pi.ub)
  stopifnot(max(abs(out$upper - (out$mu + qnorm(0.975) * sqrt(out$tau^2 + out$se_mu^2)))) < 1e-10,
            max(abs(out$lower - (out$mu - qnorm(0.975) * sqrt(out$tau^2 + out$se_mu^2)))) < 1e-10)
  mutate(out, class = factor(case_when(lower > 0 ~ "Consistent winner", upper < 0 ~ "Consistent loser",
                                       TRUE ~ "Context-dependent"), levels = class_levels))
}
m_C <- fit_c(d_c)
species_summary <- summarise_c(m_C, d_c)
stopifnot(sum(table(species_summary$class)) == 77) # classes sum to the eligible count

# sensitivity: the same model on step 02's all-years strength, same species
d_all <- readRDS("data/derived/02_strength/strength_df.rds") |>
  filter(species %in% eligible) |> mutate(species = factor(species))
stopifnot(nrow(d_all) == 552)
summary_all <- summarise_c(fit_c(d_all), d_all)

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
ranges <- d_c |> mutate(species = as.character(species)) |>
  summarise(min_s = min(strength), max_s = max(strength), .by = species) |>
  inner_join(species_summary, by = "species") |>
  arrange(mu) |> mutate(species = factor(species, levels = species))
figS6 <- ggplot(ranges, aes(y = species)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_segment(aes(x = min_s, xend = max_s, yend = species, colour = class), linewidth = 0.9) +
  geom_point(aes(x = mu, colour = class), size = 2.6) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(title = "a)", x = "Losing-winning strength across experiments (line) and pooled mean (point)", y = NULL) +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 24), legend.position = "top",
        legend.text = element_text(size = 14, face = "bold"), axis.text.y = element_text(face = "italic", size = 10))
ggsave("figures/diagnostics/04_figS6_draft.png", figS6, width = 10, height = 14, dpi = 150, bg = "white")

# ---- Report --------------------------------------------------------------
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
f4 <- function(x) sprintf("%.4f", x)
cross <- table(`years >= 3 (main)` = species_summary$class,
               `all years (step 02)` = summary_all$class[match(species_summary$species, summary_all$species)])
cross_md <- as.data.frame.matrix(cross); cross_md <- cbind(`years >= 3 (main)` = rownames(cross_md), cross_md)
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
         "(below zero) or context-dependent (spanning zero). The same model on step 02's all-years strength is the ",
         "sensitivity run."),
  "", "## Model as fitted (for reviewers)", "", "```r",
  paste0("m_C <- ", paste(trimws(deparse(m_C$call, width.cutoff = 500)), collapse = " ")),
  "# d: eligible species (>= 4 experiments); mods = ~ 0 + species, random = ~ species | effect_id, struct = \"DIAG\"",
  "p <- predict(m_C, newmods = diag(n_species), tau2.levels = species)",
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
  "", "## Sensitivity: all-years strength (step 02)", "",
  "Rows: class on the strength from year 3 (main); columns: class on the all-years strength.", "", md_table(cross_md),
  "", "## Figures (drafts)", "",
  "![Fig 2 draft](../figures/diagnostics/04_fig2_draft.png)", "",
  "Fig 2 corner panels follow once four species are named (MILESTONES).", "",
  "![Fig S6 draft](../figures/diagnostics/04_figS6_draft.png)",
  "", "## Warnings captured", "", if (length(warns)) paste("-", unique(warns)) else "None.",
  "", "## Deviations from PLAN.md", "",
  paste("None: PLAN.md Step C was revised on 30 Oct (Bills) to fixed species means with predict() intervals. Decisions",
        "approved with the plan: se_mu added to species_summary; the all-years sensitivity run here; draft figures in",
        "figures/diagnostics/ rather than data/derived/04_classification/.")
)
writeLines(report, "reports/04_classification.md")
message("Step 04 written: ", paste(sprintf("%s %d", names(table(species_summary$class)), table(species_summary$class)),
                                   collapse = "; "))
