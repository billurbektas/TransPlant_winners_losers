# Step 02: strength per species x experiment (Step A, residents).
#
# PLAN.md: Step A (model, covariates, data rules); Step D (A0, property
#   variables and their standardisation, niche n_obs_sub >= 30).
# Reads:  data/derived/01_relative_cover/cover_long.rds, a0.rds;
#         data/derived/prep/experiment_climate.rds, traits.rds,
#         functional_groups.rds, phylogeny_axes.rds, niche_optima.rds
# Writes: data/derived/02_strength/models/<experiment_id>.rds (one fit per
#         experiment, `::` written `__`), by_year.rds, strength_df.rds,
#         strength_df.csv; figures/diagnostics/02_strength.png;
#         reports/02_strength.md
# Checks: 108,884 resident rows in origin controls and warmed plots, 43
#         experiments, 1,655 species x experiment combinations, 528 species;
#         by_year has one row per species x survey year; strength_df one row
#         per species x experiment; every strength_se finite and positive; W,
#         T, P standardised across the 43 experiments, A0 across the 1,655
#         effects; every join keeps the 1,655 rows.
#
# Decided 30 Oct (Bills): one residual variance per species (glmmTMB,
# dispformula = ~ species, in R/fit_experiment.R), so each species gets its own
# strength_se; degrees of freedom = origin-control plots + warmed plots - 2
# per experiment; the plot x species variance stays one per experiment;
# PEV1 to PEV3 only; species properties standardised across the 528 residents
# that have each one (as legacy/R/models/predictors.R:170-184), traits logged
# first (same lines); graminoid and woody as 0/1, not standardised.

suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
source("R/fit_experiment.R")
source("R/md_table.R")

warns <- character()
globalCallingHandlers(warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
})
out_dir <- "data/derived/02_strength"
model_dir <- file.path(out_dir, "models")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("figures/diagnostics", recursive = TRUE, showWarnings = FALSE)
z <- function(x) (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)

# ---- Resident rows ---------------------------------------------------------
cover_long <- readRDS("data/derived/01_relative_cover/cover_long.rds")
resident <- cover_long |>
  filter(species_pool == "resident", treatment %in% c("origin_control", "warmed")) |>
  mutate(treatment = factor(treatment, levels = c("origin_control", "warmed")), species = factor(species))
combos <- distinct(resident, experiment_id, gradient_id, species)
stopifnot(nrow(resident) == 108884, n_distinct(resident$experiment_id) == 43) # residents, step 01
stopifnot(nrow(combos) == 1655, n_distinct(combos$species) == 528)

# ---- One fit per experiment, cached ---------------------------------------
experiment_ids <- sort(unique(resident$experiment_id))
fits <- list()
for (e in experiment_ids) {
  f <- file.path(model_dir, paste0(gsub("::", "__", e), ".rds"))
  if (!file.exists(f)) {
    message("fitting ", e)
    saveRDS(fit_experiment(filter(resident, experiment_id == e), control_level = "origin_control"), f)
  }
  fits[[e]] <- readRDS(f)
}

by_year <- bind_rows(lapply(experiment_ids, function(e) mutate(fits[[e]]$by_year, experiment_id = e))) |>
  transmute(species = as.character(species), experiment_id, year, estimate, SE, p.value)
expected_by_year <- resident |> distinct(experiment_id, species, year) |> nrow()
stopifnot(nrow(by_year) == expected_by_year) # one contrast per species x survey year
stopifnot(all(is.finite(by_year$estimate)), all(is.finite(by_year$SE)))

strength <- bind_rows(lapply(experiment_ids, function(e) mutate(fits[[e]]$strength, experiment_id = e))) |>
  transmute(species = as.character(species), experiment_id, strength = estimate, strength_se = SE)
stopifnot(nrow(strength) == 1655, !anyDuplicated(strength[c("species", "experiment_id")])) # one per species x experiment
stopifnot(all(is.finite(strength$strength_se)), all(strength$strength_se > 0))

# ---- Climate covariates, standardised across experiments -------------------
climate <- readRDS("data/derived/prep/experiment_climate.rds") |>
  filter(experiment_id %in% experiment_ids) |>
  transmute(experiment_id, W_raw = W, T_raw = T, P_raw = P_summer, W = z(W), T = z(T), P = z(P_summer), TP = T * P)
stopifnot(nrow(climate) == 43) # every experiment has climate
stopifnot(all(abs(colMeans(climate[c("W", "T", "P")])) < 1e-8), all(abs(sapply(climate[c("W", "T", "P")], stats::sd) - 1) < 1e-8))

# ---- A0, standardised across effects -------------------------------------
a0 <- readRDS("data/derived/01_relative_cover/a0.rds") |>
  semi_join(strength, by = c("experiment_id", "species"))
stopifnot(nrow(a0) == 1655) # every effect has both A0 values
a0 <- mutate(a0, A0_origin_raw = A0_origin, A0_dest_raw = A0_dest, A0_origin = z(A0_origin), A0_dest = z(A0_dest))

# ---- Species properties, standardised across the residents that have them --
species_res <- sort(unique(strength$species))
traits <- readRDS("data/derived/prep/traits.rds") |>
  filter(species %in% species_res) |>
  transmute(species, ldmc_z = z(log(LDMC)), lma_z = z(log(LMA)), height_z = z(log(height)))
groups <- readRDS("data/derived/prep/functional_groups.rds") |>
  filter(species %in% species_res) |>
  transmute(species, graminoid = as.integer(graminoid), woody = as.integer(woody))
pev <- readRDS("data/derived/prep/phylogeny_axes.rds")$scores |>
  filter(species %in% species_res) |>
  transmute(species, PEV1 = z(PEV1), PEV2 = z(PEV2), PEV3 = z(PEV3))
niche_all <- readRDS("data/derived/prep/niche_optima.rds") |> filter(species %in% species_res)
niche <- niche_all |>
  filter(n_obs_sub >= 30) |> # PLAN.md Step D, niche paragraph
  transmute(species, T_opt = z(T_opt), P_opt = z(P_opt))
stopifnot(nrow(groups) == 528, nrow(pev) == 528) # growth form and phylogeny for every resident

strength_df <- strength |>
  left_join(distinct(combos, experiment_id, gradient_id) |> distinct(), by = "experiment_id") |>
  left_join(select(climate, experiment_id, W, T, P, TP), by = "experiment_id") |>
  left_join(select(a0, experiment_id, species, A0_origin, A0_dest), by = c("experiment_id", "species")) |>
  left_join(traits, by = "species") |>
  left_join(groups, by = "species") |>
  left_join(pev, by = "species") |>
  left_join(niche, by = "species") |>
  mutate(effect_id = paste(experiment_id, species, sep = "::")) |>
  select(effect_id, species, experiment_id, gradient_id, strength, strength_se, W, T, P, TP, A0_origin, A0_dest,
         ldmc_z, lma_z, height_z, graminoid, woody, PEV1, PEV2, PEV3, T_opt, P_opt) |>
  arrange(experiment_id, species)
stopifnot(nrow(strength_df) == 1655, !anyDuplicated(strength_df$effect_id)) # joins kept every effect
stopifnot(!anyNA(strength_df[c("gradient_id", "W", "T", "P", "TP", "A0_origin", "A0_dest", "graminoid", "woody",
                               "PEV1", "PEV2", "PEV3")]))
stopifnot(n_distinct(strength_df$gradient_id) == 21)

# ---- Write -------------------------------------------------------------------
saveRDS(by_year, file.path(out_dir, "by_year.rds"))
saveRDS(strength_df, file.path(out_dir, "strength_df.rds"))
utils::write.csv(strength_df, file.path(out_dir, "strength_df.csv"), row.names = FALSE)

long <- tidyr::pivot_longer(select(strength_df, strength, strength_se), everything())
p1 <- ggplot(long, aes(value)) +
  geom_histogram(bins = 60, fill = "grey60", colour = "white") +
  facet_wrap(~ name, scales = "free", labeller = as_labeller(c(strength = "strength (warmed - origin control)",
                                                               strength_se = "strength_se"))) +
  labs(x = "Relative cover (proportion)", y = "Species x experiment") +
  theme_classic()
p2 <- ggplot(strength_df, aes(strength_se, strength)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_point(alpha = 0.4, size = 1) +
  labs(x = "strength_se", y = "strength") +
  theme_classic()
ggsave("figures/diagnostics/02_strength.png", patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1, 1.2)),
       width = 7, height = 7, dpi = 150)

# ---- Report --------------------------------------------------------------
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
fit_table <- tibble::tibble(
  experiment_id = experiment_ids,
  rows = as.integer(table(resident$experiment_id)[experiment_ids]),
  species = sapply(experiment_ids, function(e) nlevels(droplevels(filter(resident, experiment_id == e)$species))),
  years = sapply(experiment_ids, function(e) n_distinct(filter(resident, experiment_id == e)$year)),
  df = sapply(fits, `[[`, "df_plots"),
  convergence = sapply(fits, `[[`, "convergence"),
  pdHess = sapply(fits, `[[`, "pdHess"),
  sd_plot_species = signif(sapply(fits, `[[`, "sd_plot_species"), 3),
  se_min = signif(sapply(experiment_ids, function(e) min(strength$strength_se[strength$experiment_id == e])), 3),
  se_median = signif(sapply(experiment_ids, function(e) stats::median(strength$strength_se[strength$experiment_id == e])), 3),
  se_max = signif(sapply(experiment_ids, function(e) max(strength$strength_se[strength$experiment_id == e])), 3),
  warnings = sapply(fits, function(x) length(x$warnings)))
mean_cover <- summarise(resident, mean_cover = mean(relative_cover), .by = c(experiment_id, species)) |>
  mutate(species = as.character(species))
se_cover <- inner_join(strength, mean_cover, by = c("experiment_id", "species"))

# significance at each experiment's own t cut-off (df = control + warmed plots - 2)
df_exp <- sapply(fits, `[[`, "df_plots")
stopifnot(all(df_exp >= 4))
crit <- stats::qt(0.975, df_exp[strength_df$experiment_id])
n_sig <- sum(abs(strength_df$strength) > crit * strength_df$strength_se)

# share of each effect's SE variance that comes from the shared plot x species
# variance (balanced form: sd_u^2 / (sd_u^2 + sigma_s^2 / n_years))
plot_share <- bind_rows(lapply(experiment_ids, function(e) {
  f <- fits[[e]]
  n_years <- n_distinct(filter(resident, experiment_id == e)$year)
  tibble::tibble(experiment_id = e, species = names(f$sigma_species), sigma = f$sigma_species,
                 share = f$sd_plot_species^2 / (f$sd_plot_species^2 + f$sigma_species^2 / n_years))
}))
# a species' own plot-to-plot spread: SD over plots of its mean deviation from
# the treatment x year cell mean
plot_dev <- resident |>
  mutate(dev = relative_cover - mean(relative_cover), .by = c(experiment_id, species, treatment, year)) |>
  summarise(dev = mean(dev), .by = c(experiment_id, species, plot_id)) |>
  summarise(plot_sd = stats::sd(dev), .by = c(experiment_id, species)) |>
  mutate(species = as.character(species)) |>
  inner_join(plot_share, by = c("experiment_id", "species"))
fit_warnings <- unlist(lapply(experiment_ids, function(e) if (length(fits[[e]]$warnings)) paste0(e, ": ", fits[[e]]$warnings)))
q <- function(x) paste(sprintf("%.4f", stats::quantile(x, c(0, 0.25, 0.5, 0.75, 1))), collapse = " / ")
has <- function(d, col) d$species[!is.na(d[[col]])]
trait_raw <- readRDS("data/derived/prep/traits.rds") |> filter(species %in% species_res)
coverage <- tibble::tibble(
  property = c("LDMC", "LMA", "height", "any of the three traits", "growth form (graminoid, woody)", "PEV1 to PEV3",
               "niche optima (any GBIF cell)", "niche optima (n_obs_sub >= 30, used)"),
  species = c(length(has(traits, "ldmc_z")), length(has(traits, "lma_z")), length(has(traits, "height_z")),
              n_distinct(trait_raw$species[!is.na(trait_raw$LDMC) | !is.na(trait_raw$LMA) | !is.na(trait_raw$height)]),
              nrow(groups), nrow(pev), nrow(niche_all), nrow(niche)),
  effects = c(sum(!is.na(strength_df$ldmc_z)), sum(!is.na(strength_df$lma_z)), sum(!is.na(strength_df$height_z)),
              sum(!is.na(strength_df$ldmc_z) | !is.na(strength_df$lma_z) | !is.na(strength_df$height_z)),
              sum(!is.na(strength_df$graminoid)), sum(!is.na(strength_df$PEV1)),
              sum(strength_df$species %in% niche_all$species), sum(!is.na(strength_df$T_opt))))
no_trait <- setdiff(species_res, trait_raw$species[!is.na(trait_raw$LDMC) | !is.na(trait_raw$LMA) | !is.na(trait_raw$height)])
counts <- tibble::tibble(
  rule = c("resident rows, origin controls and warmed plots (step 01)", "experiments fitted", "gradients",
           "species", "species x experiment (effects)", "by_year contrasts (species x survey year)", "strength_df rows"),
  n = fmt(c(nrow(resident), length(experiment_ids), n_distinct(strength_df$gradient_id), length(species_res),
            nrow(strength), nrow(by_year), nrow(strength_df))))

report <- c(
  "# Step 02: Strength (Step A, residents)", "",
  paste0("One Gaussian `glmmTMB(relative_cover ~ treatment * year * species + (1 | plot_id:species), ",
         "dispformula = ~ species, REML = TRUE)` (one residual variance per species) was fitted per experiment on ",
         "the resident species' rows in origin controls and warmed plots (", fmt(nrow(resident)), " rows, ",
         length(experiment_ids), " experiments), year as a factor. `emmeans` gave the warmed minus origin-control contrast ",
         "per species and survey year (`by_year`, for step 03) and averaged over each experiment's years (`strength`), ",
         "with degrees of freedom = origin-control plots + warmed plots - 2 of the experiment. Each effect carries the experiment's W, T and P ",
         "(summer, log) standardised across the ", nrow(climate), " experiments, TP = T x P, A0 at origin and destination ",
         "standardised across the ", fmt(nrow(strength_df)), " effects, and the species properties standardised across ",
         "the residents that have them. Strength is in relative-cover units (proportion of plot cover)."),
  "", "## Model as fitted (for reviewers)", "",
  "The call, as stored in the fitted model objects, and the contrasts (`R/fit_experiment.R`), run once per experiment:", "",
  "```r",
  paste0("m <- ", paste(trimws(deparse(fits[[1]]$model$call, width.cutoff = 500)), collapse = " ")),
  "# d: one experiment's resident rows in origin controls and warmed plots;",
  "# treatment (origin_control, warmed), year (factor), species (factor), plot_id",
  "df_plots <- n_distinct(plot_id[origin_control]) + n_distinct(plot_id[warmed]) - 2",
  "by_year  <- contrast(emmeans(m, ~ treatment | species * year, df = df_plots), method = \"revpairwise\")",
  "strength <- contrast(emmeans(m, ~ treatment | species, df = df_plots), method = \"revpairwise\")",
  "# keep \"warmed - origin_control\"; emmeans averages the years with equal weights",
  "```", "",
  paste0("- Fixed part: one mean per treatment x year x species cell. Residual variance: one per species ",
         "(`dispformula = ~ species`, decided 30 Oct). REML; Gaussian on relative cover (bounded, zero-heavy: accepted ",
         "approximation, PLAN.md Known weaknesses)."),
  paste0("- Degrees of freedom: origin-control plots + warmed plots - 2 of the experiment (decided 30 Oct), because ",
         "treatment is applied to whole plots and the yearly surveys of a plot are repeated measures; ",
         min(df_exp), " to ", max(df_exp), " (median ", stats::median(df_exp), "), so the two-sided 5% cut-off ",
         "is t = ", sprintf("%.2f", stats::qt(0.975, max(df_exp))), " to ", sprintf("%.2f", stats::qt(0.975, min(df_exp))),
         " instead of z = 1.96. They change p-values only, not `strength` or `strength_se`."),
  paste0("- Random part, a known limitation kept by decision (30 Oct): `(1 | plot_id:species)` has one plot x species ",
         "variance shared by all species of an experiment, while the residual variance is per species. This shared ",
         "term is more than half of the SE variance in ", sum(plot_share$share > 0.5), " of ", nrow(plot_share),
         " effects (median share ", sprintf("%.2f", stats::median(plot_share$share)), "), and a species' own ",
         "plot-to-plot spread grows with its residual SD (Spearman ",
         sprintf("%.2f", cor(plot_dev$plot_sd, plot_dev$sigma, method = "spearman", use = "complete.obs")),
         "): the shared term overstates the SE of rare species and understates that of dominant ones."),
  "", "## Counts", "", md_table(counts),
  "", "## Strength", "",
  paste0("Quantiles (min / 25% / median / 75% / max). strength: ", q(strength_df$strength), "; strength_se: ",
         q(strength_df$strength_se), ". Effects significant at p < 0.05 (|strength| / SE above the experiment's ",
         "t cut-off): ", fmt(n_sig), " of ", fmt(nrow(strength_df)), "."),
  "", "![strength and strength_se](../figures/diagnostics/02_strength.png)", "",
  paste0("With one residual variance per species, the SE follows each species' own scatter: within an experiment ",
         "the largest SE is a median of ", round(stats::median(fit_table$se_max / fit_table$se_min)), " times the smallest ",
         "(range ", round(min(fit_table$se_max / fit_table$se_min)), " to ", round(max(fit_table$se_max / fit_table$se_min)),
         "), and the Spearman correlation of SE with the species' mean relative cover in the experiment is ",
         sprintf("%.2f", cor(se_cover$strength_se, se_cover$mean_cover, method = "spearman")), ". Rare species therefore ",
         "get the smallest SEs on the absolute cover scale, and the largest weights in the later `rma.mv` calls."),
  "", "## Fits per experiment", "", md_table(fit_table),
  "", "## Property coverage on the residents (added 26 Oct)", "",
  "Species and effects with each property, after steps 01 and 02:", "", md_table(coverage),
  "", sprintf("Residents without any of the three traits (%d): %s.", length(no_trait), paste(no_trait, collapse = ", ")),
  "", "## Covariates", "",
  "Raw values before standardisation:", "",
  md_table(tibble::tibble(variable = c("W (destination - origin JJA T, degC)", "T (destination JJA T, degC)",
                                       "P (log destination JJA precipitation, mm)", "A0_origin (log(% + 1))", "A0_dest (log(% + 1))"),
                          range = c(q(climate$W_raw), q(climate$T_raw), q(climate$P_raw), q(a0$A0_origin_raw), q(a0$A0_dest_raw)))),
  "", "## Warnings captured", "",
  if (length(fit_warnings) + length(warns)) paste("-", c(fit_warnings, unique(warns))) else "None.",
  "", "## Deviations from PLAN.md", "",
  paste("None in the model: PLAN.md Step A was revised on 30 Oct (Bills) to one residual variance per species",
        "(`glmmTMB`, `dispformula = ~ species`). Decisions taken where PLAN.md is silent (30 Oct, Bills): degrees",
        "of freedom = origin-control plots + warmed plots - 2 per experiment; the plot x species variance kept shared",
        "(see Model as fitted); PEV1 to PEV3 carried; properties standardised across the residents that have them, traits",
        "logged first; graminoid and woody as 0/1. MILESTONES' expected count (about 1,596 effects in 41 experiments)",
        "is 1,655 in 43, after FR_AlpeHuez (26 Oct) and US_Arizona (28 Oct) were included.")
)
writeLines(report, "reports/02_strength.md")
message("Step 02 written: ", length(experiment_ids), " experiments, ", fmt(nrow(strength_df)), " effects")
