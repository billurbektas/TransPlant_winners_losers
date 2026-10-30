# Step 03b: strength from year 3 (Step A anchored at year 3, residents).
#
# PLAN.md: "Strength anchored at year 3" (after Step B, added 30 Oct).
# Reads:  data/derived/02_strength/models/<experiment_id>.rds (the step 02
#         fits), strength_df.rds, by_year.rds
# Writes: data/derived/03b_strength_y3/strength_df.rds, strength_df.csv;
#         figures/diagnostics/03b_strength_y3.png; reports/03b_strength_y3.md
# Checks: 1,655 effects, 43 experiments, 21 gradients; every experiment has a
#         survey in year 3 or later; every strength_se finite and positive;
#         strength equals the mean of step 02's yearly contrasts over years >= 3;
#         every column other than strength and strength_se identical to step 02.
#
# Decided 30 Oct (Bills): step 03 found that the yearly contrasts settle from
# year 3, so the main strength of steps 04 to 09 averages years 3 and later,
# from the same fits (no refit, same df rule, all-zero years kept as in step
# 02); step 02's all-years strength is the sensitivity run.

suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
invisible(loadNamespace("glmmTMB")) # registers the emmeans methods for the stored fits
source("R/md_table.R")

warns <- character()
globalCallingHandlers(warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
})
out_dir <- "data/derived/03b_strength_y3"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

strength_all <- readRDS("data/derived/02_strength/strength_df.rds")
by_year <- as.data.frame(readRDS("data/derived/02_strength/by_year.rds")) |>
  mutate(year = as.integer(as.character(year)))
stopifnot(nrow(strength_all) == 1655) # step 02

# ---- Contrast averaged over years >= 3, per experiment ----------------------
experiment_ids <- sort(unique(strength_all$experiment_id))
notes <- character(); years_used <- list(); y3 <- list()
for (e in experiment_ids) {
  fit <- readRDS(file.path("data/derived/02_strength/models", paste0(gsub("::", "__", e), ".rds")))
  years <- levels(fit$model$frame$year)
  keep <- years[as.integer(years) >= 3]
  stopifnot(length(keep) >= 1) # every experiment has a survey in year 3 or later
  st <- withCallingHandlers(
    emmeans::emmeans(fit$model, ~ treatment | species, at = list(year = keep), df = fit$df_plots) |>
      emmeans::contrast(method = "revpairwise") |>
      as.data.frame(),
    message = function(m) { notes <<- c(notes, paste0(e, ": ", trimws(conditionMessage(m)))); invokeRestart("muffleMessage") })
  y3[[e]] <- st |>
    filter(contrast == "warmed - origin_control") |>
    transmute(species = as.character(species), experiment_id = e, strength = estimate, strength_se = SE,
              df = fit$df_plots)
  years_used[[e]] <- tibble::tibble(experiment_id = e, years = paste(keep, collapse = ", "), n_years = length(keep))
}
y3 <- bind_rows(y3); years_used <- bind_rows(years_used)
stopifnot(nrow(y3) == 1655, !anyDuplicated(y3[c("species", "experiment_id")])) # one per species x experiment
stopifnot(all(is.finite(y3$strength_se)), all(y3$strength_se > 0))

# independent check: the year-3 strength is the mean of step 02's yearly
# contrasts over years >= 3 (emmeans averages the years with equal weights)
check <- by_year |>
  filter(year >= 3) |>
  summarise(mean_contrast = mean(estimate), .by = c(experiment_id, species)) |>
  inner_join(y3, by = c("experiment_id", "species"))
stopifnot(nrow(check) == 1655, max(abs(check$mean_contrast - check$strength)) < 1e-10)

# ---- strength_df: step 02's table with strength and its SE replaced ----------
strength_df <- strength_all |>
  select(-strength, -strength_se) |>
  left_join(select(y3, experiment_id, species, strength, strength_se), by = c("experiment_id", "species")) |>
  select(all_of(names(strength_all)))
other <- setdiff(names(strength_all), c("strength", "strength_se"))
stopifnot(identical(strength_df[other], strength_all[other])) # covariates and properties unchanged
stopifnot(n_distinct(strength_df$experiment_id) == 43, n_distinct(strength_df$gradient_id) == 21)

saveRDS(strength_df, file.path(out_dir, "strength_df.rds"))
utils::write.csv(strength_df, file.path(out_dir, "strength_df.csv"), row.names = FALSE)

# ---- Figure ------------------------------------------------------------------
cmp <- inner_join(select(strength_all, effect_id, all_years = strength, se_all = strength_se),
                  select(strength_df, effect_id, from_year_3 = strength, se_y3 = strength_se), by = "effect_id")
lim <- range(c(cmp$all_years, cmp$from_year_3))
pa <- ggplot(cmp, aes(all_years, from_year_3)) +
  geom_abline(colour = "grey60") + geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(alpha = 0.4, size = 1.4) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(title = "a)", x = "Losing-winning strength, all years (step 02)", y = "Losing-winning strength, years 3 and later") +
  theme_classic(base_size = 15) + theme(plot.title = element_text(face = "bold", size = 22))
long <- bind_rows(tibble::tibble(panel = "Losing-winning strength", value = strength_df$strength),
                  tibble::tibble(panel = "Standard error", value = strength_df$strength_se))
pb <- ggplot(long, aes(value)) +
  geom_histogram(bins = 60, fill = "grey55", colour = "white") +
  facet_wrap(~ panel, scales = "free", ncol = 1) +
  labs(title = "b)", x = "Relative cover (proportion), years 3 and later", y = "Species x experiment") +
  theme_classic(base_size = 15) + theme(plot.title = element_text(face = "bold", size = 22),
                                        strip.text = element_text(face = "bold", size = 14), strip.background = element_blank())
ggsave("figures/diagnostics/03b_strength_y3.png", patchwork::wrap_plots(pa, pb, widths = c(1.1, 1)),
       width = 13, height = 6.5, dpi = 150, bg = "white")

# ---- Report --------------------------------------------------------------
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
sig <- function(s, se, df) sum(abs(s) > stats::qt(0.975, df) * se)
df_effect <- y3$df[match(paste(strength_df$experiment_id, strength_df$species), paste(y3$experiment_id, y3$species))]
single <- filter(years_used, n_years == 1)
n_single <- sum(strength_df$experiment_id %in% single$experiment_id)
agreement <- tibble::tibble(
  measure = c("Pearson correlation", "Concordance correlation coefficient", "Sign agreement",
              "Median SE, years >= 3 / all years", "Significant at p < 0.05, all years (step 02)",
              "Significant at p < 0.05, years 3 and later"),
  value = c(sprintf("%.3f", cor(cmp$all_years, cmp$from_year_3)),
            sprintf("%.3f", DescTools::CCC(cmp$all_years, cmp$from_year_3)$rho.c$est),
            sprintf("%.3f", mean(sign(cmp$all_years) == sign(cmp$from_year_3))),
            sprintf("%.2f", stats::median(cmp$se_y3 / cmp$se_all)),
            fmt(sig(strength_all$strength, strength_all$strength_se, df_effect)),
            fmt(sig(strength_df$strength, strength_df$strength_se, df_effect))))
q <- function(x) paste(sprintf("%.4f", stats::quantile(x, c(0, 0.25, 0.5, 0.75, 1))), collapse = " / ")
report <- c(
  "# Step 03b: Strength from year 3 (residents)", "",
  paste0("Losing-winning strength was recomputed from the ", length(experiment_ids), " fits of step 02, without ",
         "refitting, as the warmed minus origin-control contrast averaged over each experiment's years 3 and later ",
         "(`emmeans(m, ~ treatment | species, at = list(year = years >= 3), df = df_plots)`), because step 03 found ",
         "that the yearly contrasts settle from year 3 (decided 30 Oct). The degrees of freedom, the residual ",
         "variance per species and the all-zero years are as in step 02. The covariates and species properties ",
         "of step 02's `strength_df` are carried over unchanged; only `strength` and `strength_se` differ. This is ",
         "the main strength of steps 04 to 09; step 02's all-years strength is the sensitivity run."),
  "", "## Counts", "",
  md_table(tibble::tibble(
    rule = c("effects (species x experiment)", "experiments", "gradients", "species",
             "effects resting on a single year >= 3", "experiments with a single year >= 3"),
    n = fmt(c(nrow(strength_df), n_distinct(strength_df$experiment_id), n_distinct(strength_df$gradient_id),
              n_distinct(strength_df$species), n_single, nrow(single))))),
  "", "Each matches the expectation checked in the script (1,655 effects, 43 experiments, 21 gradients).",
  "", "## Years averaged per experiment", "", md_table(years_used),
  "", "## Against the all-years strength of step 02", "", md_table(agreement),
  "", paste0("Quantiles (min / 25% / median / 75% / max). Strength from year 3: ", q(strength_df$strength),
             "; its SE: ", q(strength_df$strength_se), ". Significance uses each experiment's t cut-off ",
             "(origin-control plots + warmed plots - 2 df), as in step 02."),
  "", "![Strength from year 3](../figures/diagnostics/03b_strength_y3.png)",
  "", "## Warnings captured", "",
  if (length(warns) + length(notes)) paste("-", unique(c(notes, warns))) else "None.",
  "", "## Deviations from PLAN.md", "", "None."
)
writeLines(report, "reports/03b_strength_y3.md")
message("Step 03b written: ", fmt(nrow(strength_df)), " effects; ", n_single, " on a single year")
