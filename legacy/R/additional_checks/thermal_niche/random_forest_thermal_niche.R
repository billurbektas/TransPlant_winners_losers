#!/usr/bin/env Rscript

# Climate, first-year abundance and climatic niche as predictors of winners and
# losers in an unseen gradient: a random forest for each predictor block and
# for every combination of blocks.
# Run: Rscript R/additional_checks/thermal_niche/random_forest_thermal_niche.R
#
# Response: the warmed - origin relative-abundance contrast of
# random_forest_year_group_contrasts_all_species.R at its reference year: the
# emmeans default for a continuous covariate, the mean experimental year of the
# whole model frame (3.13). Read from that script's saved inputs; run it first.
# US_Arizona is excluded, as there.
#
# Blocks; each piece of information enters once:
#   Climate   T_z, P_z, W_z: destination temperature, log precipitation, warming.
#   Abundance relative cover in the first recorded year of the origin controls
#             and of the destination controls, and presence at the destination.
#   Niche     where the species' niche centre sits relative to the destination,
#             dMAT_dest = niche MAT - destination MAT (Lynn's mat_meandiff) and
#             dAP_dest_log = log(niche AP / destination AP), on the log scale of
#             P_z; and the change in Gaussian thermal suitability from origin to
#             destination. With S(T) = -(T - niche MAT)^2 / (2 SD^2),
#             S(destination) - S(origin) = (destination MAT - origin MAT) *
#             (niche MAT - transplant midpoint MAT) / SD^2, where SD = (95th - 5th
#             percentile of the species' MAT) / (2 * 1.645). Warming therefore
#             helps species whose niche is warmer than the midpoint of their
#             transplant, and more so for narrow niches.
# A niche delta at the origin is not used: it equals the destination delta plus
# the warming already in Climate.
#
# Niche data: MacroContextTransplantDat_2.csv, built by J. Lynn (11.2024) with
# DataProcessing_MacroContext.R in TransPlant-Updates_oct2021 (method of Lynn et
# al. 2021, Global Change Biology). Only its columns are used; none of Lynn's
# processing is rerun. Niche centre and percentiles: CHELSA 1981-2010 MAT
# (bio1) and AP (bio12) over GBIF occurrences thinned to one per 30-arc-second
# cell (the "_sub" columns). Site MAT and AP come from the same CHELSA layers at
# each site's coordinates. Read from data/ (gitignored); TRANSPLANT_NICHE_CSV
# overrides. Lynn's rows hold a species only at sites where it was recorded, so
# the deltas are rebuilt from his niche and site columns and checked against
# his mat_meandiff / ap_meandiff.
#
# Cohort: effects whose species has a niche from >= 30 thinned occurrences (Lynn's
# threshold); every specification is fitted and scored on these same rows.
# Forests carry no species ID and no weights: species ID lowered transfer in every
# earlier check, and contrast_se is constant within experiments. min.node.size and
# mtry are tuned by leave-one-gradient-out inside each outer training set, so the
# held-out gradient never informs its own tuning. Each forest is reported beside
# the training mean and an oracle that knows each held-out experiment's true mean
# (the ceiling for any predictor that is constant within an experiment).

library(targets)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ranger)
library(parallel)

out_dir <- "output/additional_checks/thermal_niche"
dir.create(file.path(out_dir, "checkpoints"), recursive = TRUE, showWarnings = FALSE)
base_seed <- 20261023L
n_cores <- max(1L, detectCores() - 2L)
min_occurrences <- 30L

# 1. Effects, contrast and first-year abundances -----------------------------
effects <- tar_read(resident_auc_data) |>
  filter(as.character(gradient_id) != "US_Arizona") |>
  mutate(across(c(experiment_id, gradient_id, species_id), as.character)) |>
  select(-auc_effect, -auc_se)
experiments <- tar_read(experiments) |>
  filter(included, experiment_id %in% effects$experiment_id)
community <- tar_read(analysis_community)

contrast_path <- "output/additional_checks/random_forest_year_group_contrasts_all_species/analysis_inputs.rds"
if (!file.exists(contrast_path)) stop(contrast_path, " is missing: run the year-group contrast script first.")
contrast <- readRDS(contrast_path)$contrasts |>
  transmute(effect_id, reference_year, response = abundance_change)
stopifnot(setequal(contrast$effect_id, effects$effect_id), n_distinct(contrast$reference_year) == 1L)
reference_year <- contrast$reference_year[[1]]

# Mean relative cover of each focal species over the control plots of `site` in
# their first recorded year, with surveyed absences as zero (the definition of
# the all-species abundance check). Origin controls share plots with the
# contrast; destination controls never enter it.
first_year_abundance <- function(experiment, site) {
  controls <- community |>
    filter(region == experiment$region, origin_site == site,
           destination_site == site, treatment == "LocalControl")
  controls <- controls |> filter(year == min(year))
  plots <- distinct(controls, plot_year_id)
  stopifnot(nrow(plots) > 0L)
  cover <- controls |>
    group_by(plot_year_id, source_species = species) |>
    summarise(relative_abundance = sum(relative_cover), .groups = "drop")
  effects |>
    filter(experiment_id == experiment$experiment_id) |>
    select(effect_id, source_species) |>
    cross_join(plots) |>
    left_join(cover, by = c("plot_year_id", "source_species"), relationship = "many-to-one") |>
    mutate(relative_abundance = coalesce(relative_abundance, 0)) |>
    group_by(effect_id) |>
    summarise(abundance = mean(relative_abundance), .groups = "drop")
}
abundance <- list()
for (i in seq_len(nrow(experiments))) {
  experiment <- experiments[i, ]
  abundance[[i]] <- first_year_abundance(experiment, experiment$origin_site) |>
    rename(first_origin_abundance = abundance) |>
    left_join(first_year_abundance(experiment, experiment$destination_site) |>
                rename(dest_abundance = abundance),
              by = "effect_id", relationship = "one-to-one")
}
abundance <- bind_rows(abundance) |> mutate(dest_present = as.numeric(dest_abundance > 0))
stopifnot(!anyDuplicated(abundance$effect_id), nrow(abundance) == nrow(effects))

# 2. Niche centres and site climate -----------------------------------------
niche_csv <- Sys.getenv("TRANSPLANT_NICHE_CSV", "data/MacroContextTransplantDat_2.csv")
macro <- read.csv(niche_csv, stringsAsFactors = FALSE)
niche <- macro |>
  filter(!is.na(species), !is.na(mean_mat_sub), !is.na(mean_ap_sub)) |>
  distinct(niche_species = species, niche_mat = mean_mat_sub, niche_ap = mean_ap_sub,
           niche_mat_p05 = min5_mat_sub, niche_mat_p95 = max95_mat_sub, n_obs_sub)
stopifnot(!anyDuplicated(niche$niche_species))
site_climate <- macro |>
  distinct(region = Region, site = destSiteID, site_mat = mat_site, site_ap = ap_site)
stopifnot(!anyDuplicated(site_climate[c("region", "site")]))
# Lynn's own raw-name match, kept only where a raw name has a single niche species.
raw_names <- macro |>
  filter(species %in% niche$niche_species) |>
  distinct(region = Region, source_species = SpeciesName, raw_niche_species = species) |>
  group_by(region, source_species) |>
  filter(n() == 1L) |>
  ungroup()

# Match on the harmonised name first, then on the raw name within its region,
# then on the raw name alone.
effects <- effects |>
  left_join(abundance, by = "effect_id", relationship = "one-to-one") |>
  left_join(contrast, by = "effect_id", relationship = "one-to-one") |>
  left_join(select(experiments, experiment_id, region, origin_site, destination_site),
            by = "experiment_id", relationship = "many-to-one") |>
  left_join(raw_names, by = c("region", "source_species"), relationship = "many-to-one") |>
  mutate(niche_match = case_when(
      species_id %in% niche$niche_species ~ "Harmonised name",
      !is.na(raw_niche_species) ~ "Raw name within region",
      source_species %in% niche$niche_species ~ "Raw name",
      TRUE ~ "No niche"),
    niche_species = case_when(
      niche_match == "Harmonised name" ~ species_id,
      niche_match == "Raw name within region" ~ raw_niche_species,
      niche_match == "Raw name" ~ source_species)) |>
  left_join(niche, by = "niche_species", relationship = "many-to-one") |>
  left_join(rename(site_climate, origin_site = site, origin_mat = site_mat, origin_ap = site_ap),
            by = c("region", "origin_site"), relationship = "many-to-one") |>
  left_join(rename(site_climate, destination_site = site, dest_mat = site_mat, dest_ap = site_ap),
            by = c("region", "destination_site"), relationship = "many-to-one") |>
  mutate(dMAT_origin = niche_mat - origin_mat, dMAT_dest = niche_mat - dest_mat,
         dAP_origin = (niche_ap - origin_ap) / 1000, dAP_dest = (niche_ap - dest_ap) / 1000,
         dAP_dest_log = log(niche_ap / dest_ap),
         niche_sd = (niche_mat_p95 - niche_mat_p05) / (2 * qnorm(0.95)),
         niche_suitability_change = (dest_mat - origin_mat) *
           (niche_mat - (origin_mat + dest_mat) / 2) / niche_sd^2)
stopifnot(all(is.finite(c(effects$origin_mat, effects$dest_mat, effects$origin_ap, effects$dest_ap))))
lynn_differences <- macro |>
  filter(!is.na(mat_meandiff)) |>
  distinct(region = Region, site = destSiteID, niche_species = species, mat_meandiff, ap_meandiff)
for (end in c("origin", "dest")) {
  site_column <- if (end == "origin") "origin_site" else "destination_site"
  checked <- effects |>
    filter(!is.na(niche_species)) |>
    select(region, site = all_of(site_column), niche_species,
           dMAT = paste0("dMAT_", end), dAP = paste0("dAP_", end)) |>
    inner_join(lynn_differences, by = c("region", "site", "niche_species"),
               relationship = "many-to-many")
  # The CSV stores small differences to six significant digits.
  stopifnot(nrow(checked) > 0L, max(abs(checked$dMAT - checked$mat_meandiff)) < 1e-6,
            max(abs(checked$dAP - checked$ap_meandiff)) < 1e-6)
}

niche_matching <- effects |>
  group_by(species_id, niche_species, niche_match, n_obs_sub, niche_mat, niche_ap, niche_sd) |>
  summarise(source_species = paste(sort(unique(source_species)), collapse = "; "),
            n_effects = n(), n_gradients = n_distinct(gradient_id), .groups = "drop") |>
  mutate(in_cohort = !is.na(n_obs_sub) & n_obs_sub >= min_occurrences)
write.csv(niche_matching, file.path(out_dir, "niche_matching.csv"), row.names = FALSE)

cohort <- effects |>
  filter(!is.na(niche_mat), n_obs_sub >= min_occurrences) |>
  as.data.frame()
stopifnot(all(is.finite(cohort$niche_suitability_change)), all(cohort$niche_sd > 0),
          all(is.finite(cohort$response)))
gradients <- sort(unique(cohort$gradient_id))
cohort_counts <- bind_rows(
  effects |> summarise(cohort = "All non-Arizona effects", n_effects = n(),
    n_species = n_distinct(species_id), n_experiments = n_distinct(experiment_id),
    n_gradients = n_distinct(gradient_id)),
  cohort |> summarise(cohort = paste0("Niche from >= ", min_occurrences, " thinned occurrences"),
    n_effects = n(), n_species = n_distinct(species_id),
    n_experiments = n_distinct(experiment_id), n_gradients = n_distinct(gradient_id))) |>
  mutate(reference_year = reference_year)
write.csv(cohort_counts, file.path(out_dir, "cohort_counts.csv"), row.names = FALSE)
saveRDS(cohort, file.path(out_dir, "analysis_inputs.rds"))

# 3. Blocks and every combination -------------------------------------------
blocks <- list(
  Climate = c("T_z", "P_z", "W_z"),
  Abundance = c("first_origin_abundance", "dest_abundance", "dest_present"),
  Niche = c("dMAT_dest", "dAP_dest_log", "niche_suitability_change"))
specifications <- list()
for (size in seq_along(blocks)) {
  for (combination in combn(names(blocks), size, simplify = FALSE)) {
    specifications[[paste(combination, collapse = " + ")]] <- unlist(blocks[combination], use.names = FALSE)
  }
}

# 4. Nested leave-one-gradient-out forests ----------------------------------
# Half-samples without replacement decorrelate trees and stop each tree from
# memorising repeated rows; tuning searches node sizes from the earlier checks'
# 5 up to leaves of several hundred effects.
node_grid <- c(5L, 25L, 75L, 200L, 500L)
mtry_grid <- function(p) sort(unique(pmin(p, c(1L, floor(sqrt(p)), ceiling(p / 2), p))))
inner_trees <- 200L
outer_trees <- 1000L
n_permutations <- 5L

# The scripts' equal-gradient R2: every held-out gradient gets equal weight.
equal_gradient_R2 <- function(observed, predicted, gradient) {
  mse <- tapply((predicted - observed)^2, gradient, mean)
  mean_y <- tapply(observed, gradient, mean)
  mean_y2 <- tapply(observed^2, gradient, mean)
  1 - mean(mse) / (mean(mean_y2) - mean(mean_y)^2)
}
grow <- function(train, predictors, node, mtry, trees, seed, threads = 1L) {
  ranger(x = train[, predictors, drop = FALSE], y = train$response,
         num.trees = trees, mtry = mtry, min.node.size = node,
         sample.fraction = 0.5, replace = FALSE, splitrule = "variance",
         importance = "none", seed = seed, num.threads = threads)
}
forecast <- function(fit, data, predictors, threads = 1L) {
  predict(fit, data = data[, predictors, drop = FALSE], num.threads = threads)$predictions
}
# Score every grid cell by leave-one-gradient-out within `train` alone. All cells
# share one seed per inner fold, so they differ only by their settings.
tune <- function(train, predictors, seed) {
  grid <- expand.grid(min_node_size = node_grid, mtry = mtry_grid(length(predictors)))
  inner <- sort(unique(train$gradient_id))
  grid$inner_R2 <- vapply(seq_len(nrow(grid)), function(k) {
    predicted <- numeric(nrow(train))
    for (j in seq_along(inner)) {
      held <- train$gradient_id == inner[j]
      fit <- grow(train[!held, ], predictors, grid$min_node_size[k], grid$mtry[k],
                  inner_trees, seed + j)
      predicted[held] <- forecast(fit, train[held, ], predictors)
    }
    equal_gradient_R2(train$response, predicted, train$gradient_id)
  }, numeric(1))
  grid
}

ranger_version <- as.character(packageVersion("ranger"))
fold_results <- list()
for (s in seq_along(specifications)) {
  specification <- names(specifications)[[s]]
  predictors <- specifications[[s]]
  checkpoint <- file.path(out_dir, "checkpoints",
                          paste0(gsub("[^A-Za-z0-9]+", "_", specification), ".rds"))
  signature <- digest::digest(list(
    data = cohort[c("effect_id", "gradient_id", "response", predictors)],
    predictors = predictors, node_grid = node_grid, mtry = "1, sqrt, p/2, p",
    inner_trees = inner_trees, outer_trees = outer_trees, permutations = n_permutations,
    sample_fraction = 0.5, replace = FALSE, seed = base_seed + 1000L * s,
    ranger_version = ranger_version))
  if (file.exists(checkpoint) && identical(readRDS(checkpoint)$signature, signature)) {
    fold_results[[specification]] <- readRDS(checkpoint)
    next
  }
  message("Nested forests: ", specification)
  folds <- mclapply(seq_along(gradients), function(g) {
    held <- cohort$gradient_id == gradients[[g]]
    train <- cohort[!held, ]
    test <- cohort[held, ]
    seed <- base_seed + 1000L * s + 10L * g
    grid <- tune(train, predictors, seed)
    # Ties go to the larger node size, then the smaller mtry (the smoother forest).
    best <- grid[order(-grid$inner_R2, -grid$min_node_size, grid$mtry), ][1, ]
    fit <- grow(train, predictors, best$min_node_size, best$mtry, outer_trees, seed)
    # Transfer importance: shuffle one predictor among the held-out gradient's rows.
    set.seed(seed + 7L)
    permuted <- list()
    for (predictor in predictors) {
      for (r in seq_len(n_permutations)) {
        shuffled <- test
        shuffled[[predictor]] <- sample(test[[predictor]])
        permuted[[paste(predictor, r)]] <- data.frame(predictor = predictor, permutation = r,
          effect_id = test$effect_id, predicted = forecast(fit, shuffled, predictors))
      }
    }
    list(predictions = data.frame(effect_id = test$effect_id,
                                  predicted = forecast(fit, test, predictors)),
         grid = mutate(grid, heldout_gradient = gradients[[g]]),
         tuning = data.frame(heldout_gradient = gradients[[g]], min_node_size = best$min_node_size,
                             mtry = best$mtry, inner_R2 = best$inner_R2),
         permuted = bind_rows(permuted))
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  failed <- vapply(folds, inherits, logical(1), "try-error")
  if (any(failed)) stop(specification, ": ", folds[failed][[1]])
  saved <- list(signature = signature,
    predictions = bind_rows(lapply(folds, `[[`, "predictions")),
    grid = bind_rows(lapply(folds, `[[`, "grid")),
    tuning = bind_rows(lapply(folds, `[[`, "tuning")),
    permuted = bind_rows(lapply(folds, `[[`, "permuted")))
  stopifnot(nrow(saved$predictions) == nrow(cohort),
            setequal(saved$predictions$effect_id, cohort$effect_id))
  saveRDS(saved, paste0(checkpoint, ".tmp"))
  stopifnot(file.rename(paste0(checkpoint, ".tmp"), checkpoint))
  fold_results[[specification]] <- saved
}
stack_results <- function(part) {
  bind_rows(lapply(names(fold_results), function(specification) {
    mutate(fold_results[[specification]][[part]], specification = specification)
  }))
}

keys <- cohort |>
  transmute(effect_id, gradient_id, experiment_id, species_id, observed = response)
# References: the training mean of the other gradients, and the oracle.
training_mean <- vapply(cohort$gradient_id, function(g) mean(cohort$response[cohort$gradient_id != g]),
                        numeric(1))
predictions <- bind_rows(
  mutate(stack_results("predictions"), method = "Random forest"),
  data.frame(effect_id = cohort$effect_id, predicted = training_mean,
             specification = "Training mean", method = "Reference"),
  data.frame(effect_id = cohort$effect_id, predicted = ave(cohort$response, cohort$experiment_id),
             specification = "Oracle: held-out experiment mean", method = "Reference")) |>
  left_join(keys, by = "effect_id", relationship = "many-to-one")
saveRDS(predictions, file.path(out_dir, "predictions.rds"))
tuning <- stack_results("tuning")
write.csv(tuning, file.path(out_dir, "tuning.csv"), row.names = FALSE)
write.csv(stack_results("grid"), file.path(out_dir, "tuning_grid.csv"), row.names = FALSE)

# 5. Scores -----------------------------------------------------------------
# R2 on every species x experiment pair of a held-out gradient; Spearman and
# winner/loser skill on species means within that gradient (as in earlier checks).
by_gradient <- predictions |>
  group_by(specification, method, gradient_id) |>
  summarise(n_effects = n(), MSE = mean((predicted - observed)^2),
            mean_y = mean(observed), mean_y2 = mean(observed^2),
            R2pred = 1 - sum((predicted - observed)^2) / sum((observed - mean(observed))^2),
            .groups = "drop") |>
  left_join(predictions |>
      group_by(specification, method, gradient_id, species_id) |>
      summarise(observed = mean(observed), predicted = mean(predicted), .groups = "drop") |>
      group_by(specification, method, gradient_id) |>
      summarise(n_species = n(),
        Spearman = if (n() >= 2L && n_distinct(observed) > 1L && n_distinct(predicted) > 1L)
          cor(observed, predicted, method = "spearman") else NA_real_,
        winner_recall = if (any(observed > 0)) mean(predicted[observed > 0] > 0) else NA_real_,
        loser_recall = if (any(observed < 0)) mean(predicted[observed < 0] < 0) else NA_real_,
        .groups = "drop") |>
      mutate(winner_loser_skill = winner_recall + loser_recall - 1),
    by = c("specification", "method", "gradient_id"), relationship = "one-to-one")
write.csv(by_gradient, file.path(out_dir, "metrics_by_gradient.csv"), row.names = FALSE)

# Whole-gradient bootstrap: resample the 19 gradients and recompute the
# equal-gradient R2 of every specification on the same draws.
n_boot <- 2000L
set.seed(base_seed)
draws <- replicate(n_boot, sample(gradients, replace = TRUE))
boot <- by_gradient |>
  group_by(specification, method) |>
  group_modify(~ data.frame(draw = seq_len(n_boot), R2 = apply(draws, 2, function(draw) {
    s <- .x[match(draw, .x$gradient_id), ]
    1 - mean(s$MSE) / (mean(s$mean_y2) - mean(s$mean_y)^2)
  }))) |>
  ungroup()
overall <- by_gradient |>
  group_by(specification, method) |>
  summarise(n_gradients = n(), n_effects = sum(n_effects),
    R2pred = 1 - mean(MSE) / (mean(mean_y2) - mean(mean_y)^2),
    Spearman = mean(Spearman, na.rm = TRUE), n_Spearman_defined = sum(is.finite(Spearman)),
    winner_loser_skill = mean(winner_loser_skill, na.rm = TRUE), .groups = "drop") |>
  left_join(boot |> group_by(specification, method) |>
              summarise(R2_low = quantile(R2, 0.025), R2_high = quantile(R2, 0.975), .groups = "drop"),
            by = c("specification", "method"), relationship = "one-to-one")
# Training R2 of forests refitted to all gradients with each specification's most
# frequently chosen settings; the gap to held-out R2 shows how much is memorised.
modal_settings <- tuning |>
  count(specification, min_node_size, mtry) |>
  group_by(specification) |>
  slice_max(n, n = 1L, with_ties = FALSE) |>
  ungroup()
full_fits <- list()
training <- list()
for (s in seq_along(specifications)) {
  specification <- names(specifications)[[s]]
  setting <- modal_settings[modal_settings$specification == specification, ]
  fit <- grow(cohort, specifications[[s]], setting$min_node_size, setting$mtry, outer_trees,
              base_seed + s, threads = n_cores)
  full_fits[[specification]] <- fit
  training[[s]] <- data.frame(specification, method = "Random forest",
    training_R2 = equal_gradient_R2(cohort$response, forecast(fit, cohort, specifications[[s]], n_cores),
                                    cohort$gradient_id))
}
overall <- overall |>
  left_join(bind_rows(training), by = c("specification", "method"), relationship = "one-to-one")
write.csv(overall, file.path(out_dir, "equal_gradient_R2.csv"), row.names = FALSE)

# Paired comparisons: what adding one block does, on the same held-out rows,
# the same gradients and the same bootstrap draws.
comparisons <- tribble(
  ~comparison,                          ~with,                         ~without,
  "Abundance added to climate",         "Climate + Abundance",         "Climate",
  "Niche added to climate",             "Climate + Niche",             "Climate",
  "Niche added to abundance",           "Abundance + Niche",           "Abundance",
  "Niche added to climate + abundance", "Climate + Abundance + Niche", "Climate + Abundance",
  "Climate added to abundance + niche", "Climate + Abundance + Niche", "Abundance + Niche")
paired_gradients <- function(with, without) {
  by_gradient |>
    filter(method == "Random forest", specification == with) |>
    select(gradient_id, MSE_with = MSE, Spearman_with = Spearman, skill_with = winner_loser_skill) |>
    left_join(by_gradient |>
                filter(method == "Random forest", specification == without) |>
                select(gradient_id, MSE_without = MSE, Spearman_without = Spearman,
                       skill_without = winner_loser_skill),
              by = "gradient_id", relationship = "one-to-one")
}
comparison_table <- bind_rows(lapply(seq_len(nrow(comparisons)), function(i) {
  pair <- paired_gradients(comparisons$with[i], comparisons$without[i])
  R2 <- function(name) overall$R2pred[overall$specification == name & overall$method == "Random forest"]
  boot_change <- boot$R2[boot$specification == comparisons$with[i]] -
    boot$R2[boot$specification == comparisons$without[i]]
  data.frame(comparisons[i, ],
    R2_change = R2(comparisons$with[i]) - R2(comparisons$without[i]),
    R2_change_low = quantile(boot_change, 0.025), R2_change_high = quantile(boot_change, 0.975),
    n_gradients_lower_error = sum(pair$MSE_with < pair$MSE_without),
    Spearman_change = mean(pair$Spearman_with - pair$Spearman_without, na.rm = TRUE),
    n_gradients_higher_Spearman = sum(pair$Spearman_with > pair$Spearman_without, na.rm = TRUE),
    n_gradients_Spearman_defined = sum(is.finite(pair$Spearman_with - pair$Spearman_without)),
    skill_change = mean(pair$skill_with - pair$skill_without, na.rm = TRUE), row.names = NULL)
}))
write.csv(comparison_table, file.path(out_dir, "comparisons.csv"), row.names = FALSE)

# Drop in held-out equal-gradient R2 when one predictor is shuffled.
heldout_R2 <- overall |>
  filter(method == "Random forest") |>
  select(specification, heldout_R2 = R2pred)
importance <- stack_results("permuted") |>
  left_join(keys, by = "effect_id", relationship = "many-to-one") |>
  group_by(specification, predictor, permutation) |>
  summarise(R2 = equal_gradient_R2(observed, predicted, gradient_id), .groups = "drop") |>
  left_join(heldout_R2, by = "specification", relationship = "many-to-one") |>
  group_by(specification, predictor, heldout_R2) |>
  summarise(R2_drop = mean(heldout_R2 - R2), R2_drop_min = min(heldout_R2 - R2),
            R2_drop_max = max(heldout_R2 - R2), .groups = "drop")
write.csv(importance, file.path(out_dir, "transfer_importance.csv"), row.names = FALSE)

# 6. Direction of the niche effects ------------------------------------------
# Partial dependence from the all-gradient forests: the mean prediction with one
# predictor set to each value and the others left as observed.
niche_predictors <- blocks$Niche
direction_specifications <- c("Climate + Niche", "Climate + Abundance + Niche")
pd <- list()
for (specification in direction_specifications) {
  X <- cohort[, specifications[[specification]], drop = FALSE]
  for (variable in niche_predictors) {
    values <- unique(quantile(cohort[[variable]], seq(0.025, 0.975, length.out = 30L), names = FALSE))
    pd[[paste(specification, variable)]] <- data.frame(specification, variable, value = values,
      predicted = vapply(values, function(v) {
        X[[variable]] <- v
        mean(forecast(full_fits[[specification]], X, specifications[[specification]], n_cores))
      }, numeric(1)))
  }
}
pd <- bind_rows(pd)
write.csv(pd, file.path(out_dir, "partial_dependence.csv"), row.names = FALSE)
binned <- bind_rows(lapply(niche_predictors, function(variable) {
  cohort |>
    mutate(bin = ntile(.data[[variable]], 10L)) |>
    group_by(bin) |>
    summarise(value = median(.data[[variable]]), observed = mean(response),
              se = sd(response) / sqrt(n()), n_effects = n(), .groups = "drop") |>
    mutate(variable = variable)
}))
write.csv(binned, file.path(out_dir, "observed_by_niche_decile.csv"), row.names = FALSE)

# Niche position x warming, held out: experiments in warming terciles, effects in
# quintiles of dMAT_dest within each tercile; observed means beside the
# leave-one-gradient-out predictions of climate alone and climate + niche.
warming_terciles <- cohort |>
  distinct(experiment_id, W) |>
  mutate(warming_tercile = ntile(W, 3L))
tercile_labels <- warming_terciles |>
  group_by(warming_tercile) |>
  summarise(label = sprintf("Warming %.1f–%.1f °C (%d experiments)", min(W), max(W), n()),
            .groups = "drop")
interaction_bins <- predictions |>
  filter(specification %in% c("Climate", "Climate + Niche")) |>
  left_join(cohort |> select(effect_id, dMAT_dest), by = "effect_id", relationship = "many-to-one") |>
  left_join(warming_terciles |> select(experiment_id, warming_tercile), by = "experiment_id",
            relationship = "many-to-one") |>
  group_by(specification, warming_tercile) |>
  mutate(niche_bin = ntile(dMAT_dest, 5L)) |>
  group_by(specification, warming_tercile, niche_bin) |>
  summarise(dMAT_dest = median(dMAT_dest), se = sd(observed) / sqrt(n()),
            observed = mean(observed), predicted = mean(predicted), n_effects = n(),
            .groups = "drop") |>
  left_join(tercile_labels, by = "warming_tercile", relationship = "many-to-one")
write.csv(interaction_bins, file.path(out_dir, "warming_by_niche_position.csv"), row.names = FALSE)

# 7. Comparison figures -----------------------------------------------------
ink <- "#0b0b0b"
ink_2 <- "#52514e"
grid_colour <- "#e4e3df"
specification_colours <- c("Climate" = "#2a78d6", "Abundance" = "#eb6834", "Niche" = "#1baf7a",
  "Climate + Abundance" = "#eda100", "Climate + Niche" = "#e87ba4",
  "Abundance + Niche" = "#008300", "Climate + Abundance + Niche" = "#6250d6")
reference_colour <- "#8a8984"
theme_check <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.3),
        axis.text = element_text(colour = ink_2), axis.title = element_text(colour = ink_2),
        plot.title = element_text(colour = ink, face = "bold"),
        plot.subtitle = element_text(colour = ink_2), legend.position = "top",
        legend.title = element_blank(), strip.text = element_text(colour = ink, hjust = 0))
response_label <- "Warmed − origin relative abundance"
contrast_label <- sprintf("mean-year contrast (year %.2f)", reference_year)
oracle_R2 <- overall$R2pred[overall$specification == "Oracle: held-out experiment mean"]

r2_plot <- overall |>
  mutate(specification = factor(specification, levels = rev(c("Training mean",
           "Oracle: held-out experiment mean", names(specifications)))),
         colour = if_else(method == "Reference", "Reference", as.character(specification))) |>
  ggplot(aes(R2pred, specification, colour = colour)) +
  geom_vline(xintercept = 0, colour = ink_2, linewidth = 0.4) +
  geom_vline(xintercept = oracle_R2, colour = ink_2, linewidth = 0.4, linetype = "dashed") +
  geom_linerange(aes(xmin = R2_low, xmax = R2_high), linewidth = 0.7) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.3f", R2pred)), colour = ink_2, size = 3, vjust = -1) +
  scale_colour_manual(values = c(specification_colours, Reference = reference_colour), guide = "none") +
  labs(x = "Equal-gradient predictive R², leave one gradient out", y = NULL,
       title = "Climate, abundance and niche: predicting the contrast in an unseen gradient",
       subtitle = paste0("Random forests, ", contrast_label, ". ", nrow(cohort), " effects, ",
                         n_distinct(cohort$species_id), " species, ", length(gradients), " gradients.",
                         "\nBars: 95% whole-gradient bootstrap. Dashed line: oracle ceiling for experiment-level predictors.")) +
  theme_check

comparison_plot <- comparison_table |>
  mutate(comparison = factor(comparison, levels = rev(comparisons$comparison)),
         label = sprintf("%+.3f  (%d/19 gradients lower error; Spearman %+.2f)",
                         R2_change, n_gradients_lower_error, Spearman_change)) |>
  ggplot(aes(R2_change, comparison)) +
  geom_vline(xintercept = 0, colour = ink_2, linewidth = 0.4) +
  geom_linerange(aes(xmin = R2_change_low, xmax = R2_change_high), colour = ink_2, linewidth = 0.7) +
  geom_point(size = 2.6, colour = ink) +
  geom_text(aes(x = R2_change_high, label = label), colour = ink_2, size = 3, hjust = -0.05) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.9))) +
  labs(x = "Change in held-out R² when the block is added (95% whole-gradient bootstrap)", y = NULL,
       title = "What each block adds", subtitle = paste0("Random forests, ", contrast_label,
       "; Spearman = mean change in within-gradient rank correlation")) +
  theme_check

importance_plot <- importance |>
  filter(specification %in% direction_specifications) |>
  mutate(specification = factor(specification, levels = direction_specifications)) |>
  ggplot(aes(R2_drop, reorder(predictor, R2_drop))) +
  geom_vline(xintercept = 0, colour = ink_2, linewidth = 0.4) +
  geom_linerange(aes(xmin = R2_drop_min, xmax = R2_drop_max, colour = specification), linewidth = 0.7) +
  geom_point(aes(colour = specification), size = 2.6) +
  scale_colour_manual(values = specification_colours, guide = "none") +
  facet_wrap(~ specification, ncol = 1, scales = "free_y") +
  labs(x = "Drop in held-out R² when the predictor is shuffled within the held-out gradient",
       y = NULL, title = "Which predictors carry the transfer?",
       subtitle = "Mean and range over five shuffles of the nested forests") +
  theme_check

variable_labels <- c(dMAT_dest = "Niche MAT − destination MAT (°C)",
                     dAP_dest_log = "log(niche AP / destination AP)",
                     niche_suitability_change = "Gaussian suitability change, origin → destination")
pd_plot <- ggplot() +
  geom_hline(yintercept = 0, colour = ink_2, linewidth = 0.4) +
  geom_pointrange(data = mutate(binned, variable = factor(variable_labels[variable], levels = variable_labels)),
                  aes(value, observed, ymin = observed - 2 * se, ymax = observed + 2 * se),
                  colour = reference_colour, size = 0.25, linewidth = 0.4) +
  geom_line(data = mutate(pd, variable = factor(variable_labels[variable], levels = variable_labels),
                          specification = factor(specification, levels = direction_specifications)),
            aes(value, predicted, colour = specification), linewidth = 0.8) +
  scale_colour_manual(values = specification_colours) +
  facet_wrap(~ variable, scales = "free_x") +
  labs(x = NULL, y = response_label, title = "Direction of the niche effects",
       subtitle = "Lines: partial dependence of the all-gradient forests. Grey: observed mean ± 2 SE by decile.") +
  theme_check + theme(panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.3))

interaction_plot <- interaction_bins |>
  mutate(label = factor(label, levels = tercile_labels$label)) |>
  ggplot(aes(dMAT_dest)) +
  geom_hline(yintercept = 0, colour = ink_2, linewidth = 0.4) +
  geom_pointrange(data = function(x) filter(x, specification == "Climate"),
                  aes(y = observed, ymin = observed - 2 * se, ymax = observed + 2 * se),
                  colour = reference_colour, size = 0.3, linewidth = 0.4) +
  geom_line(aes(y = predicted, colour = specification), linewidth = 0.8) +
  geom_point(aes(y = predicted, colour = specification), size = 1.8) +
  scale_colour_manual(values = specification_colours) +
  facet_wrap(~ label, nrow = 1) +
  labs(x = "Niche MAT − destination MAT (°C), quintile medians", y = response_label,
       title = "Niche position matters more under stronger warming?",
       subtitle = "Grey: observed mean ± 2 SE. Lines: leave-one-gradient-out forest predictions.") +
  theme_check + theme(panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.3))

cairo_pdf(file.path(out_dir, "thermal_niche.pdf"), width = 10, height = 6.5, onefile = TRUE)
print(r2_plot)
print(comparison_plot)
print(importance_plot)
print(pd_plot)
print(interaction_plot)
invisible(dev.off())

# 8. Per-gradient score densities, laid out as the earlier checks'
# predictive_metric_densities.pdf: rows are comparisons, columns are metrics,
# and each panel overlays three specifications. A specification keeps its
# colour across rows.
density_sets <- list(
  "Each block alone" = c("Climate", "Abundance", "Niche"),
  "Added to climate" = c("Climate", "Climate + Abundance", "Climate + Niche"),
  "Niche added to abundance" = c("Climate + Abundance", "Abundance + Niche", "Climate + Abundance + Niche"))
metric_labels <- c(R2pred = "Predictive R²", Spearman = "Spearman rank correlation",
                   winner_loser_skill = "Winner-loser prediction skill")
density_metrics <- bind_rows(lapply(names(density_sets), function(set_name) {
  by_gradient |>
    filter(method == "Random forest", specification %in% density_sets[[set_name]]) |>
    mutate(set = set_name, position = match(specification, density_sets[[set_name]]))
})) |>
  select(set, position, specification, gradient_id, R2pred, Spearman, winner_loser_skill) |>
  pivot_longer(c(R2pred, Spearman, winner_loser_skill), names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric_labels[metric], levels = metric_labels),
         set = factor(set, levels = names(density_sets)),
         specification = factor(specification, levels = names(specifications)))
density_counts <- density_metrics |>
  group_by(set, metric, position) |>
  summarise(n_defined = sum(is.finite(value)), .groups = "drop") |>
  arrange(set, metric, position) |>
  group_by(set, metric) |>
  summarise(label = paste0("n = ", paste(n_defined, collapse = " | ")), .groups = "drop")
density_metrics <- density_metrics |> filter(is.finite(value))
density_means <- density_metrics |>
  group_by(set, metric, specification, position) |>
  summarise(mean_score = mean(value), .groups = "drop") |>
  left_join(overall |>
      filter(method == "Random forest") |>
      transmute(specification = factor(specification, levels = names(specifications)),
                metric = factor(metric_labels[["R2pred"]], levels = metric_labels), overall_R2 = R2pred),
    by = c("specification", "metric"), relationship = "many-to-one") |>
  mutate(label = if_else(is.na(overall_R2), sprintf("Mean = %.3f", mean_score),
                         sprintf("Mean = %.3f | Overall = %.3f", mean_score, overall_R2)),
         label_line = 1.2 + (position - 1L) * 1.4)
write.csv(density_means |> mutate(across(c(set, metric, specification), as.character)),
          file.path(out_dir, "density_means.csv"), row.names = FALSE)
density_ranges <- density_metrics |>
  group_by(set, metric) |>
  summarise(low = if (first(metric) == metric_labels[["R2pred"]])
              min(value) - 0.1 * diff(range(value)) else -1,
            high = 1, .groups = "drop") |>
  pivot_longer(c(low, high), values_to = "value")

cairo_pdf(file.path(out_dir, "predictive_metric_densities.pdf"), width = 20, height = 13, onefile = TRUE)
# Training R2 predicts the rows used to fit each all-gradient forest; it is not
# held-out or out-of-bag. Held-out R2 is the leave-one-gradient-out score.
tiles <- overall |>
  filter(method == "Random forest") |>
  select(specification, "Training R²" = training_R2, "Held-out R²" = R2pred) |>
  pivot_longer(-specification, names_to = "score", values_to = "value") |>
  mutate(specification = factor(specification, levels = rev(names(specifications))),
         score = factor(score, levels = c("Training R²", "Held-out R²")))
print(ggplot(tiles, aes(score, specification)) +
  geom_tile(aes(fill = specification), alpha = 0.18, colour = "white", linewidth = 2) +
  geom_text(aes(label = sprintf("%.3f", value)), size = 10) +
  scale_fill_manual(values = specification_colours, guide = "none") +
  scale_x_discrete(position = "top", expand = expansion(add = 0)) +
  scale_y_discrete(expand = expansion(add = 0)) +
  labs(x = NULL, y = NULL, title = "Full-training and held-out R² of the random forests",
       subtitle = paste("Relative-abundance", contrast_label), tag = "A",
       caption = paste(
         "Training R2: predict the same observations used to fit each all-gradient forest (each specification's most often tuned settings). Not held-out or out-of-bag.",
         "Held-out R2: leave one gradient out, forests tuned inside each training set. Both weight each gradient equally and use individual species x experiment effects.",
         "Cell colour identifies the specification, not R2 magnitude. US_Arizona is excluded.",
         sep = "\n")) +
  theme_minimal(base_size = 25) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24), plot.title = element_text(size = 30, face = "bold"),
        plot.subtitle = element_text(size = 24), plot.tag = element_text(size = 30, face = "bold"),
        plot.tag.position = "topleft", plot.caption = element_text(size = 15, hjust = 0),
        plot.margin = margin(20, 25, 20, 25)))
print(ggplot(density_metrics, aes(value, fill = specification, colour = specification)) +
  geom_blank(data = density_ranges, aes(value, 0), inherit.aes = FALSE) +
  geom_density(data = filter(density_metrics, metric == metric_labels[["R2pred"]]),
               bounds = c(-Inf, 1), n = 4096, alpha = 0.25, linewidth = 0.9) +
  geom_density(data = filter(density_metrics, metric != metric_labels[["R2pred"]]),
               bounds = c(-1, 1), n = 4096, alpha = 0.25, linewidth = 0.9) +
  geom_rug(sides = "b", alpha = 0.5, linewidth = 0.4, show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey20", linewidth = 0.6) +
  geom_text(data = density_counts, aes(label = label), x = Inf, y = Inf,
            hjust = 1.1, vjust = 1.3, size = 5, inherit.aes = FALSE) +
  geom_text(data = density_means, aes(label = label, colour = specification, vjust = label_line),
            x = -Inf, y = Inf, hjust = -0.08, size = 4.8, inherit.aes = FALSE, show.legend = FALSE) +
  facet_wrap(vars(set, metric), ncol = 3, scales = "free",
             labeller = label_wrap_gen(width = 50, multi_line = FALSE)) +
  scale_fill_manual(values = specification_colours, name = NULL) +
  scale_colour_manual(values = specification_colours, name = NULL) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.5), nrow = 1)) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 5), expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 4), expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = "Density", title = "Random forests: climate, abundance and niche",
       subtitle = paste("Relative-abundance", contrast_label, "— every held-out species"), tag = "B",
       caption = paste(
         sprintf("Response: fitted warmed minus origin-control relative abundance at experimental year %.4f (overall model-frame mean; emmeans default for continuous year).", reference_year),
         "Leave one gradient out; forests tuned by an inner leave-one-gradient-out loop, no species ID, no weights.",
         "Within each held-out gradient, R2 compares every species x experiment contrast with its prediction. Rank and category scores use one mean contrast per species per gradient.",
         "Mean = average of valid gradient scores. Overall = equal-gradient R2 across all effect pairs, not the mean of gradient R2 values.",
         "Prediction skill = winner detection rate + loser detection rate - 1; winner/loser uses the contrast sign, not significance. Zero is the reference for all three metrics.",
         "Climate = destination temperature + log precipitation + warming. Abundance = first recorded year of origin and destination controls, and presence at the destination.",
         "Niche = Lynn's niche MAT minus destination MAT, log niche AP / destination AP, and the Gaussian thermal suitability change from origin to destination.",
         "Automatic density smoothing; ticks show gradient scores; n gives valid counts in legend order. US_Arizona and species with < 30 thinned GBIF occurrences excluded.",
         sep = "\n")) +
  theme_bw(base_size = 20) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92"), panel.spacing = grid::unit(1, "lines"),
        strip.text = element_text(size = 17, face = "bold"),
        strip.background = element_rect(fill = "grey96", colour = "grey75"),
        axis.title.y = element_text(size = 22), axis.text = element_text(size = 16),
        plot.title = element_text(size = 24, face = "bold"), plot.subtitle = element_text(size = 20),
        plot.tag = element_text(size = 26, face = "bold"), plot.tag.position = "topleft",
        plot.caption = element_text(size = 13, hjust = 0), legend.position = "top",
        legend.text = element_text(size = 17), legend.key.width = grid::unit(1.1, "cm"),
        plot.margin = margin(15, 20, 15, 20)))
invisible(dev.off())

print(overall |>
  arrange(method, desc(R2pred)) |>
  select(specification, method, R2pred, R2_low, R2_high, Spearman, winner_loser_skill, training_R2) |>
  as.data.frame(), digits = 3)
print(comparison_table, digits = 3)
