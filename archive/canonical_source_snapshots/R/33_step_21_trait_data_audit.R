#!/usr/bin/env Rscript

# Complete pre-model audit of LDMC, LMA, and vegetative height for the retained
# resident species universe. This script does not fit an AUC model, alter source
# traits, impute missing values, remove outliers, or change taxonomy.

options(lifecycle_verbosity = "quiet")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
output_dir <- file.path(root, "output", "step_by_step_modeling", "trait_data_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  analysis = file.path(root, "output", "tables", "analysis_dataset.csv"),
  phylogeny_species = file.path(
    root, "output", "step_by_step_modeling", "resident_phylogeny_pcoa",
    "species_pev_scores.csv"
  ),
  functional_groups = file.path(root, "output", "tables", "functional_groups.csv"),
  vegetation_taxonomy = file.path(root, "output", "tables", "included_vegetation_species.csv"),
  try_taxonomy = file.path(root, "output", "tables", "tnrs_try_review.csv"),
  trait_summary = file.path(root, "output", "tables", "try_species_traits.csv"),
  trait_individual = file.path(root, "output", "intermediate", "try_traits_individual.rds"),
  trait_map = file.path(root, "config", "try_trait_map.csv"),
  raw_try = file.path(root, "data", "try_47163.txt"),
  prior_scaling = file.path(
    root, "output", "step_by_step_modeling", "species_both_traits",
    "trait_row_covariance_scaling_audit.csv"
  )
)

missing_inputs <- names(paths)[!file.exists(unlist(paths))]
if (length(missing_inputs)) {
  stop("Missing required trait-audit input(s): ", paste(missing_inputs, collapse = ", "))
}

excluded_species <- c(
  "Gymnocarpium dryopteris", "Athyrium spinulosum", "Botrychium lunaria",
  "Equisetum arvense", "Selaginella selaginoides", "Juniperus communis"
)
traits <- c("LDMC", "LMA", "vegetative_height")
trait_labels <- c(LDMC = "LDMC", LMA = "LMA", vegetative_height = "Vegetative height")

collapse_values <- function(x) {
  x <- iconv(as.character(x), from = "", to = "UTF-8", sub = "?")
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (!length(x)) return(NA_character_)
  paste(sort(unique(as.character(x))), collapse = "; ")
}

format_values <- function(x) {
  paste(format(x, digits = 15, trim = TRUE, scientific = FALSE), collapse = "; ")
}

skewness_moment <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L || stats::sd(x) == 0) return(NA_real_)
  mean((x - mean(x))^3) / stats::sd(x)^3
}

table_lines <- function(x) {
  capture.output(print(tibble::as_tibble(x), n = Inf, width = Inf))
}

read_cleaned_try_observations <- function() {
  if (requireNamespace("targets", quietly = TRUE) && dir.exists(file.path(root, "_targets"))) {
    target <- tryCatch(
      targets::tar_read_raw("try_climate_filtered", store = file.path(root, "_targets")),
      error = function(e) NULL
    )
    if (!is.null(target) && is.list(target) && "observations" %in% names(target)) {
      return(list(data = target$observations, source = "targets::tar_read(try_climate_filtered)$observations"))
    }
  }
  x <- readRDS(paths$trait_individual)
  list(data = x, source = "output/intermediate/try_traits_individual.rds (limited provenance fallback)")
}

# -----------------------------------------------------------------------------
# Retained species universe and response-independent species/experiment keys
# -----------------------------------------------------------------------------

phy <- read_csv(paths$phylogeny_species, show_col_types = FALSE) |>
  transmute(
    species_id = as.character(.data$species_id),
    tnrs_accepted_name = as.character(.data$tnrs_accepted_name),
    family = as.character(.data$tnrs_family)
  ) |>
  distinct()

if (nrow(phy) != 526L || n_distinct(phy$species_id) != 526L) {
  stop("Expected exactly 526 unique retained species in the current phylogeny species set.")
}
if (any(phy$species_id %in% excluded_species)) {
  stop("A declared downstream-excluded species occurs in the retained species universe.")
}

functional <- read_csv(paths$functional_groups, show_col_types = FALSE) |>
  transmute(
    species_id = as.character(.data$accepted_name),
    graminoid_status = as.character(.data$graminoid_status),
    woody_status = as.character(.data$woody_status)
  ) |>
  distinct()

if (anyDuplicated(functional$species_id)) {
  stop("Functional-status table is not one row per accepted species.")
}

analysis <- read_csv(paths$analysis, show_col_types = FALSE)
effect_keys <- analysis |>
  filter(
    .data$transform == "raw", as.character(.data$horizon) == "max",
    .data$contrast == "to_origin", .data$pool == "resident",
    !.data$accepted_name %in% excluded_species
  ) |>
  transmute(
    species_id = as.character(.data$accepted_name),
    experiment_id = as.character(.data$experiment_id),
    T = as.numeric(.data$destination_T_JJA),
    P = log(as.numeric(.data$destination_P_JJA)),
    W = as.numeric(.data$diff_T_JJA)
  ) |>
  distinct()

if (nrow(effect_keys) != 1604L || n_distinct(effect_keys$species_id) != 526L ||
    n_distinct(effect_keys$experiment_id) != 42L) {
  stop("Current downstream resident key set is not 1,604 effects / 526 species / 42 experiments.")
}
if (!setequal(effect_keys$species_id, phy$species_id)) {
  stop("Phylogeny and downstream resident species universes do not agree.")
}
if (anyDuplicated(effect_keys[c("species_id", "experiment_id")])) {
  stop("Downstream resident species-experiment keys are duplicated.")
}

experiment_climate <- effect_keys |>
  select(.data$experiment_id, .data$T, .data$P, .data$W) |>
  distinct()
if (nrow(experiment_climate) != 42L || anyDuplicated(experiment_climate$experiment_id)) {
  stop("Experiment-level climate values are inconsistent.")
}
for (v in c("T", "P", "W")) {
  experiment_climate[[paste0(v, "_z")]] <-
    (experiment_climate[[v]] - mean(experiment_climate[[v]])) / stats::sd(experiment_climate[[v]])
}

# -----------------------------------------------------------------------------
# Eligible individual TRY source records and current species medians
# -----------------------------------------------------------------------------

try_loaded <- read_cleaned_try_observations()
try_source <- try_loaded$source
try_obs_source <- try_loaded$data

full_provenance <- all(c(
  "SpeciesName", "AccSpeciesName", "TraitName", "StdValue", "UnitName",
  "DatasetID", "Dataset", "ObservationID", "accepted_name", "climate_code"
) %in% names(try_obs_source))

if (full_provenance) {
  trait_observations <- try_obs_source |>
    filter(.data$TraitName %in% c("LDMC", "SLA", "vegetative_height")) |>
    transmute(
      species_id = as.character(.data$accepted_name),
      submitted_species_name = as.character(.data$SpeciesName),
      try_accepted_species_name = as.character(.data$AccSpeciesName),
      trait = if_else(.data$TraitName == "SLA", "LMA", as.character(.data$TraitName)),
      source_trait = as.character(.data$TraitName),
      source_value = as.numeric(.data$StdValue),
      value = if_else(.data$TraitName == "SLA", 1 / as.numeric(.data$StdValue), as.numeric(.data$StdValue)),
      source_unit = as.character(.data$UnitName),
      unit = if_else(.data$TraitName == "SLA", paste0("1/(", .data$UnitName, ")"), as.character(.data$UnitName)),
      dataset_id = as.character(.data$DatasetID),
      dataset = as.character(.data$Dataset),
      observation_id = as.character(.data$ObservationID),
      climate_code = as.character(.data$climate_code)
    )
} else {
  trait_observations <- try_obs_source |>
    transmute(
      species_id = as.character(.data$accepted_name),
      submitted_species_name = NA_character_,
      try_accepted_species_name = as.character(.data$accepted_name),
      trait = as.character(.data$trait),
      source_trait = if_else(.data$trait == "LMA", "SLA", as.character(.data$trait)),
      source_value = if_else(.data$trait == "LMA", 1 / as.numeric(.data$value), as.numeric(.data$value)),
      value = as.numeric(.data$value),
      source_unit = as.character(.data$source_unit),
      unit = as.character(.data$unit),
      dataset_id = as.character(.data$dataset_id),
      dataset = paste0("TRY DatasetID ", .data$dataset_id),
      observation_id = NA_character_,
      climate_code = as.character(.data$climate_code)
    )
}

trait_observations_retained <- trait_observations |>
  filter(.data$species_id %in% phy$species_id)

saved_summary <- read_csv(paths$trait_summary, show_col_types = FALSE) |>
  filter(.data$accepted_name %in% phy$species_id, .data$trait %in% traits) |>
  transmute(
    species_id = as.character(.data$accepted_name),
    trait = as.character(.data$trait),
    unit = as.character(.data$unit),
    value = as.numeric(.data$median),
    saved_minimum = as.numeric(.data$minimum),
    saved_maximum = as.numeric(.data$maximum),
    saved_n_observations = as.integer(.data$n_observations),
    saved_n_datasets = as.integer(.data$n_datasets),
    saved_datasets = as.character(.data$datasets)
  )

if (anyDuplicated(saved_summary[c("species_id", "trait")])) {
  stop("Saved trait summary has duplicate retained species-trait rows.")
}

recalculated_summary <- trait_observations_retained |>
  summarise(
    recalculated_median = median(.data$value),
    recalculated_minimum = min(.data$value),
    recalculated_maximum = max(.data$value),
    recalculated_n = n(),
    .by = c("species_id", "trait")
  )

summary_validation <- saved_summary |>
  left_join(recalculated_summary, by = c("species_id", "trait")) |>
  mutate(
    median_equal = abs(.data$value - .data$recalculated_median) <= 1e-12,
    minimum_equal = abs(.data$saved_minimum - .data$recalculated_minimum) <= 1e-12,
    maximum_equal = abs(.data$saved_maximum - .data$recalculated_maximum) <= 1e-12,
    n_equal = .data$saved_n_observations == .data$recalculated_n
  )

source_metadata <- trait_observations_retained |>
  summarise(
    raw_records = n(),
    source_dataset_ids = collapse_values(.data$dataset_id),
    source_datasets = collapse_values(.data$dataset),
    source_units = collapse_values(.data$source_unit),
    unit = collapse_values(.data$unit),
    submitted_names = collapse_values(.data$submitted_species_name),
    try_accepted_names = collapse_values(.data$try_accepted_species_name),
    aggregation_rule = "median of eligible cleaned individual TRY records",
    .by = c("species_id", "trait")
  )

trait_wide <- saved_summary |>
  select(.data$species_id, .data$trait, .data$value) |>
  pivot_wider(names_from = .data$trait, values_from = .data$value)

meta_wide <- source_metadata |>
  pivot_wider(
    names_from = .data$trait,
    values_from = c(
      .data$raw_records, .data$source_dataset_ids, .data$source_datasets,
      .data$source_units, .data$unit, .data$submitted_names, .data$try_accepted_names,
      .data$aggregation_rule
    ),
    names_glue = "{trait}_{.value}"
  )

species_audit <- phy |>
  left_join(functional, by = "species_id") |>
  left_join(trait_wide, by = "species_id") |>
  left_join(meta_wide, by = "species_id") |>
  mutate(
    across(ends_with("_raw_records"), ~ replace_na(as.integer(.x), 0L))
  ) |>
  select(
    .data$species_id, .data$tnrs_accepted_name, .data$family,
    .data$graminoid_status, .data$woody_status,
    .data$LDMC, .data$LMA, .data$vegetative_height,
    .data$LDMC_raw_records, .data$LMA_raw_records, .data$vegetative_height_raw_records,
    .data$LDMC_source_dataset_ids, .data$LMA_source_dataset_ids,
    .data$vegetative_height_source_dataset_ids,
    .data$LDMC_source_datasets, .data$LMA_source_datasets,
    .data$vegetative_height_source_datasets,
    .data$LDMC_unit, .data$LMA_unit, .data$vegetative_height_unit,
    .data$LDMC_source_units, .data$LMA_source_units,
    .data$vegetative_height_source_units,
    .data$LDMC_aggregation_rule, .data$LMA_aggregation_rule,
    .data$vegetative_height_aggregation_rule,
    .data$LDMC_submitted_names, .data$LMA_submitted_names,
    .data$vegetative_height_submitted_names,
    .data$LDMC_try_accepted_names, .data$LMA_try_accepted_names,
    .data$vegetative_height_try_accepted_names
  ) |>
  arrange(.data$species_id)

# -----------------------------------------------------------------------------
# Duplicate/conflict and exact-repeat audits
# -----------------------------------------------------------------------------

record_groups <- trait_observations_retained |>
  summarise(
    n_records = n(),
    n_distinct_values = n_distinct(.data$value),
    n_units = n_distinct(.data$unit),
    n_sources = n_distinct(.data$dataset_id),
    .by = c("species_id", "trait")
  )

conflicts <- trait_observations_retained |>
  arrange(.data$species_id, .data$trait, .data$value) |>
  summarise(
    accepted_name = first(.data$species_id),
    all_raw_values = format_values(.data$source_value),
    units = collapse_values(.data$source_unit),
    all_current_trait_values = format_values(.data$value),
    current_trait_units = collapse_values(.data$unit),
    sources = collapse_values(paste0("TRY ", .data$dataset_id, ": ", .data$dataset)),
    current_aggregation_rule = "median of eligible cleaned individual TRY records",
    range = max(.data$value) - min(.data$value),
    mean = mean(.data$value),
    median = median(.data$value),
    n_raw_records = n(),
    n_distinct_values = n_distinct(.data$value),
    .by = c("species_id", "trait")
  ) |>
  filter(.data$n_distinct_values > 1L) |>
  select(
    .data$species_id, .data$accepted_name, .data$trait, .data$all_raw_values,
    .data$units, .data$all_current_trait_values, .data$current_trait_units,
    .data$sources, .data$current_aggregation_rule, .data$range,
    .data$mean, .data$median, .data$n_raw_records, .data$n_distinct_values
  ) |>
  arrange(.data$trait, .data$species_id)

species_values <- species_audit |>
  select(.data$species_id, accepted_name = .data$tnrs_accepted_name, .data$family,
         .data$graminoid_status, .data$woody_status, all_of(traits)) |>
  pivot_longer(all_of(traits), names_to = "trait", values_to = "value") |>
  left_join(
    source_metadata |>
      select(.data$species_id, .data$trait, .data$raw_records,
             source = .data$source_datasets, unit = .data$unit),
    by = c("species_id", "trait")
  )

repeat_frequency <- species_values |>
  filter(is.finite(.data$value)) |>
  count(.data$trait, .data$value, name = "n_species") |>
  group_by(.data$trait) |>
  arrange(desc(.data$n_species), .data$value, .by_group = TRUE) |>
  mutate(
    rank = row_number(),
    proportion = .data$n_species / sum(.data$n_species)
  ) |>
  ungroup()

# -----------------------------------------------------------------------------
# Value integrity, descriptive extremes, and transformations
# -----------------------------------------------------------------------------

value_summary <- bind_rows(lapply(traits, function(tr) {
  x_all <- species_values$value[species_values$trait == tr]
  x <- x_all[is.finite(x_all)]
  tibble(
    trait = tr,
    n_with_data = length(x),
    n_missing = sum(is.na(x_all)),
    minimum = min(x),
    p01 = quantile(x, 0.01, names = FALSE),
    p05 = quantile(x, 0.05, names = FALSE),
    p25 = quantile(x, 0.25, names = FALSE),
    median = median(x),
    mean = mean(x),
    p75 = quantile(x, 0.75, names = FALSE),
    p95 = quantile(x, 0.95, names = FALSE),
    p99 = quantile(x, 0.99, names = FALSE),
    maximum = max(x),
    sd = stats::sd(x),
    coefficient_of_variation = stats::sd(x) / mean(x),
    skewness = skewness_moment(x),
    zero = sum(x == 0),
    negative = sum(x < 0),
    non_finite_provided = sum(!is.na(x_all) & !is.finite(x_all)),
    unit = collapse_values(species_values$unit[species_values$trait == tr])
  )
}))

transformation_summary <- value_summary |>
  select(.data$trait, raw_skewness = .data$skewness) |>
  left_join(
    species_values |>
      filter(is.finite(.data$value), .data$value > 0) |>
      summarise(
        log10_skewness = skewness_moment(log10(.data$value)),
        skewness_absolute_change = abs(skewness_moment(log10(.data$value))) -
          abs(skewness_moment(.data$value)),
        .by = .data$trait
      ),
    by = "trait"
  )

extreme_rows <- lapply(traits, function(tr) {
  x <- species_values |>
    filter(.data$trait == tr, is.finite(.data$value)) |>
    arrange(.data$value, .data$species_id)
  q <- quantile(x$value, c(0.25, 0.75), names = FALSE)
  iqr <- diff(q)
  med <- median(x$value)
  mad_value <- stats::mad(x$value)
  x |>
    mutate(
      low10 = row_number() <= 10L,
      high10 = row_number() > n() - 10L,
      iqr_outlier = .data$value < q[[1]] - 1.5 * iqr | .data$value > q[[2]] + 1.5 * iqr,
      robust_z = if (is.finite(mad_value) && mad_value > 0) (.data$value - med) / mad_value else NA_real_,
      mad_outlier = is.finite(.data$robust_z) & abs(.data$robust_z) > 3.5,
      outlier_criteria = paste0(
        if_else(.data$low10, "lowest 10; ", ""),
        if_else(.data$high10, "highest 10; ", ""),
        if_else(.data$iqr_outlier, "1.5x IQR; ", ""),
        if_else(.data$mad_outlier, "|robust z| > 3.5; ", "")
      ),
      outlier_criteria = sub("; $", "", .data$outlier_criteria)
    ) |>
    filter(nzchar(.data$outlier_criteria)) |>
    transmute(
      .data$species_id, .data$accepted_name, .data$family, .data$trait,
      .data$value, .data$unit, .data$outlier_criteria,
      robust_standardized_deviation = .data$robust_z,
      number_of_raw_observations = .data$raw_records,
      .data$source
    )
})
extreme_values <- bind_rows(extreme_rows) |>
  arrange(.data$trait, .data$value, .data$species_id)

extreme_concentration_family <- extreme_values |>
  count(.data$trait, .data$family, name = "n_extreme_species") |>
  group_by(.data$trait) |>
  slice_max(.data$n_extreme_species, n = 3L, with_ties = FALSE) |>
  ungroup()

extreme_concentration_source <- trait_observations_retained |>
  semi_join(extreme_values, by = c("species_id", "trait")) |>
  distinct(.data$species_id, .data$trait, .data$dataset_id, .data$dataset) |>
  count(.data$trait, .data$dataset_id, .data$dataset, name = "n_extreme_species") |>
  group_by(.data$trait) |>
  slice_max(.data$n_extreme_species, n = 3L, with_ties = FALSE) |>
  ungroup()

extreme_concentration_climate <- trait_observations_retained |>
  semi_join(extreme_values, by = c("species_id", "trait")) |>
  distinct(.data$species_id, .data$trait, .data$climate_code) |>
  count(.data$trait, .data$climate_code, name = "n_extreme_species") |>
  group_by(.data$trait) |>
  slice_max(.data$n_extreme_species, n = 3L, with_ties = FALSE) |>
  ungroup()

extreme_record_support <- extreme_values |>
  select(.data$species_id, .data$trait) |>
  distinct() |>
  left_join(record_groups, by = c("species_id", "trait")) |>
  summarise(
    candidate_species = n(),
    based_on_one_raw_record = sum(.data$n_records == 1L),
    based_on_multiple_raw_records = sum(.data$n_records > 1L),
    median_raw_records = median(.data$n_records),
    maximum_raw_records = max(.data$n_records),
    .by = .data$trait
  )

# -----------------------------------------------------------------------------
# Standardization audit
# -----------------------------------------------------------------------------

current_scaling <- species_values |>
  filter(is.finite(.data$value), .data$value > 0) |>
  mutate(log_value = log(.data$value)) |>
  summarise(
    unique_species_n = n(),
    unique_species_log_mean = mean(.data$log_value),
    unique_species_log_sd = stats::sd(.data$log_value),
    .by = .data$trait
  )

row_weighted_scaling <- effect_keys |>
  inner_join(
    species_values |> select(.data$species_id, .data$trait, .data$value),
    by = "species_id", relationship = "many-to-many"
  ) |>
  filter(is.finite(.data$value), .data$value > 0) |>
  mutate(log_value = log(.data$value)) |>
  summarise(
    expanded_rows_n = n(),
    expanded_rows_log_mean = mean(.data$log_value),
    expanded_rows_log_sd = stats::sd(.data$log_value),
    .by = .data$trait
  )

prior_scaling <- read_csv(paths$prior_scaling, show_col_types = FALSE) |>
  transmute(
    trait = as.character(.data$trait),
    prior_trait_transform = as.character(.data$trait_transform),
    prior_center = as.numeric(.data$trait_center_log_scale),
    prior_scale = as.numeric(.data$trait_scale_log_scale),
    prior_n_species = as.integer(.data$trait_scaling_n_species),
    prior_weighting = as.character(.data$trait_scaling_weighting)
  )

standardization_audit <- current_scaling |>
  left_join(row_weighted_scaling, by = "trait") |>
  left_join(prior_scaling, by = "trait") |>
  mutate(
    prior_method = if_else(
      .data$prior_weighting == "distinct species",
      "A: unique-species standardization",
      "B: expanded species x experiment row-weighting"
    )
  )

# Verify that the legacy z-trait columns in analysis_dataset.csv are row-weighted.
legacy_z_audit <- bind_rows(lapply(traits, function(tr) {
  z_col <- paste0("z_", tr)
  log_x <- log(analysis[[tr]])
  z_row <- (log_x - mean(log_x, na.rm = TRUE)) / stats::sd(log_x, na.rm = TRUE)
  unique_values <- analysis |>
    distinct(.data$accepted_name, value = .data[[tr]]) |>
    filter(is.finite(.data$value), .data$value > 0)
  z_unique <- (log_x - mean(log(unique_values$value))) / stats::sd(log(unique_values$value))
  tibble(
    trait = tr,
    max_abs_difference_from_expanded_row_z = max(abs(analysis[[z_col]] - z_row), na.rm = TRUE),
    max_abs_difference_from_unique_species_z = max(abs(analysis[[z_col]] - z_unique), na.rm = TRUE)
  )
}))

# -----------------------------------------------------------------------------
# Coverage, missingness, and cross-trait consistency
# -----------------------------------------------------------------------------

coverage_summary <- bind_rows(lapply(traits, function(tr) {
  species_with <- species_values |>
    filter(.data$trait == tr, is.finite(.data$value)) |>
    pull(.data$species_id)
  matched_effects <- effect_keys |> filter(.data$species_id %in% species_with)
  matched_experiments <- experiment_climate |>
    filter(.data$experiment_id %in% matched_effects$experiment_id)
  species_per_experiment <- matched_effects |>
    count(.data$experiment_id, name = "n_species")
  tibble(
    trait = tr,
    retained_species_with_trait = length(species_with),
    retained_species_without_trait = 526L - length(species_with),
    auc_effects_represented = nrow(matched_effects),
    experiments_represented = n_distinct(matched_effects$experiment_id),
    unique_experiments = nrow(matched_experiments),
    T_z_min = min(matched_experiments$T_z), T_z_max = max(matched_experiments$T_z),
    P_z_min = min(matched_experiments$P_z), P_z_max = max(matched_experiments$P_z),
    W_z_min = min(matched_experiments$W_z), W_z_max = max(matched_experiments$W_z),
    species_per_experiment_min = min(species_per_experiment$n_species),
    species_per_experiment_median = median(species_per_experiment$n_species),
    species_per_experiment_max = max(species_per_experiment$n_species)
  )
}))

full_coverage <- tibble(
  trait = "FULL_526_UNIVERSE",
  retained_species_with_trait = 526L,
  retained_species_without_trait = 0L,
  auc_effects_represented = nrow(effect_keys),
  experiments_represented = n_distinct(effect_keys$experiment_id),
  unique_experiments = nrow(experiment_climate),
  T_z_min = min(experiment_climate$T_z), T_z_max = max(experiment_climate$T_z),
  P_z_min = min(experiment_climate$P_z), P_z_max = max(experiment_climate$P_z),
  W_z_min = min(experiment_climate$W_z), W_z_max = max(experiment_climate$W_z),
  species_per_experiment_min = min((effect_keys |> count(.data$experiment_id))$n),
  species_per_experiment_median = median((effect_keys |> count(.data$experiment_id))$n),
  species_per_experiment_max = max((effect_keys |> count(.data$experiment_id))$n)
)

missing_species_rows <- species_values |>
  filter(is.na(.data$value) | !is.finite(.data$value)) |>
  transmute(
    record_type = "missing_species",
    .data$trait,
    group_variable = NA_character_, group_value = NA_character_,
    .data$species_id, tnrs_accepted_name = .data$accepted_name,
    .data$family, .data$graminoid_status, .data$woody_status,
    n_species_total = NA_integer_, n_species_with_trait = NA_integer_,
    n_species_missing = NA_integer_, coverage_percent = NA_real_,
    auc_effects_represented = NA_integer_, experiments_represented = NA_integer_,
    T_z_min = NA_real_, T_z_max = NA_real_, P_z_min = NA_real_, P_z_max = NA_real_,
    W_z_min = NA_real_, W_z_max = NA_real_,
    species_per_experiment_min = NA_real_, species_per_experiment_median = NA_real_,
    species_per_experiment_max = NA_real_
  )

group_missingness_one <- function(group_var) {
  species_values |>
    mutate(group_value_internal = as.character(.data[[group_var]])) |>
    summarise(
      n_species_total = n(),
      n_species_with_trait = sum(is.finite(.data$value)),
      n_species_missing = sum(!is.finite(.data$value) | is.na(.data$value)),
      coverage_percent = 100 * .data$n_species_with_trait / .data$n_species_total,
      .by = c("trait", "group_value_internal")
    ) |>
    transmute(
      record_type = "group_missingness", .data$trait,
      group_variable = group_var, group_value = .data$group_value_internal,
      species_id = NA_character_, tnrs_accepted_name = NA_character_, family = NA_character_,
      graminoid_status = NA_character_, woody_status = NA_character_,
      .data$n_species_total, .data$n_species_with_trait, .data$n_species_missing,
      .data$coverage_percent,
      auc_effects_represented = NA_integer_, experiments_represented = NA_integer_,
      T_z_min = NA_real_, T_z_max = NA_real_, P_z_min = NA_real_, P_z_max = NA_real_,
      W_z_min = NA_real_, W_z_max = NA_real_,
      species_per_experiment_min = NA_real_, species_per_experiment_median = NA_real_,
      species_per_experiment_max = NA_real_
    )
}

coverage_rows <- bind_rows(full_coverage, coverage_summary) |>
  transmute(
    record_type = "coverage_summary", .data$trait,
    group_variable = NA_character_, group_value = NA_character_,
    species_id = NA_character_, tnrs_accepted_name = NA_character_, family = NA_character_,
    graminoid_status = NA_character_, woody_status = NA_character_,
    n_species_total = 526L,
    n_species_with_trait = .data$retained_species_with_trait,
    n_species_missing = .data$retained_species_without_trait,
    coverage_percent = 100 * .data$retained_species_with_trait / 526,
    .data$auc_effects_represented, .data$experiments_represented,
    .data$T_z_min, .data$T_z_max, .data$P_z_min, .data$P_z_max,
    .data$W_z_min, .data$W_z_max,
    .data$species_per_experiment_min, .data$species_per_experiment_median,
    .data$species_per_experiment_max
  )

experiment_rows <- bind_rows(lapply(traits, function(tr) {
  species_with <- species_values |>
    filter(.data$trait == tr, is.finite(.data$value)) |>
    pull(.data$species_id)
  effect_keys |>
    mutate(has_trait = .data$species_id %in% species_with) |>
    summarise(
      n_species_total = n(),
      n_species_with_trait = sum(.data$has_trait),
      n_species_missing = sum(!.data$has_trait),
      coverage_percent = 100 * .data$n_species_with_trait / .data$n_species_total,
      .by = .data$experiment_id
    ) |>
    transmute(
      record_type = "experiment_coverage", trait = tr,
      group_variable = "experiment_id", group_value = .data$experiment_id,
      species_id = NA_character_, tnrs_accepted_name = NA_character_, family = NA_character_,
      graminoid_status = NA_character_, woody_status = NA_character_,
      .data$n_species_total, .data$n_species_with_trait, .data$n_species_missing,
      .data$coverage_percent,
      auc_effects_represented = NA_integer_, experiments_represented = NA_integer_,
      T_z_min = NA_real_, T_z_max = NA_real_, P_z_min = NA_real_, P_z_max = NA_real_,
      W_z_min = NA_real_, W_z_max = NA_real_,
      species_per_experiment_min = NA_real_, species_per_experiment_median = NA_real_,
      species_per_experiment_max = NA_real_
    )
}))

missingness_coverage <- bind_rows(
  coverage_rows, missing_species_rows,
  group_missingness_one("family"),
  group_missingness_one("graminoid_status"),
  group_missingness_one("woody_status"),
  experiment_rows
)

trait_presence <- species_audit |>
  transmute(
    .data$species_id,
    LDMC = is.finite(.data$LDMC),
    LMA = is.finite(.data$LMA),
    vegetative_height = is.finite(.data$vegetative_height),
    n_traits = rowSums(across(all_of(traits)))
  )
trait_count_summary <- trait_presence |>
  count(.data$n_traits, name = "n_species") |>
  complete(n_traits = 0:3, fill = list(n_species = 0L)) |>
  mutate(category = case_when(
    .data$n_traits == 3 ~ "all three",
    .data$n_traits == 2 ~ "exactly two",
    .data$n_traits == 1 ~ "exactly one",
    TRUE ~ "none"
  )) |>
  arrange(desc(.data$n_traits))

pairs <- list(c("LDMC", "LMA"), c("LDMC", "vegetative_height"), c("LMA", "vegetative_height"))
pairwise_data <- bind_rows(lapply(pairs, function(pair) {
  dat <- species_audit |>
    select(.data$species_id, x = all_of(pair[[1]]), y = all_of(pair[[2]])) |>
    filter(is.finite(.data$x), is.finite(.data$y))
  tibble(
    pair = paste(trait_labels[[pair[[1]]]], "vs", trait_labels[[pair[[2]]]]),
    n = nrow(dat),
    pearson = cor(dat$x, dat$y, method = "pearson"),
    spearman = cor(dat$x, dat$y, method = "spearman")
  )
}))

pair_scatter <- bind_rows(lapply(pairs, function(pair) {
  x_name <- pair[[1]]
  y_name <- pair[[2]]
  pair_name <- paste(unname(trait_labels[x_name]), "vs", unname(trait_labels[y_name]))
  species_audit |>
    transmute(
      .data$species_id,
      pair = .env$pair_name,
      x_trait = unname(.env$trait_labels[.env$x_name]),
      y_trait = unname(.env$trait_labels[.env$y_name]),
      x = .data[[x_name]], y = .data[[y_name]]
    ) |>
    filter(is.finite(.data$x), is.finite(.data$y))
})) |>
  left_join(pairwise_data, by = "pair") |>
  mutate(
    facet = sprintf("%s\nN = %d; Pearson r = %.2f; Spearman rho = %.2f",
                    .data$pair, .data$n, .data$pearson, .data$spearman)
  )

# -----------------------------------------------------------------------------
# Taxonomic match audit and readiness checks
# -----------------------------------------------------------------------------

try_taxonomy <- read_csv(paths$try_taxonomy, show_col_types = FALSE)
vegetation_taxonomy <- read_csv(paths$vegetation_taxonomy, show_col_types = FALSE)

multiple_try_names <- try_taxonomy |>
  filter(.data$accepted_name_reviewed %in% phy$species_id) |>
  distinct(.data$submitted_name, .data$accepted_name_reviewed) |>
  count(.data$accepted_name_reviewed, name = "n_submitted_names") |>
  filter(.data$n_submitted_names > 1L)

multiple_ecological_names <- vegetation_taxonomy |>
  filter(.data$accepted_name_reviewed %in% phy$species_id) |>
  distinct(.data$name_original, .data$accepted_name_reviewed) |>
  count(.data$accepted_name_reviewed, name = "n_original_names") |>
  filter(.data$n_original_names > 1L)

record_to_multiple_species <- if (full_provenance) {
  trait_observations_retained |>
    filter(!is.na(.data$observation_id)) |>
    summarise(n_species = n_distinct(.data$species_id), .by = c("observation_id", "source_trait")) |>
    filter(.data$n_species > 1L)
} else {
  tibble(observation_id = character(), source_trait = character(), n_species = integer())
}

impossible_species_values <- species_values |>
  filter(
    !is.na(.data$value) & (!is.finite(.data$value) | .data$value <= 0 |
      (.data$trait == "LDMC" & .data$value > 1))
  )
unexpected_units <- trait_observations_retained |>
  distinct(.data$trait, .data$source_unit, .data$unit) |>
  filter(
    (.data$trait == "LDMC" & .data$source_unit != "g/g") |
      (.data$trait == "LMA" & .data$source_unit != "mm2 mg-1") |
      (.data$trait == "vegetative_height" & .data$source_unit != "m")
  )
major_group_absence <- missingness_coverage |>
  filter(
    .data$record_type == "group_missingness", .data$n_species_total >= 10L,
    .data$n_species_with_trait == 0L
  )

readiness_checks <- c(
  retained_universe_is_526 = nrow(species_audit) == 526L && n_distinct(species_audit$species_id) == 526L,
  downstream_sets_match = setequal(effect_keys$species_id, species_audit$species_id),
  excluded_absent = !any(species_audit$species_id %in% excluded_species) &&
    !any(trait_observations_retained$species_id %in% excluded_species),
  canonical_names_one_to_one = !anyDuplicated(phy$species_id) && !anyDuplicated(phy$tnrs_accepted_name) &&
    all(phy$species_id == phy$tnrs_accepted_name),
  summaries_are_unique = !anyDuplicated(saved_summary[c("species_id", "trait")]),
  saved_medians_reproduced = nrow(summary_validation) == nrow(saved_summary) &&
    all(summary_validation$median_equal & summary_validation$minimum_equal &
          summary_validation$maximum_equal & summary_validation$n_equal),
  no_impossible_values = nrow(impossible_species_values) == 0L,
  standardized_units_consistent = nrow(unexpected_units) == 0L,
  no_record_maps_multiple_species = nrow(record_to_multiple_species) == 0L,
  no_major_group_completely_absent = nrow(major_group_absence) == 0L,
  prior_isolated_models_used_unique_species_scaling = all(prior_scaling$prior_weighting == "distinct species")
)
ready <- all(readiness_checks)

# -----------------------------------------------------------------------------
# Requested static figures
# -----------------------------------------------------------------------------

distribution_data <- species_values |>
  filter(is.finite(.data$value), .data$value > 0) |>
  transmute(
    trait = factor(trait_labels[.data$trait], levels = unname(trait_labels)),
    raw = .data$value,
    log10 = log10(.data$value)
  ) |>
  pivot_longer(c("raw", "log10"), names_to = "scale", values_to = "display_value") |>
  mutate(scale = factor(.data$scale, levels = c("raw", "log10"), labels = c("Raw species median", "log10(species median)")))

p_distribution <- ggplot(distribution_data, aes(x = .data$display_value)) +
  geom_histogram(bins = 30, color = "white", fill = "#2C7FB8", linewidth = 0.25) +
  facet_grid(rows = vars(.data$trait), cols = vars(.data$scale), scales = "free") +
  labs(
    x = "Trait value (panel-specific scale and documented unit)",
    y = "Number of retained species"
  ) +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey94"), panel.grid.minor = element_blank())

ggsave(
  file.path(output_dir, "trait_distribution_audit.pdf"), p_distribution,
  width = 10, height = 8, device = cairo_pdf
)

p_pairwise <- ggplot(pair_scatter, aes(x = .data$x, y = .data$y)) +
  geom_point(alpha = 0.65, size = 1.6, color = "#2C7FB8") +
  facet_wrap(vars(.data$facet), scales = "free", ncol = 3) +
  labs(x = "First trait (raw species median)", y = "Second trait (raw species median)") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey94"), panel.grid.minor = element_blank())

ggsave(
  file.path(output_dir, "trait_pairwise_audit.pdf"), p_pairwise,
  width = 12, height = 4.5, device = cairo_pdf
)

# -----------------------------------------------------------------------------
# Requested CSVs
# -----------------------------------------------------------------------------

write_csv(species_audit, file.path(output_dir, "trait_species_level_audit.csv"), na = "")
write_csv(conflicts, file.path(output_dir, "trait_duplicate_conflict_audit.csv"), na = "")
write_csv(extreme_values, file.path(output_dir, "trait_extreme_values.csv"), na = "")
write_csv(missingness_coverage, file.path(output_dir, "trait_missingness_and_coverage.csv"), na = "")

# -----------------------------------------------------------------------------
# Text audit report
# -----------------------------------------------------------------------------

source_table <- tibble::tribble(
  ~trait, ~raw_source_column, ~species_name_columns, ~documented_unit, ~current_aggregation, ~previous_model_transform, ~previous_model_standardization,
  "LDMC", "TraitName aliases -> StdValue", "SpeciesName; AccSpeciesName; TNRS-reviewed accepted_name", "g/g", "species median", "natural log", "mean/SD over distinct trait-complete species",
  "LMA", "SLA TraitName aliases -> 1 / StdValue", "SpeciesName; AccSpeciesName; TNRS-reviewed accepted_name", "1/(mm2 mg-1) (the pipeline's explicit reciprocal-SLA label)", "species median after record-level reciprocal", "natural log", "mean/SD over distinct trait-complete species",
  "vegetative_height", "Plant height vegetative aliases -> StdValue", "SpeciesName; AccSpeciesName; TNRS-reviewed accepted_name", "m", "species median", "natural log", "mean/SD over distinct trait-complete species"
)

duplicate_summary <- record_groups |>
  summarise(
    species_trait_groups = n(),
    groups_with_multiple_records = sum(.data$n_records > 1L),
    multiple_records_all_identical = sum(.data$n_records > 1L & .data$n_distinct_values == 1L),
    conflicting_distinct_values = sum(.data$n_distinct_values > 1L),
    groups_with_multiple_units = sum(.data$n_units > 1L),
    .by = .data$trait
  )

common_values <- repeat_frequency |>
  filter(.data$rank <= 10L) |>
  mutate(value = signif(.data$value, 8), proportion = round(100 * .data$proportion, 2)) |>
  select(.data$trait, .data$rank, .data$value, .data$n_species, percent_of_trait_complete_species = .data$proportion)

possible_rounding_spikes <- repeat_frequency |>
  filter(.data$n_species >= 5L, .data$proportion >= 0.05) |>
  transmute(
    .data$trait, value = signif(.data$value, 8), .data$n_species,
    percent_of_trait_complete_species = round(100 * .data$proportion, 2),
    interpretation = "possible coarse-rounding spike; not evidence of imputation/default by itself"
  )

candidate_extreme_summary <- extreme_values |>
  summarise(
    candidate_species = n_distinct(.data$species_id),
    IQR_candidates = sum(grepl("IQR", .data$outlier_criteria)),
    MAD_candidates = sum(grepl("robust z", .data$outlier_criteria)),
    .by = .data$trait
  )

source_stage <- read_csv(file.path(root, "output", "tables", "try_filter_audit.csv"), show_col_types = FALSE)
range_gate <- read_csv(file.path(root, "output", "tables", "try_range_validation.csv"), show_col_types = FALSE)

report <- c(
  "COMPLETE PRE-MODEL TRAIT-DATA AUDIT",
  "===================================",
  "",
  "No AUC model was fitted. No trait value, outlier, missing value, or taxonomic name was changed.",
  "All trait summaries below are species-level and precede expansion to species x experiment rows.",
  "",
  "1. DOWNSTREAM SPECIES UNIVERSE",
  "------------------------------",
  sprintf("Retained downstream resident species: %d", nrow(species_audit)),
  sprintf("Resident species-experiment keys used only for coverage accounting: %d", nrow(effect_keys)),
  sprintf("Experiments: %d", n_distinct(effect_keys$experiment_id)),
  paste("Excluded downstream but retained upstream:", paste(excluded_species, collapse = "; ")),
  sprintf("All six excluded taxa absent from species and matched-trait audit sets: %s", readiness_checks[["excluded_absent"]]),
  "",
  "2. EXACT TRAIT SOURCES AND CURRENT PROCESSING",
  "----------------------------------------------",
  paste("Original TRY export:", paths$raw_try),
  paste("Current cleaned individual-record source:", try_source),
  paste("Saved individual trait table:", paths$trait_individual),
  paste("Saved species summary:", paths$trait_summary),
  paste("Trait alias dictionary:", paths$trait_map),
  "The source pipeline removes experimental TRY datasets, removes TRY duplicates, applies its quantitative ErrorRisk gate, retains mapped focal traits, and excludes arid/tropical/missing climate zones before species medians are calculated.",
  "Multiple eligible individual records per species are allowed and are summarized by the median; records are not averaged before that median.",
  table_lines(source_table),
  "Pipeline stage counts:",
  table_lines(source_stage),
  "Documented unit/range gate from the current pipeline:",
  table_lines(range_gate),
  "Country and site fields are not present in the TRY export columns used by this pipeline or in the cleaned target; Dataset, DatasetID, and Koeppen climate_code are available and audited. No country/site was inferred from dataset titles.",
  "",
  "3. SPECIES-LEVEL TABLE",
  "----------------------",
  "trait_species_level_audit.csv has exactly one row for every retained species and includes canonical TNRS identity, family, statuses, current medians, record counts, source datasets, units, aggregation rules, and source-name provenance.",
  "Saved species medians/minima/maxima/counts exactly reproduced from eligible records:",
  as.character(readiness_checks[["saved_medians_reproduced"]]),
  "",
  "4. DUPLICATE AND CONFLICT AUDIT",
  "--------------------------------",
  table_lines(duplicate_summary),
  sprintf("Conflicting species-trait groups written to trait_duplicate_conflict_audit.csv: %d", nrow(conflicts)),
  sprintf("Eligible source observation/trait keys mapping to multiple ecological species: %d", nrow(record_to_multiple_species)),
  sprintf("Retained accepted names reached from >1 TRY submitted name: %d", nrow(multiple_try_names)),
  sprintf("Retained accepted names reached from >1 original ecological name: %d", nrow(multiple_ecological_names)),
  "These many-to-one name routes terminate at one canonical ecological species_id and do not create duplicate species-level trait rows.",
  "Distinct source values are preserved in the conflict audit without choosing a preferred record; the currently declared median aggregation is reported, not reconsidered here.",
  "",
  "5. VALUE INTEGRITY",
  "------------------",
  table_lines(value_summary),
  sprintf("Provided species values that are non-finite/non-positive, or LDMC > 1 g/g: %d", nrow(impossible_species_values)),
  sprintf("Unexpected standardized units among retained eligible records: %d", nrow(unexpected_units)),
  "No arbitrary ecological outlier cutoff was imposed. LMA is checked only for finite positivity because the explicit project unit is the reciprocal of standardized SLA; height is checked for finite positivity in metres.",
  "",
  "6. CANDIDATE EXTREMES (NOT REMOVED)",
  "------------------------------------",
  table_lines(candidate_extreme_summary),
  "Top family concentrations among candidate extreme species:",
  table_lines(extreme_concentration_family),
  "Top TRY dataset concentrations among candidate extreme species:",
  table_lines(extreme_concentration_source),
  "Top climate-code concentrations among records contributing to candidate extreme species:",
  table_lines(extreme_concentration_climate),
  "Raw-record support for candidate extreme species:",
  table_lines(extreme_record_support),
  "All lowest/highest-10, IQR, and robust-MAD candidates are in trait_extreme_values.csv. Shared family/source/climate representation is descriptive and is not evidence of error by itself.",
  "",
  "7. EXACT REPEATED SPECIES MEDIANS / ROUNDING",
  "---------------------------------------------",
  "Ten most frequent exact species medians per trait (percent denominator is trait-complete retained species):",
  table_lines(common_values),
  "Repeated-value spikes meeting the descriptive >=5 species and >=5% screen:",
  if (nrow(possible_rounding_spikes)) table_lines(possible_rounding_spikes) else "None.",
  "Exact repetitions are reported as possible rounding/coarseness signals only. No repeated value was relabeled as imputation or corrected without source evidence.",
  "",
  "8. TRANSFORMATION DIAGNOSTIC",
  "----------------------------",
  table_lines(transformation_summary),
  "trait_distribution_audit.pdf shows raw and log10 species-level distributions. The log10 values are inspection-only.",
  "The previous isolated M_both trait models used the natural logarithm for LDMC, LMA, and vegetative height. This audit makes no new transformation choice.",
  "",
  "9. STANDARDIZATION AUDIT",
  "-------------------------",
  table_lines(standardization_audit),
  "The previous isolated M_both trait models (R/28_step_16_both_trait_models.R via trait_complete_data()) used method A: natural-log trait centers and SDs calculated over distinct trait-complete species, then joined to repeated effects.",
  "The generic z_LDMC, z_LMA, and z_vegetative_height columns saved in analysis_dataset.csv were instead calculated after row expansion. Numerical check (near-zero difference identifies the method):",
  table_lines(legacy_z_audit),
  "Therefore the legacy analysis_dataset.csv z-trait columns must not be reused in the rerun. The intended unique-species scaling is already implemented by trait_complete_data(); current 526-species centers/SDs above are audit calculations only and did not overwrite data.",
  "",
  "10. MATCHING AND CLIMATE/EXPERIMENT COVERAGE",
  "--------------------------------------------",
  table_lines(bind_rows(full_coverage, coverage_summary)),
  "T_z, P_z, and W_z were calculated for this coverage comparison over the 42 distinct experiments (P = log destination JJA precipitation), so species-rich experiments do not control these ranges.",
  "All missing retained species are listed by trait in trait_missingness_and_coverage.csv, which also includes experiment-specific coverage.",
  "",
  "11. FUNCTIONAL AND TAXONOMIC MISSINGNESS",
  "-----------------------------------------",
  "Coverage by graminoid status:",
  table_lines(missingness_coverage |> filter(.data$record_type == "group_missingness", .data$group_variable == "graminoid_status") |> select(.data$trait, .data$group_value, .data$n_species_total, .data$n_species_with_trait, .data$n_species_missing, .data$coverage_percent)),
  "Coverage by woody status:",
  table_lines(missingness_coverage |> filter(.data$record_type == "group_missingness", .data$group_variable == "woody_status") |> select(.data$trait, .data$group_value, .data$n_species_total, .data$n_species_with_trait, .data$n_species_missing, .data$coverage_percent)),
  sprintf("Major family/status groups (>=10 retained species) with zero coverage: %d", nrow(major_group_absence)),
  "Family-level missingness is included in trait_missingness_and_coverage.csv.",
  "",
  "12. CROSS-TRAIT CONSISTENCY",
  "---------------------------",
  table_lines(pairwise_data),
  table_lines(trait_count_summary |> select(.data$category, .data$n_species)),
  "trait_pairwise_audit.pdf shows raw species medians only; no multivariate ecological model was fitted.",
  "",
  "13. TAXONOMIC MATCH AUDIT",
  "-------------------------",
  sprintf("species_id is identical to the canonical TNRS accepted name for all 526 taxa: %s", all(phy$species_id == phy$tnrs_accepted_name)),
  sprintf("Duplicate canonical species_id values: %d", sum(duplicated(phy$species_id))),
  sprintf("Duplicate TNRS accepted names: %d", sum(duplicated(phy$tnrs_accepted_name))),
  sprintf("Trait records mapping to more than one ecological species: %d", nrow(record_to_multiple_species)),
  sprintf("Excluded downstream taxa retained in final matched trait records: %d", n_distinct(intersect(trait_observations_retained$species_id, excluded_species))),
  "No retained species was renamed by this audit. Existing many-to-one synonym/original-name mappings are reported above and retain the established canonical TNRS species_id.",
  "",
  "14. READINESS CHECKS",
  "--------------------",
  table_lines(tibble(check = names(readiness_checks), pass = unname(readiness_checks))),
  "Candidate biological extremes, multiple source values, and incomplete coverage remain visible audit characteristics; they were not treated as automatic integrity failures because units are consistent, values are physically admissible under the documented gates, and the median aggregation rule is explicit.",
  "",
  paste0("TRAIT DATA READY FOR MODELING: ", if (ready) "TRUE" else "FALSE"),
  if (ready) {
    "No unresolved data-integrity problem was found. Use freshly calculated unique-species log-trait z-scores in the rerun; do not reuse the legacy row-weighted z-trait columns in analysis_dataset.csv."
  } else {
    paste("Issues requiring a decision:", paste(names(readiness_checks)[!readiness_checks], collapse = "; "))
  }
)

writeLines(report, file.path(output_dir, "trait_data_audit_report.txt"), useBytes = TRUE)

message("Trait-data audit complete: ", output_dir)
message("TRAIT DATA READY FOR MODELING: ", ready)
