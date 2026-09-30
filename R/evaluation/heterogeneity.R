# Module: R/evaluation/heterogeneity
# Responsibility: Implement the named scientific workflow stage or shared infrastructure.
# Inputs: Function arguments and explicitly declared target dependencies.
# Outputs: Validated in-memory objects and declared workflow artifacts.
# Side effects: Limited to functions whose contracts explicitly write, copy, or render files.

# Response-scale species heterogeneity summaries for model synthesis.

#' Response Heterogeneity Metrics.
#'
#' @param fit A fitted model object with the structure checked by this function.
#' @param data A data frame whose required columns are checked before use.
#' @return An R object described by the function name; its structure is validated before later stages use it.
response_heterogeneity_metrics <- function(fit, data) {
  climates <- data |>
    dplyr::distinct(.data$experiment_id, .data$T_z, .data$P_z, .data$W_z) |>
    dplyr::mutate(TP_z = .data$T_z * .data$P_z)
  G <- context_covariance(fit)
  X_context <- as.matrix(climates[c("T_z", "P_z", "TP_z")])
  H_overall <- as.numeric(fit$sigma2[[1]])
  H_W <- mean(climates$W_z^2 * as.numeric(fit$tau2[[1]]))
  H_context <- mean(rowSums((X_context %*% G) * X_context))
  c(
    H_overall = H_overall,
    H_W = H_W,
    H_context = H_context,
    H_total = H_overall + H_W + H_context
  )
}

#' Heterogeneity Resolution.
#'
#' @param baseline A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param predictor A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
heterogeneity_resolution <- function(baseline, predictor) {
  if (!identical(names(baseline), names(predictor))) {
    stop("Heterogeneity components do not match.", call. = FALSE)
  }
  100 * (baseline - predictor) / baseline
}

#' Assemble Heterogeneity Components.
#'
#' @param label A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param data A data frame whose required columns are checked before use.
#' @param baseline A value supplied by the calling workflow; admissible values are checked in the function body.
#' @param predictor A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
assemble_heterogeneity_components <- function(label, data, baseline, predictor) {
  common <- list(
    Predictor = label,
    n_effects = nrow(data),
    n_species = dplyr::n_distinct(data$species_id),
    n_experiments = dplyr::n_distinct(data$experiment_id)
  )
  dplyr::bind_rows(
    c(common, list(model_role = "matched baseline"), as.list(baseline)),
    c(common, list(model_role = "full predictor"), as.list(predictor))
  ) |>
    dplyr::select(
      "Predictor", "model_role", "n_effects", "n_species",
      "n_experiments", dplyr::everything()
    )
}

#' Heterogeneity Comparison Text.
#'
#' @param row A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
heterogeneity_comparison_text <- function(row) {
  sprintf(
    "LR=%.5g; df=%d; p=%.5g",
    row$likelihood_ratio, row$df, row$p_value
  )
}

#' Heterogeneity from matched predictor models.
#'
#' @param model_pair Matched climate-only and full predictor models.
#' @return An R object described by the function name; its structure is validated before later stages use it.
heterogeneity_from_predictor_models <- function(model_pair) {
  baseline <- response_heterogeneity_metrics(
    model_pair$fits$climate_REML$fit, model_pair$data
  )
  predictor <- response_heterogeneity_metrics(
    model_pair$fits$full_REML$fit, model_pair$data
  )
  components <- assemble_heterogeneity_components(
    model_pair$label, model_pair$data, baseline, predictor
  )
  climate_statistics <- model_pair$statistics |>
    dplyr::filter(.data$model == "climate only")
  full_statistics <- model_pair$statistics |>
    dplyr::filter(.data$model == "full predictor-climate")
  resolved <- heterogeneity_resolution(baseline, predictor)
  synthesis <- tibble::tibble(
    Predictor = model_pair$label,
    `N effects` = nrow(model_pair$data),
    `N species` = dplyr::n_distinct(model_pair$data$species_id),
    `N experiments` = dplyr::n_distinct(model_pair$data$experiment_id),
    `H_overall baseline` = baseline[["H_overall"]],
    `H_overall predictor` = predictor[["H_overall"]],
    `Overall heterogeneity resolved (%)` = resolved[["H_overall"]],
    `H_W baseline` = baseline[["H_W"]],
    `H_W predictor` = predictor[["H_W"]],
    `Warming heterogeneity resolved (%)` = resolved[["H_W"]],
    `H_context baseline` = baseline[["H_context"]],
    `H_context predictor` = predictor[["H_context"]],
    `Context heterogeneity resolved (%)` = resolved[["H_context"]],
    `H_total baseline` = baseline[["H_total"]],
    `H_total predictor` = predictor[["H_total"]],
    `Total heterogeneity resolved (%)` = resolved[["H_total"]],
    `ML climate-only vs full result` = heterogeneity_comparison_text(
      model_pair$comparisons[1, ]
    ),
    `Delta AIC full minus climate only` =
      full_statistics$AIC - climate_statistics$AIC,
    `Delta AICc full minus climate only` =
      full_statistics$AICc - climate_statistics$AICc
  )
  list(components = components, synthesis = synthesis)
}

#' Assemble Heterogeneity Synthesis.
#'
#' @param ... Additional arguments forwarded to the delegated operation.
#' @return An R object described by the function name; its structure is validated before later stages use it.
assemble_heterogeneity_synthesis <- function(...) {
  model_pairs <- list(...)
  individual <- lapply(model_pairs, heterogeneity_from_predictor_models)
  result <- list(
    components = dplyr::bind_rows(lapply(individual, `[[`, "components")),
    synthesis = dplyr::bind_rows(lapply(individual, `[[`, "synthesis"))
  )
  if (any(grepl("Family", result$synthesis$Predictor, fixed = TRUE))) {
    stop("Family was incorrectly included in species-information heterogeneity synthesis.",
         call. = FALSE)
  }
  result
}
