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

#' Heterogeneity from Predictor Sequence.
#'
#' @param sequence A value supplied by the calling workflow; admissible values are checked in the function body.
#' @return An R object described by the function name; its structure is validated before later stages use it.
heterogeneity_from_predictor_sequence <- function(sequence) {
  baseline <- response_heterogeneity_metrics(
    sequence$fits$M0_REML$fit, sequence$data
  )
  predictor <- response_heterogeneity_metrics(
    sequence$fits$full_REML$fit, sequence$data
  )
  components <- assemble_heterogeneity_components(
    sequence$label, sequence$data, baseline, predictor
  )
  M0_statistics <- sequence$statistics |>
    dplyr::filter(.data$model == "M0")
  full_name <- sequence$full_name
  full_statistics <- sequence$statistics |>
    dplyr::filter(.data$model == .env$full_name)
  resolved <- heterogeneity_resolution(baseline, predictor)
  synthesis <- tibble::tibble(
    Predictor = sequence$label,
    `N effects` = nrow(sequence$data),
    `N species` = dplyr::n_distinct(sequence$data$species_id),
    `N experiments` = dplyr::n_distinct(sequence$data$experiment_id),
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
    `ML M0 vs M1 result` = heterogeneity_comparison_text(
      sequence$comparisons[1, ]
    ),
    `ML M1 vs full/context result` = heterogeneity_comparison_text(
      sequence$comparisons[2, ]
    ),
    `Delta AIC` = full_statistics$AIC - M0_statistics$AIC,
    `Delta AICc` = full_statistics$AICc - M0_statistics$AICc
  )
  list(components = components, synthesis = synthesis)
}

#' Assemble Heterogeneity Synthesis.
#'
#' @param ... Additional arguments forwarded to the delegated operation.
#' @return An R object described by the function name; its structure is validated before later stages use it.
assemble_heterogeneity_synthesis <- function(...) {
  sequences <- list(...)
  individual <- lapply(sequences, heterogeneity_from_predictor_sequence)
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
