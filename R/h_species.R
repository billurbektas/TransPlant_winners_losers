# Heterogeneity pieces read at every ladder model (PLAN.md Step D, "Route 1");
# used by steps 05, 06 and 07.

tau2       <- function(m) m$sigma2[m$s.names == "effect_id"]       # species x experiment piece
sigma2_exp <- function(m) m$sigma2[m$s.names == "experiment_id"]   # shared site shift

# Variance of the species-specific climate response at the observed experiment
# climates: the W slope variance and the T, P, TP slope covariance, evaluated at
# each experiment and averaged over experiments (legacy
# calculate_response_heterogeneity(), docs/model_analysis_review.R:575)
species_slope_variance <- function(m, data) {
  e <- dplyr::distinct(data, experiment_id, W, T, P, TP)
  X <- as.matrix(e[c("T", "P", "TP")])
  mean(m$tau2 * e$W^2) + mean(rowSums((X %*% m$H) * X))
}

# PLAN.md tests is.null(m$tau2); rma.mv returns tau2 = 0, not NULL, for a model
# without an ~ inner | outer term, so the test here is withG
H_species <- function(m, data) {
  intercept <- m$sigma2[m$s.names == "species"]
  slope     <- if (isTRUE(m$withG)) species_slope_variance(m, data) else 0
  c(intercept = intercept, slope = slope)
}

# The same slope variance per effect, at that effect's experiment climate (for
# the effective weights)
slope_variance_effect <- function(m, data) {
  if (!isTRUE(m$withG)) return(rep(0, nrow(data)))
  X <- as.matrix(data[c("T", "P", "TP")])
  m$tau2 * data$W^2 + rowSums((X %*% m$H) * X)
}
