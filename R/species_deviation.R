# Species deviations for prediction (PLAN.md Step D, "Route 2", P-species); used
# by steps 06 and 07. PLAN.md writes species_deviation(ranef(m), test), but
# metafor 5.0.1 gives no BLUPs for a GEN or GDIAG slope block, so the BLUPs are
# computed by hand, jointly for every random term, from the fitted variance
# components, V and beta (Bills, 2 Nov): b = D Z' M^-1 (y - X beta), with M the
# fit's marginal covariance.
blup <- function(m, data) {
  r <- as.vector(Matrix::solve(m$M, m$yi - m$X %*% m$beta))
  out <- list(species    = tapply(r, data$species, sum) * m$sigma2[m$s.names == "species"],
              experiment = tapply(r, data$experiment_id, sum) * m$sigma2[m$s.names == "experiment_id"])
  if (isTRUE(m$withG)) out$W <- tapply(r * data$W, data$species, sum) * m$tau2
  if (isTRUE(m$withH))
    out$TPTP <- sapply(c("T", "P", "TP"), function(v) tapply(r * data[[v]], data$species, sum)) %*% m$H
  out
}

# Deviation of each test row from the fixed-effect prediction: the species
# intercept plus its W, T, P and TP slopes at the row's climate. The experiment
# intercept is left out (a new site's shift is unknown, expected 0), and a species
# not in `train` gets 0.
species_deviation <- function(m, train, test) {
  b <- blup(m, train)
  k <- match(test$species, names(b$species))
  dev <- b$species[k]
  if (!is.null(b$W)) dev <- dev + test$W * b$W[k]
  if (!is.null(b$TPTP)) dev <- dev + test$T * b$TPTP[k, 1] + test$P * b$TPTP[k, 2] + test$TP * b$TPTP[k, 3]
  unname(ifelse(is.na(k), 0, dev))
}
