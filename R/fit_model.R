# The ladder's meta-regression (PLAN.md Step D, "The six models"); used by steps
# 05, 06 and 07.

re_base   <- list(~ 1 | species, ~ 1 | experiment_id, ~ 1 | effect_id)
re_slopes <- c(re_base, list(~ 0 + W | species,
                             ~ 0 + T + P + TP | species))

# As PLAN.md writes it, plus `control`. Every ladder fit starts from metafor's
# default (Bills, 1 Nov): carried-over starting values put the boundary P and TP
# slope variances at 1e-10, where nlminb stalled on two fits. `control` carries
# only the step 05 timing run's start from the gate fit.
fit_model <- function(mods, random, data, struct = "GEN", method = "REML", control = list()) {
  metafor::rma.mv(yi = strength, V = strength_se^2, mods = mods, random = random,
                  struct = struct, data = data, method = method, sparse = TRUE, control = control)
}

# Species-dependent and property models: the T, P, TP slope block is diagonal
# ladder-wide (Bills, 1 Nov). The step 05 timing run, GEN started from the gate
# fit's own optimum (rho(T, P) = 0.989), did not converge, so the correlations are
# taken as unidentified, the fallback of decision 25 Oct. metafor spells a
# diagonal block of several inner variables "GDIAG"; the W block has one variable.
slope_struct <- c("GEN", "GDIAG")
