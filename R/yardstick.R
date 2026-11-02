# The R² yardstick (PLAN.md Step D, "Route 2"); used by steps 06 and 07.
# yardstick: gradient-weighted variance of observed strengths = null-model error.
# R2fit: V from all rows. R2pred: V from the fold's training gradients only.
yardstick <- function(d) {
  g <- d |>
    dplyr::group_by(gradient_id) |>
    dplyr::summarise(m = mean(strength), q = mean(strength^2), .groups = "drop")
  mean(g$q) - mean(g$m)^2
}
