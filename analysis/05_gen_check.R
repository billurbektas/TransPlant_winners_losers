# Step 05 side check: the 12 slope fits under GEN, against the diagonal ladder.
#
# Asked by Bills (1 Nov): after the diagonal ladder of analysis/05_ladder_fit.R,
# refit every slope fit with the T, P, TP block GEN, store it apart, and compare
# the slope variances, the fixed effects and what changes in the ladder's reading.
# Not a step of MILESTONES.md; it does not feed any later step. Kept here for
# now, to move to an additional-checks section later (Bills, 2 Nov). Paused on
# 2 Nov with 6 of 12 GEN fits cached; rerunning resumes from them.
# Reads:  data/derived/03b_strength_y3/strength_df.rds;
#         data/derived/05_ladder_fit/models/<rowset>/<slug>.rds (the diagonal fits)
# Writes: data/derived/05_ladder_fit/gen_check/models/<rowset>/<slug>.rds,
#         gen_check/variance_components.csv, gen_check/fixed_effects.csv;
#         figures/diagnostics/05_gen_vs_diag.png; reports/05_gen_vs_diag.md
# Checks: the same rowset counts as step 05; the GEN fits carry the same fixed
#         part as their diagonal twins.

suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(metafor)})
source("R/fit_model.R"); source("R/h_species.R"); source("R/md_table.R")

t_start <- Sys.time()
diag_dir <- "data/derived/05_ladder_fit/models"
out_dir  <- "data/derived/05_ladder_fit/gen_check"
cores <- 8L

strength_df <- readRDS("data/derived/03b_strength_y3/strength_df.rds")
stopifnot(nrow(strength_df) == 1655) # step 03b

# ---- Rowsets and formulas, as in analysis/05_ladder_fit.R ---------------------
rows <- list(all         = strength_df,
             ldmc        = filter(strength_df, !is.na(ldmc_z)),
             lma         = filter(strength_df, !is.na(lma_z)),
             height      = filter(strength_df, !is.na(height_z)),
             growth_form = filter(strength_df, !is.na(graminoid), !is.na(woody)),
             phylogeny   = filter(strength_df, !is.na(PEV1), !is.na(PEV2), !is.na(PEV3)),
             niche       = filter(strength_df, !is.na(T_opt), !is.na(P_opt)))
stopifnot(identical(unname(sapply(rows, nrow)), c(1655L, 1350L, 1386L, 1308L, 1655L, 1655L, 1618L))) # step 05 rowsets
f_both <- ~ W * T * P + A0_origin + A0_dest + A0_origin:W + A0_dest:W
prop_vars <- list(ldmc = "ldmc_z", lma = "lma_z", height = "height_z", graminoid = "graminoid", woody = "woody",
                  phylogeny = c("PEV1", "PEV2", "PEV3"), niche = c("T_opt", "P_opt"))
prop_rowset <- c(ldmc = "ldmc", lma = "lma", height = "height", graminoid = "growth_form", woody = "growth_form",
                 phylogeny = "phylogeny", niche = "niche")
prop_label <- c(ldmc = "LDMC", lma = "LMA", height = "Height", graminoid = "Graminoid", woody = "Woody",
                phylogeny = "Phylogeny", niche = "Niche")
prop_terms <- function(v) unlist(lapply(v, function(x) paste0(x, c("", ":W", ":T", ":P", ":TP"))))
f_prop <- lapply(prop_vars, function(v) reformulate(c(attr(terms(f_both), "term.labels"), prop_terms(v))))
tasks <- c(lapply(c("all", "ldmc", "lma", "height", "niche"), function(r)
             list(rowset = r, slug = "species", model = "Species-dependent", mods = f_both)),
           lapply(names(prop_vars), function(p)
             list(rowset = prop_rowset[[p]], slug = p, model = prop_label[[p]], mods = f_prop[[p]])))

# ---- GEN fits, from metafor's default start, cached; a failure is kept ----------
run_gen <- function(task) {
  path <- file.path(out_dir, "models", task$rowset, paste0(task$slug, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  w <- character(); t0 <- Sys.time()
  res <- withCallingHandlers(
    tryCatch(list(m = fit_model(task$mods, re_slopes, rows[[task$rowset]], struct = c("GEN", "GEN"))),
             error = function(e) list(error = conditionMessage(e))),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  res <- c(res, task[c("rowset", "slug", "model")],
           list(minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")), warnings = list(w)))
  saveRDS(res, path)
  res
}
gen <- parallel::mclapply(tasks, run_gen, mc.cores = cores, mc.preschedule = FALSE)
diag <- lapply(tasks, function(t) { p <- file.path(diag_dir, t$rowset, paste0(t$slug, ".rds"))
  if (file.exists(p)) readRDS(p) else list(error = "diagonal fit missing") })
ok <- function(f) is.list(f) && is.null(f$error) && !is.null(f$m)

# ---- Variance components, side by side ----------------------------------------
vc_row <- function(f, t, struct) {
  base <- tibble::tibble(rowset = t$rowset, model = t$model, struct = struct,
                         status = if (ok(f)) "converged" else if (is.list(f)) f$error else as.character(f),
                         minutes = if (is.list(f) && !is.null(f$minutes)) f$minutes else NA_real_)
  if (!ok(f)) return(base)
  m <- f$m; rho <- if (!is.null(m$phi) && length(m$phi) == 3) m$phi else rep(0, 3)
  mutate(base, logLik_REML = as.numeric(logLik(m)), n_parameters = m$parms,
         sigma2_species = m$sigma2[1], sigma2_experiment = sigma2_exp(m), tau2 = tau2(m), tau2_W = m$tau2,
         var_T = m$gamma2[1], var_P = m$gamma2[2], var_TP = m$gamma2[3],
         rho_T_P = rho[1], rho_T_TP = rho[2], rho_P_TP = rho[3],
         H_slope = H_species(m, rows[[t$rowset]])[["slope"]])
}
vc <- bind_rows(lapply(seq_along(tasks), function(i) bind_rows(vc_row(diag[[i]], tasks[[i]], "DIAG"),
                                                                vc_row(gen[[i]], tasks[[i]], "GEN"))))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(vc, file.path(out_dir, "variance_components.csv"), row.names = FALSE)
pair <- vc |> select(rowset, model, struct, any_of(c("logLik_REML", "var_T", "var_P", "var_TP", "tau2_W", "tau2",
                                                    "sigma2_experiment", "sigma2_species", "H_slope", "rho_T_P", "rho_T_TP", "rho_P_TP"))) |>
  tidyr::pivot_wider(names_from = struct, values_from = -c(rowset, model, struct))
both_ok <- vapply(seq_along(tasks), function(i) ok(gen[[i]]) && ok(diag[[i]]), logical(1))
gen_ok <- vapply(gen, ok, logical(1))
rho_one <- vapply(gen, function(f) ok(f) && max(abs(f$m$phi)) >= 0.999, logical(1)) # rho = 1 (Bills, 1 Nov)

# ---- Fixed effects with CR1, side by side ---------------------------------------
cr1 <- function(m, d) { r <- robust(m, cluster = d$experiment_id, adjust = TRUE)
  tibble::tibble(term = rownames(r$beta), estimate = as.vector(r$beta), se_CR1 = r$se, p_CR1 = r$pval) }
fe <- bind_rows(lapply(which(both_ok), function(i) {
  t <- tasks[[i]]; d <- rows[[t$rowset]]
  stopifnot(identical(colnames(gen[[i]]$m$X), colnames(diag[[i]]$m$X))) # same fixed part
  inner_join(cr1(diag[[i]]$m, d), cr1(gen[[i]]$m, d), by = "term", suffix = c("_DIAG", "_GEN")) |>
    mutate(rowset = t$rowset, model = t$model, .before = 1)
})) |> mutate(shift_in_SE = (estimate_GEN - estimate_DIAG) / se_CR1_DIAG, se_ratio = se_CR1_GEN / se_CR1_DIAG,
              sign_flip = sign(estimate_GEN) != sign(estimate_DIAG), crosses_0.05 = (p_CR1_DIAG < 0.05) != (p_CR1_GEN < 0.05))
write.csv(fe, file.path(out_dir, "fixed_effects.csv"), row.names = FALSE)

# ---- The ladder's readings under each structure --------------------------------
read_ladder <- function(s) {
  v <- filter(vc, struct == s, status == "converged")
  sp <- filter(v, model == "Species-dependent") |> select(rowset, tau2_sp = tau2, sig_sp = sigma2_experiment, Hs_sp = H_slope)
  # growth_form and phylogeny rows equal all rows: their Species-dependent is the all-rows fit (step 05)
  filter(v, model != "Species-dependent") |>
    mutate(same_rows = ifelse(rowset %in% c("growth_form", "phylogeny"), "all", rowset)) |>
    inner_join(sp, by = c(same_rows = "rowset")) |>
    transmute(model, rowset, struct = s, tau2_move = (tau2 - tau2_sp) / (tau2_sp + sig_sp),
              sigma2_experiment_move = (sigma2_experiment - sig_sp) / (tau2_sp + sig_sp),
              Q3d_H_slope_change = H_slope - Hs_sp)
}
ladder_read <- bind_rows(read_ladder("DIAG"), read_ladder("GEN"))

# ---- Figure -----------------------------------------------------------------------
th <- theme_bw(base_size = 18) + theme(plot.title = element_text(face = "bold", size = 24), panel.grid.minor = element_blank(),
                                       strip.text = element_text(face = "bold", size = 15), strip.background = element_blank(),
                                       legend.position = "bottom")
vc_long <- vc |> filter(status == "converged") |>
  select(rowset, model, struct, var_T, var_P, var_TP, tau2_W, tau2, sigma2_experiment, sigma2_species, H_slope) |>
  tidyr::pivot_longer(-c(rowset, model, struct), names_to = "component") |>
  mutate(fit = paste0(model, " (", rowset, ")"),
         component = factor(component, levels = c("var_T", "var_P", "var_TP", "tau2_W", "H_slope", "tau2", "sigma2_experiment", "sigma2_species"),
                            labels = c("T slope variance", "P slope variance", "TP slope variance", "W slope variance",
                                       "H_species slope piece", "Species x site variance (tau2)", "Experiment variance", "Species intercept variance")))
pa <- ggplot(vc_long, aes(1e4 * value, fit, colour = struct, shape = struct)) + geom_point(size = 4, stroke = 1.3) +
  facet_wrap(~ component, ncol = 4, scales = "free_x") + scale_shape_manual(values = c(DIAG = 16, GEN = 4)) +
  scale_colour_manual(values = c(DIAG = "grey30", GEN = "#C8102E")) +
  labs(title = "a)", subtitle = "Variance components of each slope fit, diagonal against GEN", x = "Variance (x 10^-4)", y = NULL,
       colour = "T, P, TP block", shape = "T, P, TP block") + th
pb <- if (nrow(fe)) ggplot(filter(fe, term != "intrcpt"), aes(estimate_DIAG, estimate_GEN)) + geom_abline(colour = "grey60") +
  geom_point(aes(colour = crosses_0.05), size = 3.5, alpha = 0.8) +
  scale_colour_manual(values = c(`FALSE` = "grey30", `TRUE` = "#C8102E"), labels = c(`FALSE` = "same side of p = 0.05", `TRUE` = "crosses p = 0.05"), name = NULL) +
  labs(title = "b)", subtitle = "Fixed effects of every slope fit (CR1)", x = "Estimate, diagonal", y = "Estimate, GEN") + th else
  ggplot() + labs(title = "b)", subtitle = "no pair of fits converged") + theme_void()
pc <- if (nrow(fe)) ggplot(filter(fe, term != "intrcpt"), aes(shift_in_SE, paste0(model, " (", rowset, ")"))) +
  geom_vline(xintercept = 0, colour = "grey60") + geom_point(size = 3, alpha = 0.7) +
  labs(title = "c)", subtitle = "Shift of each fixed effect, GEN minus diagonal, in diagonal CR1 standard errors", x = "Shift (standard errors)", y = NULL) + th else
  ggplot() + theme_void()
ggsave("figures/diagnostics/05_gen_vs_diag.png", patchwork::wrap_plots(pa, patchwork::wrap_plots(pb, pc, widths = c(1, 1.3)), ncol = 1, heights = c(1.3, 1)),
       width = 26, height = 26, dpi = 100, bg = "white")

# ---- Report -----------------------------------------------------------------------
e4 <- function(x) ifelse(is.na(x), "", sprintf("%.3f", 1e4 * x))
f3 <- function(x) ifelse(is.na(x), "", sprintf("%.3f", x))
vc_tab <- vc |> transmute(fit = paste0(model, " (", rowset, ")"), struct, status = substr(status, 1, 60), minutes = sprintf("%.0f", minutes),
                          logLik = ifelse(is.na(logLik_REML), "", sprintf("%.2f", logLik_REML)),
                          `T` = e4(var_T), `P` = e4(var_P), `TP` = e4(var_TP), `ρ(T,P)` = f3(rho_T_P), `ρ(T,TP)` = f3(rho_T_TP), `ρ(P,TP)` = f3(rho_P_TP),
                          W = e4(tau2_W), `τ²` = e4(tau2), `σ²_exp` = e4(sigma2_experiment), `σ²_species` = e4(sigma2_species), `H slope` = e4(H_slope))
lr <- pair |> filter(!is.na(logLik_REML_GEN), !is.na(logLik_REML_DIAG)) |>
  transmute(fit = paste0(model, " (", rowset, ")"), `LR (GEN vs DIAG, 3 df)` = sprintf("%.2f", 2 * (logLik_REML_GEN - logLik_REML_DIAG)),
            `p, naive χ²₃` = sprintf("%.3f", pchisq(pmax(2 * (logLik_REML_GEN - logLik_REML_DIAG), 0), 3, lower.tail = FALSE)),
            `H slope change, GEN − DIAG (×10⁻⁴)` = e4(H_slope_GEN - H_slope_DIAG))
fe_sum <- fe |> filter(term != "intrcpt") |> group_by(fit = paste0(model, " (", rowset, ")")) |>
  summarise(terms = n(), `largest shift (SE)` = sprintf("%.2f", shift_in_SE[which.max(abs(shift_in_SE))]),
            `term with largest shift` = term[which.max(abs(shift_in_SE))],
            `SE ratio GEN/DIAG, range` = sprintf("%.2f to %.2f", min(se_ratio), max(se_ratio)),
            `sign flips` = sum(sign_flip), `cross p = 0.05` = sum(crosses_0.05), .groups = "drop")
crossers <- fe |> filter(crosses_0.05) |> transmute(fit = paste0(model, " (", rowset, ")"), term,
  `estimate DIAG` = sprintf("%.4f", estimate_DIAG), `p DIAG` = sprintf("%.3f", p_CR1_DIAG),
  `estimate GEN` = sprintf("%.4f", estimate_GEN), `p GEN` = sprintf("%.3f", p_CR1_GEN))
pct <- function(x) ifelse(is.na(x), "", sprintf("%+.2f%%", 100 * x))
lad_tab <- if (!nrow(ladder_read)) tibble::tibble() else ladder_read |> tidyr::complete(tidyr::nesting(model, rowset), struct = c("DIAG", "GEN")) |>
  tidyr::pivot_wider(names_from = struct, values_from = c(tau2_move, sigma2_experiment_move, Q3d_H_slope_change)) |>
  transmute(model, rowset, `τ² move DIAG` = pct(tau2_move_DIAG), `τ² move GEN` = pct(tau2_move_GEN),
            `σ²_exp move DIAG` = pct(sigma2_experiment_move_DIAG), `σ²_exp move GEN` = pct(sigma2_experiment_move_GEN),
            `Q3d: H slope, property − Species-dependent, DIAG (×10⁻⁴)` = e4(Q3d_H_slope_change_DIAG),
            `same, GEN (×10⁻⁴)` = e4(Q3d_H_slope_change_GEN))
ft <- filter(fe, term != "intrcpt")
report <- c(
  "# Step 05 side check: GEN against the diagonal T, P, TP block", "",
  paste0("Every slope fit of the step 05 ladder (Species-dependent on the five rowsets and the seven property models) was refitted with ",
         "the T, P, TP slope block GEN, from metafor's default start, and set beside its diagonal fit (asked by Bills, 1 Nov). ",
         "GEN converged in ", sum(gen_ok), " of ", length(gen), " fits; ", sum(rho_one), " of those have a correlation at |ρ| ≥ 0.999 ",
         "(the ρ = 1 rule). Over the ", nrow(ft), " fixed effects of the ", sum(both_ok), " pairs where both converged, the median absolute ",
         "shift is ", if (nrow(ft)) sprintf("%.2f", median(abs(ft$shift_in_SE))) else "NA", " diagonal CR1 standard errors, the largest ",
         if (nrow(ft)) sprintf("%.2f", max(abs(ft$shift_in_SE))) else "NA", "; ", sum(ft$sign_flip), " change sign and ", sum(ft$crosses_0.05),
         " cross p = 0.05. GEN wall time ", sprintf("%.0f", as.numeric(difftime(Sys.time(), t_start, units = "mins"))), " min on ", cores, " cores."),
  "", "## Variance components (×10⁻⁴; ρ under DIAG is 0 by construction)", "", md_table(vc_tab),
  "", "## Likelihood and the slope piece", "",
  "Same fixed part, so the REML logLik are comparable. GEN adds the three correlations; the χ² p is naive (one variance near its bound).", "",
  md_table(lr),
  "", "## Fixed effects (CR1), per fit", "", md_table(fe_sum),
  "", "Terms whose CR1 p crosses 0.05 between the two structures:", "",
  if (nrow(crossers)) md_table(crossers) else "None.",
  "", "## The ladder's reading under each structure", "",
  "The 5% property check (move against τ² + σ²_experiment of Species-dependent on the same rows) and the Q3d change in the H_species slope piece:", "",
  if (nrow(lad_tab)) md_table(lad_tab) else "No property model converged next to its Species-dependent fit.",
  "", "## Figure", "", "![GEN against diagonal](../figures/diagnostics/05_gen_vs_diag.png)",
  "", "## Warnings captured", "",
  { w <- unique(unlist(lapply(gen, function(f) if (is.list(f) && length(f$warnings[[1]])) paste0(f$rowset, " / ", f$model, ": ", f$warnings[[1]]))))
    if (length(w)) paste("-", w) else "None." },
  "", "## Files", "",
  paste0("GEN fits in `", out_dir, "/models/`; tables `variance_components.csv`, `fixed_effects.csv` in `", out_dir, "/`."))
writeLines(report, "reports/05_gen_vs_diag.md")
message("GEN check written: ", sum(gen_ok), "/", length(gen), " GEN fits converged")
