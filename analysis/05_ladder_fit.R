# Step 05: the ladder, fitted (Step D, route 1, residents).
#
# PLAN.md: Step D, "The six models", "The formulas, written out", "Route 1".
# Reads:  data/derived/03b_strength_y3/strength_df.rds (strength from year 3);
#         data/derived/04_classification/species_summary.csv (Step C spread)
# Writes: data/derived/05_ladder_fit/models/<rowset>/<model>.rds,
#         heterogeneity.csv, weights.csv, coefficients.csv, tableS4.csv,
#         tableS5.csv; figures/appendix/figS11.png; reports/05_ladder_fit.md
# Checks: 1,004 effects in, 512 ladder rows (71 species, 39 experiments, 20
#         gradients): residents in at least 3 gradients that carry every
#         property; every property rowset identical to all;
#         random structure exactly re_base / re_slopes; property terms exactly
#         X, X:W, X:T, X:P, X:TP; A0 interactions only A0:W; tau2 at Climate
#         within 5% of Null; tau2 and sigma2_experiment at each property model
#         within 5% of tau2 + sigma2_experiment of Species-dependent on the same
#         rows; sigma2_species at Null off the boundary; all 8 slope fits with
#         the diagonal T, P, TP block.
#
# Decided 1 Nov (Bills): the T, P, TP slope block is diagonal in every slope fit
# (the 1 Nov GEN timing run did not converge from the gate fit's own optimum;
# MILESTONES step 05 notes), every fit from metafor's default start; no ML refits
# (AICc dropped), no profile CIs, no all-years sensitivity. Decided 2 Nov
# (Bills): the ladder runs on residents in at least 3 gradients, and the species
# lacking a property (Carex norvegica, Veronica serpyllifolia, Bistorta
# macrophylla) leave the whole ladder, so every model is fitted on the same rows;
# step 02's standardisation is kept; the GEN timing run is not repeated; 7 cores.
# Fig 4 is drawn by step 06 and Fig S8 is not drawn here (Bills, 2 Nov: both deleted).

suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(metafor)})
source("R/fit_model.R"); source("R/h_species.R"); source("R/md_table.R")

warns <- character()
globalCallingHandlers(warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
})
t_start <- Sys.time()
out_dir <- "data/derived/05_ladder_fit"
cores <- 7L

strength_df <- readRDS("data/derived/03b_strength_y3/strength_df.rds")
stopifnot(nrow(strength_df) == 1004) # step 03b

# ---- Ladder rows: residents in at least 3 gradients with every property (2 Nov)
props <- c("ldmc_z", "lma_z", "height_z", "graminoid", "woody", "PEV1", "PEV2", "PEV3", "T_opt", "P_opt")
g3 <- summarise(strength_df, n_g = n_distinct(gradient_id), .by = species) |> filter(n_g >= 3)
strength_df <- filter(strength_df, species %in% g3$species) |> filter(if_all(all_of(props), ~ !is.na(.x)))
stopifnot(nrow(strength_df) == 512, n_distinct(strength_df$species) == 71,
          n_distinct(strength_df$experiment_id) == 39, n_distinct(strength_df$gradient_id) == 20) # ladder rows

# ---- Rowsets: each property's complete-case rows, all equal to the ladder rows -
rows <- list(all         = strength_df,
             ldmc        = filter(strength_df, !is.na(ldmc_z)),
             lma         = filter(strength_df, !is.na(lma_z)),
             height      = filter(strength_df, !is.na(height_z)),
             growth_form = filter(strength_df, !is.na(graminoid), !is.na(woody)),
             phylogeny   = filter(strength_df, !is.na(PEV1), !is.na(PEV2), !is.na(PEV3)),
             niche       = filter(strength_df, !is.na(T_opt), !is.na(P_opt)))
counts <- tibble::tibble(rowset = names(rows), effects = sapply(rows, nrow),
                         species = sapply(rows, function(x) n_distinct(x$species)),
                         experiments = sapply(rows, function(x) n_distinct(x$experiment_id)))
stopifnot(all(counts$effects == 512), all(counts$experiments == 39)) # complete cases per property
stopifnot(all(sapply(rows, function(r) identical(r$effect_id, rows$all$effect_id)))) # one set of rows: reuse the all-rows ladder
# subspecies niche check (26 Oct): neither species is a resident at CH_Lavey
stopifnot(!any(grepl("^CH_Lavey", strength_df$experiment_id) &
               strength_df$species %in% c("Solidago virgaurea", "Ranunculus acris")))

# ---- Formulas, verbatim from PLAN.md ---------------------------------------
f_env   <- ~ W * T * P
f_start <- ~ A0_origin + A0_dest + A0_origin:W + A0_dest:W
f_both  <- ~ W * T * P + A0_origin + A0_dest + A0_origin:W + A0_dest:W
prop_vars <- list(ldmc = "ldmc_z", lma = "lma_z", height = "height_z", graminoid = "graminoid", woody = "woody",
                  phylogeny = c("PEV1", "PEV2", "PEV3"), niche = c("T_opt", "P_opt"))
prop_rowset <- c(ldmc = "ldmc", lma = "lma", height = "height", graminoid = "growth_form", woody = "growth_form",
                 phylogeny = "phylogeny", niche = "niche")
prop_label <- c(ldmc = "LDMC", lma = "LMA", height = "Height", graminoid = "Graminoid", woody = "Woody",
                phylogeny = "Phylogeny", niche = "Niche")
# X + X:W + X:T + X:P + X:TP per variable of the property
prop_terms <- function(v) unlist(lapply(v, function(x) paste0(x, c("", ":W", ":T", ":P", ":TP"))))
f_prop <- lapply(prop_vars, function(v) reformulate(c(attr(terms(f_both), "term.labels"), prop_terms(v))))
floor_models <- list(null = list("Null", ~ 1), start = list("Start", f_start),
                     climate = list("Climate", f_env), climate_start = list("Climate + start", f_both))
fit_rowsets <- "all" # every rowset equals the ladder rows (2 Nov)

# ---- One fit, cached: skipped when its file exists ----------------------------
run_fit <- function(task) {
  path <- file.path(out_dir, "models", task$rowset, paste0(task$slug, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  w <- character(); t0 <- Sys.time()
  struct <- if (task$kind == "floor") "GEN" else slope_struct
  res <- withCallingHandlers(
    tryCatch(list(m = fit_model(task$mods, if (task$kind == "floor") re_base else re_slopes, rows[[task$rowset]],
                                struct = struct, control = task$control)),
             error = function(e) list(error = conditionMessage(e))),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  res <- c(res, list(rowset = task$rowset, slug = task$slug, model = task$model, slopes = task$kind != "floor",
                     struct = if (task$kind == "floor") "none (no slopes)" else "DIAG",
                     minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")), warnings = list(w)))
  if (is.null(res$error)) saveRDS(res, path)
  res
}
run_wave <- function(tasks) {
  out <- parallel::mclapply(tasks, run_fit, mc.cores = cores, mc.preschedule = FALSE)
  failed <- vapply(seq_along(out), function(i) inherits(out[[i]], "try-error") ||
                     !is.null(out[[i]]$error), logical(1))
  if (any(failed)) stop("fits failed: ", paste(vapply(tasks[failed], function(t) paste(t$rowset, t$slug), ""), collapse = ", "),
                        "\n", paste(vapply(out[failed], function(r) as.character(if (is.list(r)) r$error else r), ""), collapse = "\n"))
  out
}
task <- function(rowset, slug, model, mods, kind, control = list())
  list(rowset = rowset, slug = slug, model = model, mods = mods, kind = kind, control = control)

# ---- All fits in one wave: 4 floor fits and 8 diagonal slope fits, each from
# metafor's default start --------------------------------------------------------
wave <- c(lapply(fit_rowsets, function(r) task(r, "species", "Species-dependent", f_both, "slopes")),
          lapply(names(prop_vars), function(p) task(prop_rowset[[p]], p, prop_label[[p]], f_prop[[p]], "slopes")),
          unlist(lapply(fit_rowsets, function(r) lapply(names(floor_models), function(s)
            task(r, s, floor_models[[s]][[1]], floor_models[[s]][[2]], "floor"))), recursive = FALSE))
fits <- run_wave(wave)
stopifnot(length(fits) == 12) # 4 floor, Species-dependent, 7 property models
names(fits) <- vapply(fits, function(f) paste(f$rowset, f$slug, sep = "/"), "")
fit_of <- function(rowset, slug) fits[[paste(if (!slug %in% names(prop_vars)) "all" else rowset, slug, sep = "/")]]

# ---- The ladder table: one row per rowset x model -------------------------
ladder <- bind_rows(lapply(names(rows), function(r) {
  slugs <- c(names(floor_models), "species", names(prop_rowset)[prop_rowset == r])
  bind_rows(lapply(slugs, function(s) {
    f <- fit_of(r, s); m <- f$m; H <- H_species(m, rows[[r]])
    rho <- if (isTRUE(m$withH) && !is.null(m$phi)) m$phi else rep(NA_real_, 3)
    tibble::tibble(rowset = r, model = f$model, slug = s, fitted_on = f$rowset,
                   n_effects = m$k, n_species = n_distinct(rows[[r]]$species), n_experiments = n_distinct(rows[[r]]$experiment_id),
                   tau2 = tau2(m), sigma2_experiment = sigma2_exp(m), sigma2_species = H[["intercept"]],
                   H_intercept = H[["intercept"]], H_slope = H[["slope"]],
                   tau2_W = if (isTRUE(m$withG)) m$tau2 else NA_real_,
                   var_T = if (isTRUE(m$withH)) m$gamma2[1] else NA_real_, var_P = if (isTRUE(m$withH)) m$gamma2[2] else NA_real_,
                   var_TP = if (isTRUE(m$withH)) m$gamma2[3] else NA_real_,
                   rho_T_P = rho[1], rho_T_TP = rho[2], rho_P_TP = rho[3],
                   logLik_REML = as.numeric(logLik(m)), n_parameters = m$parms,
                   convergence = "converged", struct = if (f$slopes) f$struct else "none (no slopes)",
                   minutes = f$minutes)
  }))
}))
stopifnot(nrow(ladder) == 7 * 5 + 7) # five ladder models per rowset plus the seven property models
write.csv(select(ladder, rowset, model, fitted_on, n_effects, n_species, n_experiments, tau2, sigma2_experiment, sigma2_species,
                 H_intercept, H_slope, convergence, struct, minutes),
          file.path(out_dir, "heterogeneity.csv"), row.names = FALSE)
write.csv(select(ladder, rowset, model, fitted_on, n_effects, n_species, n_experiments, n_parameters, logLik_REML,
                 sigma2_species, sigma2_experiment, tau2, tau2_W, var_T, var_P, var_TP, rho_T_P, rho_T_TP, rho_P_TP,
                 H_intercept, H_slope, struct),
          file.path(out_dir, "tableS4.csv"), row.names = FALSE)

# ---- Checks -----------------------------------------------------------------
fitted <- filter(ladder, fitted_on == rowset)
for (f in fits) { # random structure exactly re_base / re_slopes
  stopifnot(identical(unname(f$m$s.names), c("species", "experiment_id", "effect_id")))
  if (f$slopes) stopifnot(identical(unname(f$m$g.names), c("W", "species")), identical(unname(f$m$h.names), c("T", "P", "TP", "species")))
  else stopifnot(!isTRUE(f$m$withG), !isTRUE(f$m$withH))
}
canon <- function(x) vapply(strsplit(x, ":"), function(p) paste(sort(p), collapse = ":"), "")
for (p in names(prop_vars)) { # property terms exactly X, X:W, X:T, X:P, X:TP; no three-way property term
  added <- setdiff(colnames(fit_of(prop_rowset[[p]], p)$m$X), colnames(fit_of(prop_rowset[[p]], "species")$m$X))
  stopifnot(setequal(canon(added), canon(prop_terms(prop_vars[[p]]))))
}
a0 <- grep("A0", colnames(fit_of("all", "climate_start")$m$X), value = TRUE) # A0 interactions limited to A0:W
stopifnot(setequal(canon(a0), canon(c("A0_origin", "A0_dest", "A0_origin:W", "A0_dest:W"))))
climate_check <- fitted |> filter(model %in% c("Null", "Climate")) |>
  select(rowset, model, tau2) |> tidyr::pivot_wider(names_from = model, values_from = tau2) |>
  mutate(change = Climate / Null - 1)
print(climate_check)
stopifnot(all(abs(climate_check$change) <= 0.05)) # climate's main effects cannot reach tau2
prop_check <- bind_rows(lapply(names(prop_vars), function(p) {
  r <- prop_rowset[[p]]; m4 <- fit_of(r, p)$m; m3 <- fit_of(r, "species")$m
  denom <- tau2(m3) + sigma2_exp(m3)
  tibble::tibble(property = prop_label[[p]], rowset = r, base = denom,
                 move_tau2 = (tau2(m4) - tau2(m3)) / denom, move_sigma2_experiment = (sigma2_exp(m4) - sigma2_exp(m3)) / denom)
}))
print(prop_check)
stopifnot(all(abs(prop_check$move_tau2) <= 0.05), all(abs(prop_check$move_sigma2_experiment) <= 0.05)) # a property reaches only H_species
species_summary <- read.csv("data/derived/04_classification/species_summary.csv")
var_mu <- var(species_summary$mu)
stopifnot(var_mu > 0, all(filter(fitted, model == "Null")$sigma2_species > 1e-8)) # no collapse of the species intercept
slope_rows <- filter(fitted, struct == "DIAG")
stopifnot(nrow(slope_rows) == 8, all(vapply(fits, function(f) !f$slopes || identical(f$m$struct, slope_struct), logical(1))))

# ---- Fixed effects, model-based (REML Wald) intervals (Table S5, Fig S4) -------
# model-based rather than CR1 clustered on experiment (Bills, 2 Nov)
coefs <- bind_rows(lapply(fits, function(f) {
  if (f$slug == "null") return(NULL)
  s <- coef(summary(f$m))
  tibble::tibble(rowset = f$rowset, model = f$model, term = rownames(s), estimate = s$estimate, se = s$se,
                 ci_lower = s$ci.lb, ci_upper = s$ci.ub, zval = s$zval, pval = s$pval)
}))
write.csv(select(coefs, rowset, model, term, estimate, se, ci_lower, ci_upper), file.path(out_dir, "coefficients.csv"), row.names = FALSE)
write.csv(coefs, file.path(out_dir, "tableS5.csv"), row.names = FALSE)

# ---- Effective weights 1/(V + tau2 + sigma2_experiment + sigma2_species [+ slopes]) --
weights <- bind_rows(lapply(fits, function(f) {
  d <- rows[[f$rowset]]; m <- f$m
  tibble::tibble(rowset = f$rowset, model = f$model, effect_id = d$effect_id, experiment_id = d$experiment_id,
                 A0_origin = d$A0_origin, V = d$strength_se^2,
                 weight = 1 / (d$strength_se^2 + sum(m$sigma2) + slope_variance_effect(m, d)))
}))
write.csv(weights, file.path(out_dir, "weights.csv"), row.names = FALSE)
weight_summary <- function(w, a0, exp_id, label) {
  low <- a0 <= median(a0)
  tibble::tibble(weights = label, `max/min within experiment (median)` = sprintf("%.1f", median(tapply(w, exp_id, function(x) max(x) / min(x)))),
                 `share, low-A0_origin half` = sprintf("%.2f", sum(w[low]) / sum(w)),
                 `share, top 10% of effects` = sprintf("%.2f", sum(sort(w, decreasing = TRUE)[seq_len(round(0.1 * length(w)))]) / sum(w)))
}
w_all <- filter(weights, rowset == "all")
weights_table <- bind_rows(
  weight_summary(1 / rows$all$strength_se^2, rows$all$A0_origin, rows$all$experiment_id, "raw 1/V"),
  bind_rows(lapply(c("Null", "Start", "Climate", "Climate + start", "Species-dependent"), function(mo) {
    x <- filter(w_all, model == mo); weight_summary(x$weight, x$A0_origin, x$experiment_id, mo) })))

# ---- Fig S4: coefficients with model-based intervals ------------------------
coef_panel <- function(x, tag, ttl) {
  x <- mutate(x, excl = ci_lower > 0 | ci_upper < 0, term = factor(term, levels = rev(unique(term))))
  ggplot(x, aes(estimate, term)) + geom_vline(xintercept = 0, colour = "grey60") +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0, linewidth = 0.9) +
    geom_point(aes(fill = excl), shape = 21, size = 4, stroke = 1) +
    scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "white"), labels = c(`TRUE` = "95% interval excludes zero", `FALSE` = "includes zero"), name = NULL) +
    labs(title = tag, subtitle = ttl, x = "Effect on losing-winning strength (per SD)", y = NULL) +
    theme_classic(base_size = 17) + theme(plot.title = element_text(face = "bold", size = 22), legend.position = "bottom",
                                          strip.text = element_text(face = "bold", size = 15), strip.background = element_blank())
}
s4a <- coef_panel(filter(coefs, rowset == "all", model == "Climate + start", term != "intrcpt"), "a)", "Residents, Climate + start (model-based intervals)")
prop_coefs <- bind_rows(lapply(names(prop_vars), function(p) filter(coefs, model == prop_label[[p]]) |>
  filter(canon(term) %in% canon(prop_terms(prop_vars[[p]]))) |> mutate(property = prop_label[[p]])))
s4c <- coef_panel(prop_coefs, "b)", "Property terms, each from its own model (Climate + start + property, species slopes kept)") +
  facet_wrap(~ property, scales = "free", ncol = 4)
ggsave("figures/appendix/figS11.png", patchwork::wrap_plots(patchwork::wrap_plots(s4a, patchwork::plot_spacer(), widths = c(1.4, 1)), s4c, ncol = 1, heights = c(1, 2)), # colonizers: Fig S15
       width = 20, height = 22, dpi = 100, bg = "white")


# ---- Report ----------------------------------------------------------------
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
e4 <- function(x) ifelse(is.na(x), "", sprintf("%.3f", 1e4 * x))
# step 02's standardisation, read on the ladder rows (kept, Bills 2 Nov)
exp_cov <- distinct(strength_df, experiment_id, W, T, P)
sp_prop <- distinct(strength_df, species, across(all_of(c("ldmc_z", "lma_z", "height_z", "PEV1", "PEV2", "PEV3", "T_opt", "P_opt"))))
drift <- bind_rows(
  tibble::tibble(variable = c("W", "T", "P"), across = "the 39 experiments",
                 mean = sapply(exp_cov[c("W", "T", "P")], mean), sd = sapply(exp_cov[c("W", "T", "P")], sd)),
  tibble::tibble(variable = c("A0_origin", "A0_dest"), across = "the 512 effects",
                 mean = sapply(strength_df[c("A0_origin", "A0_dest")], mean), sd = sapply(strength_df[c("A0_origin", "A0_dest")], sd)),
  tibble::tibble(variable = names(sp_prop)[-1], across = "the 71 species",
                 mean = sapply(sp_prop[-1], mean), sd = sapply(sp_prop[-1], sd))) |>
  mutate(mean = sprintf("%.2f", mean), sd = sprintf("%.2f", sd))
het_tab <- ladder |> transmute(rowset, model, fitted_on, effects = fmt(n_effects), species = fmt(n_species),
  `τ²` = e4(tau2), `σ²_experiment` = e4(sigma2_experiment), `H intercept` = e4(H_intercept), `H slope` = e4(H_slope),
  struct, minutes = ifelse(fitted_on == rowset, sprintf("%.1f", minutes), "reused"))
fit_code <- c("```r",
  "re_base   <- list(~ 1 | species, ~ 1 | experiment_id, ~ 1 | effect_id)",
  "re_slopes <- c(re_base, list(~ 0 + W | species, ~ 0 + T + P + TP | species))",
  "rma.mv(yi = strength, V = strength_se^2, mods = <formula>, random = <re_base or re_slopes>,",
  "       struct = c(\"GEN\", \"GDIAG\"), data = <rowset>, method = \"REML\", sparse = TRUE)",
  "# struct only acts on re_slopes: W slope; T, P, TP slopes diagonal. No starting values: metafor's default for every fit",
  paste0("# Null ~ 1; Start ", deparse(f_start), "; Climate ", deparse(f_env)),
  paste0("# Climate + start and Species-dependent ", paste(deparse(f_both), collapse = "")),
  "# property models: Climate + start + X + X:W + X:T + X:P + X:TP per variable, with re_slopes",
  "coef(summary(m))   # model-based (REML Wald) intervals for Table S5 and Fig S4", "```")
report <- c(
  "# Step 05: The ladder, fitted (Step D, route 1, residents)", "",
  paste0("The model ladder was fitted on the strength from year 3 (step 03b), on one set of rows (2 Nov): the ",
         fmt(nrow(strength_df)), " effects of the ", n_distinct(strength_df$species), " residents found in at least 3 gradients ",
         "that carry every property (", n_distinct(strength_df$experiment_id), " experiments, ", n_distinct(strength_df$gradient_id),
         " gradients). Null, Start, Climate, Climate + start, Species-dependent and one property model per property (LDMC, LMA, ",
         "height, graminoid, woody, phylogeny as PEV1 to PEV3, niche as T_opt and P_opt) with the species slopes kept, ",
         nrow(fitted), " REML fits, all on the same rows. Each model reports τ² (species × site), σ²_experiment (shared site ",
         "shift) and the two H_species pieces as raw values. The T, P, TP slope block is diagonal in every slope fit (Bills, 1 Nov). ",
         "Fixed effects carry model-based (REML Wald) standard errors (Bills, 2 Nov; CR1 before). Fits ran in parallel on ", cores, " cores, ",
         sprintf("%.1f", sum(fitted$minutes)), " fit-minutes in all; this run took ",
         sprintf("%.0f", as.numeric(difftime(Sys.time(), t_start, units = "mins"))), " min."),
  "", "## Slope structure", "",
  paste0("The 1 Nov GEN timing run (Species-dependent on the 1,655 rows of the 28 Oct pools, started from the gate fit's own ",
         "optimum) did not converge, so the T, P, TP slope block is diagonal in every slope fit; the case is in MILESTONES step 05 ",
         "notes. It is not repeated on the new rows (Bills, 2 Nov). The experiments' T and P are correlated at r = ",
         sprintf("%.2f", cor(exp_cov$T, exp_cov$P)), " on the ladder rows."),
  "", "## Model as fitted (for reviewers)", "", fit_code,
  "", "## Counts per rowset", "", md_table(mutate(counts, across(-rowset, fmt))),
  "", "Each matches the expectation checked in the script. Every property rowset equals the ladder rows, so Null to Species-dependent are fitted once, on `all` (marked `reused` elsewhere).",
  "", "Step 02's standardisation, kept (Bills, 2 Nov), read on the ladder rows:", "", md_table(drift),
  "", "## Heterogeneity per model (×10⁻⁴)", "", md_table(het_tab),
  "", paste0("All ", nrow(slope_rows), " slope fits converged with the diagonal T, P, TP block, each from metafor's default start ",
             "(Bills, 1 Nov). Across the ", nrow(slope_rows), " fits the T slope variance ranges ",
             e4(min(slope_rows$var_T)), " to ", e4(max(slope_rows$var_T)), " ×10⁻⁴; the P and TP slope variances are at most ",
             sprintf("%.1e", max(c(slope_rows$var_P, slope_rows$var_TP))), ", at the boundary."),
  "", "## Checks", "",
  "τ² at Climate against Null (climate's main effects cannot reach the species × site piece; stop beyond 5%):", "",
  md_table(transmute(climate_check, rowset, `τ² Null` = e4(Null), `τ² Climate` = e4(Climate), change = sprintf("%+.1f%%", 100 * change))),
  "", "Each property model against Species-dependent on the same rows, as a share of that model's τ² + σ²_experiment (a property reaches only H_species; stop beyond 5%, Bills 1 Nov):", "",
  md_table(transmute(prop_check, property, rowset, `τ² + σ²_experiment (×10⁻⁴)` = e4(base),
                     `τ² move` = sprintf("%+.2f%%", 100 * move_tau2), `σ²_experiment move` = sprintf("%+.2f%%", 100 * move_sigma2_experiment))),
  "", paste0("σ²_species at Null on all rows is ", e4(filter(fitted, rowset == "all", model == "Null")$sigma2_species),
             " ×10⁻⁴, off the boundary; the Step C spread of μ_s (step 04, ", nrow(species_summary), " species) is ", e4(var_mu),
             " ×10⁻⁴. Random structure, property terms (X, X:W, X:T, X:P, X:TP) and A0 interactions (A0:W only) are as PLAN.md writes them. ",
             "Subspecies niche check (26 Oct): neither Solidago virgaurea nor Ranunculus acris is a resident at CH_Lavey under the original pools, so there is nothing to refit."),
  "", "## Effective weights (all rows)", "", md_table(weights_table),
  "", "The low-A0_origin half is the effects at or below the median A0_origin. Species-dependent weights include the slope variance at each effect's climate.",
  "", "## Figures", "",
  "![Fig S4, coefficients](../figures/appendix/figS11.png)", "", "Fig 4 is drawn by step 06; Fig S8 is not drawn (Bills, 2 Nov).",
  "", "## Warnings captured", "",
  { w <- unique(c(unlist(lapply(fits, function(f) paste0(f$rowset, " / ", f$model, ": ", f$warnings[[1]])[length(f$warnings[[1]]) > 0])), warns))
    if (length(w)) paste("-", w) else "None." },
  "", "## Deviations from PLAN.md", "",
  paste0("- H_species tests `m$withG` where PLAN.md tests `is.null(m$tau2)`: rma.mv returns tau2 = 0, not NULL, for a model without slopes. ",
         "Same quantity."),
  "- Fig S4 is written to `figures/appendix/figS11.png`, not to `data/derived/05_ladder_fit/`, so that the report can link it.",
  "- Decided 1 and 2 Nov (Bills) and in PLAN.md, not deviations: the T, P, TP slope block diagonal; no ML refits and no AICc column; no profile CIs; the 5% property check measured against τ² + σ²_experiment; Fig 4 with twelve columns, drawn by step 06; the ladder on residents in at least 3 gradients that carry every property, one set of rows for every model; step 02's standardisation kept.")
writeLines(report, "reports/05_ladder_fit.md")
message("Step 05 written: ", nrow(fitted), " fits in ", sprintf("%.1f", sum(fitted$minutes)), " fit-minutes")
