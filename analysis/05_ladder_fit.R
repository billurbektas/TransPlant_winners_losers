# Step 05: the ladder, fitted (Step D, route 1, residents).
#
# PLAN.md: Step D, "The six models", "The formulas, written out", "Route 1".
# Reads:  data/derived/03b_strength_y3/strength_df.rds (strength from year 3);
#         data/derived/05_ladder_fit/gate_start_values.csv (variance components
#         of the 1 Nov gate fit of Species-dependent, the timing run's start);
#         data/derived/04_classification/species_summary.csv (Step C spread)
# Writes: data/derived/05_ladder_fit/models/<rowset>/<model>.rds,
#         heterogeneity.csv, weights.csv, coefficients.csv, tableS4.csv,
#         tableS5.csv; figures/diagnostics/05_fig4_ladder.png, 05_figS4.png,
#         05_figS8.png; reports/05_ladder_fit.md
# Checks: rowsets all 1,655 / ldmc 1,350 / lma 1,386 / height 1,308 /
#         growth_form 1,655 / phylogeny 1,655 / niche 1,618 effects, 43
#         experiments each; growth_form and phylogeny rows identical to all;
#         random structure exactly re_base / re_slopes; property terms exactly
#         X, X:W, X:T, X:P, X:TP; A0 interactions only A0:W; tau2 at Climate
#         within 5% of Null; tau2 and sigma2_experiment at each property model
#         within 5% of tau2 + sigma2_experiment of Species-dependent on the same
#         rows; sigma2_species at Null off the boundary; all 12 slope fits with
#         the diagonal T, P, TP block.
#
# Decided 1 Nov (Bills): the step opens with a timing run, Species-dependent on
# all rows under GEN, started from the gate fit, reporting time, convergence and
# rho ("stability across refits"). It did not converge (nlminb, 2.5 min) from the
# gate fit's own optimum, so the T, P, TP correlations are taken as unidentified
# and the slope block is diagonal in every fit (decision 25 Oct fallback, applied
# ladder-wide; the GEN attempt is kept and reported). Every other fit starts
# from metafor's default: two diagonal fits started from the all-rows fit did
# not converge, and both did from the default (Bills, 1 Nov). No ML refits (AICc
# dropped, decision-log row 7), no profile CIs, no all-years sensitivity.
# growth_form and phylogeny rows equal all rows, so those rowsets reuse the
# all-rows Null to Species-dependent fits. Fits run in parallel on 8 cores.

suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(metafor)})
source("R/fit_model.R"); source("R/h_species.R"); source("R/fig4_ladder.R"); source("R/md_table.R")

warns <- character()
globalCallingHandlers(warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
})
t_start <- Sys.time()
out_dir <- "data/derived/05_ladder_fit"
cores <- 8L

strength_df <- readRDS("data/derived/03b_strength_y3/strength_df.rds")
stopifnot(nrow(strength_df) == 1655) # step 03b

# ---- Rowsets: each property's complete-case rows (PLAN.md Step D) ----------
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
stopifnot(identical(unname(counts$effects), c(1655L, 1350L, 1386L, 1308L, 1655L, 1655L, 1618L)),
          all(counts$experiments == 43)) # complete cases per property
stopifnot(identical(rows$growth_form$effect_id, rows$all$effect_id),
          identical(rows$phylogeny$effect_id, rows$all$effect_id)) # reuse the all-rows ladder
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
fit_rowsets <- c("all", "ldmc", "lma", "height", "niche")

# ---- One fit, cached: skipped when its file exists ----------------------------
run_fit <- function(task) {
  path <- file.path(out_dir, "models", task$rowset, paste0(task$slug, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  w <- character(); t0 <- Sys.time()
  struct <- switch(task$kind, floor = "GEN", slopes = slope_struct, gen_timing = c("GEN", "GEN"))
  res <- withCallingHandlers(
    tryCatch(list(m = fit_model(task$mods, if (task$kind == "floor") re_base else re_slopes, rows[[task$rowset]],
                                struct = struct, control = task$control)),
             error = function(e) list(error = conditionMessage(e))),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  res <- c(res, list(rowset = task$rowset, slug = task$slug, model = task$model, slopes = task$kind != "floor",
                     struct = switch(task$kind, floor = "none (no slopes)", slopes = "DIAG", gen_timing = "GEN"),
                     minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")), warnings = list(w)))
  # the GEN timing run is kept whether or not it converged: it is the reason for DIAG
  if (is.null(res$error) || task$kind == "gen_timing") saveRDS(res, path)
  res
}
run_wave <- function(tasks) {
  out <- parallel::mclapply(tasks, run_fit, mc.cores = cores, mc.preschedule = FALSE)
  failed <- vapply(seq_along(out), function(i) inherits(out[[i]], "try-error") ||
                     (!is.null(out[[i]]$error) && tasks[[i]]$kind != "gen_timing"), logical(1))
  if (any(failed)) stop("fits failed: ", paste(vapply(tasks[failed], function(t) paste(t$rowset, t$slug), ""), collapse = ", "),
                        "\n", paste(vapply(out[failed], function(r) as.character(if (is.list(r)) r$error else r), ""), collapse = "\n"))
  out
}
task <- function(rowset, slug, model, mods, kind, control = list())
  list(rowset = rowset, slug = slug, model = model, mods = mods, kind = kind, control = control)

# ---- All fits in one wave: the GEN timing run (from the gate fit), the 20 floor
# fits and the 12 diagonal slope fits, the last two from metafor's default start --
gate <- read.csv(file.path(out_dir, "gate_start_values.csv"))
gv <- setNames(gate$value, gate$component)
start_gate <- list(sigma2.init = unname(gv[c("sigma2_species", "sigma2_experiment_id", "sigma2_effect_id")]),
                   tau2.init = unname(gv["tau2_W"]), gamma2.init = unname(gv[c("gamma2_T", "gamma2_P", "gamma2_TP")]),
                   phi.init = unname(gv[c("phi_T_P", "phi_T_TP", "phi_P_TP")]))
wave <- c(list(task("all", "species_gen_timing", "Species-dependent (GEN timing run)", f_both, "gen_timing", start_gate)),
          lapply(fit_rowsets, function(r) task(r, "species", "Species-dependent", f_both, "slopes")),
          lapply(names(prop_vars), function(p) task(prop_rowset[[p]], p, prop_label[[p]], f_prop[[p]], "slopes")),
          unlist(lapply(fit_rowsets, function(r) lapply(names(floor_models), function(s)
            task(r, s, floor_models[[s]][[1]], floor_models[[s]][[2]], "floor"))), recursive = FALSE))
out <- run_wave(wave)
timing <- out[[1]]; fits <- out[-1]
stopifnot(!is.null(timing$error) || !is.null(timing$m)) # the GEN attempt either converged or left its error
names(fits) <- vapply(fits, function(f) paste(f$rowset, f$slug, sep = "/"), "")
fit_of <- function(rowset, slug) fits[[paste(if (rowset %in% c("growth_form", "phylogeny") &&
                                               !slug %in% names(prop_vars)) "all" else rowset, slug, sep = "/")]]

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
stopifnot(nrow(slope_rows) == 12, all(vapply(fits, function(f) !f$slopes || identical(f$m$struct, slope_struct), logical(1))))

# ---- Fixed effects with CR1 clustered on experiment (Table S5) ---------------
coefs <- bind_rows(lapply(fits, function(f) {
  if (f$slug == "null") return(NULL)
  r <- robust(f$m, cluster = rows[[f$rowset]]$experiment_id, adjust = TRUE)
  s <- coef(summary(r))
  tibble::tibble(rowset = f$rowset, model = f$model, term = rownames(s), estimate = s$estimate, se_CR1 = s$se,
                 ci_lower = s$ci.lb, ci_upper = s$ci.ub, df = s$df, pval = s$pval)
}))
write.csv(select(coefs, rowset, model, term, estimate, se_CR1, ci_lower, ci_upper), file.path(out_dir, "coefficients.csv"), row.names = FALSE)
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

# ---- Fig 4, the ladder: row a) from the table above, b) to d) wait for step 06 --
het_fig <- bind_rows(
  filter(ladder, rowset == "all", model %in% ladder_columns[1:5]) |> mutate(column = model, ref_cd = NA_real_, ref_H = NA_real_),
  bind_rows(lapply(names(prop_vars), function(p) {
    r <- prop_rowset[[p]]; sp <- filter(ladder, rowset == r, slug == "species")
    filter(ladder, rowset == r, slug == p) |>
      mutate(column = prop_label[[p]], ref_cd = sp$tau2 + sp$sigma2_experiment, ref_H = sp$H_intercept + sp$H_slope) })))
stopifnot(nrow(het_fig) == 12, setequal(het_fig$column, ladder_columns))
fig4_ladder(het_fig, "figures/diagnostics/05_fig4_ladder.png",
            caption = "Row a): REML fits of step 05 on the strength from year 3; property columns on their own complete rows, the black mark is Species-dependent on those rows. Rows b) to d): step 06.")

# ---- Fig S4: coefficients with CR1 intervals ---------------------------------
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
s4a <- coef_panel(filter(coefs, rowset == "all", model == "Climate + start", term != "intrcpt"), "a)", "Residents, Climate + start (CR1)")
s4b <- ggplot() + annotate("text", x = 0, y = 0, label = "Colonizers: step 07", size = 7, colour = "grey50") +
  labs(title = "b)") + theme_void(base_size = 17) + theme(plot.title = element_text(face = "bold", size = 22))
prop_coefs <- bind_rows(lapply(names(prop_vars), function(p) filter(coefs, model == prop_label[[p]]) |>
  filter(canon(term) %in% canon(prop_terms(prop_vars[[p]]))) |> mutate(property = prop_label[[p]])))
s4c <- coef_panel(prop_coefs, "c)", "Property terms, each from its own model (Climate + start + property, species slopes kept)") +
  facet_wrap(~ property, scales = "free", ncol = 4)
ggsave("figures/diagnostics/05_figS4.png", patchwork::wrap_plots(patchwork::wrap_plots(s4a, s4b, widths = c(1.4, 1)), s4c, ncol = 1, heights = c(1, 2)),
       width = 20, height = 22, dpi = 100, bg = "white")

# ---- Fig S8: fitted surfaces over T and P (fixed effects, A0 at its mean) ----
pred_fixed <- function(m, grid) { b <- coef(m)
  X <- sapply(names(b), function(nm) if (nm == "intrcpt") rep(1, nrow(grid)) else Reduce(`*`, grid[strsplit(nm, ":")[[1]]]))
  as.vector(X %*% b) }
exps <- distinct(strength_df, experiment_id, T, P)
grid0 <- expand.grid(T = seq(min(exps$T), max(exps$T), length.out = 60), P = seq(min(exps$P), max(exps$P), length.out = 60)) |>
  mutate(TP = T * P, A0_origin = 0, A0_dest = 0)
surf <- bind_rows(
  bind_rows(lapply(c(-1, 1), function(w) mutate(grid0, W = w, strength = pred_fixed(fit_of("all", "climate_start")$m, cbind(grid0, W = w)),
                                                  panel = sprintf("Climate + start, W %+d SD", w)))),
  bind_rows(lapply(names(prop_vars), function(p) bind_rows(lapply(prop_vars[[p]], function(v) {
    lv <- if (v %in% c("graminoid", "woody")) c(0, 1) else c(-1, 1)
    bind_rows(lapply(lv, function(x) { g <- cbind(grid0, W = 0); for (u in prop_vars[[p]]) g[[u]] <- 0; g[[v]] <- x
      mutate(g, strength = pred_fixed(fit_of(prop_rowset[[p]], p)$m, g),
             panel = if (v %in% c("graminoid", "woody")) sprintf("%s = %d, W at mean", v, x) else sprintf("%s %+d SD, W at mean", v, x)) }))
  }))))) |> mutate(panel = factor(panel, levels = unique(panel)))
lim <- max(abs(surf$strength))
s8 <- ggplot(surf, aes(T, P, fill = strength)) + geom_raster() +
  geom_point(data = exps, aes(T, P), inherit.aes = FALSE, shape = 21, fill = "white", size = 2) +
  facet_wrap(~ panel, ncol = 4) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", limits = c(-lim, lim), name = "Losing-winning\nstrength") +
  labs(x = "Destination summer temperature T (SD)", y = "Destination summer precipitation P (SD)") +
  theme_classic(base_size = 16) + theme(strip.text = element_text(face = "bold", size = 13), strip.background = element_blank())
ggsave("figures/diagnostics/05_figS8.png", s8, width = 20, height = 5 * ceiling(nlevels(surf$panel) / 4), dpi = 100, bg = "white", limitsize = FALSE)

# ---- Report ----------------------------------------------------------------
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
e4 <- function(x) ifelse(is.na(x), "", sprintf("%.3f", 1e4 * x))
md_all <- fits[["all/species"]]$m
gate_comp <- tibble::tibble(component = gate$component[1:7], `gate fit, GEN` = signif(gate$value[1:7], 4),
                            `step 05, diagonal` = signif(c(md_all$sigma2, md_all$tau2, md_all$gamma2), 4))
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
  "robust(m, cluster = experiment_id, adjust = TRUE)   # CR1 for Table S5", "```")
report <- c(
  "# Step 05: The ladder, fitted (Step D, route 1, residents)", "",
  paste0("The model ladder was fitted on the strength from year 3 (step 03b): Null, Start, Climate, Climate + start and ",
         "Species-dependent on all ", fmt(nrow(strength_df)), " effects and again on each property's complete rows, plus one ",
         "property model per property (LDMC, LMA, height, graminoid, woody, phylogeny as PEV1 to PEV3, niche as T_opt and P_opt) ",
         "with the species slopes kept, ", nrow(fitted), " REML fits in all. Each model reports τ² (species × site), ",
         "σ²_experiment (shared site shift) and the two H_species pieces as raw values. The step opened with a timing run of ",
         "Species-dependent under GEN from the gate fit's variance components (Bills, 1 Nov); it did not converge, so the T, P, TP ",
         "slope block is diagonal in every slope fit (Bills, 1 Nov). Fixed effects carry CR1 standard errors clustered ",
         "on experiment. Fits ran in parallel on ", cores, " cores, ", sprintf("%.0f", sum(fitted$minutes)),
         " fit-minutes in all (each fit's own time in the heterogeneity table; this run, reusing cached fits, took ",
         sprintf("%.0f", as.numeric(difftime(Sys.time(), t_start, units = "mins"))), " min)."),
  "", "## Stability across refits: the timing run", "",
  paste0("Species-dependent on all rows, GEN, `sparse = TRUE`, started from `gate_start_values.csv` (the variance components ",
         "of the 1 Nov gate fit of the same model on the same rows, ρ(T, P) = ", sprintf("%.3f", gate$value[gate$component == "phi_T_P"]),
         ", fitted from metafor's default start in 54 min, scratchpad run before this step): ",
         if (is.null(timing$error)) paste0("converged in ", sprintf("%.1f", timing$minutes), " min, ρ(T, P) = ", sprintf("%.4f", timing$m$phi[1]))
         else paste0("did not converge, stopping after ", sprintf("%.1f", timing$minutes), " min with \"", timing$error, "\""),
         ". A fit that converges from metafor's default start but not from its own optimum sits on a flat ridge: with the ",
         "experiments' T and P correlated at r = ", sprintf("%.2f", cor(distinct(strength_df, experiment_id, T, P)$T, distinct(strength_df, experiment_id, T, P)$P)),
         ", the T and P slope correlation is not identified. Decided 1 Nov (Bills): the T, P, TP slope block is diagonal in all ",
         "12 slope fits (the 25 Oct fallback, applied ladder-wide so that one structure runs through the ladder). The variance ",
         "components of the all-rows Species-dependent fit, diagonal, against the gate fit:"), "", md_table(gate_comp),
  "", "## Model as fitted (for reviewers)", "", fit_code,
  "", "## Counts per rowset", "", md_table(mutate(counts, across(-rowset, fmt))),
  "", "Each matches the expectation checked in the script. growth_form and phylogeny have the same rows as all, so their Null to Species-dependent are the all-rows fits (marked `reused`).",
  "", "## Heterogeneity per model (×10⁻⁴)", "", md_table(het_tab),
  "", paste0("All ", nrow(slope_rows), " slope fits converged with the diagonal T, P, TP block, each from metafor's default start. ",
             "Started instead from the all-rows fit's variance components, Species-dependent on the niche rows and the Woody model did ",
             "not converge (nlminb): the P and TP slope variances, at the boundary, started at 1e-10, where the optimiser stalled. ",
             "So no fit carries starting values from another (Bills, 1 Nov). Across the ", nrow(slope_rows), " fits the T slope variance ranges ",
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
  "", "## Figures", "", "![Fig 4, the ladder (row a)](../figures/diagnostics/05_fig4_ladder.png)", "",
  "![Fig S4, coefficients](../figures/diagnostics/05_figS4.png)", "", "![Fig S8, surfaces](../figures/diagnostics/05_figS8.png)",
  "", "## Warnings captured", "",
  { w <- unique(c(unlist(lapply(fits, function(f) paste0(f$rowset, " / ", f$model, ": ", f$warnings[[1]])[length(f$warnings[[1]]) > 0])), warns))
    if (length(w)) paste("-", w) else "None." },
  "", "## Deviations from PLAN.md", "",
  paste0("- H_species tests `m$withG` where PLAN.md tests `is.null(m$tau2)`: rma.mv returns tau2 = 0, not NULL, for a model without slopes. ",
         "Same quantity."),
  "- Fig S4 and Fig S8 are written to `figures/diagnostics/` (05_figS4.png, 05_figS8.png), not to `data/derived/05_ladder_fit/`, so that the report can link them; Fig S8's property surfaces are at W = its mean.",
  "- The T, P, TP slope block is diagonal (GDIAG) in every slope fit, where PLAN.md writes GEN with DIAG as the fallback if the correlations are unidentified (decision 25 Oct); applied ladder-wide after the GEN timing run did not converge (Bills, 1 Nov).",
  "- Decided 1 Nov (Bills), not deviations: no ML refits and no AICc column; no profile CIs; the 5% property check is measured against τ² + σ²_experiment; Fig 4 drawn here with twelve columns.")
writeLines(report, "reports/05_ladder_fit.md")
message("Step 05 written: ", nrow(fitted), " fits; GEN timing run ", sprintf("%.1f", timing$minutes), " min, ", if (is.null(timing$error)) "converged" else "did not converge")
