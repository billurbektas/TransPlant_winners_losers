# Step 07c: the colonizer ladder, fitted (Step D route 1, colonizers).
#
# PLAN.md: The colonizer story, "The model ladder for arrival"; Step D as step 05.
# Reads:  data/derived/07a_colonizers/strength_df_y3.rds (warmed cover, years >= 3);
#         data/derived/02_strength/strength_df.rds (the residents, for the property
#         scale); data/derived/prep/traits.rds, functional_groups.rds,
#         phylogeny_axes.rds, niche_optima.rds
# Writes: data/derived/07c_ladder_fit/ladder_rows.rds, models/<slug>.rds,
#         heterogeneity.csv, coefficients.csv; figures/diagnostics/07c_figS4b.png;
#         reports/07c_ladder_fit.md
# Checks: 379 colonizer effects in; ladder rows 110 effects, 27 species, 24
#         experiments, 18 gradients; no A0_origin in any model; 12 fits; random
#         structure exactly re_base / re_slopes; property terms exactly X, X:W, X:T,
#         X:P, X:TP; A0 interactions only A0_dest:W; slope fits with the diagonal
#         T, P, TP block. Step 05's tau2 and property checks are computed and reported.
#
# Decided 2 Nov (Bills): the residents' ladder (step 05) on colonizers, the same
# reduction (species in at least 3 gradients that carry every property) and the
# same 12 models with A0_origin dropped (PLAN.md); properties on the residents'
# step 02 scale (means and SDs across the 346 residents, traits logged), so both
# pools share one scale; A0_dest (raw destination-control cover averaged over
# years, step 07a) standardised across these effects as the residents' A0;
# model-based coefficient intervals as step 05; run 07c to 07f without stopping.

suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(metafor)})
source("R/fit_model.R"); source("R/h_species.R"); source("R/md_table.R")

warns <- character()
globalCallingHandlers(warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
})
out_dir <- "data/derived/07c_ladder_fit"
dir.create(file.path(out_dir, "models"), recursive = TRUE, showWarnings = FALSE)
cores <- 7L

d <- readRDS("data/derived/07a_colonizers/strength_df_y3.rds")
stopifnot(nrow(d) == 379, !"A0_origin" %in% names(d)) # step 07a

# ---- Properties on the residents' scale (step 02: across the 346 residents) ----------
residents <- unique(readRDS("data/derived/02_strength/strength_df.rds")$species)
stopifnot(length(residents) == 346)
z_res <- function(x, species, ref = residents) (x - mean(x[species %in% ref], na.rm = TRUE)) / sd(x[species %in% ref], na.rm = TRUE)
traits <- readRDS("data/derived/prep/traits.rds") |>
  transmute(species, ldmc_z = z_res(log(LDMC), species), lma_z = z_res(log(LMA), species), height_z = z_res(log(height), species))
groups <- readRDS("data/derived/prep/functional_groups.rds") |> transmute(species, graminoid = as.integer(graminoid), woody = as.integer(woody))
pev <- readRDS("data/derived/prep/phylogeny_axes.rds")$scores |>
  transmute(species, PEV1 = z_res(PEV1, species), PEV2 = z_res(PEV2, species), PEV3 = z_res(PEV3, species))
niche <- readRDS("data/derived/prep/niche_optima.rds") |> filter(n_obs_sub >= 30) |> # PLAN.md niche paragraph
  transmute(species, T_opt = z_res(T_opt, species), P_opt = z_res(P_opt, species))
# check: the residents' step 02 values are reproduced exactly on this scale
res02 <- distinct(readRDS("data/derived/02_strength/strength_df.rds"), species, ldmc_z, PEV1, T_opt)
chk <- inner_join(res02, select(traits, species, l2 = ldmc_z), by = "species") |> inner_join(select(pev, species, p2 = PEV1), by = "species") |>
  left_join(select(niche, species, t2 = T_opt), by = "species")
stopifnot(isTRUE(all.equal(chk$ldmc_z, chk$l2)), isTRUE(all.equal(chk$PEV1, chk$p2)), isTRUE(all.equal(chk$T_opt, chk$t2)))

d <- d |> mutate(A0_dest_raw = A0_dest, A0_dest = (A0_dest - mean(A0_dest)) / sd(A0_dest)) |>
  left_join(traits, by = "species") |> left_join(groups, by = "species") |> left_join(pev, by = "species") |> left_join(niche, by = "species")
props <- c("ldmc_z", "lma_z", "height_z", "graminoid", "woody", "PEV1", "PEV2", "PEV3", "T_opt", "P_opt")
g3 <- summarise(d, n_g = n_distinct(gradient_id), .by = species) |> filter(n_g >= 3)
lacking <- setdiff(g3$species, filter(d, species %in% g3$species) |> filter(if_all(all_of(props), ~ !is.na(.x))) |> pull(species))
rows <- filter(d, species %in% g3$species) |> filter(if_all(all_of(props), ~ !is.na(.x)))
stopifnot(nrow(g3) == 29, nrow(rows) == 110, n_distinct(rows$species) == 27, n_distinct(rows$experiment_id) == 24,
          n_distinct(rows$gradient_id) == 18) # at least 3 gradients, every property
saveRDS(rows, file.path(out_dir, "ladder_rows.rds"))

# ---- Formulas: step 05's, A0_origin dropped (PLAN.md) -----------------------------
f_env   <- ~ W * T * P
f_start <- ~ A0_dest + A0_dest:W
f_both  <- ~ W * T * P + A0_dest + A0_dest:W
prop_vars <- list(ldmc = "ldmc_z", lma = "lma_z", height = "height_z", graminoid = "graminoid", woody = "woody",
                  phylogeny = c("PEV1", "PEV2", "PEV3"), niche = c("T_opt", "P_opt"))
prop_label <- c(ldmc = "LDMC", lma = "LMA", height = "Height", graminoid = "Graminoid", woody = "Woody", phylogeny = "Phylogeny", niche = "Niche")
prop_terms <- function(v) unlist(lapply(v, function(x) paste0(x, c("", ":W", ":T", ":P", ":TP"))))
f_prop <- lapply(prop_vars, function(v) reformulate(c(attr(terms(f_both), "term.labels"), prop_terms(v))))
levels_fit <- tibble::tibble(slug = c("null", "start", "climate", "climate_start", "species", names(prop_vars)),
                             model = c("Null", "Start", "Climate", "Climate + start", "Species-dependent", prop_label),
                             slopes = c(rep(FALSE, 4), rep(TRUE, 8)))
mods <- c(list(null = ~ 1, start = f_start, climate = f_env, climate_start = f_both, species = f_both), f_prop)

run_fit <- function(i) {
  lv <- levels_fit[i, ]; path <- file.path(out_dir, "models", paste0(lv$slug, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  w <- character(); t0 <- Sys.time(); optimizer <- "nlminb"
  fit <- function(control) withCallingHandlers(
    tryCatch(fit_model(mods[[lv$slug]], if (lv$slopes) re_slopes else re_base, rows, struct = if (lv$slopes) slope_struct else "GEN",
                       control = control), error = function(e) conditionMessage(e)),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  m <- fit(list()); if (is.character(m)) { optimizer <- "bobyqa"; m <- fit(list(optimizer = "bobyqa")) }
  if (is.character(m)) stop(lv$model, ": ", m)
  res <- list(m = m, slug = lv$slug, model = lv$model, slopes = lv$slopes, optimizer = optimizer,
              minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")), warnings = unique(w))
  saveRDS(res, path); res
}
fits <- parallel::mclapply(seq_len(nrow(levels_fit)), run_fit, mc.cores = cores, mc.preschedule = FALSE)
if (any(sapply(fits, inherits, "try-error"))) stop("fits failed: ", paste(unlist(fits[sapply(fits, inherits, "try-error")]), collapse = "\n"))
names(fits) <- levels_fit$slug
stopifnot(length(fits) == 12) # 4 floor, Species-dependent, 7 property models

# ---- Heterogeneity per model --------------------------------------------------------
ladder <- bind_rows(lapply(fits, function(f) { m <- f$m; H <- H_species(m, rows)
  tibble::tibble(model = f$model, n_effects = m$k, tau2 = tau2(m), sigma2_experiment = sigma2_exp(m), sigma2_species = H[["intercept"]],
                 H_intercept = H[["intercept"]], H_slope = H[["slope"]], var_T = if (f$slopes) m$gamma2[1] else NA_real_,
                 var_P = if (f$slopes) m$gamma2[2] else NA_real_, var_TP = if (f$slopes) m$gamma2[3] else NA_real_,
                 logLik_REML = as.numeric(logLik(m)), optimizer = f$optimizer, minutes = f$minutes) }))
write.csv(ladder, file.path(out_dir, "heterogeneity.csv"), row.names = FALSE)

# ---- Checks, as step 05 -------------------------------------------------------------
for (f in fits) { # random structure exactly re_base / re_slopes
  stopifnot(identical(unname(f$m$s.names), c("species", "experiment_id", "effect_id")))
  if (f$slopes) stopifnot(identical(unname(f$m$g.names), c("W", "species")), identical(unname(f$m$h.names), c("T", "P", "TP", "species")),
                          identical(f$m$struct, slope_struct))
  else stopifnot(!isTRUE(f$m$withG), !isTRUE(f$m$withH))
  stopifnot(!any(grepl("A0_origin", colnames(f$m$X)))) # A0_origin absent (PLAN.md)
}
canon <- function(x) vapply(strsplit(x, ":"), function(p) paste(sort(p), collapse = ":"), "")
# property terms X, X:W, X:T, X:P, X:TP; metafor drops a term that is redundant on these rows (listed)
dropped <- unlist(lapply(names(prop_vars), function(p) { added <- canon(setdiff(colnames(fits[[p]]$m$X), colnames(fits$species$m$X)))
  stopifnot(all(added %in% canon(prop_terms(prop_vars[[p]])))); setdiff(canon(prop_terms(prop_vars[[p]])), added) }))
stopifnot(identical(dropped, "TP:woody")) # one woody species among the 27: woody:TP is redundant
stopifnot(setequal(canon(grep("A0", colnames(fits$climate_start$m$X), value = TRUE)), canon(c("A0_dest", "A0_dest:W"))))
# step 05 stops on these two; for colonizers they are reported (Bills, 2 Nov: run without stopping)
climate_change <- tau2(fits$climate$m) / tau2(fits$null$m) - 1
prop_check <- bind_rows(lapply(names(prop_vars), function(p) { m4 <- fits[[p]]$m; m3 <- fits$species$m; denom <- tau2(m3) + sigma2_exp(m3)
  tibble::tibble(property = prop_label[[p]], move_tau2 = (tau2(m4) - tau2(m3)) / denom, move_sigma2_experiment = (sigma2_exp(m4) - sigma2_exp(m3)) / denom) }))

# ---- Fixed effects, model-based intervals (Fig S4 b, as step 05's a) ------------------
coefs <- bind_rows(lapply(fits[-1], function(f) { s <- coef(summary(f$m))
  tibble::tibble(model = f$model, term = rownames(s), estimate = s$estimate, se = s$se, ci_lower = s$ci.lb, ci_upper = s$ci.ub, pval = s$pval) }))
write.csv(coefs, file.path(out_dir, "coefficients.csv"), row.names = FALSE)
coef_panel <- function(x, tag, ttl) { # as analysis/05_ladder_fit.R
  x <- mutate(x, excl = ci_lower > 0 | ci_upper < 0, term = factor(term, levels = rev(unique(term))))
  ggplot(x, aes(estimate, term)) + geom_vline(xintercept = 0, colour = "grey60") +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0, linewidth = 0.9) +
    geom_point(aes(fill = excl), shape = 21, size = 4, stroke = 1) +
    scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "white"), labels = c(`TRUE` = "95% interval excludes zero", `FALSE` = "includes zero"), name = NULL) +
    labs(title = tag, subtitle = ttl, x = "Effect on warmed cover (per SD)", y = NULL) +
    theme_classic(base_size = 17) + theme(plot.title = element_text(face = "bold", size = 22), legend.position = "bottom",
                                          strip.text = element_text(face = "bold", size = 15), strip.background = element_blank())
}
s4b <- coef_panel(filter(coefs, model == "Climate + start", term != "intrcpt"), "b)", "Colonizers, Climate + destination abundance (model-based intervals)")
prop_coefs <- bind_rows(lapply(names(prop_vars), function(p) filter(coefs, model == prop_label[[p]], canon(term) %in% canon(prop_terms(prop_vars[[p]]))) |>
  mutate(property = prop_label[[p]])))
s4c <- coef_panel(prop_coefs, "c)", "Colonizers: property terms, each from its own model") + facet_wrap(~ property, scales = "free", ncol = 4)
ggsave("figures/diagnostics/07c_figS4b.png", patchwork::wrap_plots(s4b, s4c, ncol = 1, heights = c(1, 2)), width = 20, height = 20, dpi = 100, bg = "white")

# ---- Report --------------------------------------------------------------------------
e4 <- function(x) ifelse(is.na(x), "", sprintf("%.3f", 1e4 * x))
excl <- filter(coefs, ci_lower > 0 | ci_upper < 0, term != "intrcpt")
report <- c(
  "# Step 07c: The colonizer ladder, fitted (Step D route 1, colonizers)", "",
  paste0("The residents' ladder (step 05) was fitted to the colonizers (Bills, 2 Nov): warmed cover from years 3 and later (step 07a) for the ",
         nrow(rows), " effects of the ", n_distinct(rows$species), " colonizer species found in at least 3 gradients that carry every property (",
         n_distinct(rows$experiment_id), " experiments, ", n_distinct(rows$gradient_id), " gradients). The same 12 models, with A0_origin dropped ",
         "(a colonizer has none) so that Start is the destination abundance A0_dest and its interaction with W. Properties are on the residents' ",
         "step 02 scale (checked), A0_dest standardised across these effects. ", sum(rows$strength_se == 0), " effects are measured absences (SE 0)."),
  "", "## Counts", "",
  md_table(tibble::tibble(rule = c("colonizer effects (07a, years 3 and later)", "species in >= 3 gradients", "of which lacking a property",
                                   "ladder effects", "species", "experiments", "gradients", "measured absences among them"),
                          n = c(379, nrow(g3), length(lacking), nrow(rows), n_distinct(rows$species), n_distinct(rows$experiment_id),
                                n_distinct(rows$gradient_id), sum(rows$strength_se == 0)))),
  "", paste0("Species lacking a property: ", paste(lacking, collapse = ", "), "."),
  "", paste0("Growth form on these rows: ", sum(tapply(rows$woody, rows$species, max)), " woody and ", sum(tapply(rows$graminoid, rows$species, max)), " graminoid species of ", n_distinct(rows$species), "; with one woody species metafor drops the redundant term ", paste(dropped, collapse = ", "), " from the Woody model, which rests on that species alone."),
  "", "## Heterogeneity per model (×10⁻⁴)", "",
  md_table(transmute(ladder, model, `τ²` = e4(tau2), `σ²_experiment` = e4(sigma2_experiment), `H intercept` = e4(H_intercept), `H slope` = e4(H_slope),
                     optimizer, minutes = sprintf("%.1f", minutes))),
  "", "## Step 05's checks, reported", "",
  paste0("τ² at Climate against Null: ", sprintf("%+.1f%%", 100 * climate_change), " (step 05 stops beyond 5%)."), "",
  md_table(transmute(prop_check, property, `τ² move` = sprintf("%+.1f%%", 100 * move_tau2), `σ²_experiment move` = sprintf("%+.1f%%", 100 * move_sigma2_experiment))),
  "", paste0("Moves beyond step 05's 5% limit: ", sum(abs(c(prop_check$move_tau2, prop_check$move_sigma2_experiment)) > 0.05), " of ", 2 * nrow(prop_check),
             ". On ", nrow(rows), " effects the variance components are not pinned down well enough for that check to hold; they are reported, not enforced (Bills, 2 Nov)."),
  "", "## Coefficients whose 95% interval excludes zero", "",
  if (nrow(excl)) md_table(transmute(excl, model, term, estimate = sprintf("%.4f", estimate), ci = sprintf("%.4f to %.4f", ci_lower, ci_upper))) else "None.",
  "", "## Figures", "", "![Fig S4 b and c, colonizers](../figures/diagnostics/07c_figS4b.png)",
  "", "## Warnings captured", "", { w <- unique(c(unlist(lapply(fits, function(f) if (length(f$warnings)) paste0(f$model, ": ", f$warnings))), warns))
    if (length(w)) paste("-", w) else "None." },
  "", "## Deviations from PLAN.md", "",
  "- The colonizer ladder on the residents' reduction and models (Bills, 2 Nov); properties on the residents' scale.",
  "- Step 05's tau2 and property checks are reported, not enforced (Bills, 2 Nov: run without stopping).",
  "- Fig S4 panels b and c drawn here as their own file; step 10 assembles Fig S4.")
writeLines(report, "reports/07c_ladder_fit.md")
message("Step 07c written: ", nrow(rows), " effects, ", length(fits), " fits")
