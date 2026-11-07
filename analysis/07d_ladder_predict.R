# Step 07d: the colonizer ladder, predicted (Step D route 2, colonizers).
#
# PLAN.md: The colonizer story, "The model ladder for arrival"; Step D Route 2 as
#   step 06; Fig S18 (colonizer performance) and Fig S21 (per gradient).
# Reads:  data/derived/07c_ladder_fit/ladder_rows.rds, models/<slug>.rds
# Writes: data/derived/07d_ladder_predict/folds/<level>_<gradient_id>.rds,
#         pred_metrics.csv, fit_metrics.csv, summary_metrics.csv;
#         figures/appendix/figS18.png, figS21.png; reports/07d_ladder_predict.md
# Checks: 110 ladder rows (07c); 8 scored gradients (at least 5 effects), one fold
#         each per fitted level; the held-out gradient absent from every training
#         set; species deviations only for species seen in training; Null R²
#         exactly 0 on every scored gradient, in sample and held out.
#
# Decided 2 Nov (Bills): step 06 on the colonizer ladder rows (every gradient
# trains, those with at least 5 effects are held out and scored; fixed-effects
# prediction, hand-BLUP species deviations, supplement <property>+species levels;
# R² against the Null on the same gradient and Spearman, medians); no loser/winner
# skill (colonizers have no loser class, PLAN.md); an R² below -100 is taken out of
# the R² row and the medians and listed.
# Without abundance (Bills, 7 Nov): step 06's new levels; Experimental warming and
# Background climate have no Spearman (one prediction per experiment).
# Warm start (Bills, 7 Nov), as step 06: every property fold starts its variance
# components at those of the same fold's Species-specific responses fit (stored with
# that fold) and is fitted with optim BFGS; a failing one falls back to nlminb from
# metafor's default start, then bobyqa, then NA.

suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(metafor)})
source("R/fit_model.R"); source("R/metrics_gradient.R"); source("R/species_deviation.R"); source("R/ladder_figures.R"); source("R/md_table.R")
source("R/cache.R"); cache_version <- "07d-v1" # raise when the fold code changes

warns <- character()
globalCallingHandlers(warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
})
out_dir <- "data/derived/07d_ladder_predict"
fold_dir <- file.path(out_dir, "folds")
dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
cores <- 4L # of 11; 7 choked the machine (Bills, 7 Nov); results do not depend on it

rows <- readRDS("data/derived/07c_ladder_fit/ladder_rows.rds")
stopifnot(nrow(rows) == 110, n_distinct(rows$gradient_id) == 18) # step 07c
n_eff <- count(rows, gradient_id)
scored <- sort(n_eff$gradient_id[n_eff$n >= 5])
stopifnot(length(scored) == 8) # held out and scored: at least 5 effects

# ---- Levels and formulas, as 07c ------------------------------------------------------
f_warm <- ~ W; f_env <- ~ W * T * P
prop_vars <- list(ldmc = "ldmc_z", lma = "lma_z", height = "height_z", graminoid = "graminoid", woody = "woody",
                  phylogeny = c("PEV1", "PEV2", "PEV3"), niche = c("T_opt", "P_opt"))
prop_label <- c(ldmc = "LDMC", lma = "LMA", height = "Height", graminoid = "Graminoid", woody = "Woody", phylogeny = "Phylogeny", niche = "Niche")
prop_terms <- function(v) unlist(lapply(v, function(x) paste0(x, c("", ":W", ":T", ":P", ":TP"))))
levels_fit <- tibble::tibble(level = c("null", "warming", "climate", "species", names(prop_vars)),
                             column = c("Null", "Experimental warming", "Background climate", "Species-specific responses", prop_label),
                             slopes = c(rep(FALSE, 3), rep(TRUE, 8)))
stopifnot(identical(unname(levels_fit$column), ladder_columns))
mods <- c(list(null = ~ 1, warming = f_warm, climate = f_env, species = f_env),
          lapply(prop_vars, function(v) reformulate(c(attr(terms(f_env), "term.labels"), prop_terms(v)))))
supp <- paste0(names(prop_vars), "+species")
level_label <- c(setNames(levels_fit$column, levels_fit$level), setNames(paste(prop_label, "+ species"), supp))
level_display <- level_label # the model names are the display names (Bills, 7 Nov)
one_per_experiment <- c("warming", "climate") # predictions vary only between experiments
fixed_part <- function(m, f, d) { X <- model.matrix(f, d); colnames(X)[colnames(X) == "(Intercept)"] <- "intrcpt"
  as.numeric(X[, names(coef(m)), drop = FALSE] %*% coef(m)) }

# ---- One fold, cached ---------------------------------------------------------------
fit_fold <- function(task) {
  path <- file.path(fold_dir, paste0(task$level, "_", task$gradient_id, ".rds"))
  lv <- filter(levels_fit, level == task$level); g <- task$gradient_id
  start <- list()
  if (task$level %in% names(prop_vars)) { # warm start from this fold's Species-specific responses fit (Bills, 7 Nov)
    sp <- readRDS(file.path(fold_dir, paste0("species_", g, ".rds")))
    stopifnot(sp$summary$status == "converged", !is.null(sp$start))
    start <- c(list(optimizer = "optim", optmethod = "BFGS"), sp$start)
  }
  # slope folds store their variance components as the start of the property folds (7 Nov)
  key <- cache_key(rows[fit_columns(mods[[lv$level]])], deparse(mods[[lv$level]]), lv$slopes, g, cache_version,
                   if (lv$slopes) list("stores start", start))
  hit <- read_cache(path, key); if (!is.null(hit)) return(hit)
  train <- filter(rows, gradient_id != g); test <- filter(rows, gradient_id == g)
  stopifnot(!g %in% train$gradient_id, nrow(test) > 0)
  w <- character()
  fit <- function(control) withCallingHandlers(
    tryCatch(fit_model(mods[[lv$level]], if (lv$slopes) re_slopes else re_base, train, struct = if (lv$slopes) slope_struct else "GEN",
                       control = control), error = function(e) conditionMessage(e)),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  if (length(start)) { optimizer <- "optim BFGS, warm start"; m <- fit(start)
    if (is.character(m)) { optimizer <- "nlminb"; m <- fit(list()) } } else { optimizer <- "nlminb"; m <- fit(list()) }
  if (is.character(m)) { optimizer <- "bobyqa"; m <- fit(list(optimizer = "bobyqa")) }
  ok <- !is.character(m)
  pred <- tibble::tibble(level = lv$level, gradient_id = g, effect_id = test$effect_id, species = test$species, observed = test$strength,
                         fixed = if (ok) fixed_part(m, mods[[lv$level]], test) else NA_real_,
                         deviation = if (ok && lv$slopes) species_deviation(m, train, test) else if (ok) 0 else NA_real_,
                         species_known = test$species %in% train$species)
  stopifnot(all(pred$deviation[!pred$species_known] %in% c(0, NA))) # deviations only for species seen in training
  res <- list(pred = pred, summary = tibble::tibble(level = lv$level, gradient_id = g, status = if (ok) "converged" else m, optimizer = optimizer),
              warnings = unique(w), # boundary values lifted to 1e-8 so that the optimiser can leave them
              start = if (ok && lv$slopes) list(sigma2.init = pmax(m$sigma2, 1e-8), tau2.init = max(m$tau2, 1e-8), gamma2.init = pmax(m$gamma2, 1e-8)))
  write_cache(res, path, key)
}
# two phases: the property folds start from the Species-specific responses folds
run_folds <- function(levels) {
  tasks <- tidyr::expand_grid(level = levels, gradient_id = scored)
  out <- parallel::mclapply(lapply(seq_len(nrow(tasks)), function(i) as.list(tasks[i, ])), fit_fold, mc.cores = cores, mc.preschedule = FALSE)
  if (any(sapply(out, inherits, "try-error"))) stop("worker crashed; rerun resumes from the cached folds")
  out
}
out <- run_folds(setdiff(levels_fit$level, names(prop_vars)))
out <- c(out, run_folds(names(prop_vars)))
folds <- bind_rows(lapply(out, `[[`, "pred")); fold_summary <- bind_rows(lapply(out, `[[`, "summary"))
stopifnot(nrow(fold_summary) == 11 * 8, !anyDuplicated(fold_summary[c("level", "gradient_id")])) # one fold per gradient per level

pred_long <- bind_rows(mutate(folds, predicted = fixed + ifelse(level == "species", deviation, 0)),
                       filter(folds, level %in% names(prop_vars)) |> mutate(level = paste0(level, "+species"), predicted = fixed + deviation))
fit_long <- bind_rows(lapply(seq_len(nrow(levels_fit)), function(i) {
  lv <- levels_fit[i, ]; m <- readRDS(file.path("data/derived/07c_ladder_fit/models", paste0(lv$level, ".rds")))$m
  fx <- fixed_part(m, mods[[lv$level]], rows); dev <- if (lv$slopes) species_deviation(m, rows, rows) else 0
  base <- tibble::tibble(gradient_id = rows$gradient_id, effect_id = rows$effect_id, observed = rows$strength)
  bind_rows(mutate(base, level = lv$level, predicted = fx + if (lv$level == "species") dev else 0),
            if (lv$level %in% names(prop_vars)) mutate(base, level = paste0(lv$level, "+species"), predicted = fx + dev))
}))

# ---- Metrics: R² against the Null on the same gradient, Spearman; no skill ---------------
score <- function(d) {
  d <- inner_join(d, filter(d, level == "null") |> select(effect_id, null_pred = predicted), by = "effect_id")
  bind_rows(lapply(split(d, list(d$level, d$gradient_id), drop = TRUE), function(x) {
    met <- if (anyNA(x$predicted)) tibble::tibble(r2 = NA_real_, spearman = NA_real_)
           else metrics_gradient(x$observed, x$predicted, mean((x$observed - x$null_pred)^2)) |> select(r2, spearman)
    if (x$level[1] %in% one_per_experiment) met$spearman <- NA_real_ # cannot order species within a gradient (Bills, 2 and 7 Nov)
    tibble::tibble(level = x$level[1], gradient_id = x$gradient_id[1], met, n = nrow(x))
  }))
}
pred_metrics <- score(pred_long); fit_metrics <- score(filter(fit_long, gradient_id %in% scored))
stopifnot(setequal(pred_metrics$gradient_id, scored), all(filter(pred_metrics, level == "null")$r2 == 0), all(filter(fit_metrics, level == "null")$r2 == 0))
write.csv(pred_metrics, file.path(out_dir, "pred_metrics.csv"), row.names = FALSE)
write.csv(fit_metrics, file.path(out_dir, "fit_metrics.csv"), row.names = FALSE)
# an R² below -100 is taken out of the R² row and the medians, and listed (Bills, 2 Nov)
outliers <- bind_rows(mutate(filter(pred_metrics, r2 < -100), version = "held out"), mutate(filter(fit_metrics, r2 < -100), version = "in sample"))
trim <- function(m) mutate(m, r2 = ifelse(r2 < -100, NA_real_, r2))
all_levels <- c(levels_fit$level, supp)
summary_metrics <- bind_rows(lapply(c("r2", "spearman"), function(k) {
  f <- trim(fit_metrics) |> summarise(fit = median(.data[[k]], na.rm = TRUE), n_gradients_fit = sum(!is.na(.data[[k]])), .by = level)
  p <- trim(pred_metrics) |> summarise(pred = median(.data[[k]], na.rm = TRUE), n_gradients_pred = sum(!is.na(.data[[k]])),
                                       n_better_than_null = if (k == "r2") sum(.data[[k]] > 0, na.rm = TRUE) else NA_integer_, .by = level)
  full_join(f, p, by = "level") |> mutate(metric = k, .after = level)
})) |> arrange(match(level, all_levels), metric)
write.csv(summary_metrics, file.path(out_dir, "summary_metrics.csv"), row.names = FALSE)

# ---- Fig S18: colonizer performance (as step 06's Fig 6, without the skill row) ----------
grads <- bind_rows(mutate(trim(fit_metrics), type = "fit"), mutate(trim(pred_metrics), type = "pred")) |>
  filter(level %in% levels_fit$level, level != "null") |> tidyr::pivot_longer(c(r2, spearman), names_to = "metric") |>
  transmute(column = level_label[level], metric = c(r2 = "R2", spearman = "Spearman")[metric], type, value) |> filter(!is.na(value))
ladder_performance(grads, "figures/appendix/figS18.png", skill = FALSE)

# ---- Fig S21 (colonizers): every metric x level x gradient, as step 06 ----------------------
grad_info <- filter(rows, gradient_id %in% scored) |> summarise(T_dest = mean(T), .by = gradient_id) |> arrange(T_dest)
s7 <- bind_rows(mutate(fit_metrics, version = "In sample"), mutate(pred_metrics, version = "Held out")) |>
  tidyr::pivot_longer(c(r2, spearman), names_to = "metric") |>
  mutate(gradient_id = factor(gradient_id, levels = rev(grad_info$gradient_id)), level = factor(level_display[level], levels = level_display[all_levels]),
         version = factor(version, levels = c("In sample", "Held out")))
s7_panel <- function(k, title, tags) lapply(1:2, function(j) {
  v <- levels(s7$version)[j]; d <- filter(s7, metric == k, version == v)
  ggplot(d, aes(level, gradient_id, fill = pmax(pmin(value, 1), -1))) + geom_tile(colour = "white") +
    geom_text(aes(label = ifelse(is.na(value), "", ifelse(value < -10, "< -10", sprintf("%.2f", value)))), size = 7) + # readable at print size (Bills, 2 Nov)
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", limits = c(-1, 1), na.value = "grey85", name = title) +
    labs(title = tags[j], subtitle = paste0(title, ", ", tolower(v)), x = NULL, y = NULL) + theme_minimal(base_size = 30) +
    theme(plot.title = element_text(face = "bold", size = 40), plot.subtitle = element_text(size = 32), axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(size = 26), legend.key.height = unit(2, "cm"), legend.key.width = unit(1, "cm"),
          panel.grid = element_blank(), plot.background = element_rect(fill = "white", colour = NA))
})
ggsave("figures/appendix/figS21.png", patchwork::wrap_plots(c(s7_panel("r2", "R² against Null", c("a)", "b)")),
                                                                   s7_panel("spearman", "Spearman ρ", c("c)", "d)"))), ncol = 2),
       width = 48, height = 30, dpi = 100, bg = "white")

# ---- Report ---------------------------------------------------------------------------
f3 <- function(x) ifelse(is.na(x), "", sprintf("%.3f", x))
sum_tab <- summary_metrics |> transmute(level = level_display[level], metric = c(r2 = "R² against Null", spearman = "Spearman ρ")[metric],
                                        fit = f3(fit), pred = f3(pred), `gradients (fit / pred)` = paste(n_gradients_fit, "/", n_gradients_pred),
                                        `held out better than Null` = ifelse(is.na(n_better_than_null), "", paste(n_better_than_null, "of", n_gradients_pred)))
failed <- filter(fold_summary, status != "converged")
fold_warn <- unique(unlist(lapply(out, function(x) if (length(x$warnings)) paste0(x$summary$level, " / ", x$summary$gradient_id, ": ", x$warnings))))
# metafor's "struct disregarded" note is listed once with its count (CLAUDE.md: warnings are reported, never dropped)
n_struct <- sum(grepl("struct' argument is disregaded", fold_warn))
fold_warn <- c(fold_warn[!grepl("struct' argument is disregaded", fold_warn)],
               if (n_struct) paste0("metafor: \"'struct' argument is disregaded since no ~ inner | outer term\" (", n_struct, " fold fits without slopes)"))
report <- c(
  "# Step 07d: The colonizer ladder, predicted (Step D route 2, colonizers)", "",
  paste0("Step 06 on the ", nrow(rows), " colonizer ladder rows of step 07c: every level was refitted ", length(scored), " times, each time without one of the ",
         length(scored), " gradients with at least 5 effects, and asked to predict its warmed cover; fixed effects with the experiment term at zero, ",
         "Species-specific responses and the supplement levels adding the species' deviations learnt in training. Scored per gradient by R² against the ",
         "Null on the same gradient and Spearman ρ, in sample and held out, summarised by medians; no loser/winner skill, since colonizers have no ",
         "loser class. ", nrow(fold_summary), " fold fits (", sum(fold_summary$optimizer == "optim BFGS, warm start"), " property folds warm-started), ", sum(fold_summary$optimizer == "bobyqa"), " needing the bobyqa retry, ", nrow(failed), " failed."),
  "", if (nrow(failed)) c("Folds that did not converge with nlminb nor with the bobyqa retry, left NA (Bills, 2 Nov):", "",
      md_table(transmute(failed, level = level_display[level], gradient_id, `bobyqa error` = status))) else "Every fold fit converged.",
  "", "## Counts", "",
  md_table(transmute(n_eff, gradient_id, effects = n, held_out = ifelse(gradient_id %in% scored, "yes", "no (fewer than 5 effects; trains only)"))),
  "", "## Prediction metrics, median over gradients", "", md_table(sum_tab),
  "", if (nrow(outliers)) paste0("Taken out of the R² row and the medians (R² below −100, Bills, 2 Nov): ",
                                  paste(sprintf("%s, %s at %s (%s)", level_display[outliers$level], outliers$version, outliers$gradient_id, sprintf("%.1f", outliers$r2)), collapse = "; "), ".")
      else "No R² below −100.",
  "", paste0("Every held-out species had training data: ", all(folds$species_known), "."),
  "", "## Figures", "", "![Fig S18, colonizer performance](../figures/appendix/figS18.png)", "",
  paste0("Fig S18 as Fig 6: half violins of the ", length(scored), " scored gradients, explanatory (grey) beside predictive (blue), medians written in; ",
         "R² against the Null on that gradient. No skill row. Experimental warming and Background climate have no Spearman (one prediction per experiment)."), "",
  "![Fig S21, colonizers per gradient](../figures/appendix/figS21.png)",
  "", "## Warnings captured", "", { w <- unique(c(fold_warn, warns)); if (length(w)) paste("-", w) else "None." },
  "", "## Deviations from PLAN.md", "",
  "- The colonizer ladder predicted on the residents' reduction (Bills, 2 Nov).",
  "- Property folds warm-started from the same fold's Species-specific responses fit with optim BFGS (Bills, 7 Nov, as steps 05 and 06).",
  "- An R² below −100 is taken out of the R² row and the medians and listed (Bills, 2 Nov).",
  "- Without abundance (Bills, 7 Nov): Fig S18 shows the performance rows only; the colonizer variance components are Figs S19 and S20 (step 07c).")
writeLines(report, "reports/07d_ladder_predict.md")
message("Step 07d written: ", nrow(fold_summary), " folds, ", nrow(failed), " failed")
