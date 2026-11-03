# Step 06: the ladder, predicted (Step D, route 2, residents).
#
# PLAN.md: Step D, "Route 2. Prediction of an unseen gradient"; "Predicting a
#          gradient with the experiment term in the model".
# Reads:  data/derived/03b_strength_y3/strength_df.rds (strength from year 3);
#         data/derived/05_ladder_fit/models/<rowset>/<model>.rds (fit metrics);
#         data/derived/05_ladder_fit/heterogeneity.csv (Fig 6 row a)
# Writes: data/derived/06_ladder_predict/folds/<level>_<gradient_id>.rds,
#         pred_metrics.csv, fit_metrics.csv, summary_metrics.csv,
#         fold_summary.csv; figures/main/fig6.png (the ladder, main Fig 6),
#         06_figS7.png; reports/06_ladder_predict.md
# Checks: 1,004 effects in, 512 ladder rows (71 species, 39 experiments, 20
#         gradients) as in step 05, every rowset equal to them; 16 scored
#         gradients (at least 5 effects), one fold each per fitted level; every
#         held-out species seen in training; the held-out
#         gradient absent from every training set; intercept renamed before
#         matching coef(m); species deviations only for species seen in
#         training; hand BLUPs equal ranef() to 1e-8; Null R2 exactly 0 on
#         every scored gradient, in sample and held out; every slope fold with
#         the diagonal T, P, TP block.
#
# Decided 2 Nov (Bills): the ladder rows of step 05 (residents in at least 3
# gradients that carry every property); every gradient trains, but only the 16
# with at least 5 effects are held out and scored; one fold each; species deviations
# from hand BLUPs (R/species_deviation.R), metafor giving none for a GDIAG slope
# block; a fold that fails is refitted once with bobyqa, then NA for that level x
# gradient, named in the report; Climate's Spearman and skill NA in every
# gradient (climate cannot order species within a gradient); supplement levels
# <property>+species from the property fold models; no all-years sensitivity;
# per fold the predictions and a one-row summary, no model objects.
# Reopened 2 Nov (Bills): R2 per gradient against the Null model on the same
# gradient, 1 - MSE(level) / MSE(Null), replacing the shared yardstick V; every
# metric summarised by its median over the scored gradients; Fig 6 rows b) to d)
# as half violins; display labels "Initial abundances", "Loser/winner skill".

suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(metafor)})
source("R/fit_model.R"); source("R/h_species.R"); source("R/metrics_gradient.R")
source("R/species_deviation.R"); source("R/fig4_ladder.R"); source("R/md_table.R")

warns <- character()
globalCallingHandlers(warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
})
t_start <- Sys.time()
out_dir <- "data/derived/06_ladder_predict"
fold_dir <- file.path(out_dir, "folds")
dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
cores <- 7L # of 11; 8 strained the machine (Bills, 2 Nov)

strength_df <- readRDS("data/derived/03b_strength_y3/strength_df.rds")
stopifnot(nrow(strength_df) == 1004) # step 03b

# ---- Ladder rows, as in step 05: residents in at least 3 gradients with every property
props <- c("ldmc_z", "lma_z", "height_z", "graminoid", "woody", "PEV1", "PEV2", "PEV3", "T_opt", "P_opt")
g3 <- summarise(strength_df, n_g = n_distinct(gradient_id), .by = species) |> filter(n_g >= 3)
strength_df <- filter(strength_df, species %in% g3$species) |> filter(if_all(all_of(props), ~ !is.na(.x)))
stopifnot(nrow(strength_df) == 512, n_distinct(strength_df$species) == 71,
          n_distinct(strength_df$experiment_id) == 39, n_distinct(strength_df$gradient_id) == 20) # step 05 rows

# ---- Rowsets and formulas, as in step 05 (every rowset equals the ladder rows) ---
rows <- list(all         = strength_df,
             ldmc        = filter(strength_df, !is.na(ldmc_z)),
             lma         = filter(strength_df, !is.na(lma_z)),
             height      = filter(strength_df, !is.na(height_z)),
             growth_form = filter(strength_df, !is.na(graminoid), !is.na(woody)),
             phylogeny   = filter(strength_df, !is.na(PEV1), !is.na(PEV2), !is.na(PEV3)),
             niche       = filter(strength_df, !is.na(T_opt), !is.na(P_opt)))
stopifnot(all(sapply(rows, function(r) identical(r$effect_id, rows$all$effect_id)))) # one set of rows (step 05)
gradients <- sort(unique(strength_df$gradient_id)) # all 20 train
# held out and scored: gradients with at least 5 effects (Bills, 2 Nov); the others train only
n_eff <- count(strength_df, gradient_id)
scored <- sort(n_eff$gradient_id[n_eff$n >= 5]); unscored <- filter(n_eff, n < 5)
stopifnot(length(gradients) == 20, length(scored) == 16,
          setequal(unscored$gradient_id, c("CN_Damxung", "US_Montana", "US_Colorado", "CN_Gongga")))

f_env   <- ~ W * T * P
f_start <- ~ A0_origin + A0_dest + A0_origin:W + A0_dest:W
f_both  <- ~ W * T * P + A0_origin + A0_dest + A0_origin:W + A0_dest:W
prop_vars <- list(ldmc = "ldmc_z", lma = "lma_z", height = "height_z", graminoid = "graminoid", woody = "woody",
                  phylogeny = c("PEV1", "PEV2", "PEV3"), niche = c("T_opt", "P_opt"))
prop_rowset <- c(ldmc = "ldmc", lma = "lma", height = "height", graminoid = "growth_form", woody = "growth_form",
                 phylogeny = "phylogeny", niche = "niche")
prop_label <- c(ldmc = "LDMC", lma = "LMA", height = "Height", graminoid = "Graminoid", woody = "Woody",
                phylogeny = "Phylogeny", niche = "Niche")
prop_terms <- function(v) unlist(lapply(v, function(x) paste0(x, c("", ":W", ":T", ":P", ":TP"))))
f_prop <- lapply(prop_vars, function(v) reformulate(c(attr(terms(f_both), "term.labels"), prop_terms(v))))

# ---- Prediction levels (PLAN.md level_formula, one level per property) --------
levels_fit <- tibble::tibble(
  level  = c("null", "start", "climate", "climate_start", "species", names(prop_vars)),
  rowset = c(rep("all", 5), prop_rowset[names(prop_vars)]),
  slopes = c(rep(FALSE, 4), rep(TRUE, 8)),
  column = c("Null", "Start", "Climate", "Climate + start", "Species-dependent", prop_label[names(prop_vars)]))
level_mods <- c(list(null = ~ 1, start = f_start, climate = f_env, climate_start = f_both, species = f_both), f_prop)
supp <- paste0(names(prop_vars), "+species") # supplement: property model plus species deviations
level_label <- c(setNames(levels_fit$column, levels_fit$level), setNames(paste(prop_label, "+ species"), supp))
# what readers see (Bills, 2 Nov); level_label stays the key of Fig 6's columns
level_display <- replace(level_label, c("start", "climate_start"), c("Initial abundances", "Climate + initial abundances"))

fixed_part <- function(m, mods, d) { # experiment term silent at prediction: its expected value, 0
  X <- model.matrix(mods, d)
  colnames(X)[colnames(X) == "(Intercept)"] <- "intrcpt" # metafor's name for the intercept
  stopifnot(all(names(coef(m)) %in% colnames(X)))
  as.numeric(X[, names(coef(m)), drop = FALSE] %*% coef(m))
}

# ---- Hand BLUPs against ranef() on the step 05 all-rows Species-dependent fit --
m05 <- readRDS("data/derived/05_ladder_fit/models/all/species.rds")$m
b05 <- blup(m05, rows$all)
m05_int <- m05; m05_int$withG <- FALSE; m05_int$withH <- FALSE # ranef() then returns the two intercept terms from the full M
rf05 <- ranef(m05_int)
blup_diff <- c(species = max(abs(b05$species[rownames(rf05$species)] - rf05$species$intrcpt)),
               experiment = max(abs(b05$experiment[rownames(rf05$experiment_id)] - rf05$experiment_id$intrcpt)))
stopifnot(all(blup_diff < 1e-8)) # Bills, 2 Nov
blup_sd <- c(intercept = sd(b05$species), W = sd(b05$W), T = sd(b05$TPTP[, 1]), P = sd(b05$TPTP[, 2]), TP = sd(b05$TPTP[, 3]))

# ---- One fold, cached: fit on the other gradients, predict the held-out one ----
fit_fold <- function(task) {
  path <- file.path(fold_dir, paste0(task$level, "_", task$gradient_id, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  lv <- filter(levels_fit, level == task$level); mods <- level_mods[[task$level]]; g <- task$gradient_id
  train <- filter(rows[[lv$rowset]], gradient_id != g); test <- filter(rows[[lv$rowset]], gradient_id == g)
  stopifnot(!g %in% train$gradient_id, nrow(test) > 0) # held-out gradient absent from training
  random <- if (lv$slopes) re_slopes else re_base
  struct <- if (lv$slopes) slope_struct else "GEN"
  w <- character(); t0 <- Sys.time()
  try_fit <- function(control) withCallingHandlers(
    tryCatch(fit_model(mods, random, train, struct = struct, control = control), error = function(e) conditionMessage(e)),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  m <- try_fit(list()); optimizer <- "nlminb"; first_error <- NA_character_
  if (is.character(m)) { # retry once from the default start with bobyqa (Bills, 2 Nov)
    first_error <- m; m <- try_fit(list(optimizer = "bobyqa")); optimizer <- "bobyqa" }
  ok <- !is.character(m)
  if (ok && lv$slopes) stopifnot(identical(m$struct, slope_struct)) # T, P, TP block diagonal
  pred <- tibble::tibble(level = task$level, gradient_id = g, effect_id = test$effect_id, species = test$species,
                         experiment_id = test$experiment_id, observed = test$strength,
                         fixed = if (ok) fixed_part(m, mods, test) else NA_real_,
                         deviation = if (ok && lv$slopes) species_deviation(m, train, test) else if (ok) 0 else NA_real_,
                         species_known = test$species %in% train$species)
  stopifnot(all(pred$deviation[!pred$species_known] %in% c(0, NA))) # deviations only for species seen in training
  summary <- tibble::tibble(level = task$level, gradient_id = g, status = if (ok) "converged" else m,
                            optimizer = optimizer, first_error = first_error,
                            struct = if (!lv$slopes) "none (no slopes)" else "DIAG",
                            minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")),
                            tau2 = if (ok) tau2(m) else NA_real_, sigma2_experiment = if (ok) sigma2_exp(m) else NA_real_,
                            sigma2_species = if (ok) m$sigma2[m$s.names == "species"] else NA_real_,
                            tau2_W = if (ok && lv$slopes) m$tau2 else NA_real_,
                            var_T = if (ok && lv$slopes) m$gamma2[1] else NA_real_, var_P = if (ok && lv$slopes) m$gamma2[2] else NA_real_,
                            var_TP = if (ok && lv$slopes) m$gamma2[3] else NA_real_,
                            rho_T_P = NA_real_, rho_T_TP = NA_real_, rho_P_TP = NA_real_) # diagonal: no correlations
  res <- list(pred = pred, summary = summary, warnings = unique(w))
  saveRDS(res, path)
  res
}

# Species-dependent first, then the four levels without species slopes, then the
# properties, so that the floor and ceiling of the ladder are done first (Bills, 2 Nov)
queue <- c("species", "null", "start", "climate", "climate_start", names(prop_vars))
tasks <- tidyr::expand_grid(level = queue, gradient_id = scored)
tasks <- lapply(seq_len(nrow(tasks)), function(i) as.list(tasks[i, ]))
out <- parallel::mclapply(tasks, fit_fold, mc.cores = cores, mc.preschedule = FALSE)
crashed <- vapply(out, function(x) inherits(x, "try-error") || !is.list(x), logical(1))
if (any(crashed)) stop("worker crashed on ", sum(crashed), " folds; rerun resumes from the cached folds:\n",
                       paste(unique(vapply(out[crashed], as.character, "")), collapse = "\n"))
folds <- bind_rows(lapply(out, `[[`, "pred"))
fold_summary <- bind_rows(lapply(out, `[[`, "summary"))
stopifnot(nrow(fold_summary) == nrow(levels_fit) * length(scored), all(folds$species_known),
          !anyDuplicated(fold_summary[c("level", "gradient_id")])) # one fold per gradient per level
write.csv(fold_summary, file.path(out_dir, "fold_summary.csv"), row.names = FALSE)

# ---- Predictions per level: fitted levels and the supplement levels -------------
pred_long <- bind_rows(
  mutate(folds, predicted = fixed + ifelse(level == "species", deviation, 0)),
  filter(folds, level %in% names(prop_vars)) |> mutate(level = paste0(level, "+species"), predicted = fixed + deviation))

# ---- In sample: the same rule on the step 05 full-data fits ---------------------
fit_long <- bind_rows(lapply(seq_len(nrow(levels_fit)), function(i) {
  lv <- levels_fit[i, ]; d <- rows[[lv$rowset]]
  dir <- if (lv$level %in% names(prop_vars)) lv$rowset else "all" # step 05 keeps the ladder of reused rowsets under all
  m <- readRDS(file.path("data/derived/05_ladder_fit/models", dir, paste0(lv$level, ".rds")))$m
  fx <- fixed_part(m, level_mods[[lv$level]], d)
  dev <- if (lv$slopes) species_deviation(m, d, d) else 0
  base <- tibble::tibble(gradient_id = d$gradient_id, effect_id = d$effect_id, observed = d$strength)
  bind_rows(mutate(base, level = lv$level, predicted = fx + if (lv$level == "species") dev else 0),
            if (lv$level %in% names(prop_vars)) mutate(base, level = paste0(lv$level, "+species"), predicted = fx + dev))
}))

# ---- Metrics per level x gradient ---------------------------------------------------
# R2 against the Null on the same gradient (Bills, 2 Nov): V is the Null's error on
# the same effects, its held-out fold for pred and its full-data fit for fit
score <- function(d) {
  d <- inner_join(d, filter(d, level == "null") |> select(effect_id, null_pred = predicted), by = "effect_id")
  bind_rows(lapply(split(d, list(d$level, d$gradient_id), drop = TRUE), function(x) {
    if (anyNA(x$predicted)) met <- tibble::tibble(r2 = NA_real_, spearman = NA_real_, skill = NA_real_)
    else met <- metrics_gradient(x$observed, x$predicted, mean((x$observed - x$null_pred)^2))
    if (x$level[1] == "climate") met$spearman <- met$skill <- NA_real_ # cannot order species within a gradient (Bills, 2 Nov)
    tibble::tibble(level = x$level[1], gradient_id = x$gradient_id[1], met, n = nrow(x))
  }))
}
pred_metrics <- score(pred_long)
fit_metrics  <- score(filter(fit_long, gradient_id %in% scored)) # the same 16 gradients as held out (Bills, 2 Nov)
stopifnot(setequal(pred_metrics$gradient_id, scored), setequal(fit_metrics$gradient_id, scored),
          all(filter(pred_metrics, level == "null")$r2 == 0), all(filter(fit_metrics, level == "null")$r2 == 0)) # Null is the reference
write.csv(pred_metrics, file.path(out_dir, "pred_metrics.csv"), row.names = FALSE)
write.csv(fit_metrics, file.path(out_dir, "fit_metrics.csv"), row.names = FALSE)

all_levels <- c(levels_fit$level, supp)
summary_metrics <- bind_rows(lapply(c("r2", "spearman", "skill"), function(k) {
  f <- fit_metrics |> group_by(level) |> summarise(fit = median(.data[[k]], na.rm = TRUE), n_gradients_fit = sum(!is.na(.data[[k]])))
  p <- pred_metrics |> group_by(level) |>
    summarise(pred = median(.data[[k]], na.rm = TRUE), n_gradients_pred = sum(!is.na(.data[[k]])),
              n_better_than_null = if (k == "r2") sum(.data[[k]] > 0, na.rm = TRUE) else NA_integer_)
  full_join(f, p, by = "level") |> mutate(metric = k, .after = level)
})) |>
  arrange(match(level, all_levels), match(metric, c("r2", "spearman", "skill")))
write.csv(summary_metrics, file.path(out_dir, "summary_metrics.csv"), row.names = FALSE)
# every held-out species was seen in training (checked above with the fold counts)

# Gradients missing per level and metric, with the reason
failed <- filter(fold_summary, status != "converged")
missing <- pred_metrics |> tidyr::pivot_longer(c(r2, spearman, skill), names_to = "metric") |> filter(is.na(value)) |>
  mutate(reason = case_when(
    paste(sub("\\+species$", "", level), gradient_id) %in% paste(failed$level, failed$gradient_id) ~ "fold fit failed",
    level == "climate" & metric != "r2" ~ "Climate: undefined within a gradient",
    metric == "spearman" & n < 5 ~ "fewer than 5 effects",
    metric == "spearman" ~ "one prediction for every effect",
    metric == "skill" ~ "one observed sign")) |>
  group_by(level, metric, reason) |> summarise(gradients = paste(gradient_id, collapse = ", "), n = n(), .groups = "drop") |>
  arrange(match(level, all_levels)) |> mutate(level = level_label[level])

# ---- Figures ------------------------------------------------------------------------
met_name <- c(r2 = "R2", spearman = "Spearman", skill = "Skill")
het05 <- read.csv("data/derived/05_ladder_fit/heterogeneity.csv")
het_fig <- bind_rows(
  filter(het05, rowset == "all", model %in% ladder_columns[1:5]) |> mutate(column = model, ref_cd = NA_real_, ref_H = NA_real_),
  bind_rows(lapply(names(prop_vars), function(p) {
    r <- prop_rowset[[p]]; sp <- filter(het05, rowset == r, model == "Species-dependent")
    filter(het05, rowset == r, model == prop_label[[p]]) |>
      mutate(column = prop_label[[p]], ref_cd = sp$tau2 + sp$sigma2_experiment, ref_H = sp$H_intercept + sp$H_slope) })))
stopifnot(nrow(het_fig) == 12)
# rows b) to d): every scored gradient, in sample (fit) and held out (pred); Null is the R² reference and not drawn
fig4_grad <- bind_rows(mutate(fit_metrics, type = "fit"), mutate(pred_metrics, type = "pred")) |>
  filter(level %in% levels_fit$level, level != "null") |>
  tidyr::pivot_longer(c(r2, spearman, skill), names_to = "metric") |>
  transmute(column = level_label[level], metric = met_name[metric], type, value) |> filter(!is.na(value))
# no caption on the figure (Bills, 2 Nov): the report says what each row shows
fig4_ladder(het_fig, "figures/main/fig6.png", gradients = fig4_grad) # main Fig 6, anchored (Bills, 2 Nov)

# Fig S12 rows: gradients ordered by destination summer temperature (Fig 5 dropped, Bills, 2 Nov)
grad_info <- filter(strength_df, gradient_id %in% scored) |> group_by(gradient_id) |>
  summarise(T_dest = mean(T), spread = sd(strength), .groups = "drop") |> arrange(T_dest)

# Fig S12: the complete per-gradient view, every metric x level, in sample and held out
s7 <- bind_rows(mutate(fit_metrics, version = "In sample"), mutate(pred_metrics, version = "Held out")) |>
  tidyr::pivot_longer(c(r2, spearman, skill), names_to = "metric") |>
  mutate(gradient_id = factor(gradient_id, levels = rev(grad_info$gradient_id)),
         level = factor(level_display[level], levels = level_display[all_levels]),
         version = factor(version, levels = c("In sample", "Held out")))
s7_panel <- function(k, title, tags) lapply(seq_along(levels(s7$version)), function(j) {
  v <- levels(s7$version)[j]; d <- filter(s7, metric == k, version == v)
  ggplot(d, aes(level, gradient_id, fill = pmax(pmin(value, 1), -1))) + geom_tile(colour = "white") +
    geom_text(aes(label = ifelse(is.na(value), "", ifelse(value < -10, "< -10", sprintf("%.2f", value)))), size = 7) + # readable at print size (Bills, 2 Nov)
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", limits = c(-1, 1), na.value = "grey85", name = title) +
    labs(title = tags[j], subtitle = paste0(title, ", ", tolower(v)), x = NULL, y = NULL) +
    theme_minimal(base_size = 30) +
    theme(plot.title = element_text(face = "bold", size = 40), plot.subtitle = element_text(size = 32), axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(size = 26), legend.key.height = unit(2, "cm"), legend.key.width = unit(1, "cm"),
          panel.grid = element_blank(), plot.background = element_rect(fill = "white", colour = NA))
})
figS7 <- patchwork::wrap_plots(c(s7_panel("r2", "R² against Null", c("a)", "b)")), s7_panel("spearman", "Spearman ρ", c("c)", "d)")),
                                 s7_panel("skill", "Loser/winner skill", c("e)", "f)"))), ncol = 2)
ggsave("figures/appendix/figS12.png", figS7, width = 46, height = 48, dpi = 100, bg = "white")

# ---- Report ---------------------------------------------------------------------------
f3 <- function(x) ifelse(is.na(x), "", sprintf("%.3f", x))
sum_tab <- summary_metrics |> transmute(level = level_display[level], metric = c(r2 = "R² against Null", spearman = "Spearman ρ", skill = "loser/winner skill")[metric],
                                        fit = f3(fit), pred = f3(pred), `gradients (fit / pred)` = paste(n_gradients_fit, "/", n_gradients_pred),
                                        `held out better than Null` = ifelse(is.na(n_better_than_null), "", paste(n_better_than_null, "of", n_gradients_pred)))
fold_tab <- fold_summary |> group_by(level = level_display[level]) |>
  summarise(folds = n(), converged = sum(status == "converged"), `bobyqa retries` = sum(optimizer == "bobyqa"),
            `median minutes` = sprintf("%.1f", median(minutes)), `max minutes` = sprintf("%.1f", max(minutes)), .groups = "drop") |>
  arrange(match(level, level_display))
grad_counts <- n_eff |> transmute(gradient_id, effects = n, species = sapply(gradient_id, function(g) n_distinct(strength_df$species[strength_df$gradient_id == g])),
  held_out = ifelse(gradient_id %in% scored, "yes", "no (fewer than 5 effects; trains only)"))
fold_warn <- unique(unlist(lapply(out, function(x) if (length(x$warnings)) paste0(x$summary$level, " / ", x$summary$gradient_id, ": ", x$warnings))))
# metafor's "struct disregarded" note is listed once with its count (CLAUDE.md: warnings are reported, never dropped)
n_struct <- sum(grepl("struct' argument is disregaded", fold_warn))
fold_warn <- c(fold_warn[!grepl("struct' argument is disregaded", fold_warn)],
               if (n_struct) paste0("metafor: \"'struct' argument is disregaded since no ~ inner | outer term\" (", n_struct, " fold fits without slopes)"))
report <- c(
  "# Step 06: The ladder, predicted (Step D, route 2, residents)", "",
  paste0("Every level of the ladder was refitted ", length(scored), " times on the ", nrow(strength_df), " ladder rows of step 05, each time without one of the ", length(scored), " gradients that have at least 5 effects, and asked to predict ",
         "the strengths of the gradient left out: Null, Initial abundances, Climate, Climate + initial abundances, Species-dependent and one level per property, ",
         nrow(fold_summary), " fold fits on ", cores, " cores (", sprintf("%.0f", sum(fold_summary$minutes)), " fit-minutes). ",
         "Predictions are fixed effects only, the experiment term at its expected value zero; Species-dependent adds each species' ",
         "intercept and slope deviations learnt from the training gradients (zero for a species not seen there), and the seven ",
         "supplement levels add the same to each property model. Each level is scored per gradient by R² against the Null model on the same gradient ",
         "(1 − its error / the Null's error, so Null is 0), Spearman ρ within the gradient and loser/winner skill, held out and in sample ",
         "(step 05 full-data fits, the same rule), then summarised by the median over gradients (Bills, 2 Nov)."),
  "", "## Counts", "", paste0(length(gradients), " gradients train; the ", length(scored), " with at least 5 effects are held out and scored, each in one fold (Bills, 2 Nov). Effects and species per gradient:"), "", md_table(grad_counts),
  "", "## Fold fits", "", md_table(fold_tab),
  "", if (nrow(failed)) md_table(transmute(failed, level, gradient_id, `nlminb error` = first_error, `bobyqa error` = status))
      else "Every fold fit converged.",
  if (any(fold_summary$optimizer == "bobyqa")) c("", "Folds that needed the bobyqa retry:", "",
      md_table(filter(fold_summary, optimizer == "bobyqa") |> transmute(level, gradient_id, `nlminb error` = first_error, status))),
  "", "## Prediction metrics, median over gradients", "", md_table(sum_tab),
  "", "Gradients without a value, per level and metric:", "",
  if (nrow(missing)) md_table(missing) else "None.",
  "", "Null's R² is 0 on every scored gradient, in sample and held out (checked): it is the reference of every R².",
  "", "## Species deviations", "",
  paste0("metafor 5.0.1 gives no BLUPs for a GEN or GDIAG slope block, so the species deviations are computed by hand, jointly for every ",
         "random term, as D Z′ M⁻¹ (y − Xβ) from each fit's variance components, V and β (`R/species_deviation.R`). On the step 05 ",
         "all-rows Species-dependent fit they equal `ranef()` for the species and experiment intercepts to ",
         sprintf("%.1e", blup_diff[["species"]]), " and ", sprintf("%.1e", blup_diff[["experiment"]]), " (checked to 1e-8). ",
         "Standard deviation of the BLUPs across species on that fit: intercept ", sprintf("%.2e", blup_sd[["intercept"]]),
         ", W slope ", sprintf("%.2e", blup_sd[["W"]]), ", T slope ", sprintf("%.2e", blup_sd[["T"]]), ", P slope ",
         sprintf("%.2e", blup_sd[["P"]]), ", TP slope ", sprintf("%.2e", blup_sd[["TP"]]),
         ". The P BLUPs are about zero because the P slope variance sits at the boundary (step 05): a BLUP is its variance times ",
         "the species' weighted residual, so in practice the deviation is the species' intercept and its W, T and TP slopes."),
  "", paste0("Every held-out species had training data: each ladder species occurs in at least 3 gradients, so at least 2 remain in training (checked)."),
  "", "## Figures", "", "![Fig 6, the ladder](../figures/main/fig6.png)", "",
  paste0("Fig 6. Row a): step 05, every column on the same ", nrow(strength_df), " rows; dashed lines carry Species-dependent across the property columns. ",
         "Rows b) to d): each half violin is the ", length(scored), " scored gradients (at least 5 effects); grey = explanatory (step 05 fit on all ", length(gradients),
         " gradients), blue = predictive (the fold that left the gradient out); points and numbers = medians. R² = 1 − error(level) / error(Null) on that gradient, ",
         "so Null is 0 and not drawn. Dashed lines from the Species-dependent column on: its explanatory (grey) and predictive (blue) medians. ",
         "Climate has no Spearman or skill: one climate per experiment cannot order species within a gradient."), "",
  "![Fig S12, per-gradient view](../figures/appendix/figS12.png)",
  "", "## Warnings captured", "",
  { w <- unique(c(fold_warn, warns)); if (length(w)) paste("-", w) else "None." },
  "", "## Deviations from PLAN.md", "",
  "- Species deviations from hand BLUPs (`R/species_deviation.R`), not `ranef(m)`: metafor 5.0.1 gives none for the GDIAG slope block (Bills, 2 Nov).",
  "- Climate's Spearman and skill are NA in every gradient (Bills, 2 Nov): with several experiments in a gradient they would rank experiments, not species.",
  "- A fold that does not converge is refitted once with bobyqa, then NA for that level and gradient (Bills, 2 Nov); PLAN.md names no fallback.",
  "- Spearman is also NA where every prediction in the gradient is the same (MILESTONES), beside PLAN.md's fewer than 5 effects.",
  "- The ladder figure is written to `figures/main/fig6.png` (main Fig 6), the per-gradient view to `figures/appendix/`; the ladder's prediction rows widen their y range to the data.")
writeLines(report, "reports/06_ladder_predict.md")
message("Step 06 written: ", nrow(fold_summary), " folds, ", nrow(failed), " failed; ",
        sprintf("%.0f", as.numeric(difftime(Sys.time(), t_start, units = "mins"))), " min")
