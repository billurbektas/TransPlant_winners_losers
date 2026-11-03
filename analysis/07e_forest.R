# Step 07e: sensitivity, random forest, colonizers (held out).
#
# PLAN.md: Sensitivity analyses, "Random forest" (as step 06b), on the colonizer ladder.
# Reads:  data/derived/07c_ladder_fit/ladder_rows.rds; data/derived/07d_ladder_predict/
#         folds/null_<gradient_id>.rds (the metafor Null held out) and pred_metrics.csv
# Writes: data/derived/07e_forest/folds/<setting>/<level>_<gradient_id>.rds, tuning.csv,
#         forest_metrics.csv, forest_vs_metafor.csv; figures/appendix/figS19.png;
#         reports/07e_forest.md
# Checks: 110 colonizer ladder rows (07c); 8 scored gradients, the same as 07d; test rows
#         of every fold equal 07d's; the held-out gradient absent from training; 8
#         gradients per level.
#
# Decided 2 Nov (Bills): step 06b on the colonizer ladder rows (forests without case
# weights, one setting chosen by out-of-bag R² over the levels, 3000 trees, a
# second-seed check, species also in the property forests); A0_origin drops; no
# loser/winner skill (no loser class); an R² below -100 is taken out of the R² row and
# the medians and listed, as 07d; the half-violin code copied from step 06b.

suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
source("R/metrics_gradient.R"); source("R/md_table.R"); source("R/cache.R")
cache_version <- "07e-v2" # v2 (3 Nov): folds return their warnings

warns <- character()
globalCallingHandlers(warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
})
t_start <- Sys.time()
out_dir <- "data/derived/07e_forest"
fold_dir <- file.path(out_dir, "folds")
dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
cores <- 7L # as step 06 (Bills, 2 Nov)

strength_df <- readRDS("data/derived/07c_ladder_fit/ladder_rows.rds") |> mutate(species = factor(species))
stopifnot(nrow(strength_df) == 110, n_distinct(strength_df$species) == 27, n_distinct(strength_df$experiment_id) == 24,
          n_distinct(strength_df$gradient_id) == 18) # step 07c ladder rows
n_eff <- count(strength_df, gradient_id)
scored <- sort(n_eff$gradient_id[n_eff$n >= 5])
pred06 <- read.csv("data/derived/07d_ladder_predict/pred_metrics.csv") # the metafor held-out values (07d)
stopifnot(length(scored) == 8, setequal(scored, pred06$gradient_id)) # 07d's scored gradients

# ---- Predictors per level: the variables of step 06's formulas, no interaction
# terms (the forest builds them); TP left out, it is T x P
v_start <- c("W", "A0_dest"); v_both <- c("W", "T", "P", "A0_dest") # A0_origin drops (PLAN.md)
prop_vars <- list(ldmc = "ldmc_z", lma = "lma_z", height = "height_z", graminoid = "graminoid", woody = "woody",
                  phylogeny = c("PEV1", "PEV2", "PEV3"), niche = c("T_opt", "P_opt"))
level_vars <- c(list(start = v_start, climate = c("W", "T", "P"), climate_start = v_both, species = c(v_both, "species")),
                lapply(prop_vars, function(v) c(v_both, v)),
                setNames(lapply(prop_vars, function(v) c(v_both, v, "species")), paste0(names(prop_vars), "+species")))
stopifnot(length(level_vars) == 18)
prop_label <- c(ldmc = "LDMC", lma = "LMA", height = "Height", graminoid = "Graminoid", woody = "Woody",
                phylogeny = "Phylogeny", niche = "Niche")
level_display <- c(null = "Null", start = "Destination abundance", climate = "Climate", climate_start = "Climate + destination abundance",
                   species = "Species-dependent", prop_label, setNames(paste(prop_label, "+ species"), paste0(names(prop_vars), "+species")))
all_levels <- c("null", names(level_vars))
seeds <- c(main = 20261102L, check = 77L)
n_trees <- 3000L # medians steadier between seeds than with 1000 (scratchpad check, 2 Nov)
# one parameter set for every level and fold, chosen from this grid (Bills, 2 Nov)
grid <- tidyr::expand_grid(node = c(5L, 10L, 20L, 40L, 80L), mtry = c("sqrt", "half")) |>
  mutate(setting = paste0("node", node, "_", mtry))

# ---- One fold: train on the other gradients, predict the held-out one ---------
fit_fold <- function(task, set, seed, cache) {
  path <- file.path(fold_dir, set$setting, paste0(task$level, "_", task$gradient_id, ".rds"))
  g <- task$gradient_id; vars <- level_vars[[task$level]]
  train <- filter(strength_df, gradient_id != g); test <- filter(strength_df, gradient_id == g)
  null06 <- readRDS(file.path("data/derived/07d_ladder_predict/folds", paste0("null_", g, ".rds")))$pred
  stopifnot(!g %in% train$gradient_id, setequal(test$effect_id, null06$effect_id)) # same fold as 07d
  key <- cache_key(train[c(vars, "strength")], test[c("effect_id", vars, "strength")], null06[c("effect_id", "fixed")], set, seed, n_trees, cache_version)
  if (cache) { hit <- read_cache(path, key); if (!is.null(hit)) return(hit) }
  # warnings are collected here and returned: a forked worker cannot write to the parent's list (3 Nov)
  w <- character()
  catch <- function(expr) withCallingHandlers(expr, warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  mtry <- if (set$mtry == "sqrt") floor(sqrt(length(vars))) else ceiling(length(vars) / 2)
  rf <- catch(ranger::ranger(x = train[vars], y = train$strength, num.trees = n_trees, # no case weights (Bills, 2 Nov)
                       mtry = max(1L, mtry), min.node.size = set$node, replace = TRUE,
                       respect.unordered.factors = "order", seed = seed + match(g, scored), num.threads = 1L))
  pred <- tibble::tibble(level = task$level, gradient_id = g, effect_id = test$effect_id, species = as.character(test$species),
                         observed = test$strength, predicted = catch(predict(rf, test[vars], num.threads = 1L))$predictions,
                         null_metafor = null06$fixed[match(test$effect_id, null06$effect_id)],
                         null_forest = mean(train$strength)) # what a forest without predictors predicts
  res <- list(pred = pred, warnings = unique(w))
  if (cache) res <- write_cache(res, path, key)
  res
}
tasks <- tidyr::expand_grid(level = names(level_vars), gradient_id = scored)
tasks <- lapply(seq_len(nrow(tasks)), function(i) as.list(tasks[i, ]))
run <- function(set, seed, cache) {
  out <- parallel::mclapply(tasks, fit_fold, set = set, seed = seed, cache = cache, mc.cores = cores)
  if (any(sapply(out, inherits, "try-error"))) stop("forest worker failed")
  warns <<- c(warns, unlist(lapply(out, `[[`, "warnings")))
  bind_rows(lapply(out, `[[`, "pred"))
}

# ---- Metrics: R² against the metafor Null on the same gradient, as step 06 ------
score <- function(d) bind_rows(lapply(split(d, list(d$level, d$gradient_id), drop = TRUE), function(x) {
  met <- select(metrics_gradient(x$observed, x$predicted, mean((x$observed - x$null_metafor)^2)), r2, spearman) # no skill: no loser class
  if (x$level[1] == "climate") met$spearman <- NA_real_ # cannot order species within a gradient (Bills, 2 Nov)
  tibble::tibble(level = x$level[1], gradient_id = x$gradient_id[1], met, n = nrow(x))
}))

# ---- Tuning: one setting for every level and fold, the one with the highest median
# out-of-bag R² over the 18 levels, each forest fitted once on all 110 rows (Bills,
# 2 Nov). Out-of-bag uses none of the held-out scores, so they stay unbiased.
oob_out <- parallel::mclapply(seq_len(nrow(grid)), function(i) { w <- character(); tab <- withCallingHandlers(bind_rows(lapply(names(level_vars), function(l) {
  vars <- level_vars[[l]]
  mtry <- if (grid$mtry[i] == "sqrt") floor(sqrt(length(vars))) else ceiling(length(vars) / 2)
  rf <- ranger::ranger(x = strength_df[vars], y = strength_df$strength, num.trees = n_trees, mtry = max(1L, mtry),
                       min.node.size = grid$node[i], replace = TRUE, respect.unordered.factors = "order",
                       seed = seeds[["main"]], num.threads = 1L)
  tibble::tibble(setting = grid$setting[i], level = l, oob_r2 = rf$r.squared)
})), warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  list(tab = tab, warnings = w) }, mc.cores = cores)
if (any(sapply(oob_out, inherits, "try-error"))) stop("tuning worker failed")
warns <- c(warns, unlist(lapply(oob_out, `[[`, "warnings")))
oob <- bind_rows(lapply(oob_out, `[[`, "tab"))
tuning <- left_join(grid, summarise(oob, oob_min = min(oob_r2), oob_max = max(oob_r2), oob_r2 = median(oob_r2), .by = setting), by = "setting")
best <- slice_max(tuning, oob_r2, n = 1, with_ties = FALSE)
stopifnot(nrow(tuning) == nrow(grid), nrow(oob) == nrow(grid) * 18, nrow(best) == 1)
write.csv(tuning, file.path(out_dir, "tuning.csv"), row.names = FALSE)

folds <- run(best, seeds[["main"]], cache = TRUE)
stopifnot(nrow(distinct(folds, level, gradient_id)) == 18 * 8)
# the forest Null: the mean of the training effects, one per fold
folds <- bind_rows(folds, distinct(filter(folds, level == "start"), gradient_id, effect_id, .keep_all = TRUE) |>
                     mutate(level = "null", predicted = null_forest))
forest_metrics <- score(folds) |> arrange(match(level, all_levels), gradient_id)
stopifnot(all(count(forest_metrics, level)$n == 8), n_distinct(forest_metrics$level) == 19)
write.csv(forest_metrics, file.path(out_dir, "forest_metrics.csv"), row.names = FALSE)

# seed check: every forest of the chosen setting again with a second seed, not cached
check <- run(best, seeds[["check"]], cache = FALSE)
check_metrics <- score(check)
seed_diff <- inner_join(tidyr::pivot_longer(filter(forest_metrics, level != "null"), c(r2, spearman), names_to = "metric"),
                        tidyr::pivot_longer(check_metrics, c(r2, spearman), names_to = "metric", values_to = "value2"),
                        by = c("level", "gradient_id", "metric")) |> filter(!is.na(value))
seed_med <- summarise(seed_diff, d = abs(median(value) - median(value2)), .by = c(level, metric))

# ---- Forest beside metafor, per gradient and as medians ---------------------------
long <- function(d, model) tidyr::pivot_longer(d, c(r2, spearman), names_to = "metric", values_to = model) |>
  select(level, gradient_id, metric, all_of(model))
# an R² below -100 is taken out of the R² row and the medians and listed (Bills, 2 Nov; as 07d)
outliers <- bind_rows(transmute(filter(pred06, r2 < -100), model = "meta-regression", level, gradient_id, r2),
                      transmute(filter(forest_metrics, r2 < -100), model = "forest", level, gradient_id, r2))
# a gradient taken out on either side is taken out of both, so the two medians of a level are
# over the same gradients (Bills, 3 Nov)
out_pairs <- paste(outliers$level, outliers$gradient_id)
trim <- function(m) mutate(m, r2 = ifelse(paste(level, gradient_id) %in% out_pairs, NA_real_, r2))
vs <- inner_join(long(trim(pred06), "metafor"), long(trim(forest_metrics), "forest"), by = c("level", "gradient_id", "metric")) |>
  mutate(difference = forest - metafor) |> arrange(match(level, all_levels), match(metric, c("r2", "spearman")), gradient_id)
stopifnot(nrow(vs) == 19 * 8 * 2) # every level, gradient and metric
write.csv(vs, file.path(out_dir, "forest_vs_metafor.csv"), row.names = FALSE)
med <- vs |> group_by(level, metric) |>
  summarise(metafor = median(metafor, na.rm = TRUE), forest = median(forest, na.rm = TRUE),
            `median difference` = median(difference, na.rm = TRUE), gradients = sum(!is.na(difference)),
            `forest higher` = sum(difference > 0, na.rm = TRUE), .groups = "drop") |>
  arrange(match(level, all_levels), match(metric, c("r2", "spearman")))

# how far the forest predictions spread within a held-out gradient, against the observed
spread <- filter(folds, level != "null") |>
  summarise(one_sign = all(predicted > 0) | all(predicted < 0), sd_ratio = sd(predicted) / sd(observed), .by = c(level, gradient_id))

# ---- Fig S19: metafor (left) beside forest (right), held out ---------------------
cols <- c(metafor = "#2F80D8", forest = "#E69F00") # blue as Fig 6; Okabe-Ito orange
fills <- c(metafor = "#8DB8EA", forest = "#F3C766")
col_levels <- setdiff(all_levels, "null"); n <- length(col_levels) # Null is the R² reference, not drawn (Bills, 2 Nov)
xlabs <- sub(" \\+ species$", "\n+ species", sub("Climate \\+ destination abundance", "Climate +\ndest.\nabundance",
             sub("Destination abundance", "Destination\nabundance", sub("Species-dependent", "Species-\ndependent", level_display[col_levels]))))
xs <- scale_x_continuous(limits = c(0.4, n + 0.6), breaks = seq_len(n), labels = xlabs, expand = c(0, 0))
base <- theme_classic(base_size = 38) +
  theme(plot.title = element_text(face = "bold", size = 46), axis.title.y = element_text(face = "bold", size = 38),
        axis.text = element_text(size = 32), panel.grid.major.y = element_line(colour = "grey90"),
        plot.background = element_rect(fill = "white", colour = NA))
no_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
half <- function(v, x0, side) { # one half violin, as R/fig4_ladder.R
  if (length(unique(v)) < 3) return(NULL)
  d <- density(v, from = min(v), to = max(v), n = 128); w <- 0.4 * d$y / max(d$y)
  data.frame(x = x0 + side * c(0, w, 0), y = c(d$x[1], d$x, d$x[128]))
}
metric_row <- function(key, tag, ylab, lims, last) {
  g <- filter(vs, metric == key, level != "null") |> tidyr::pivot_longer(c(metafor, forest), names_to = "model") |>
    filter(!is.na(value)) |> mutate(x = match(level, col_levels))
  s <- summarise(g, med = median(value), .by = c(x, model))
  lims <- range(c(lims, g$value)); lims[2] <- lims[2] + 0.14 * diff(lims)
  vio <- bind_rows(lapply(split(g, list(g$x, g$model), drop = TRUE), function(z) {
    h <- half(z$value, z$x[1], if (z$model[1] == "metafor") -1 else 1)
    if (is.null(h)) NULL else cbind(h, model = z$model[1], id = paste(z$x[1], z$model[1]))
  }))
  set.seed(1)
  g$xj <- g$x + ifelse(g$model == "metafor", -1, 1) * runif(nrow(g), 0.03, 0.13)
  s$xm <- s$x + ifelse(s$model == "metafor", -0.06, 0.06)
  sw <- tidyr::pivot_wider(s, id_cols = x, names_from = model, values_from = med)
  lab <- left_join(s, summarise(g, top = max(value), .by = x), by = "x") |>
    mutate(xl = x + ifelse(model == "metafor", -0.24, 0.24), yl = top + 0.07 * diff(lims))
  # empty columns stay empty, no "NA" (as Fig 4, Bills, 2 Nov)
  ggplot() + geom_hline(yintercept = 0, colour = "grey60") + geom_polygon(data = vio, aes(x, y, group = id, fill = model), colour = NA, alpha = 0.4) +
    geom_point(data = g, aes(xj, value, colour = model), size = 2.4, alpha = 0.25, show.legend = FALSE) +
    geom_segment(data = sw, aes(x = x - 0.06, xend = x + 0.06, y = metafor, yend = forest), colour = "grey40", linewidth = 1.2) +
    geom_point(data = s, aes(xm, med, colour = model), size = 5, show.legend = FALSE) +
    geom_text(data = lab, aes(xl, yl, label = sprintf("%.2f", round(med, 2) + 0), colour = model), size = 9, fontface = "bold", show.legend = FALSE) +
    scale_colour_manual(values = cols, guide = "none") +
    scale_fill_manual(values = fills, labels = c(metafor = "Meta-regression\n(held out)", forest = "Random forest\n(held out)"), name = NULL,
                      breaks = c("metafor", "forest")) +
    xs + scale_y_continuous(limits = lims) + labs(title = tag, x = NULL, y = ylab) + base +
    theme(legend.position = if (key == "r2") "right" else "none", legend.justification = c(0, 0.5),
          legend.key.spacing.y = unit(0.4, "cm"), legend.text = element_text(size = 33), legend.key.size = unit(1.6, "cm")) +
    if (last) theme(axis.text.x = element_text(size = 26, face = "bold", colour = "black")) else no_x
}
br <- data.frame(x = c(0.7, 3.75, 4.75, 11.75), xend = c(3.3, 4.25, 11.25, 18.3), y = 1,
                 lab = c("floor: without the species", "ceiling", "one property", "one property + species"))
p_top <- ggplot(br) + geom_segment(aes(x = x, xend = xend, y = y, yend = y), colour = "grey40", linewidth = 0.9) +
  geom_text(aes(x = (x + xend) / 2, y = y + 0.5, label = lab), size = 13.5, colour = "grey25") +
  xs + scale_y_continuous(limits = c(0.5, 2)) + theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))
figS12 <- patchwork::wrap_plots(p_top, metric_row("r2", "a)", "R² against Null", c(-0.2, 0.6), FALSE),
                                metric_row("spearman", "b)", "Spearman ρ", c(-0.2, 0.6), TRUE), ncol = 1, heights = c(0.35, 1, 1)) +
  patchwork::plot_annotation(caption = paste0(
    "Each half violin is the ", length(scored), " scored gradients, held out (leave one gradient out); left (blue) the meta-regression of step 07d, right (orange) the random forest ",
    "on the same rows, predictors and folds, no case weights; points and numbers = medians.\nR² = 1 - error / error of the metafor Null on that gradient, so Null is 0 and not drawn (the forest Null, ",
    "the mean of the training effects, is in the report). Colonizers, no skill (no loser class); Climate has no Spearman. R² below -100 taken out (listed in the report)."),
    theme = theme(plot.caption = element_text(size = 24, colour = "grey40"), plot.background = element_rect(fill = "white", colour = NA)))
ggsave("figures/appendix/figS19.png", figS12, width = 48, height = 22, dpi = 100, bg = "white")

# ---- Report ---------------------------------------------------------------------------
f3 <- function(x) ifelse(is.na(x), "", sprintf("%.3f", x))
med_tab <- med |> transmute(level = level_display[level], metric = c(r2 = "R² against Null", spearman = "Spearman ρ")[metric],
                            metafor = f3(metafor), forest = f3(forest), `median difference` = f3(`median difference`),
                            `forest higher` = ifelse(gradients > 0, paste(`forest higher`, "of", gradients), ""))
count_tab <- tibble::tibble(rule = c("colonizer ladder rows (step 07c)", "species", "experiments",
                                     "gradients training", "gradients held out and scored (at least 5 effects)", "forest levels", "forests (levels x scored gradients)"),
                            n = c(nrow(strength_df), n_distinct(strength_df$species), n_distinct(strength_df$experiment_id),
                                  n_distinct(strength_df$gradient_id), length(scored), length(level_vars), length(tasks)))
r2_main <- filter(med, metric == "r2")
report <- c(
  "# Step 07e: Sensitivity, random forest, colonizers (held out)", "",
  paste0("Every level of the ladder except Null was refitted as a random forest (ranger ", as.character(packageVersion("ranger")), ") on the ", nrow(strength_df),
         " colonizer ladder rows of step 07c, with the variables of the step 07d formulas (A0_origin dropped) as predictors and species as a factor in the Species-dependent and the seven ",
         "property + species forests. Each forest was fitted ", length(scored), " times, once without each scored gradient, and asked to predict it, ",
         "with ", n_trees, " trees, no case weights and one parameter set for every level and fold: minimum node size ", best$node, " and mtry ",
         c(sqrt = "the square root of", half = "half")[[best$mtry]], " the number of predictors, chosen from ", nrow(grid), " settings as the one with the highest ",
         "median out-of-bag R² over the ", length(level_vars), " levels (Bills, 2 Nov). Only held-out predictions are read. ",
         "They are scored as in step 07d: R² against the metafor Null on the same gradient and Spearman ρ within the gradient (no skill: colonizers have no loser class), ",
         "and set beside 07d's held-out values per gradient and as medians. The forest Null is the mean of the training effects."),
  "", "## Counts", "", md_table(count_tab),
  "", "## Tuning: one parameter set for all forests", "",
  paste0("Each setting was fitted once per level on all ", nrow(strength_df), " rows, and scored by its out-of-bag R² (each effect predicted by the trees ",
         "that did not draw it). The setting with the highest median over the ", length(level_vars), " levels (", best$setting, ") is used for every level and fold. ",
         "Out-of-bag R² uses none of the held-out scores, so these stay unbiased; it measures prediction within known gradients (the other effects of a ",
         "gradient stay in training), so its values say nothing about transfer: the forests with species score highest out of bag and lowest held out. ",
         "In-sample R² would always pick the most flexible forest."), "",
  md_table(transmute(tuning, setting, `minimum node size` = node, mtry = c(sqrt = "√p", half = "p / 2")[mtry], `out-of-bag R², median over levels` = f3(oob_r2),
                     `lowest level` = f3(oob_min), `highest level` = f3(oob_max), chosen = ifelse(setting == best$setting, "yes", ""))),
  "", paste0("Spread of the forest predictions: within a held-out gradient their standard deviation is a median ",
             sprintf("%.0f%%", 100 * median(spread$sd_ratio)), " of the observed one, and ", sum(spread$one_sign), " of ", nrow(spread),
             " forests predict one sign for every effect of the gradient."),
  "", "## Forest beside meta-regression, held out, median over gradients", "", md_table(med_tab),
  "", if (nrow(outliers)) paste0("Taken out of the R² row and the medians (R² below −100): ", paste(sprintf("%s, %s at %s (%.1f)", outliers$model, level_display[outliers$level], outliers$gradient_id, outliers$r2), collapse = "; "), ".") else "No R² below −100.",
  "", paste0("R² against Null, median held out: the forest is above the meta-regression at ", sum(r2_main$`median difference` > 0 & r2_main$level != "null"),
             " of ", sum(r2_main$level != "null"), " levels."),
  "", "## Seed check", "",
  paste0("Every forest refitted with a second seed: the largest change in a per-gradient value is ", sprintf("%.3f", max(abs(seed_diff$value - seed_diff$value2))),
         " and the largest change in a median ", sprintf("%.3f", max(seed_med$d)), ". Largest change in a median per metric: ",
         paste0(c("R² ", "Spearman ρ "), sprintf("%.3f", sapply(c("r2", "spearman"), function(k) max(seed_med$d[seed_med$metric == k]))), collapse = ", "),
         "."),
  "", "## Figures", "", "![Fig S19, colonizers: forest beside meta-regression](../figures/appendix/figS19.png)",
  "", "## Warnings captured", "", if (length(unique(warns))) paste("-", unique(warns)) else "None.",
  "", "## Deviations from PLAN.md", "",
  "- Step 06b on the colonizer ladder rows (Bills, 2 Nov): no case weights, held-out values only, out-of-bag tuning, property + species forests, as step 06b.",
  "- No loser/winner skill (no loser class); R² below −100 taken out and listed (Bills, 2 Nov), from both meta-regression and forest for that level and gradient (3 Nov).")
writeLines(report, "reports/07e_forest.md")
message("Step 07e written: ", length(tasks), " forests (and ", length(tasks), " for the seed check); ",
        sprintf("%.1f", as.numeric(difftime(Sys.time(), t_start, units = "mins"))), " min")
