# Step 07a: colonizer strength and presence flows (Steps A and B, colonizers).
#
# PLAN.md: The colonizer story, "Strength: arrival, not change" and "Temporal
#   consistency"; Figure map and Fig 1 caption (panel b).
# Reads:  data/derived/01_relative_cover/cover_long.rds (destination pool);
#         data/derived/02_strength/strength_df.rds (W, T, P, TP per experiment);
#         data/derived/03_temporal/classes_by_anchor.rds, stability.csv,
#         flows.rds (Fig 1 panel a, redrawn unchanged)
# Writes: data/derived/07a_colonizers/models/<experiment_id>.rds, by_year.rds,
#         strength_df.rds/.csv (all years), strength_df_y3.rds/.csv (years >= 3),
#         presence_classes.csv, stability.csv, never_present.csv;
#         figures/main/fig1.png (panel a unchanged, panel b below);
#         figures/diagnostics/07a_stability.png; reports/07a_colonizers.md
# Checks: 402 colonizer effects, 236 species, 39 experiments, 20 gradients;
#         every fit converged with a positive-definite Hessian (R/fit_experiment.R
#         stops otherwise); one strength per species x experiment, SEs finite and
#         positive; absence read from raw warmed cover; no A0_origin anywhere;
#         class counts per anchor sum to the combinations surveyed there.
#
# Decided 2 Nov (Bills): step 07 split into 07a, 07b, 07c; the residents' model
# (glmmTMB, one residual variance per species, R/fit_experiment.R) with the
# destination control as reference; strength averaged over all years and over
# years 3 and later, both kept until Bills reads the stability; Fig 1 panel b
# below panel a, three classes (absent, present (ns), present); the 47 species
# never present in a warmed plot reported as a count.

suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
source("R/fit_experiment.R"); source("R/md_table.R")
invisible(loadNamespace("glmmTMB")) # emmeans reads cached glmmTMB fits only with it loaded (as 03b)

warns <- character()
globalCallingHandlers(warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
}, message = function(m) { # emmeans notes, for the report
  warns <<- c(warns, trimws(conditionMessage(m))); invokeRestart("muffleMessage")
})
out_dir <- "data/derived/07a_colonizers"
dir.create(file.path(out_dir, "models"), recursive = TRUE, showWarnings = FALSE)
anchors <- c(2L, 3L, 4L, 6L)

cover_long <- readRDS("data/derived/01_relative_cover/cover_long.rds")
colonizers <- filter(cover_long, species_pool == "destination_pool", treatment %in% c("destination_control", "warmed"))
effects <- distinct(colonizers, experiment_id, gradient_id, species)
stopifnot(nrow(effects) == 402, n_distinct(effects$species) == 236,
          n_distinct(effects$experiment_id) == 39, n_distinct(effects$gradient_id) == 20) # step 01, destination pool

# ---- Step A: one model per experiment, destination control as reference --------
experiment_ids <- sort(unique(colonizers$experiment_id))
fits <- lapply(experiment_ids, function(e) {
  f <- file.path(out_dir, "models", paste0(gsub("::", "__", e), ".rds"))
  if (!file.exists(f)) saveRDS(fit_experiment(filter(colonizers, experiment_id == e), control_level = "destination_control"), f)
  readRDS(f)
})
names(fits) <- experiment_ids

# marginal means of both groups from the same fit: warmed = strength, destination
# control = A0_dest (PLAN.md); averaged over all years, or over years >= 3 (as 03b)
# an experiment with one colonizer species has no species term (R/fit_experiment.R): its name is added back
with_species <- function(fit, x) if (fit$one_species) cbind(species = names(fit$sigma_species), x) else x
means_of <- function(fit, years) {
  em <- emmeans::emmeans(fit$model, if (fit$one_species) ~ treatment else ~ treatment | species,
                         at = list(year = years), df = fit$df_plots) |> as.data.frame() |> with_species(fit = fit)
  w <- filter(em, treatment == "warmed"); d <- filter(em, treatment == "destination_control")
  tibble::tibble(species = as.character(w$species), strength = w$emmean, strength_se = w$SE,
                 A0_dest = d$emmean[match(w$species, d$species)], A0_dest_se = d$SE[match(w$species, d$species)], df = fit$df_plots)
}
strength_of <- function(from_year) bind_rows(lapply(experiment_ids, function(e) {
  years <- levels(fits[[e]]$model$frame$year)
  keep <- years[as.integer(years) >= from_year]
  if (!length(keep)) return(NULL) # no survey in year 3 or later: no year-3 strength for this experiment
  mutate(means_of(fits[[e]], keep), experiment_id = e, n_years = length(keep))
}))
covariates <- distinct(readRDS("data/derived/02_strength/strength_df.rds"), experiment_id, gradient_id, W, T, P, TP)
finish <- function(s) s |>
  left_join(covariates, by = "experiment_id") |>
  mutate(effect_id = paste(experiment_id, species, sep = "::")) |>
  select(effect_id, species, experiment_id, gradient_id, strength, strength_se, A0_dest, A0_dest_se, df, n_years, W, T, P, TP)
strength_df <- finish(strength_of(1L))
strength_df_y3 <- finish(strength_of(3L))
stopifnot(nrow(strength_df) == 402, !anyDuplicated(strength_df$effect_id), !anyNA(strength_df$W), # every colonizer effect once
          all(is.finite(strength_df$strength_se)), all(strength_df$strength_se > 0),
          all(strength_df_y3$effect_id %in% strength_df$effect_id), !any(grepl("A0_origin", c(names(strength_df), names(strength_df_y3)))))
n_y3_missing <- nrow(strength_df) - nrow(strength_df_y3)

# ---- Step B: warmed cover per year, presence classes at the anchors ------------
by_year <- bind_rows(lapply(experiment_ids, function(e) {
  em <- emmeans::emmeans(fits[[e]]$model, if (fits[[e]]$one_species) ~ treatment | year else ~ treatment | species * year,
                         df = fits[[e]]$df_plots)
  as.data.frame(emmeans::test(em, side = ">")) |> with_species(fit = fits[[e]]) |> # p-value for warmed cover above zero
    filter(treatment == "warmed") |>
    transmute(experiment_id = e, species = as.character(species), year = as.integer(as.character(year)), estimate = emmean, SE, p.value)
}))
stopifnot(!anyDuplicated(by_year[c("experiment_id", "species", "year")]))
present <- colonizers |> filter(treatment == "warmed") |> mutate(year = as.integer(as.character(year))) |>
  summarise(present = any(relative_cover > 0), .by = c(experiment_id, species, year))
class_levels <- c("Present", "Present (ns)", "Absent") # top down in Fig 1b (Bills, 2 Nov)
classes <- by_year |> filter(year %in% anchors) |>
  left_join(present, by = c("experiment_id", "species", "year")) |>
  mutate(class = factor(case_when(!present ~ "Absent", # read from the raw warmed cover (PLAN.md)
                                  p.value < 0.05 ~ "Present", TRUE ~ "Present (ns)"), levels = class_levels))
stopifnot(!anyNA(classes$present), !any(classes$year == 1L))
class_counts <- count(classes, year, class) |> tidyr::pivot_wider(names_from = class, values_from = n, values_fill = 0)
stopifnot(all(rowSums(class_counts[intersect(class_levels, names(class_counts))]) == as.integer(table(classes$year))))

# presence agreement and CCC of warmed cover between adjacent anchors
pair <- function(y1, y2, cl) {
  wide <- cl |> filter(year %in% c(y1, y2)) |>
    tidyr::pivot_wider(id_cols = c(species, experiment_id), names_from = year, values_from = c(estimate, class)) |>
    tidyr::drop_na()
  p1 <- wide[[paste0("class_", y1)]] != "Absent"; p2 <- wide[[paste0("class_", y2)]] != "Absent"
  tibble::tibble(pair = paste(y1, "to", y2), n = nrow(wide), presence_agreement = mean(p1 == p2),
                 ccc = DescTools::CCC(wide[[paste0("estimate_", y1)]], wide[[paste0("estimate_", y2)]])$rho.c$est,
                 flows = list(table(factor(wide[[paste0("class_", y1)]], class_levels), factor(wide[[paste0("class_", y2)]], class_levels))))
}
stability <- bind_rows(pair(2L, 3L, classes), pair(3L, 4L, classes), pair(4L, 6L, classes))
stopifnot(all(sapply(stability$flows, sum) == stability$n))

never <- colonizers |> filter(treatment == "warmed") |>
  summarise(any = any(relative_cover > 0), .by = c(experiment_id, gradient_id, species)) |>
  summarise(never = all(!any), n_gradients = n_distinct(gradient_id), gradients = paste(sort(unique(gradient_id)), collapse = ", "),
            n_experiments = n(), .by = species) |> filter(never) |> select(-never) |> arrange(desc(n_gradients), species)
stopifnot(nrow(never) == 47, all(never$n_gradients < 3)) # none reaches the Step C rule (Bills, 2 Nov)

# ---- Write ----------------------------------------------------------------------
saveRDS(by_year, file.path(out_dir, "by_year.rds"))
saveRDS(strength_df, file.path(out_dir, "strength_df.rds")); write.csv(strength_df, file.path(out_dir, "strength_df.csv"), row.names = FALSE)
saveRDS(strength_df_y3, file.path(out_dir, "strength_df_y3.rds")); write.csv(strength_df_y3, file.path(out_dir, "strength_df_y3.csv"), row.names = FALSE)
write.csv(select(classes, experiment_id, species, year, estimate, p.value, present, class), file.path(out_dir, "presence_classes.csv"), row.names = FALSE)
write.csv(select(stability, pair, n, presence_agreement, ccc), file.path(out_dir, "stability.csv"), row.names = FALSE)
write.csv(never, file.path(out_dir, "never_present.csv"), row.names = FALSE)

# ---- Fig 1: panel a as step 03 draws it (unchanged), panel b the colonizers below --
# copied from analysis/03_temporal.R:102-164 (draw_fig1), with the classes,
# colours, sides and box text as arguments so that panel b can use it
draw_panel <- function(classes, stability, class_levels, cols, side, strong, box_text, box_title, tag) {
  gap <- 0.025 * max(table(classes$year)); w <- 0.08; xpos <- setNames(seq_along(anchors), anchors)
  strata <- count(classes, year, class) |> arrange(year, class) |>
    mutate(top = -(cumsum(n) - n + (as.integer(class) - 1) * gap), bottom = top - n, x = xpos[as.character(year)], .by = year)
  ribbons <- list(); switch_labels <- list()
  for (k in seq_len(nrow(stability))) {
    y1 <- anchors[k]; y2 <- anchors[k + 1]; fl <- as.data.frame(stability$flows[[k]], stringsAsFactors = FALSE)
    names(fl) <- c("from", "to", "n")
    fl <- filter(fl, n > 0) |> mutate(from = factor(from, class_levels), to = factor(to, class_levels))
    src <- filter(strata, year == y1); dst <- filter(strata, year == y2)
    fl <- fl |> arrange(from, to) |> mutate(s0 = src$top[match(from, src$class)] - (cumsum(n) - n), .by = from) |>
      arrange(to, from) |> mutate(t0 = dst$top[match(to, dst$class)] - (cumsum(n) - n), .by = to)
    t <- seq(0, 1, length.out = 40); s <- (1 - cos(pi * t)) / 2
    x <- xpos[as.character(y1)] + w / 2 + t * (xpos[as.character(y2)] - xpos[as.character(y1)] - w)
    for (i in seq_len(nrow(fl))) {
      upper <- fl$s0[i] + (fl$t0[i] - fl$s0[i]) * s
      ribbons[[length(ribbons) + 1]] <- tibble::tibble(id = paste(k, i), class = fl$from[i],
        switch = side(fl$from[i]) != side(fl$to[i]) & (fl$from[i] %in% strong | fl$to[i] %in% strong),
        x = c(x, rev(x)), y = c(upper, rev(upper - fl$n[i])))
    }
    switch_labels[[k]] <- fl |> filter(side(from) != side(to)) |>
      mutate(share = 100 * n / src$n[match(from, src$class)]) |>
      transmute(x = xpos[as.character(y2)] - w / 2 - 0.01, y = t0 - n / 2, class = from,
                label = sprintf("%d (%s%%)", n, ifelse(share < 1, sprintf("%.1f", share), sprintf("%.0f", share))))
  }
  ribbons <- bind_rows(ribbons)
  switch_labels <- bind_rows(switch_labels) |> arrange(x, desc(y)) |>
    mutate(y = Reduce(function(prev, cur) min(cur, prev - 1.3 * gap), y, accumulate = TRUE), .by = x)
  low <- min(strata$bottom)
  boxes <- tibble::tibble(x = head(xpos, -1) + 0.5, y = low - 6.5 * gap, label = box_text)
  ggplot() +
    geom_polygon(data = filter(ribbons, !switch), aes(x, y, group = id, fill = class), alpha = 0.35) +
    geom_polygon(data = filter(ribbons, switch), aes(x, y, group = id, fill = class), alpha = 0.95) +
    geom_rect(data = strata, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = bottom, ymax = top, fill = class)) +
    geom_text(data = filter(strata, n >= 40), aes(x, (top + bottom) / 2, label = n), colour = "white", size = 3.8, fontface = "bold") +
    geom_label(data = switch_labels, aes(x, y, label = label, colour = class), hjust = 1, size = 3.6,
               fontface = "bold", fill = "white", label.size = 0, label.padding = unit(0.08, "lines")) +
    geom_text(data = filter(strata, year == 6L), aes(x + w / 2 + 0.03, (top + bottom) / 2, label = class, colour = class),
              hjust = 0, size = 5.8, fontface = "bold") +
    annotate("text", x = xpos, y = 3 * gap, label = paste("Year", anchors), size = 6, fontface = "bold", colour = "grey30") +
    geom_label(data = boxes, aes(x, y, label = label), size = 5.8, label.size = 0.3, colour = "grey20", lineheight = 1.1) +
    annotate("text", x = mean(boxes$x), y = low - 1.6 * gap, size = 5, fontface = "bold", colour = "grey30", label = box_title) +
    scale_fill_manual(values = cols, guide = "none") + scale_colour_manual(values = cols, guide = "none") +
    coord_cartesian(xlim = c(0.85, length(anchors) + 0.75), ylim = c(low - 9.5 * gap, 4 * gap), clip = "off") +
    labs(title = tag) +
    theme_void() + theme(plot.title = element_text(face = "bold", size = 24, colour = "grey20"), plot.margin = margin(10, 10, 10, 10))
}
res_levels <- c("Loser", "Loser (ns)", "Winner (ns)", "Winner")
res_classes <- readRDS("data/derived/03_temporal/classes_by_anchor.rds") |> mutate(class = factor(class, levels = res_levels))
res_stab <- read.csv("data/derived/03_temporal/stability.csv") |> mutate(flows = unname(readRDS("data/derived/03_temporal/flows.rds")))
pa <- draw_panel(res_classes, res_stab, res_levels,
                 c("Loser" = "#D95F2B", "Loser (ns)" = "#F0A07A", "Winner (ns)" = "#8FD3B6", "Winner" = "#1B9E77"),
                 side = function(cl) ifelse(cl %in% c("Loser", "Loser (ns)"), "loser side", "winner side"), strong = c("Loser", "Winner"),
                 box_text = sprintf("Sign agreement: %.2f\nConcordance correlation\ncoefficient: %.2f", res_stab$sign_agreement, res_stab$ccc),
                 box_title = "Agreement of losing-winning strength between adjacent years", tag = "a)   Resident species")
pb <- draw_panel(classes, stability, class_levels,
                 c("Absent" = "#CFC6DD", "Present (ns)" = "#A68BC8", "Present" = "#5E3C99"), # colonizers purple (PLAN.md)
                 side = function(cl) ifelse(cl == "Absent", "absent", "present"), strong = "Present",
                 box_text = sprintf("Presence agreement: %.2f\nConcordance correlation\ncoefficient: %.2f", stability$presence_agreement, stability$ccc),
                 box_title = "Agreement of colonizer warmed cover between adjacent years", tag = "b)   Colonizing species")
ggsave("figures/main/fig1.png", patchwork::wrap_plots(pa, pb, ncol = 1), width = 14, height = 17, dpi = 150, bg = "white")

# ---- Diagnostic: all years against years >= 3, and agreement by anchor pair ---------
cmp <- inner_join(select(strength_df, effect_id, all_years = strength), select(strength_df_y3, effect_id, from_year_3 = strength), by = "effect_id")
lim <- range(c(cmp$all_years, cmp$from_year_3))
d1 <- ggplot(cmp, aes(all_years, from_year_3)) + geom_abline(colour = "grey60") + geom_point(alpha = 0.5, size = 2, colour = "#5E3C99") +
  coord_equal(xlim = lim, ylim = lim) + labs(title = "a)", x = "Warmed cover, all years", y = "Warmed cover, years 3 and later") +
  theme_classic(base_size = 18) + theme(plot.title = element_text(face = "bold", size = 24))
agree <- bind_rows(transmute(stability, group = "Colonizers", pair, agreement = presence_agreement, ccc),
                   transmute(res_stab, group = "Residents", pair, agreement = sign_agreement, ccc)) |>
  tidyr::pivot_longer(c(agreement, ccc), names_to = "measure") |>
  mutate(measure = c(agreement = "Presence (colonizers) or sign (residents) agreement", ccc = "Concordance correlation coefficient")[measure])
d2 <- ggplot(agree, aes(pair, value, colour = group, group = group)) + geom_line(linewidth = 1.2) + geom_point(size = 4) +
  facet_wrap(~ measure, ncol = 1) + scale_colour_manual(values = c(Colonizers = "#5E3C99", Residents = "#1B9E77"), name = NULL) +
  scale_y_continuous(limits = c(0, 1)) + labs(title = "b)", x = "Adjacent anchors (years)", y = NULL) +
  theme_classic(base_size = 18) + theme(plot.title = element_text(face = "bold", size = 24), legend.position = "top",
                                        strip.text = element_text(face = "bold", size = 15), strip.background = element_blank())
ggsave("figures/diagnostics/07a_stability.png", patchwork::wrap_plots(d1, d2, widths = c(1, 1.2)), width = 18, height = 9, dpi = 100, bg = "white")

# ---- Report -------------------------------------------------------------------------
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
f3 <- function(x) sprintf("%.3f", x)
fit_tab <- tibble::tibble(experiment_id = experiment_ids, species = sapply(fits, function(f) nlevels(f$model$frame$species)),
                          years = sapply(fits, function(f) paste(levels(f$model$frame$year), collapse = ", ")),
                          df_plots = sapply(fits, `[[`, "df_plots"), warnings = sapply(fits, function(f) length(f$warnings)))
anchor_tab <- count(classes, year, name = "combinations") |> left_join(class_counts, by = "year")
fit_warn <- unique(unlist(lapply(experiment_ids, function(e) if (length(fits[[e]]$warnings)) paste0(e, ": ", fits[[e]]$warnings))))
report <- c(
  "# Step 07a: Colonizer strength and presence flows (Steps A and B, colonizers)", "",
  paste0("Each of the ", length(experiment_ids), " experiments with destination-pool species got one mixed model of relative cover on its ",
         "destination-control and warmed plots, the residents' model (glmmTMB, one residual variance per species, plot x species random ",
         "intercept; `R/fit_experiment.R`) with the destination control as reference. Colonizer strength is the warmed marginal mean, zeros kept, ",
         "and A0_dest the destination-control marginal mean from the same fit, each averaged over all years and, separately, over years 3 and ",
         "later; both are kept until the choice is made (Bills, 2 Nov). Each colonizer was classified at years 2, 3, 4 and 6 as absent (no cover ",
         "in any warmed plot that year, from the raw data), present (warmed mean above zero at one-sided p < 0.05) or present (ns), and ",
         "presence agreement and the concordance correlation of warmed cover were computed between adjacent anchors."),
  "", "## Counts", "",
  md_table(tibble::tibble(rule = c("colonizer effects (species x experiment)", "species", "experiments", "gradients",
                                   "effects with a strength from years 3 and later", "species never present in a warmed plot (any gradient)"),
                          n = c(nrow(strength_df), n_distinct(strength_df$species), length(experiment_ids), n_distinct(strength_df$gradient_id),
                                nrow(strength_df_y3), nrow(never)))),
  "", paste0(n_y3_missing, " effects have no strength from years 3 and later (their experiment has no survey in year 3 or later)."),
  "", "## Fits", "", md_table(fit_tab),
  "", "## Presence classes per anchor", "", md_table(anchor_tab),
  "", "## Stability between anchors (`stability.csv`)", "",
  md_table(mutate(select(stability, pair, n, presence_agreement, ccc), presence_agreement = f3(presence_agreement), ccc = f3(ccc))),
  "", paste0("Residents (step 03), for the year-3 read: ", paste(sprintf("%s sign agreement %.2f, CCC %.2f", res_stab$pair, res_stab$sign_agreement, res_stab$ccc), collapse = "; "), "."),
  "", paste0("All years against years 3 and later, on the ", nrow(cmp), " effects with both: Pearson r ", f3(cor(cmp$all_years, cmp$from_year_3)),
             ", concordance correlation ", f3(DescTools::CCC(cmp$all_years, cmp$from_year_3)$rho.c$est), "."),
  "", "![Stability](../figures/diagnostics/07a_stability.png)",
  "", "## Fig 1", "", "![Fig 1](../figures/main/fig1.png)", "",
  paste0("Panel a is step 03's, redrawn from its saved classes and flows with the same code; panel b the colonizers below it (Bills, 2 Nov). ",
         "Ribbons between absent and present are labelled where they enter the next anchor, with their count and their share of the class they started from; ",
         "those that start or end in Present are drawn darker."),
  "", "## Never present", "",
  paste0(nrow(never), " colonizer species never appear in a warmed plot of any experiment where they are in the pool; ",
         sum(never$n_gradients == 1), " of them occur in one gradient and none in three or more, so none enters Step C (`never_present.csv`)."),
  "", "## Warnings captured", "", { w <- unique(c(fit_warn, warns)); if (length(w)) paste("-", w) else "None." },
  "", "## Deviations from PLAN.md", "",
  "- The residents' glmmTMB model with one residual variance per species, not PLAN.md's `lmer` for colonizers (Bills, 2 Nov).",
  "- Strength kept twice, all years and years 3 and later, until the stability is read (Bills, 2 Nov).",
  "- Presence is tested one-sided (warmed mean above zero); PLAN.md's `test()` call leaves the side open.",
  "- Figures written to `figures/main/` and `figures/diagnostics/` so that the report can link them.")
writeLines(report, "reports/07a_colonizers.md")
message("Step 07a written: ", length(experiment_ids), " fits; ", paste(sprintf("%s CCC %.2f", stability$pair, stability$ccc), collapse = "; "))
