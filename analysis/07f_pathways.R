# Step 07f: the pathways (Q4), Fig 4.
#
# PLAN.md: The colonizer story, "Q4, the pathways figure", its fallback; the agreed
#   mock plan_figures/mock_fig7_pathways.png.
# Reads:  data/derived/04_classification/species_summary.csv (residents, Fig 2);
#         data/derived/07b_classification/species_summary.csv (colonizers, Fig S2);
#         data/derived/prep/traits.rds, functional_groups.rds, niche_optima.rds
# Writes: data/derived/07f_pathways/group_properties.csv, group_summary.csv;
#         figures/main/fig4.png; reports/07f_pathways.md
# Checks: 74 resident species (step 04) split by the sign of their net strength;
#         28 colonizer species with net warmed cover above zero (07b); every species
#         has growth form; group sizes in the labels.
#
# Decided 2 Nov (Bills): PLAN.md's fallback, since Step C found almost nothing
# consistent: residents grouped by the sign of their pooled net strength mu (Fig 2),
# losers and winners, against colonizers with pooled net warmed cover above zero
# (Fig S2), one point per species; Fig 2's loser and winner colours, Fig S2's purple.
# Placed after Fig 2 and S2 as Fig 4 (Fig 3 is the climate space). Raw property values
# in their units, on log axes except the thermal optimum; niche optima only where
# n_obs_sub >= 30, as in the models. Descriptive, no model.

suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
source("R/md_table.R")
warns <- character()
globalCallingHandlers(warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
})
out_dir <- "data/derived/07f_pathways"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Groups by the sign of the net effect (Fig 2, Fig S2) -------------------------------
res <- read.csv("data/derived/04_classification/species_summary.csv")
stopifnot(nrow(res) == 74, all(res$mu != 0)) # step 04's resident species, every one with a sign
res <- transmute(res, species, mu, group = ifelse(mu < 0, "Resident losers", "Resident winners"))
col <- read.csv("data/derived/07b_classification/species_summary.csv") |> filter(!is.na(mu), mu > 0) # net warmed cover above zero
stopifnot(nrow(col) == 28) # 07b's colonizers less Ranunculus acris, which never arrived
col <- transmute(col, species, mu, group = "Colonizers")
groups_def <- c("Resident losers", "Resident winners", "Colonizers")
cols <- c("Resident losers" = "#D95F2B", "Resident winners" = "#1B9E77", "Colonizers" = "#5E3C99") # Fig 2, Fig S2
members <- bind_rows(res, col)
sizes <- count(members, group) |> arrange(match(group, groups_def))
both <- intersect(res$species, col$species) # residents at some sites, colonizers at others

# ---- Raw properties (units as the mock) -------------------------------------------------
traits <- readRDS("data/derived/prep/traits.rds") |>
  transmute(species, `LDMC (mg g⁻¹)` = 1000 * LDMC, `LMA (g m⁻²)` = 1000 * LMA, `Vegetative height (cm)` = 100 * height)
niche <- readRDS("data/derived/prep/niche_optima.rds") |> filter(n_obs_sub >= 30) |>
  transmute(species, `T optimum (°C)` = T_opt, `P optimum (mm)` = exp(P_opt)) # P_opt is log mean bio18 (prep/R/niche.R)
growth <- readRDS("data/derived/prep/functional_groups.rds") |> transmute(species, graminoid, woody)
props <- members |> left_join(traits, by = "species") |> left_join(niche, by = "species") |> left_join(growth, by = "species")
stopifnot(nrow(props) == nrow(members), !anyNA(props$graminoid), !anyNA(props$woody)) # growth form for every species
write.csv(props, file.path(out_dir, "group_properties.csv"), row.names = FALSE)
prop_names <- c("LDMC (mg g⁻¹)", "LMA (g m⁻²)", "Vegetative height (cm)", "T optimum (°C)", "P optimum (mm)")
long <- props |> tidyr::pivot_longer(all_of(prop_names), names_to = "property") |> filter(!is.na(value)) |>
  mutate(property = factor(property, levels = prop_names), group = factor(group, levels = groups_def))
meds <- summarise(long, median = median(value), n = n(), .by = c(property, group))
shares <- summarise(props, `graminoid share` = mean(graminoid), `woody share` = mean(woody), n = n(), .by = group) |>
  mutate(group = factor(group, levels = groups_def)) |> arrange(group)
write.csv(left_join(shares, tidyr::pivot_wider(select(meds, -n), names_from = property, values_from = median), by = "group"),
          file.path(out_dir, "group_summary.csv"), row.names = FALSE)

# ---- Figure: dot strips with group medians, growth-form shares below ------------------
glab <- setNames(paste0(groups_def, "\n(n = ", sizes$n[match(groups_def, sizes$group)], ")"), groups_def)
# one panel per property, so that each axis carries its own units; log axes for the skewed properties
set.seed(1)
panels <- lapply(prop_names, function(p) {
  ggplot(filter(long, property == p), aes(group, value, colour = group)) +
    geom_vline(xintercept = 2.5, colour = "grey75", linewidth = 0.8) + # residents | colonizers
    geom_jitter(width = 0.2, height = 0, size = 4, alpha = 0.6) +
    geom_crossbar(data = filter(meds, property == p), aes(group, median, ymin = median, ymax = median), width = 0.6, linewidth = 1.6, middle.linewidth = 1.6) +
    (if (p == "T optimum (°C)") scale_y_continuous() else scale_y_log10(breaks = scales::breaks_log(6), labels = scales::label_number(big.mark = ""))) +
    scale_colour_manual(values = cols, guide = "none") + scale_x_discrete(labels = NULL) +
    labs(title = if (p == prop_names[1]) "a)" else " ", x = NULL, y = p) +
    theme_classic(base_size = 30) +
    theme(plot.title = element_text(face = "bold", size = 40), axis.title.y = element_text(face = "bold", size = 30),
          axis.ticks.x = element_blank(), axis.line.x = element_blank(), plot.background = element_rect(fill = "white", colour = NA))
})
top <- patchwork::wrap_plots(panels, nrow = 1)
share_long <- tidyr::pivot_longer(shares, c(`graminoid share`, `woody share`), names_to = "share")
bottom <- ggplot(share_long, aes(y = 1)) +
  geom_col(data = mutate(share_long, full = 1), aes(x = full), fill = "grey92", width = 0.7, orientation = "y") + # 100% background
  geom_col(aes(x = value, fill = group), width = 0.7, orientation = "y") +
  geom_text(aes(x = value + 0.03, label = sprintf("%.0f%%", 100 * value)), hjust = 0, size = 9, colour = "grey25") +
  facet_grid(share ~ group, switch = "y", labeller = labeller(group = glab)) +
  scale_fill_manual(values = cols, guide = "none") + scale_x_continuous(limits = c(0, 1.35), expand = c(0, 0)) +
  labs(title = "b)", x = NULL, y = NULL) +
  theme_void(base_size = 30) +
  theme(plot.title = element_text(face = "bold", size = 40), strip.text.y.left = element_text(size = 28, angle = 0, hjust = 1),
        strip.text.x = element_text(face = "bold", size = 25, lineheight = 0.9, margin = margin(b = 10)),
        strip.placement = "outside", plot.background = element_rect(fill = "white", colour = NA))
ggsave("figures/main/fig4.png", patchwork::wrap_plots(top, bottom, ncol = 1, heights = c(2.2, 1)), width = 32, height = 20, dpi = 100, bg = "white")

# ---- Report --------------------------------------------------------------------------
report <- c(
  "# Step 07f: The pathways (Q4), Fig 4", "",
  paste0("PLAN.md's fallback for the pathways figure (Bills, 2 Nov), since Step C found almost nothing consistent: resident species grouped by the sign of ",
         "their pooled net strength from step 04 (Fig 2), losers below zero and winners above, against colonizer species whose pooled net warmed cover ",
         "from 07b (Fig S2) is above zero; one point per species. As the agreed mock: LDMC, LMA, vegetative height, and the thermal and precipitation ",
         "optima in their units (log axes except the thermal optimum) as dot strips with the group median, and the share of graminoids and woody species ",
         "per group. ", length(both), " species are residents at some sites and colonizers at others, and appear in both pools. Descriptive; no model."),
  "", "## Groups", "", md_table(tibble::tibble(rule = c("resident species (step 04)", "net strength below zero", "net strength above zero",
                                                       "colonizer species (07b)", "net warmed cover above zero", "in both pools"),
                                              n = c(nrow(res), sum(res$mu < 0), sum(res$mu > 0), nrow(col) + 1, nrow(col), length(both)))),
  "", "## Group medians and shares", "",
  md_table(left_join(transmute(shares, group = as.character(group), species = n, `graminoid share` = sprintf("%.0f%%", 100 * `graminoid share`),
                               `woody share` = sprintf("%.0f%%", 100 * `woody share`)),
                     tidyr::pivot_wider(transmute(meds, group = as.character(group), property, median = paste0(signif(median, 3), " (", n, ")")),
                                        names_from = property, values_from = median), by = "group")),
  "", "Medians with the number of species that have the property in brackets.",
  "", paste0("In both pools: ", paste(both, collapse = ", "), "."),
  "", "## Figure", "", "![Fig 4, pathways](../figures/main/fig4.png)",
  "", "## Warnings captured", "", if (length(warns)) paste("-", unique(warns)) else "None.",
  "", "## Deviations from PLAN.md", "",
  "- PLAN.md's fallback, with resident losers added beside the winners (Bills, 2 Nov): groups by the sign of the net effect, not Step C classes.",
  "- The figure is Fig 4, after Fig 2 and Fig S2 (Bills, 2 Nov); written to `figures/main/fig4.png`.")
writeLines(report, "reports/07f_pathways.md")
message("Step 07f written: ", paste(sprintf("%s %d", sizes$group, sizes$n), collapse = "; "))
