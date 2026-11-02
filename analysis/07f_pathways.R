# Step 07f: the pathways (Q4), placed after Fig 1.
#
# PLAN.md: The colonizer story, "Q4, the pathways figure"; the agreed mock
#   plan_figures/mock_fig7_pathways.png.
# Reads:  data/derived/03_temporal/classes_by_anchor.rds (resident classes, Fig 1a);
#         data/derived/07a_colonizers/presence_classes.csv (colonizer classes, Fig 1b);
#         data/derived/prep/traits.rds, functional_groups.rds, niche_optima.rds
# Writes: data/derived/07f_pathways/group_properties.csv, group_summary.csv;
#         figures/main/fig1_pathways.png; reports/07f_pathways.md
# Checks: the year-3 class counts equal Fig 1's (residents 879, colonizers 331);
#         every species has growth form; group sizes in the labels.
#
# Decided 2 Nov (Bills): with almost nothing consistent in Step C, the groups are
# Fig 1's year-specific classes at year 3, with Fig 1's colours: residents Loser,
# Loser (ns), Winner (ns), Winner; colonizers Present, Present (ns), Absent. A class
# belongs to a species x experiment, so each point is one species in one
# experiment, as Fig 1 counts them; the figure follows Fig 1. Raw property values in
# their units, on log axes except the thermal optimum (Bills, 2 Nov: values on the
# y axis); niche optima only where n_obs_sub >= 30, as in the models. Descriptive,
# no model.

suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
source("R/md_table.R")
warns <- character()
globalCallingHandlers(warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
})
out_dir <- "data/derived/07f_pathways"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Fig 1's classes at year 3 -------------------------------------------------------------
res <- readRDS("data/derived/03_temporal/classes_by_anchor.rds") |> filter(year == 3) |>
  transmute(species, experiment_id, group = as.character(class), pool = "Resident species")
col <- read.csv("data/derived/07a_colonizers/presence_classes.csv") |> filter(year == 3) |>
  transmute(species, experiment_id, group = class, pool = "Colonizing species")
stopifnot(nrow(res) == 879, nrow(col) == 331) # Fig 1 a and b at year 3
groups_def <- c("Loser", "Loser (ns)", "Winner (ns)", "Winner", "Present", "Present (ns)", "Absent")
cols <- c("Loser" = "#D95F2B", "Loser (ns)" = "#F0A07A", "Winner (ns)" = "#8FD3B6", "Winner" = "#1B9E77", # Fig 1a
          "Present" = "#5E3C99", "Present (ns)" = "#A68BC8", "Absent" = "#CFC6DD")                       # Fig 1b
members <- bind_rows(res, col)
stopifnot(all(members$group %in% groups_def))
sizes <- count(members, group) |> arrange(match(group, groups_def))

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
shares <- summarise(props, `graminoid share` = mean(graminoid), `woody share` = mean(woody), n = n(), n_species = n_distinct(species), .by = group) |>
  mutate(group = factor(group, levels = groups_def)) |> arrange(group)
write.csv(left_join(shares, tidyr::pivot_wider(select(meds, -n), names_from = property, values_from = median), by = "group"),
          file.path(out_dir, "group_summary.csv"), row.names = FALSE)

# ---- Figure: dot strips with group medians, growth-form shares below ------------------
glab <- setNames(paste0(groups_def, "\n(n = ", sizes$n[match(groups_def, sizes$group)], ")"), groups_def)
# one panel per property, so that each axis carries its own units; log axes for the skewed properties
set.seed(1)
panels <- lapply(prop_names, function(p) {
  ggplot(filter(long, property == p), aes(group, value, colour = group)) +
    geom_vline(xintercept = 4.5, colour = "grey75", linewidth = 0.8) + # residents | colonizers
    geom_jitter(width = 0.2, height = 0, size = 2.6, alpha = 0.45) +
    geom_crossbar(data = filter(meds, property == p), aes(group, median, ymin = median, ymax = median), width = 0.6, linewidth = 1.6, middle.linewidth = 1.6) +
    (if (p == "T optimum (°C)") scale_y_continuous() else scale_y_log10(breaks = scales::breaks_log(6), labels = scales::label_number(big.mark = ""))) +
    scale_colour_manual(values = cols, guide = "none") + scale_x_discrete(labels = NULL) +
    labs(title = if (p == prop_names[1]) "a)" else " ", subtitle = p, x = NULL, y = NULL) +
    theme_classic(base_size = 30) +
    theme(plot.title = element_text(face = "bold", size = 40), plot.subtitle = element_text(face = "bold", size = 28, hjust = 0.5),
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
ggsave("figures/main/fig1_pathways.png", patchwork::wrap_plots(top, bottom, ncol = 1, heights = c(2.2, 1)), width = 38, height = 20, dpi = 100, bg = "white")

# ---- Report --------------------------------------------------------------------------
report <- c(
  "# Step 07f: The pathways (Q4), after Fig 1", "",
  paste0("Species properties for Fig 1's year-specific classes at year 3 (Bills, 2 Nov), with Fig 1's colours: residents as Loser, Loser (ns), ",
         "Winner (ns) and Winner (step 03), colonizers as Present, Present (ns) and Absent (step 07a). A class belongs to a species x experiment, so each point is ",
         "one species in one experiment, as Fig 1 counts them. As the agreed mock: LDMC, LMA and vegetative height, and the thermal and ",
         "precipitation optima in their units (log axes except the thermal optimum) as dot strips with the group median, and the share of graminoids and woody species per group. Descriptive; no model."),
  "", "## Groups (year 3)", "", md_table(transmute(shares, group = as.character(group), `species x experiment` = n, species = n_species)),
  "", "## Group medians and shares", "",
  md_table(left_join(transmute(shares, group = as.character(group), `graminoid share` = sprintf("%.0f%%", 100 * `graminoid share`),
                               `woody share` = sprintf("%.0f%%", 100 * `woody share`)),
                     tidyr::pivot_wider(transmute(meds, group = as.character(group), property, median = paste0(signif(median, 3), " (", n, ")")),
                                        names_from = property, values_from = median), by = "group")),
  "", "Medians with the number of species x experiment that have the property in brackets; shares over species x experiment.",
  "", "## Figure", "", "![Pathways](../figures/main/fig1_pathways.png)",
  "", "## Warnings captured", "", if (length(warns)) paste("-", unique(warns)) else "None.",
  "", "## Deviations from PLAN.md", "",
  "- Groups are Fig 1's year-3 classes, not Step C's consistent classes, which were nearly empty (Bills, 2 Nov); points are species x experiment.",
  "- The figure follows Fig 1 (`figures/main/fig1_pathways.png`); final numbering in step 10.")
writeLines(report, "reports/07f_pathways.md")
message("Step 07f written: ", paste(sprintf("%s %d", sizes$group, sizes$n), collapse = "; "))
