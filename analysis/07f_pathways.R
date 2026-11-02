# Step 07f: the pathways (Q4, Fig 7).
#
# PLAN.md: The colonizer story, "Q4, the pathways figure"; Fig 7 caption; the agreed
#   mock plan_figures/mock_fig7_pathways.png.
# Reads:  data/derived/04_classification/species_summary.csv (resident classes);
#         data/derived/07b_classification/species_summary.csv (colonizer classes);
#         data/derived/prep/traits.rds, functional_groups.rds, niche_optima.rds
# Writes: data/derived/07f_pathways/group_properties.csv, group_summary.csv;
#         figures/main/fig7.png; reports/07f_pathways.md
# Checks: group sizes equal the step 04 and 07b class counts (74 residents: 3 classes;
#         28 context-dependent colonizers); every species has growth form; group
#         sizes in the labels.
#
# Decided 2 Nov (Bills): the mock's groups without consistent colonizers, since 07b
# found none: consistent resident winners, context-dependent residents, consistent
# resident losers, context-dependent colonizers. Raw property values (traits logged,
# as the mock's axes); niche optima only where n_obs_sub >= 30, as in the models.
# Descriptive, no model.

suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
source("R/md_table.R")
warns <- character()
globalCallingHandlers(warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
})
out_dir <- "data/derived/07f_pathways"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

res <- read.csv("data/derived/04_classification/species_summary.csv")
col <- read.csv("data/derived/07b_classification/species_summary.csv")
stopifnot(nrow(res) == 74, nrow(col) == 29) # steps 04 and 07b
groups_def <- c("Consistent resident winners", "Context-dependent residents", "Consistent resident losers", "Context-dependent colonizers")
cols <- c("Consistent resident winners" = "#1B9E77", "Context-dependent residents" = "#3A86D6",
          "Consistent resident losers" = "#D95F2B", "Context-dependent colonizers" = "#A68BC8")
members <- bind_rows(
  transmute(res, species, group = c("Consistent winner" = "Consistent resident winners", "Context-dependent" = "Context-dependent residents",
                                    "Consistent loser" = "Consistent resident losers")[class]),
  transmute(filter(col, class == "Context-dependent"), species, group = "Context-dependent colonizers"))
stopifnot(!anyNA(members$group), sum(col$class == "Consistent colonizer") == 0) # no consistent colonizers (07b)
sizes <- count(members, group)
stopifnot(sizes$n[sizes$group == "Context-dependent colonizers"] == 28,
          sum(sizes$n[sizes$group != "Context-dependent colonizers"]) == 74) # class counts of steps 04 and 07b

# ---- Raw properties (units as the mock) -------------------------------------------------
traits <- readRDS("data/derived/prep/traits.rds") |>
  transmute(species, `LDMC (mg g⁻¹, log)` = log(1000 * LDMC), `LMA (g m⁻², log)` = log(1000 * LMA), `Vegetative height (cm, log)` = log(100 * height))
niche <- readRDS("data/derived/prep/niche_optima.rds") |> filter(n_obs_sub >= 30) |>
  transmute(species, `T optimum (°C)` = T_opt, `P optimum (mm, log)` = P_opt)
growth <- readRDS("data/derived/prep/functional_groups.rds") |> transmute(species, graminoid, woody)
props <- members |> left_join(traits, by = "species") |> left_join(niche, by = "species") |> left_join(growth, by = "species")
stopifnot(!anyNA(props$graminoid), !anyNA(props$woody)) # growth form for every species
write.csv(props, file.path(out_dir, "group_properties.csv"), row.names = FALSE)
prop_names <- c("LDMC (mg g⁻¹, log)", "LMA (g m⁻², log)", "Vegetative height (cm, log)", "T optimum (°C)", "P optimum (mm, log)")
long <- props |> tidyr::pivot_longer(all_of(prop_names), names_to = "property") |> filter(!is.na(value)) |>
  mutate(property = factor(property, levels = prop_names), group = factor(group, levels = groups_def))
meds <- summarise(long, median = median(value), n = n(), .by = c(property, group))
shares <- summarise(props, `graminoid share` = mean(graminoid), `woody share` = mean(woody), n = n(), .by = group) |>
  mutate(group = factor(group, levels = groups_def))
write.csv(left_join(shares, tidyr::pivot_wider(select(meds, -n), names_from = property, values_from = median), by = "group"),
          file.path(out_dir, "group_summary.csv"), row.names = FALSE)

# ---- Fig 7: dot strips with group medians, growth-form shares below ------------------
glab <- setNames(paste0(sub(" (residents|resident|colonizers)", "\n\\1", groups_def), " (n = ", sizes$n[match(groups_def, sizes$group)], ")"), groups_def)
set.seed(1)
top <- ggplot(long, aes(group, value, colour = group)) +
  geom_jitter(width = 0.18, height = 0, size = 4, alpha = 0.55) +
  geom_crossbar(data = meds, aes(group, median, ymin = median, ymax = median), width = 0.55, linewidth = 1.4, fatten = 1) +
  facet_wrap(~ property, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = cols, guide = "none") + scale_x_discrete(labels = NULL) +
  labs(title = "a)", x = NULL, y = NULL) +
  theme_classic(base_size = 30) +
  theme(plot.title = element_text(face = "bold", size = 40), strip.text = element_text(face = "bold", size = 28), strip.background = element_blank(),
        axis.ticks.x = element_blank(), axis.line.x = element_blank(), plot.background = element_rect(fill = "white", colour = NA))
share_long <- tidyr::pivot_longer(shares, c(`graminoid share`, `woody share`), names_to = "share")
bottom <- ggplot(share_long, aes(y = 1)) +
  geom_col(data = mutate(share_long, full = 1), aes(x = full), fill = "grey92", width = 0.7, orientation = "y") + # 100% background
  geom_col(aes(x = value, fill = group), width = 0.7, orientation = "y") +
  geom_text(aes(x = value + 0.03, label = sprintf("%.0f%%", 100 * value)), hjust = 0, size = 10, colour = "grey25") +
  facet_grid(share ~ group, switch = "y", labeller = labeller(group = glab)) +
  scale_fill_manual(values = cols, guide = "none") + scale_x_continuous(limits = c(0, 1.25), expand = c(0, 0)) +
  labs(title = "b)", x = NULL, y = NULL) +
  theme_void(base_size = 30) +
  theme(plot.title = element_text(face = "bold", size = 40), strip.text.y.left = element_text(size = 28, angle = 0, hjust = 1),
        strip.text.x = element_text(face = "bold", size = 26, lineheight = 0.9, margin = margin(b = 10)),
        strip.placement = "outside", plot.background = element_rect(fill = "white", colour = NA))
ggsave("figures/main/fig7.png", patchwork::wrap_plots(top, bottom, ncol = 1, heights = c(2.2, 1)), width = 34, height = 20, dpi = 100, bg = "white")

# ---- Report --------------------------------------------------------------------------
f2 <- function(x) sprintf("%.2f", x)
report <- c(
  "# Step 07f: The pathways (Q4, Fig 7)", "",
  paste0("Species-level properties for the classes of Fig 2 (step 04, residents in at least 3 gradients) and Fig S2 (step 07b, colonizers in at least ",
         "3 gradients), as the agreed mock: LDMC, LMA and vegetative height (logged), and the thermal and precipitation optima as dot strips with ",
         "the group median, and the share of graminoids and woody species per group. Step 07b found no consistent colonizers, so the colonizers ",
         "form one group, context-dependent colonizers (Bills, 2 Nov). Descriptive; no model."),
  "", "## Groups", "", md_table(transmute(sizes, group, species = n)),
  "", "## Group medians and shares", "",
  md_table(left_join(transmute(shares, group = as.character(group), `graminoid share` = sprintf("%.0f%%", 100 * `graminoid share`),
                               `woody share` = sprintf("%.0f%%", 100 * `woody share`)),
                     tidyr::pivot_wider(transmute(meds, group = as.character(group), property, median = paste0(f2(median), " (", n, ")")),
                                        names_from = property, values_from = median), by = "group")),
  "", "Medians with the number of species that have the property in brackets.",
  "", "## Fig 7", "", "![Fig 7](../figures/main/fig7.png)",
  "", "## Warnings captured", "", if (length(warns)) paste("-", unique(warns)) else "None.",
  "", "## Deviations from PLAN.md", "",
  "- No consistent-colonizer group: step 07b found none (Bills, 2 Nov: keep the mock without it).")
writeLines(report, "reports/07f_pathways.md")
message("Step 07f written: ", paste(sprintf("%s %d", sizes$group, sizes$n), collapse = "; "))
