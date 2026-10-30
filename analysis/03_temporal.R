# Step 03: temporal consistency (Step B, residents).
#
# PLAN.md: Step B (classes at anchors 2, 4, 6; sign agreement and Lin's CCC
#   between adjacent anchors; flows for Fig 1); Figure map and Fig 1 caption.
# Reads:  data/derived/02_strength/by_year.rds;
#         data/derived/01_relative_cover/cover_long.rds (absence per year)
# Writes: data/derived/03_temporal/classes_by_anchor.rds, stability.csv,
#         flows.rds; figures/diagnostics/03_fig1_residents_draft.png;
#         reports/03_temporal.md
# Checks: 7,571 yearly contrasts; 3,571 at the anchors (1,655, 1,150, 766),
#         of which 489 absent from both groups (232, 125, 132), leaving 1,423,
#         1,025 and 634 classified; year 1 absent from every output; class
#         counts per anchor sum to the combinations at that anchor; n per pair
#         equals the combinations present at both anchors, counted separately.
#
# Decided 30 Oct (Bills): a resident with no cover in any origin-control or
# warmed plot in an anchor year is not classified at that anchor (no stratum,
# no ribbon, not in sign agreement or CCC); absence is read from the raw cover,
# not from estimate == 0, as PLAN.md does for colonizers. Anchor labels
# "Year 2/4/6". The draft figure goes to figures/diagnostics/.

suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
source("R/md_table.R")

warns <- character()
globalCallingHandlers(warning = function(w) {
  warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
})
out_dir <- "data/derived/03_temporal"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
anchors <- c(2L, 4L, 6L)
class_levels <- c("Loser", "Loser (ns)", "Winner (ns)", "Winner")

by_year <- as.data.frame(readRDS("data/derived/02_strength/by_year.rds"))
stopifnot(nrow(by_year) == 7571) # step 02

# ---- Absence from both groups in a year, from the raw cover ---------------
present <- readRDS("data/derived/01_relative_cover/cover_long.rds") |>
  filter(species_pool == "resident", treatment %in% c("origin_control", "warmed")) |>
  summarise(present = any(relative_cover > 0), .by = c(experiment_id, species, year)) |>
  mutate(year = as.integer(as.character(year)))

at_anchor <- by_year |>
  mutate(year = as.integer(as.character(year))) |>   # emmeans returns year as a factor
  filter(year %in% anchors) |>
  left_join(present, by = c("experiment_id", "species", "year"))
stopifnot(nrow(at_anchor) == 3571, !anyNA(at_anchor$present))
stopifnot(identical(as.integer(table(at_anchor$year)), c(1655L, 1150L, 766L))) # combinations per anchor
absent_counts <- count(filter(at_anchor, !present), year, name = "absent")
stopifnot(identical(absent_counts$absent, c(232L, 125L, 132L))) # absent from both groups (decided 30 Oct)

# ---- Classes (PLAN.md Step B) ----------------------------------------------
classes <- at_anchor |>
  filter(present) |>
  mutate(class = case_when(
    estimate > 0 & p.value < 0.05 ~ "Winner",
    estimate < 0 & p.value < 0.05 ~ "Loser",
    estimate > 0 ~ "Winner (ns)",
    TRUE ~ "Loser (ns)"),
    class = factor(class, levels = class_levels))
stopifnot(identical(as.integer(table(classes$year)), c(1423L, 1025L, 634L))) # classified per anchor
stopifnot(!any(classes$year == 1L))
class_counts <- count(classes, year, class) |> tidyr::pivot_wider(names_from = class, values_from = n, values_fill = 0)
stopifnot(all(rowSums(class_counts[class_levels]) == c(1423, 1025, 634))) # classes sum to the anchor's combinations

# ---- Sign agreement, CCC and flows between adjacent anchors -----------------
pair <- function(y1, y2) {
  wide <- classes |>
    dplyr::filter(year %in% c(y1, y2)) |>
    tidyr::pivot_wider(id_cols = c(species, experiment_id),
                       names_from = year, values_from = c(estimate, class)) |>
    tidyr::drop_na()
  d1 <- wide[[paste0("estimate_", y1)]]; d2 <- wide[[paste0("estimate_", y2)]]
  tibble::tibble(
    pair = paste(y1, "to", y2),
    n = nrow(wide),
    sign_agreement = mean(sign(d1) == sign(d2)),
    ccc = DescTools::CCC(d1, d2)$rho.c$est,
    flows = list(table(wide[[paste0("class_", y1)]], wide[[paste0("class_", y2)]]))
  )
}
stability <- dplyr::bind_rows(pair(2L, 4L), pair(4L, 6L))
both <- function(y1, y2) nrow(inner_join(filter(classes, year == y1), filter(classes, year == y2),
                                         by = c("species", "experiment_id")))
stopifnot(stability$n == c(both(2L, 4L), both(4L, 6L))) # n per pair = present at both anchors
stopifnot(all(sapply(stability$flows, sum) == stability$n))

# ---- Write -------------------------------------------------------------------
saveRDS(select(classes, species, experiment_id, year, estimate, class), file.path(out_dir, "classes_by_anchor.rds"))
utils::write.csv(select(stability, pair, n, sign_agreement, ccc), file.path(out_dir, "stability.csv"), row.names = FALSE)
saveRDS(setNames(stability$flows, stability$pair), file.path(out_dir, "flows.rds"))

# ---- Fig 1 panel a (draft) ---------------------------------------------------
# Strata are stacked bars on one absolute count scale (top down: Loser to
# Winner); a ribbon joins the combinations present at two adjacent anchors.
# Drawn by hand: a stratum also holds combinations without a ribbon.
cols <- c("Loser" = "#D95F2B", "Loser (ns)" = "#F0A07A", "Winner (ns)" = "#8FD3B6", "Winner" = "#1B9E77")
gap <- 0.025 * max(table(classes$year)); w <- 0.08; xpos <- c(`2` = 1, `4` = 2, `6` = 3)
strata <- count(classes, year, class) |>
  arrange(year, class) |>
  mutate(top = -(cumsum(n) - n + (as.integer(class) - 1) * gap), bottom = top - n, x = xpos[as.character(year)], .by = year)
ribbons <- list(); labels <- list()
for (k in seq_len(nrow(stability))) {
  y1 <- anchors[k]; y2 <- anchors[k + 1]; fl <- as.data.frame(stability$flows[[k]], stringsAsFactors = FALSE)
  names(fl) <- c("from", "to", "n")
  fl <- filter(fl, n > 0) |> mutate(from = factor(from, class_levels), to = factor(to, class_levels))
  src <- filter(strata, year == y1); dst <- filter(strata, year == y2)
  fl <- fl |>
    arrange(from, to) |> mutate(s0 = src$top[match(from, src$class)] - (cumsum(n) - n), .by = from) |>
    arrange(to, from) |> mutate(t0 = dst$top[match(to, dst$class)] - (cumsum(n) - n), .by = to)
  t <- seq(0, 1, length.out = 40); s <- (1 - cos(pi * t)) / 2
  x <- xpos[as.character(y1)] + w / 2 + t * (xpos[as.character(y2)] - xpos[as.character(y1)] - w)
  for (i in seq_len(nrow(fl))) {
    upper <- fl$s0[i] + (fl$t0[i] - fl$s0[i]) * s
    ribbons[[length(ribbons) + 1]] <- tibble::tibble(id = paste(k, i), class = fl$from[i],
                                                    x = c(x, rev(x)), y = c(upper, rev(upper - fl$n[i])))
  }
  # ribbons that cross zero, summed by direction, with the share of their starting side
  side <- function(cl) ifelse(cl %in% c("Loser", "Loser (ns)"), "loser side", "winner side")
  cross <- fl |> filter(side(from) != side(to)) |> mutate(dir = side(from)) |>
    summarise(y = mean(s0 + (t0 - s0) / 2 - n / 2), n = sum(n), .by = dir)
  start_n <- src |> mutate(dir = side(class)) |> summarise(n = sum(n), .by = dir)
  labels[[k]] <- mutate(cross, x = mean(xpos[c(as.character(y1), as.character(y2))]),
                        label = sprintf("%d (%.0f%%)", n, 100 * n / start_n$n[match(dir, start_n$dir)]))
}
ribbons <- bind_rows(ribbons); labels <- bind_rows(labels)
low <- min(strata$bottom)
boxes <- tibble::tibble(x = c(1.5, 2.5), y = low - 4 * gap,
                        label = sprintf("sign agreement %.2f\nCCC %.2f", stability$sign_agreement, stability$ccc))
fig <- ggplot() +
  geom_polygon(data = ribbons, aes(x, y, group = id, fill = class), alpha = 0.35) +
  geom_rect(data = strata, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = bottom, ymax = top, fill = class)) +
  geom_text(data = filter(strata, n >= 40), aes(x, (top + bottom) / 2, label = n), colour = "white", size = 3) +
  geom_text(data = labels, aes(x, y, label = label), size = 2.8, colour = "grey25") +
  geom_text(data = filter(strata, year == 6L), aes(x + w / 2 + 0.03, (top + bottom) / 2, label = class, colour = class),
            hjust = 0, size = 3.2) +
  annotate("text", x = xpos, y = 2.5 * gap, label = paste("Year", anchors), size = 3.6, colour = "grey30") +
  geom_label(data = boxes, aes(x, y, label = label), size = 3, label.size = 0.3, colour = "grey30") +
  scale_fill_manual(values = cols, guide = "none") + scale_colour_manual(values = cols, guide = "none") +
  coord_cartesian(xlim = c(0.85, 3.55), clip = "off") +
  labs(title = "a Residents: warmed minus origin contrast (draft)") +
  theme_void() + theme(plot.title = element_text(face = "bold", size = 12, colour = "grey30"),
                       plot.margin = margin(10, 10, 10, 10))
ggsave("figures/diagnostics/03_fig1_residents_draft.png", fig, width = 9, height = 5.5, dpi = 150)

# ---- Report --------------------------------------------------------------
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
flow_md <- function(tab, title) {
  d <- as.data.frame.matrix(tab); d <- cbind(from = rownames(d), d)
  c(paste0("**", title, "** (rows: class at the first anchor; columns: class at the second)"), "", md_table(d), "")
}
anchor_table <- tibble::tibble(year = anchors, combinations = as.integer(table(at_anchor$year)),
                               absent_both_groups = absent_counts$absent) |>
  mutate(classified = combinations - absent_both_groups) |>
  left_join(class_counts, by = "year")
report <- c(
  "# Step 03: Temporal consistency (Step B, residents)", "",
  paste0("Each resident species x experiment was classified at years 2, 4 and 6 after transplanting from its ",
         "year-specific warmed minus origin-control contrast (step 02, p-values with the plot-based df): Winner or Loser ",
         "at p < 0.05 by the sign of the estimate, otherwise Winner (ns) or Loser (ns). A combination with no cover in ",
         "any origin-control or warmed plot in that year was not classified there (decided 30 Oct). Between adjacent ",
         "anchors, for the combinations present at both, the share keeping its sign and Lin's concordance correlation ",
         "of the two contrasts were computed, with the 4 x 4 flow counts that draw Fig 1a. Year 1 is excluded."),
  "", "## Counts per anchor", "", md_table(anchor_table),
  "", "## Stability between anchors (`stability.csv`)", "",
  md_table(mutate(select(stability, pair, n, sign_agreement, ccc),
                  sign_agreement = sprintf("%.3f", sign_agreement), ccc = sprintf("%.3f", ccc))),
  "", "## Flows (`flows.rds`)", "",
  flow_md(stability$flows[[1]], "Year 2 to year 4"), flow_md(stability$flows[[2]], "Year 4 to year 6"),
  "## Fig 1 panel a (draft)", "",
  "![Fig 1a draft](../figures/diagnostics/03_fig1_residents_draft.png)", "",
  paste0("Strata are drawn to one absolute scale, so columns shrink as experiments end and as species are absent from ",
         "both groups. Ribbon labels give the combinations that crossed from the loser side (Loser, Loser (ns)) to the ",
         "winner side or back, with their share of the starting side. The final figure, with panel b (colonizers, ",
         "step 07), is drawn in step 10."),
  "", "## Warnings captured", "", if (length(warns)) paste("-", unique(warns)) else "None.",
  "", "## Deviations from PLAN.md", "",
  paste("Combinations absent from both groups in an anchor year (", sum(absent_counts$absent), " of ", fmt(nrow(at_anchor)),
        " anchor contrasts) are not classified there, where PLAN.md's rule would have put their zero contrast in",
        "Loser (ns) (decided 30 Oct). Anchors are labelled Year 2, 4, 6 rather than the mock's First 2, 4, 6 years.",
        "The draft figure is written to figures/diagnostics/ rather than data/derived/03_temporal/, so that it is",
        "committed with the report.")
)
writeLines(report, "reports/03_temporal.md")
message("Step 03 written: ", paste(sprintf("%s n = %d", stability$pair, stability$n), collapse = "; "))
