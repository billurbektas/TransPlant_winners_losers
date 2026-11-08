# The ladder figures (Bills, 7 Nov: the ladder without abundance; 8 Nov: Fig 6 is R² beside
# the raw heterogeneity, Spearman and skill and the change in unexplained variation go to the appendix).
#   ladder_performance()   Fig 6 (residents, step 06: a) R², b) raw heterogeneity), Fig S12
#                          (residents, step 06: a) Spearman, b) loser/winner skill) and Fig S19
#                          (colonizers, step 07d: a) R², b) raw heterogeneity, c) Spearman)
#   ladder_change()        Fig S13 (step 05) and Fig S20 (step 07c): the change in each
#                          kind of unexplained variation at each comparison of Table 2
#
# het: one row per model, `column` (one of ladder_columns), tau2, sigma2_experiment,
#   H_intercept, H_slope (raw variances). gradients: one row per column x metric x
#   scored gradient x type, `column`, `metric` ("R2", "Spearman", "Skill"), `type`
#   ("fit" = in sample, "pred" = held out), `value`.

ladder_columns <- c("Null", "Experimental warming", "Background climate", "Species-specific responses",
                    "LDMC", "LMA", "Height", "Graminoid", "Woody", "Phylogeny", "Niche")
ladder_labels <- c("Null", "Experimental\nwarming", "Background\nclimate", "Species-\nspecific\nresponses",
                   "LDMC", "LMA", "Height", "Graminoid", "Woody", "Phylogeny", "Niche")
# each model is compared with the one it adds to (Table 2)
ladder_compared_with <- c(`Experimental warming` = "Null", `Background climate` = "Experimental warming",
                          `Species-specific responses` = "Background climate",
                          setNames(rep("Species-specific responses", 7), ladder_columns[5:11]))
ladder_fills <- c(total = "#2B2E35", sigma2_experiment = "#B9BCC6", tau2 = "#6E7380", H_intercept = "#DDB54F", H_slope = "#8C6D1F")

ladder_base <- function() ggplot2::theme_classic(base_size = 38) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 64), # text enlarged twice (Bills, 8 Nov)
                 axis.title.y = ggplot2::element_text(face = "bold", size = 56), axis.text = ggplot2::element_text(size = 48),
                 panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
                 plot.background = ggplot2::element_rect(fill = "white", colour = NA))
ladder_no_x <- function() ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                                         axis.line.x = ggplot2::element_blank())
ladder_x <- function() ggplot2::scale_x_continuous(limits = c(0.4, length(ladder_columns) + 0.6), breaks = seq_along(ladder_columns),
                                                   labels = ladder_labels, expand = c(0, 0))

# ---- the violin panel: one metric against the Null --------------------------------
# half violins per column, left = explanatory (in sample), right = predictive (held
# out), medians drawn and written; R² is against the Null, so Null is not drawn;
# the Species-specific medians run as dashed lines to the right (Bills, 2 Nov)
ladder_metric_panel <- function(gradients, key, tag, ylab, lims, last, legend) {
  `%>%` <- dplyr::`%>%`
  library(ggplot2)
  n <- length(ladder_columns); sp <- "Species-specific responses"; sp_x <- match(sp, ladder_columns)
  g <- dplyr::filter(gradients, metric == key, !is.na(value)) %>% dplyr::mutate(x = match(column, ladder_columns))
  s <- dplyr::summarise(g, med = median(value), .by = c(column, x, type))
  lims <- range(c(lims, g$value), na.rm = TRUE); lims[2] <- lims[2] + 0.14 * diff(lims)
  half <- function(v, x0, side) {
    if (length(unique(v)) < 3) return(NULL)
    d <- density(v, from = min(v), to = max(v), n = 128)
    w <- 0.4 * d$y / max(d$y)
    data.frame(x = x0 + side * c(0, w, 0), y = c(d$x[1], d$x, d$x[128]))
  }
  vio <- dplyr::bind_rows(lapply(split(g, list(g$x, g$type), drop = TRUE), function(z) {
    side <- if (z$type[1] == "fit") -1 else 1
    h <- half(z$value, z$x[1], side); if (is.null(h)) NULL else cbind(h, type = z$type[1], id = paste(z$x[1], z$type[1]))
  }))
  set.seed(1)
  g$xj <- g$x + ifelse(g$type == "fit", -1, 1) * runif(nrow(g), 0.03, 0.13)
  s$xm <- s$x + ifelse(s$type == "fit", -0.06, 0.06)
  sw <- tidyr::pivot_wider(s, id_cols = c(column, x), names_from = type, values_from = med)
  top <- dplyr::summarise(g, top = max(value), .by = x)
  lab_xy <- dplyr::left_join(s, top, by = "x") %>% dplyr::mutate(xl = x + ifelse(type == "fit", -0.24, 0.24), yl = top + 0.08 * diff(lims))
  p <- ggplot() + geom_hline(yintercept = 0, colour = "grey60")
  ceil <- dplyr::filter(s, column == sp)
  if (nrow(ceil)) p <- p + geom_segment(data = ceil, aes(x = sp_x - 0.45, xend = n + 0.45, y = med, yend = med, colour = type),
                                        linetype = "dashed", linewidth = 1, show.legend = FALSE)
  p <- p + geom_polygon(data = vio, aes(x, y, group = id, fill = type), colour = NA, alpha = 0.4) +
    geom_point(data = g, aes(xj, value, colour = type), size = 2.4, alpha = 0.25, show.legend = FALSE) +
    geom_segment(data = sw, aes(x = x - 0.06, xend = x + 0.06, y = fit, yend = pred), colour = "grey40", linewidth = 1.2) +
    geom_point(data = dplyr::filter(s, type == "fit"), aes(xm, med), colour = "#7D8894", size = 6) +
    geom_text(data = lab_xy, aes(xl, yl, label = sprintf("%.2f", round(med, 2) + 0), colour = type), size = 16.5, fontface = "bold", show.legend = FALSE) +
    scale_colour_manual(values = c(fit = "#7D8894", pred = "#2F80D8"), guide = "none") +
    geom_point(data = dplyr::filter(s, type == "pred"), aes(xm, med), colour = "#2F80D8", size = 6) +
    scale_fill_manual(values = c(fit = "#C9D3DE", pred = "#8DB8EA"),
                      labels = c(fit = "Explanatory\n(in sample)", pred = "Predictive\n(held-out gradient)"), name = NULL,
                      guide = if (legend) "legend" else "none") # one key per figure, also when the keys are collected
  p + ladder_x() + scale_y_continuous(limits = lims) + labs(title = tag, x = NULL, y = ylab) + ladder_base() +
    # the legend goes with the first violin panel: inside Spearman, top left, where the floor has no
    # values; under R² it is collected at the bottom of the figure (Bills, 8 Nov)
    theme(legend.position = if (!legend) "none" else if (key == "Spearman") "inside" else "bottom",
          legend.position.inside = c(0.005, 0.99), legend.justification = if (key == "Spearman") c(0, 1) else "left",
          legend.direction = if (key == "Spearman") "vertical" else "horizontal",
          legend.key.spacing.y = unit(0.4, "cm"), legend.key.spacing.x = unit(1.2, "cm"),
          legend.text = element_text(size = 46), legend.key.size = unit(2.4, "cm"),
          legend.background = element_rect(fill = alpha("white", 0.9), colour = NA)) +
    if (last) theme(axis.text.x = element_text(size = 42, face = "bold", colour = "black")) else ladder_no_x()
}

# ---- the raw heterogeneity panel (the old Fig 6a) ------------------------------------
# context dependency (σ² experiment under τ², two greys) beside H species (intercept
# under slope piece, two golds); the Species-specific values run as dashed lines across
# the property columns, since every model is fitted on the same rows. The legend sits
# below the panel so the panel keeps the width of the violin panels above it.
ladder_het_panel <- function(het, tag, last) {
  `%>%` <- dplyr::`%>%`
  library(ggplot2)
  stopifnot(setequal(het$column, ladder_columns))
  n <- length(ladder_columns); sc <- 1e4
  bars <- het %>%
    dplyr::mutate(x = match(column, ladder_columns)) %>%
    tidyr::pivot_longer(c(sigma2_experiment, tau2, H_intercept, H_slope), names_to = "piece", values_to = "value") %>%
    dplyr::mutate(value = sc * value, bar = ifelse(piece %in% c("sigma2_experiment", "tau2"), "cd", "H"),
                  xb = x + ifelse(bar == "cd", -0.17, 0.17),
                  order = match(piece, c("sigma2_experiment", "tau2", "H_intercept", "H_slope"))) %>%
    dplyr::arrange(x, bar, order) %>%
    dplyr::mutate(ymax = cumsum(value), ymin = ymax - value, .by = c(x, bar)) %>%
    dplyr::mutate(piece = factor(piece, levels = c("tau2", "sigma2_experiment", "H_slope", "H_intercept")))
  totals <- dplyr::summarise(bars, value = max(ymax), .by = c(x, bar, xb))
  sp_x <- match("Species-specific responses", ladder_columns)
  refs <- dplyr::filter(bars, x == sp_x, piece %in% c("tau2", "H_intercept", "H_slope")) %>%
    dplyr::transmute(xb, value = ymax, piece = as.character(piece), xend = n + 0.4, mark = "Species-specific responses")
  ref_cols <- c(tau2 = "grey15", H_intercept = "#B08A1E", H_slope = "#5C4510")
  ymax <- max(c(totals$value, refs$value)) * 1.15
  fill_labs <- c(tau2 = "Context dependency:\nspecies × experiment (τ²)",
                 sigma2_experiment = "Context dependency:\nexperiment component (σ² experiment)",
                 H_slope = "H species: slope piece\n(climate response)", H_intercept = "H species: intercept\npiece (average)")
  ggplot() +
    geom_rect(data = bars, aes(xmin = xb - 0.15, xmax = xb + 0.15, ymin = ymin, ymax = ymax, fill = piece)) +
    geom_segment(data = refs, aes(x = xb + 0.15, xend = xend, y = value, yend = value, linetype = mark),
                 colour = ref_cols[refs$piece], linewidth = 1) +
    geom_text(data = totals, aes(xb, value, label = sprintf("%.2f", value)), vjust = -0.5, size = 14, fontface = "bold") +
    scale_fill_manual(values = ladder_fills[names(fill_labs)], labels = fill_labs, name = NULL,
                      guide = guide_legend(order = 1, nrow = 1)) +
    scale_linetype_manual(values = "dashed", guide = "none") + ladder_x() + # the dashed lines are named in the caption (Bills, 8 Nov)
    scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
    labs(title = tag, x = NULL, y = "Heterogeneity\n(unexplained variance, ×10⁻⁴)") +
    ladder_base() + theme(legend.position = "bottom", legend.justification = "center",
                          legend.key.spacing.x = unit(1, "cm"), legend.text = element_text(size = 40),
                          legend.key.size = unit(2.2, "cm"), legend.spacing.y = unit(0.2, "cm"), legend.margin = margin(30, 0, 0, 0)) +
    if (last) theme(axis.text.x = element_text(size = 42, face = "bold", colour = "black")) else ladder_no_x()
}

# ---- Fig 6 / S12 / S19: the panels in the order given ---------------------------------
# panels: from "R2", "het", "Spearman", "Skill"; tagged a), b), ... in that order, the
# column names under the last panel only; no brackets on top (Bills, 8 Nov)
ladder_performance <- function(gradients, file, panels, het = NULL, height) {
  library(ggplot2)
  stopifnot(all(gradients$column %in% ladder_columns), all(panels %in% c("R2", "het", "Spearman", "Skill")),
            !("het" %in% panels) || !is.null(het))
  spec <- list(R2 = list("R²", c(-0.2, 0.8)), Spearman = list("Spearman ρ", c(-0.2, 0.6)), Skill = list("Loser/winner skill", c(-0.2, 0.6)))
  first_violin <- panels[panels != "het"][1]
  rows <- lapply(seq_along(panels), function(i) {
    k <- panels[i]; tag <- paste0(letters[i], ")"); last <- i == length(panels)
    if (k == "het") ladder_het_panel(het, tag, last)
    else ladder_metric_panel(gradients, k, tag, spec[[k]][[1]], spec[[k]][[2]], last, legend = k == first_violin)
  })
  fig <- patchwork::wrap_plots(rows, ncol = 1, heights = ifelse(panels == "het", 1.35, 1))
  # every key at the bottom, one row each, when the first violin panel is R² (Bills, 8 Nov)
  if (first_violin != "Spearman")
    fig <- fig + patchwork::plot_layout(guides = "collect") &
      theme(legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 42))
  ggsave(file, fig, width = 48, height = height, dpi = 100, bg = "white")
  invisible(fig)
}

# ---- Fig S13 / S20: change in unexplained variation at each comparison -----------
# One row per kind of unexplained variation, each on its own scale; every bar is the
# change against the model the column is compared with (Table 2): green label when
# it falls, rust when it rises, with the percentage of the compared model's value;
# "new" when the compared value is zero, a flat tick when both are zero (below 1e-10,
# the REML boundary). Purple arrows at Species-specific responses: the one step where
# variation is moved by design, from τ² and the intercept piece into the slope piece.
# Drawn as one panel with five bands so the arrows can cross the rows.
ladder_change_table <- function(het, zero = 1e-10) {
  `%>%` <- dplyr::`%>%`
  long <- het %>% dplyr::mutate(total = tau2 + sigma2_experiment + H_intercept + H_slope) %>%
    tidyr::pivot_longer(c(total, sigma2_experiment, tau2, H_intercept, H_slope), names_to = "component", values_to = "value")
  now <- dplyr::filter(long, column != "Null") %>% dplyr::mutate(compared_with = unname(ladder_compared_with[column]))
  dplyr::left_join(now, dplyr::select(long, compared_with = column, component, before = value), by = c("compared_with", "component")) %>%
    dplyr::mutate(change = value - before,
                  state = dplyr::case_when(abs(before) < zero & abs(value) < zero ~ "zero",
                                           abs(before) < zero ~ "new", TRUE ~ "change"),
                  percent = ifelse(state == "change", 100 * change / before, NA_real_)) %>%
    dplyr::select(column, compared_with, component, before, value, change, state, percent)
}

ladder_change <- function(het, file) {
  `%>%` <- dplyr::`%>%`
  library(ggplot2)
  stopifnot(setequal(het$column, ladder_columns))
  ch <- ladder_change_table(het)
  comps <- c("total", "sigma2_experiment", "tau2", "H_intercept", "H_slope")
  row_lab <- c(total = "total unexplained variation", sigma2_experiment = "context dependency shared by an experiment (σ² experiment)",
               tau2 = "context dependency, species × experiment (τ²)", H_intercept = "species differ in average response",
               H_slope = "species differ in response to climate")
  cols <- ladder_columns[-1]
  xpos <- setNames(c(1:3, 4.5 + 0:6) * 1.8, cols) # a gap before the properties, as Table 2
  band_h <- 1.5; gap <- 0.8 # room above each row for its name
  d <- ch %>% dplyr::mutate(v = 1e4 * change, x = xpos[column], r = match(component, comps))
  bands <- dplyr::bind_rows(lapply(seq_along(comps), function(i) {
    v <- d$v[d$r == i]; lo <- min(0, v); hi <- max(0, v); span <- max(hi - lo, 1e-3)
    data.frame(r = i, component = comps[i], lo = lo - 0.45 * span, hi = hi + 0.45 * span, # room for the labels
               base = (length(comps) - i) * (band_h + gap))
  }))
  ymap <- function(v, r) { b <- bands[match(r, bands$r), ]; b$base + (v - b$lo) / (b$hi - b$lo) * band_h }
  d <- d %>% dplyr::mutate(y0 = ymap(0, r), y1 = ymap(ifelse(state == "zero", 0, v), r))
  ticks <- dplyr::bind_rows(lapply(seq_len(nrow(bands)), function(i) {
    b <- bands[i, ]; t <- pretty(c(b$lo, b$hi), n = 4); t <- t[t >= b$lo & t <= b$hi]
    data.frame(r = b$r, t = t, y = ymap(t, b$r))
  }))
  green <- "#2E7D32"; rust <- "#B5451B"; purple <- "#7B4FA6"
  # two decimals, more in a row whose changes are all small, so no bar reads +0.00 by rounding
  digits <- dplyr::summarise(d, k = max(2, ceiling(-log10(max(abs(v), 1e-6))) + 1), .by = r)
  labs_d <- dplyr::filter(d, state != "zero") %>% dplyr::left_join(digits, by = "r") %>%
    dplyr::filter(round(v, k) != 0) %>% # a change that rounds to zero keeps its bar, without a label
    dplyr::mutate(lab = paste0(sprintf(paste0("%+.", k, "f"), v), "\n",
                               ifelse(state == "new", "new", sprintf("%+.0f%%", round(percent) + 0))),
                  colour = ifelse(v < 0, green, rust), up = v >= 0,
                  yl = ifelse(up, y1 + 0.04, y1 - 0.04))
  x_left <- 0.45; x_right <- max(xpos) + 0.8; x_tag <- -0.9
  p <- ggplot() +
    geom_segment(data = ticks, aes(x = x_left, xend = x_right, y = y, yend = y), colour = "grey92", linewidth = 0.5) +
    geom_segment(data = bands, aes(x = x_left, xend = x_right, y = ymap(0, r), yend = ymap(0, r)), colour = "grey50", linewidth = 0.7) +
    geom_segment(data = bands, aes(x = x_left, xend = x_left, y = base, yend = base + band_h), colour = "black", linewidth = 0.8) +
    geom_segment(data = ticks, aes(x = x_left - 0.06, xend = x_left, y = y, yend = y), colour = "black", linewidth = 0.6) +
    geom_text(data = ticks, aes(x = x_left - 0.1, y = y, label = format(t, trim = TRUE)), hjust = 1, size = 14, colour = "grey30") +
    geom_segment(aes(x = 3.75 * 1.8, xend = 3.75 * 1.8, y = 0, yend = max(bands$base) + band_h), colour = "grey80", linewidth = 0.8) +
    geom_rect(data = dplyr::filter(d, state != "zero"),
              aes(xmin = x - 0.5, xmax = x + 0.5, ymin = pmin(y0, y1), ymax = pmax(y0, y1), fill = component)) +
    geom_segment(data = dplyr::filter(d, state == "zero"), aes(x = x - 0.5, xend = x + 0.5, y = y0, yend = y0, colour = component),
                 linewidth = 2.2, show.legend = FALSE) +
    geom_text(data = dplyr::filter(labs_d, up), aes(x, yl, label = lab), colour = dplyr::filter(labs_d, up)$colour,
              vjust = 0, size = 13.5, fontface = "bold", lineheight = 0.85) +
    geom_text(data = dplyr::filter(labs_d, !up), aes(x, yl, label = lab), colour = dplyr::filter(labs_d, !up)$colour,
              vjust = 1, size = 13.5, fontface = "bold", lineheight = 0.85) +
    geom_text(data = data.frame(x = xpos, lab = ladder_labels[-1]), aes(x, -0.3, label = lab),
              vjust = 1, size = 15.5, fontface = "bold", lineheight = 0.9)
  # the reassignment by design: τ² and the intercept piece fall, the slope piece is new
  sx <- xpos[["Species-specific responses"]]
  at <- function(comp) dplyr::filter(d, column == "Species-specific responses", component == comp)
  slope <- at("H_slope")
  if (nrow(slope) && slope$v > 0) {
    arrows <- dplyr::bind_rows(lapply(c("tau2", "H_intercept"), function(cmp) {
      a <- at(cmp); if (!nrow(a) || a$v >= 0) return(NULL)
      data.frame(x = sx + 0.54, y = a$y1, xend = sx + 0.54, yend = slope$y1, w = if (cmp == "tau2") 2.4 else 1.3)
    }))
    if (!is.null(arrows) && nrow(arrows))
      p <- p + geom_curve(data = arrows, aes(x = x, y = y, xend = xend, yend = yend, linewidth = I(w)), colour = purple,
                          curvature = -0.15, arrow = arrow(length = unit(0.55, "cm"), type = "closed")) # no "moved" label (Bills, 8 Nov)
  }
  # each row named on its tag line, above the row (Bills, 8 Nov: larger text); drawn after the
  # arrows on white so that the arrows pass behind the names
  p <- p + geom_label(data = bands, aes(x = x_tag, y = base + band_h + 0.22, label = paste0(letters[r], ") ", row_lab[component])),
                      hjust = 0, vjust = 0, size = 17, fontface = "bold", fill = "white", border.colour = NA,
                      label.padding = unit(0.1, "lines"))
  p <- p + scale_fill_manual(values = ladder_fills, guide = "none") + scale_colour_manual(values = ladder_fills, guide = "none") +
    annotate("text", x = x_right, y = max(bands$base) + band_h + 0.9, hjust = 1, size = 14, colour = "grey35",
             label = "change against the model compared with (×10⁻⁴); green: reduced, rust: increased") +
    coord_cartesian(xlim = c(x_tag - 0.05, x_right + 0.05), ylim = c(-1.45, max(bands$base) + band_h + 1.05), expand = FALSE, clip = "off") +
    theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(20, 30, 20, 20))
  ggsave(file, p, width = 48, height = 49, dpi = 100, bg = "white")
  invisible(ch)
}
