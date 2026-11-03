# Fig 4, the ladder (PLAN.md Figure map and caption; 12 columns, Bills 1 Nov).
# Step 06 draws it, rows a) to d); step 10 renders the final figure. Used by steps
# 06 and 10 (step 05 no longer draws it, 2 Nov).
#
# het: one row per column, with `column` (a label of ladder_columns),
#   sigma2_experiment, tau2, H_intercept, H_slope (raw variances); since 2 Nov
#   every column has the same rows, and the Species-dependent totals are drawn
#   as dashed lines across the property columns (ref_cd, ref_H are ignored).
# gradients: NULL (step 05: rows b) to d) NA), or one row per column x metric x
#   scored gradient x type with `column`, `metric` ("R2", "Spearman", "Skill"),
#   `type` ("fit" = in sample, "pred" = held out) and `value`. Since 2 Nov
#   (Bills) rows b) to d) are half violins, explanatory left and predictive right,
#   with the medians drawn and written; R² is against the Null, so Null is not drawn.
#   Display labels only change here ("Initial abundances"); ladder_columns stay the keys.

ladder_columns <- c("Null", "Start", "Climate", "Climate + start", "Species-dependent",
                    "LDMC", "LMA", "Height", "Graminoid", "Woody", "Phylogeny", "Niche")

# labels: display names of the 12 columns (NULL: the residents'); skill: draw the
# loser/winner skill row (FALSE for colonizers, which have no loser class; step 07d)
fig4_ladder <- function(het, file, gradients = NULL, caption = NULL, labels = NULL, skill = TRUE) {
  `%>%` <- dplyr::`%>%`
  library(ggplot2)
  n <- length(ladder_columns)
  if (is.null(labels)) labels <- c("Null", "Initial\nabundances", "Climate", "Climate +\ninitial\nabundances", "Species-\ndependent",
              "LDMC", "LMA", "Height", "Graminoid", "Woody", "Phylogeny", "Niche")
  stopifnot(all(het$column %in% ladder_columns))
  sc <- 1e4 # drawn as x 10^-4
  xs <- scale_x_continuous(limits = c(0.4, n + 0.6), breaks = seq_len(n), labels = labels, expand = c(0, 0))
  base <- theme_classic(base_size = 38) +
    theme(plot.title = element_text(face = "bold", size = 46),
          axis.title.y = element_text(face = "bold", size = 38), axis.text = element_text(size = 32),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.background = element_rect(fill = "white", colour = NA))
  no_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())

  # ---- brackets and comparison arcs (parabolas of fixed height) ----
  arc <- function(from, to, base_y, h, q) {
    t <- seq(0, 1, length.out = 60)
    data.frame(x = from + (to - from) * t, y = base_y + h * 4 * t * (1 - t), q = q, id = paste(q, from, to))
  }
  q3d <- lapply(6:n, function(k) arc(5, k, 0.3, 0.25 + 0.17 * (k - 6), "Q3d"))
  arcs <- rbind(arc(1, 2, 1.45, 0.35, "Q3a"), arc(2, 4, 1.45, 0.5, "Q3b"),
                arc(1, 3, 0.3, 0.5, "Q3b"), arc(4, 5, 0.3, 0.35, "Q3c"), do.call(rbind, q3d))
  ends <- arcs %>% dplyr::group_by(id) %>% dplyr::slice(c(1, dplyr::n())) %>% dplyr::ungroup()
  arc_lab <- data.frame(x = c(1.5, 3, 2, 4.5, 5 + (n - 5) / 2), y = c(1.8, 1.95, 0.8, 0.65, 0.3 + 0.25 + 0.17 * (n - 6) + 0.12),
                        q = c("Q3a", "Q3b", "Q3b", "Q3c", "Q3d"))
  qcol <- c(Q3a = "#B5651D", Q3b = "#2C7A93", Q3c = "#444444", Q3d = "#7442A0")
  top_br <- data.frame(x = c(0.7, 4.75, 5.7), xend = c(4.6, 5.25, n + 0.3), y = 3.35,
                       lab = c("floor: without the species", "ceiling", "gap: one property at a time"))
  sub_br <- data.frame(x = c(5.75, 8.75, 10.75, 11.75), xend = c(8.25, 10.25, 11.25, 12.25), y = 2.55,
                       lab = c("Trait", "Growth form", "Phylogeny", "Niche"))
  p_top <- ggplot() +
    geom_segment(data = top_br, aes(x = x, xend = xend, y = y, yend = y), colour = "grey40", linewidth = 0.9) +
    geom_text(data = top_br, aes(x = (x + xend) / 2, y = y + 0.4, label = lab), size = 13.5, colour = "grey25") +
    geom_segment(data = sub_br, aes(x = x, xend = xend, y = y, yend = y), colour = "grey60", linewidth = 0.7) +
    geom_text(data = sub_br, aes(x = (x + xend) / 2, y = y + 0.3, label = lab), size = 11.5, colour = "grey35") +
    geom_path(data = arcs, aes(x, y, group = id, colour = q), linewidth = 1.1) +
    geom_point(data = ends, aes(x, y, colour = q), size = 3.5) +
    geom_label(data = arc_lab, aes(x, y, label = q, colour = q), size = 11.5, fontface = "bold",
               fill = "white", linewidth = 0, label.padding = unit(0.08, "lines")) +
    scale_colour_manual(values = qcol, guide = "none") + xs + scale_y_continuous(limits = c(0, 3.9)) +
    theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))

  # ---- a) heterogeneity ----
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
  # since 2 Nov every column has the same rows: the Species-dependent totals are
  # drawn as two dashed lines from its two bars across all property columns (Bills)
  sd_x <- match("Species-dependent", ladder_columns)
  # context dependency: its total; H_species: the top of the intercept piece and the
  # top of the slope piece, each in a darker shade of its piece (Bills, 2 Nov)
  refs <- dplyr::filter(bars, x == sd_x, piece %in% c("tau2", "H_intercept", "H_slope")) %>%
    dplyr::transmute(xb, value = ymax, piece = as.character(piece), xend = n + 0.4, mark = "Species-dependent")
  ref_cols <- c(tau2 = "grey15", H_intercept = "#B08A1E", H_slope = "#5C4510")
  na_a <- setdiff(seq_len(n), unique(bars$x))
  ymax <- max(c(totals$value, refs$value)) * 1.15
  fills <- c(tau2 = "#6E7380", sigma2_experiment = "#B9BCC6", H_slope = "#8C6D1F", H_intercept = "#DDB54F")
  fill_labs <- c(tau2 = "Context dependency:\nspecies x experiment (τ²)",
                 sigma2_experiment = "Context dependency:\nshared by site (σ² experiment)",
                 H_slope = "H species: slope piece\n(climate response)",
                 H_intercept = "H species: intercept\npiece (average)")
  pa <- ggplot() +
    geom_rect(data = bars, aes(xmin = xb - 0.15, xmax = xb + 0.15, ymin = ymin, ymax = ymax, fill = piece)) +
    geom_segment(data = refs, aes(x = xb + 0.15, xend = xend, y = value, yend = value, linetype = mark),
                 colour = ref_cols[refs$piece], linewidth = 1) +
    geom_text(data = totals, aes(xb, value, label = sprintf("%.2f", value)), vjust = -0.5, size = 9.5, fontface = "bold") +
    annotate("text", x = na_a, y = ymax * 0.25, label = rep("NA", length(na_a)), size = 12.5, colour = "grey55", fontface = "bold") +
    scale_fill_manual(values = fills, labels = fill_labs, name = NULL, guide = guide_legend(order = 1)) +
    scale_linetype_manual(values = "dashed", name = NULL, guide = guide_legend(order = 2)) + xs +
    scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
    labs(title = "a)", x = NULL, y = "Heterogeneity (×10⁻⁴)") +
    base + no_x + theme(legend.position = "right", legend.justification = c(0, 0.5), legend.key.spacing.y = unit(0.4, "cm"),
                        legend.text = element_text(size = 33), legend.key.size = unit(1.6, "cm"), legend.spacing.y = unit(0, "cm"),
                        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA))

  # ---- b)-d) prediction metrics (step 06) ----
  # half violins per column, left = explanatory (in sample), right = predictive
  # (held out); without gradients (step 05) the row is empty with NA in each column
  metric_row <- function(key, tag, ylab, lims, last) {
    if (is.null(gradients))
      return(ggplot() + annotate("text", x = 2:n, y = 0, label = "NA", size = 12.5, colour = "grey55", fontface = "bold") +
               xs + labs(title = tag, x = NULL, y = ylab) + base +
               if (last) theme(axis.text.x = element_text(size = 33, face = "bold", colour = "black")) else no_x)
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
    lab_xy <- dplyr::left_join(s, top, by = "x") %>% dplyr::mutate(xl = x + ifelse(type == "fit", -0.22, 0.22), yl = top + 0.07 * diff(lims))
    vcol <- c(fit = "#C9D3DE", pred = "#8DB8EA")
    # empty columns stay empty (no "NA"); the Species-dependent medians, explanatory
    # in grey and predictive in blue, run as dashed lines from its column to the right (Bills, 2 Nov)
    p <- ggplot() + geom_hline(yintercept = 0, colour = "grey60")
    sd_x <- match("Species-dependent", ladder_columns)
    ceil <- dplyr::filter(s, column == "Species-dependent")
    if (nrow(ceil)) p <- p + geom_segment(data = ceil, aes(x = sd_x - 0.45, xend = n + 0.45, y = med, yend = med, colour = type),
                                          linetype = "dashed", linewidth = 1, show.legend = FALSE)
    p <- p + geom_polygon(data = vio, aes(x, y, group = id, fill = type), colour = NA, alpha = 0.4) +
      geom_point(data = g, aes(xj, value, colour = type), size = 2.4, alpha = 0.25, show.legend = FALSE) +
      geom_segment(data = sw, aes(x = x - 0.06, xend = x + 0.06, y = fit, yend = pred), colour = "grey40", linewidth = 1.2) +
      geom_point(data = dplyr::filter(s, type == "fit"), aes(xm, med), colour = "#7D8894", size = 5) +
      geom_text(data = lab_xy, aes(xl, yl, label = sprintf("%.2f", round(med, 2) + 0), colour = type), size = 11, fontface = "bold", show.legend = FALSE) +
      scale_colour_manual(values = c(fit = "#7D8894", pred = "#2F80D8"), guide = "none") +
      geom_point(data = dplyr::filter(s, type == "pred"), aes(xm, med), colour = "#2F80D8", size = 5) +
      scale_fill_manual(values = vcol, labels = c(fit = "Explanatory\n(in sample)", pred = "Predictive\n(held-out gradient)"), name = NULL)
    p + xs + scale_y_continuous(limits = lims) + labs(title = tag, x = NULL, y = ylab) + base +
      theme(legend.position = if (key == "R2") "right" else "none", legend.justification = c(0, 0.5), legend.key.spacing.y = unit(0.4, "cm"),
            legend.text = element_text(size = 33), legend.key.size = unit(1.6, "cm")) +
      if (last) theme(axis.text.x = element_text(size = 33, face = "bold", colour = "black")) else no_x
  }
  pb <- metric_row("R2", "b)", "R²", c(-0.2, 0.8), FALSE)
  pc <- metric_row("Spearman", "c)", "Spearman ρ", c(-0.2, 0.6), !skill)
  pd <- if (skill) metric_row("Skill", "d)", "Loser/winner skill", c(-0.2, 0.6), TRUE)

  fig <- if (skill) patchwork::wrap_plots(p_top, pa, pb, pc, pd, ncol = 1, heights = c(1.1, 1.7, 1, 1, 1)) else
    patchwork::wrap_plots(p_top, pa, pb, pc, ncol = 1, heights = c(1.1, 1.7, 1, 1))
  if (!is.null(caption))
    fig <- fig + patchwork::plot_annotation(caption = caption, theme = theme(
      plot.caption = element_text(size = 26, colour = "grey40"), plot.background = element_rect(fill = "white", colour = NA)))
  ggsave(file, fig, width = 38, height = 38, dpi = 100, bg = "white")
  invisible(fig)
}
