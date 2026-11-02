# Fig 4, the ladder (PLAN.md Figure map and caption; 12 columns, Bills 1 Nov).
# Step 05 draws row a) from heterogeneity.csv with rows b) to d) NA; step 06 adds
# the prediction metrics; step 10 renders the final figure. Used by steps 05, 06, 10.
#
# het: one row per column, with `column` (a label of ladder_columns),
#   sigma2_experiment, tau2, H_intercept, H_slope (raw variances); since 2 Nov
#   every column has the same rows, and the Species-dependent totals are drawn
#   as dashed lines across the property columns (ref_cd, ref_H are ignored).
# metrics: NULL, or one row per column x metric with `column`, `metric`
#   ("R2", "Spearman", "Skill"), `fit`, `pred`.
# gradients: NULL, or one row per column x metric x held-out gradient with
#   `column`, `metric`, `pred`.

ladder_columns <- c("Null", "Start", "Climate", "Climate + start", "Species-dependent",
                    "LDMC", "LMA", "Height", "Graminoid", "Woody", "Phylogeny", "Niche")

fig4_ladder <- function(het, file, metrics = NULL, gradients = NULL, caption = NULL) {
  `%>%` <- dplyr::`%>%`
  library(ggplot2)
  n <- length(ladder_columns)
  labels <- c("Null", "Start", "Climate", "Climate\n+ start", "Species-\ndependent",
              "LDMC", "LMA", "Height", "Graminoid", "Woody", "Phylogeny", "Niche")
  stopifnot(all(het$column %in% ladder_columns))
  sc <- 1e4 # drawn as x 10^-4
  xs <- scale_x_continuous(limits = c(0.4, n + 0.6), breaks = seq_len(n), labels = labels, expand = c(0, 0))
  base <- theme_classic(base_size = 22) +
    theme(plot.title = element_text(face = "bold", size = 28),
          axis.title.y = element_text(face = "bold", size = 22), axis.text = element_text(size = 18),
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
    geom_text(data = top_br, aes(x = (x + xend) / 2, y = y + 0.4, label = lab), size = 8, colour = "grey25") +
    geom_segment(data = sub_br, aes(x = x, xend = xend, y = y, yend = y), colour = "grey60", linewidth = 0.7) +
    geom_text(data = sub_br, aes(x = (x + xend) / 2, y = y + 0.3, label = lab), size = 6.5, colour = "grey35") +
    geom_path(data = arcs, aes(x, y, group = id, colour = q), linewidth = 1.1) +
    geom_point(data = ends, aes(x, y, colour = q), size = 2.6) +
    geom_label(data = arc_lab, aes(x, y, label = q, colour = q), size = 6.5, fontface = "bold",
               fill = "white", linewidth = 0, label.padding = unit(0.08, "lines")) +
    annotate("text", x = 0.45, y = c(1.45, 0.3), label = c("compared", "in"), hjust = 0, vjust = -0.4, size = 6, colour = "grey40") +
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
  fill_labs <- c(tau2 = "Context dependency: species x site (τ²)",
                 sigma2_experiment = "Context dependency: shared by site (σ² experiment)",
                 H_slope = "H species: slope piece (climate response)",
                 H_intercept = "H species: intercept piece (average)")
  pa <- ggplot() +
    geom_rect(data = bars, aes(xmin = xb - 0.15, xmax = xb + 0.15, ymin = ymin, ymax = ymax, fill = piece)) +
    geom_segment(data = refs, aes(x = xb + 0.15, xend = xend, y = value, yend = value, linetype = mark),
                 colour = ref_cols[refs$piece], linewidth = 1) +
    geom_text(data = totals, aes(xb, value, label = sprintf("%.2f", value)), vjust = -0.5, size = 5.2, fontface = "bold") +
    annotate("text", x = na_a, y = ymax * 0.25, label = rep("NA", length(na_a)), size = 8, colour = "grey55", fontface = "bold") +
    scale_fill_manual(values = fills, labels = fill_labs, name = NULL) +
    scale_linetype_manual(values = "dashed", name = NULL) + xs +
    scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
    labs(title = "a)", x = NULL, y = "Heterogeneity (×10⁻⁴)") +
    base + no_x + theme(legend.position = "inside", legend.position.inside = c(0.995, 0.99), legend.justification = c(1, 1),
                        legend.text = element_text(size = 17), legend.key.size = unit(1, "cm"), legend.spacing.y = unit(0, "cm"),
                        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA))

  # ---- b)-d) prediction metrics (step 06) ----
  metric_row <- function(key, tag, ylab, lims, last) {
    s <- if (is.null(metrics)) NULL else dplyr::filter(metrics, metric == key) %>% dplyr::mutate(x = match(column, ladder_columns))
    g <- if (is.null(gradients)) NULL else dplyr::filter(gradients, metric == key) %>% dplyr::mutate(x = match(column, ladder_columns))
    have <- if (is.null(s)) integer() else s$x[!is.na(s$pred) | !is.na(s$fit)]
    p <- ggplot() + geom_hline(yintercept = 0, colour = "grey60")
    if (length(setdiff(seq_len(n), have)))
      p <- p + annotate("text", x = setdiff(seq_len(n), have), y = mean(lims), label = "NA", size = 8, colour = "grey55", fontface = "bold")
    if (!is.null(g) && nrow(g)) p <- p + geom_jitter(data = g, aes(x, pred), width = 0.12, height = 0, colour = "#8DB8EA", size = 2.2, alpha = 0.7)
    if (!is.null(s) && nrow(s)) {
      ceil <- s$pred[s$column == "Species-dependent"]
      if (length(ceil) && !is.na(ceil)) p <- p + geom_hline(yintercept = ceil, linetype = "dashed", colour = "grey30")
      p <- p + geom_segment(data = s, aes(x = x, xend = x, y = fit, yend = pred), colour = "#2F80D8") +
        geom_point(data = s, aes(x, fit), shape = 21, fill = "white", colour = "#2F80D8", size = 5, stroke = 1.6) +
        geom_point(data = s, aes(x, pred), colour = "#2F80D8", size = 5)
    }
    p + xs + scale_y_continuous(limits = lims) + labs(title = tag, x = NULL, y = ylab) + base +
      if (last) theme(axis.text.x = element_text(size = 19, face = "bold", colour = "black")) else no_x
  }
  pb <- metric_row("R2", "b)", "R²", c(-0.2, 0.8), FALSE)
  pc <- metric_row("Spearman", "c)", "Spearman ρ", c(-0.2, 0.6), FALSE)
  pd <- metric_row("Skill", "d)", "Winner/loser skill", c(-0.2, 0.6), TRUE)

  fig <- patchwork::wrap_plots(p_top, pa, pb, pc, pd, ncol = 1, heights = c(1.1, 1.7, 1, 1, 1))
  if (!is.null(caption))
    fig <- fig + patchwork::plot_annotation(caption = caption, theme = theme(
      plot.caption = element_text(size = 15, colour = "grey40"), plot.background = element_rect(fill = "white", colour = NA)))
  ggsave(file, fig, width = 28, height = 24, dpi = 100, bg = "white")
  invisible(fig)
}
