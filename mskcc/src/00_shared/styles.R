# ============================================================================
# Shared helper: sourced by the analysis scripts, never run on its own.
# ============================================================================
#####################
# Author: Ethan Tse
# Date: Oct 9th, 2024
#
# Define project-wide styles
#
######################

# Colour palettes: https://r-charts.com/color-palettes/#discrete
library(paletteer)

# theme for publication figures
theme_std <- function(base_size = 7) {
  theme_cowplot(font_size = base_size) %+replace%
    theme(
      line = element_line(colour = "black", linewidth = 0.1),
      axis.line.x = element_line(colour = "black", linewidth = 0.1),
      axis.line.y = element_line(colour = "black", linewidth = 0.1),
      axis.ticks = element_line(linewidth = 0.1),
      axis.ticks.length = unit(1, "pt")
    ) +
    # format facets
    theme(
      strip.background = element_blank()
    )
}

#########################
# Forest plot directional arrows
#########################

# A pair of directional arrow annotations (+ labels) for forest plots, sitting
# above the top row and pointing away from the reference line (x = 0) to show
# which side of the estimate corresponds to which direction of effect.
# x_lim: the (symmetric) axis half-width already used for scale_x_continuous.
# y: numeric y position, in the discrete y-scale's underlying integer coords,
#    to anchor the arrows at (typically nlevels(y_var) + 0.75).
# labels: c(left, right) text; colors: c(left, right) arrow/text color.
# facet_data: NULL (default) draws via annotate(), which repeats the arrows on
#    every facet panel. For a faceted plot where the arrows should only sit
#    above one panel (e.g. the top-most facet of a stacked facet_wrap), pass a
#    one-row data frame with the facet variable(s) set to that panel's value
#    -- geom_segment()/geom_text() (unlike annotate()) only draw in panels
#    matching their own data.
# Requires the plot to also add theme(plot.margin = margin(t = <enough>)) and
# coord_cartesian(clip = "off") so the arrows aren't clipped above the panel.
# Default colors match Figure 5A's volcano-plot direction arrows
# (snvs_bmi_association_volcano.R) so all directional-arrow panels share one look.
forest_arrows <- function(x_lim, y, labels = c("Lower", "Higher"),
                           colors = c("#1f77b4", "#d62728"),
                           gap_frac = 0.10, len_frac = 0.40,
                           text_size = 2, text_y_offset = 0.55,
                           facet_data = NULL) {
  gap <- x_lim * gap_frac
  len <- x_lim * len_frac
  if (is.null(facet_data)) {
    list(
      annotate("segment",
        x = -gap, xend = -(gap + len), y = y, yend = y,
        arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
        color = colors[1], linewidth = 0.1
      ),
      annotate("segment",
        x = gap, xend = gap + len, y = y, yend = y,
        arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
        color = colors[2], linewidth = 0.1
      ),
      annotate("text",
        x = -gap, y = y + text_y_offset, label = labels[1],
        color = colors[1], size = text_size, hjust = 1
      ),
      annotate("text",
        x = gap, y = y + text_y_offset, label = labels[2],
        color = colors[2], size = text_size, hjust = 0
      )
    )
  } else {
    list(
      geom_segment(
        data = facet_data, inherit.aes = FALSE,
        aes(x = -gap, xend = -(gap + len), y = y, yend = y),
        arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
        color = colors[1], linewidth = 0.1
      ),
      geom_segment(
        data = facet_data, inherit.aes = FALSE,
        aes(x = gap, xend = gap + len, y = y, yend = y),
        arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
        color = colors[2], linewidth = 0.1
      ),
      geom_text(
        data = facet_data, inherit.aes = FALSE,
        aes(x = -gap, y = y + text_y_offset), label = labels[1],
        color = colors[1], size = text_size, hjust = 1
      ),
      geom_text(
        data = facet_data, inherit.aes = FALSE,
        aes(x = gap, y = y + text_y_offset), label = labels[2],
        color = colors[2], size = text_size, hjust = 0
      )
    )
  }
}

# Subtitle-slot alternative to forest_arrows(): the same colour-coded
# line + closed-arrowhead direction indicator (matching Figure 5A's
# volcano-plot arrows), but rendered as a custom ggplot2 theme element for
# plot.subtitle instead of annotate()d inside the panel. Use this when the
# plot's top margin can't be relied on to stay wide enough to hold
# forest_arrows()'s in-panel annotations (e.g. panels whose margin is later
# overridden when composed into a multi-panel manuscript figure), or when the
# arrows would otherwise crowd a cowplot panel letter drawn just above the
# panel -- ggplot always reserves the subtitle its own row in the plot's
# gtable regardless of plot.margin, so it can't be clipped or collide with
# neighbouring content the way an in-panel annotation can.
# Usage: labs(subtitle = "x") + theme(plot.subtitle = element_arrow_subtitle(...))
# -- the labs() string just needs to be non-empty so ggplot allocates a
# subtitle row at all; its text is ignored in favour of `labels` below.
# `zero_frac` is where x = 0 sits as a fraction of the panel width, so the
# Lower/Higher divider can be placed on zero instead of the panel midpoint.
# Defaults to 0.5 (the previous behaviour, correct only for a symmetric axis).
# Passing the true zero position lets a plot use a data-driven asymmetric x
# range without the indicator drifting away from the zero line.
element_arrow_subtitle <- function(labels = c("Lower", "Higher"),
                                    colors = c("#1f77b4", "#d62728"),
                                    text_size = 7, arrow_inset = 0.1,
                                    height = unit(0.34, "in"),
                                    zero_frac = 0.5) {
  structure(
    list(
      labels = labels, colors = colors, text_size = text_size,
      arrow_inset = arrow_inset, height = height,
      zero_frac = min(max(zero_frac, 0.05), 0.95)
    ),
    class = c("element_arrow_subtitle", "element_text", "element")
  )
}

element_grob.element_arrow_subtitle <- function(element, label = "", ...) {
  labels <- element$labels
  colors <- element$colors
  inset <- element$arrow_inset

  make_half <- function(lab, color, point_left) {
    txt <- grid::textGrob(lab, gp = grid::gpar(col = color, fontsize = element$text_size), y = unit(1, "npc"), vjust = 1)
    if (point_left) {
      x0 <- unit(1 - inset, "npc")
      x1 <- unit(inset / 2, "npc")
    } else {
      x0 <- unit(inset / 2, "npc")
      x1 <- unit(1 - inset, "npc")
    }
    seg <- grid::segmentsGrob(
      x0 = x0, x1 = x1, y0 = unit(0.12, "npc"), y1 = unit(0.12, "npc"),
      arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
      gp = grid::gpar(col = color, lwd = 0.6)
    )
    grid::gTree(children = grid::gList(txt, seg))
  }

  left <- make_half(labels[1], colors[1], point_left = TRUE)
  right <- make_half(labels[2], colors[2], point_left = FALSE)

  zf <- if (is.null(element$zero_frac)) 0.5 else element$zero_frac
  gt <- gtable::gtable(widths = unit(c(zf, 1 - zf), "null"), heights = element$height)
  gt <- gtable::gtable_add_grob(gt, left, t = 1, l = 1)
  gt <- gtable::gtable_add_grob(gt, right, t = 1, l = 2)
  gt
}

# ── Significance colouring, shared across every forest/volcano panel ──────────
# Three levels so a panel can distinguish a result that survives multiple-testing
# correction from one that is only nominally significant:
#   padj < 0.05  red    (#ad3f3f, the long-standing significance colour)
#   p < 0.05     orange (nominal only -- did NOT survive correction)
#   NS           black
# sig_level() keeps all three levels whether or not a panel draws them, and the
# scale helpers below pass drop = FALSE, so every significance legend in the
# manuscript reads padj / p / NS in that order -- a reader comparing two panels
# sees the same key in the same order even when one panel has no nominal-only
# or no non-significant result.
SIG_COLOURS <- c(
  "padj < 0.05" = "#ad3f3f",
  "p < 0.05" = "#e08214",
  "NS" = "black"
)

sig_level <- function(p_value, padj = NULL) {
  if (is.null(padj)) padj <- rep(NA_real_, length(p_value))
  out <- dplyr::case_when(
    !is.na(padj) & padj < 0.05 ~ "padj < 0.05",
    !is.na(p_value) & p_value < 0.05 ~ "p < 0.05",
    TRUE ~ "NS"
  )
  factor(out, levels = names(SIG_COLOURS))
}

# ggplot2 4.0 draws a legend key's glyph only for levels some layer actually
# has rows for: drop = FALSE above reserves the key's label but leaves its dot
# blank, so a panel with no nominal-only result printed "P < 0.05" with nothing
# beside it. Appending one invisible point per level (plotted at NA and dropped
# by na.rm) gives every key something to draw, so the padj / p / NS key is
# complete on every panel whether or not that panel draws all three.
# `size` should match the panel's own points, or the present levels get two
# superimposed glyphs of different sizes.
significance_key <- function(size = 2, shape = 16) {
  d <- data.frame(
    .sig_key_x = NA_real_,
    .sig_key_y = NA_real_,
    significance = factor(names(SIG_COLOURS), levels = names(SIG_COLOURS))
  )
  ggplot2::geom_point(
    data = d,
    mapping = ggplot2::aes(x = .sig_key_x, y = .sig_key_y, colour = significance),
    inherit.aes = FALSE, na.rm = TRUE, show.legend = TRUE,
    size = size, shape = shape
  )
}
