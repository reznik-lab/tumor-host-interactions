# ============================================================================
# Produces: Supplementary Figure 2; Supplementary Table 14
# (a figure listed without a panel letter contributes to a composite
#  panel whose letter could not be attributed unambiguously)
# ============================================================================
################
# Forest plots: one-vs-rest enrichment of myopenia, myosteatosis, and liver
# steatosis across cancer types, on the log-odds scale.
#
# Reads the adjusted logistic one-vs-rest contrasts from
# myopenia_one_vs_rest.R (cancer type vs. every other patient pooled,
# adjusted for age/sex/stage/ancestry/smoking/Yost index) for both the full
# cohort and the treatment-naive ("clean scan") sensitivity subset, and
# renders one row of panels (one per outcome), with the two cohorts
# dodged side by side per cancer type and distinguished by point color
# (colors match the "All Scans" / "Clean Scan Subset" comparison bars in
# abnormal_physiology.R). All three panels share one cancer-type ordering,
# set by the full-cohort myopenia odds ratio (see the order_key note below).
#
# Significance is marked directly on the plot as asterisks next to each
# point (not on the axis labels), following the same convention as Supp Fig
# 4A (pathway_enrichment_by_bodycomp_tx_naive_bmi.R): geom_text at the upper
# CI bound, offset with hjust, dodged alongside its point/line.
#
# @author: Ethan Tse
# @date: August 26th, 2026
################

library(tidytext)
source("./src/00_shared/load_libraries.R")

results_dir <- "./code_outputs/results/tables/body_comp/abnormal_physiology"
out_dir <- "./code_outputs/results/figs/body_comp/abnormal_physiology"
create_dir_if_needed(paste0(out_dir, "/"))

outcome_labels <- c(
  "myopenia" = "Myopenia",
  "myosteatosis" = "Myosteatosis",
  "liver_steatosis" = "Liver Steatosis"
)

# match the "All Scans" / "Clean Scan Subset" naming and colors used in
# abnormal_physiology.R's own all-vs-clean-scan comparison bar plot
cohort_colors <- c("All Scans" = "#4C72B0", "Clean Scan Subset" = "#DD8452")

full_cohort <- read_tsv(
  file.path(results_dir, "myopenia_myosteatosis_liver_steatosis_one_vs_rest.tsv"),
  show_col_types = FALSE
) %>%
  mutate(cohort = "All Scans")

clean_scan_subset <- read_tsv(
  file.path(results_dir, "myopenia_myosteatosis_liver_steatosis_one_vs_rest_clean_scan_subset.tsv"),
  show_col_types = FALSE
) %>%
  mutate(cohort = "Clean Scan Subset")

combined_results <- bind_rows(full_cohort, clean_scan_subset) %>%
  filter(!is.na(cancer_code)) %>%
  mutate(
    # Two cancer types carry a slash-separated alias in cancer_type_detailed
    # ("Uterine Serous Carcinoma/Uterine Papillary Serous Carcinoma",
    # "Uterine Carcinosarcoma/Uterine Malignant Mixed Mullerian Tumor").
    # Every other panel in Supplementary Figure 2 shows the first form only,
    # and at full length these two ran past the figure's left margin.
    cancer_type_detailed = sub("/.*$", "", cancer_type_detailed),
    outcome_label = factor(outcome_labels[outcome], levels = unname(outcome_labels)),
    cohort = factor(cohort, levels = names(cohort_colors)),
    # Significance stars, drawn on the plot rather than the axis: one star =
    # nominally significant, two = survives BH correction. Same convention as
    # the companion GLMM forest (Supplementary Figure 2f) and the pathway
    # forest in Supplementary Figure 5a.
    sig_label = case_when(
      padj < 0.05 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# One cancer-type ordering shared by all three outcome panels, set by the
# full-cohort ("All Scans") myopenia odds ratio.
#
# Previously each panel sorted by its own odds ratio via
# tidytext::reorder_within(), which meant each facet carried its own y-axis
# labels. That is affordable with short cancer codes and not with full
# cancer type names: three label columns at ~1.5in each exceeds the panel's
# whole width, and the plot areas collapse. A shared order means one label
# column serves all three facets. The cost is that cancer types are no
# longer ranked by each outcome's own effect size.
order_key <- combined_results %>%
  filter(cohort == "All Scans", outcome == "myopenia") %>%
  select(cancer_type_detailed, order_or = odds_ratio)

combined_results <- combined_results %>%
  left_join(order_key, by = "cancer_type_detailed") %>%
  # Cancer types not tested for myopenia still need a position; fall back to
  # their own estimate so they sort somewhere sensible rather than dropping.
  mutate(
    order_or = coalesce(order_or, odds_ratio),
    cancer_type_detailed = forcats::fct_reorder(cancer_type_detailed, order_or)
  )

# Stars sit in one column at the right of every facet rather than at each
# interval's own end. At the interval end they tracked conf.high, so they
# landed at a different x on every row and in the wider rows ran past the
# facet into its neighbour; a fixed column also lets a reader scan which
# cancer types are flagged without following each bar out. The three facets
# share one x scale, so a single x serves all of them; the axis expansion
# below reserves the space it needs.
#
# One column per cohort rather than one column dodged vertically: at 14
# cancer types over this panel's height each category is ~10pt, so the
# dodged halves land ~3pt apart while the asterisks themselves are ~4pt
# tall -- where both arms were flagged the two labels collided into an
# unreadable cluster. Side by side they never touch, and the star colour
# matches the arm's own points, so the cohort key already in the figure
# says which column is which.
#
# Carried as a column rather than a free variable: aes() captures the name,
# not the value, and figs.qmd re-renders this saved object in its own
# environment, where a loose `star_x` does not exist.
# Both columns are anchored to the RIGHT of their labels (hjust = 1 below)
# and placed relative to the panel's right edge, so the asterisks line up on
# a right margin rather than on a left one: "**" and "*" end together, and
# the pair of columns sits flush against the edge of every facet. 0.15 of
# the x range between the two anchors -- a "**" label is ~0.09 of the range
# wide at the size set below, so a narrower gap let the two arms' asterisks
# run together and read as one four-star label.
star_span <- diff(range(c(combined_results$conf.low, combined_results$conf.high),
                        na.rm = TRUE))
combined_results <- combined_results %>%
  mutate(star_x = max(conf.high, na.rm = TRUE) +
           star_span * (0.30 - (2 - as.integer(cohort)) * 0.15))

forest_combined <- combined_results %>%
  ggplot(aes(x = estimate, y = cancer_type_detailed)) +
  ggstats::geom_stripped_rows() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1) +
  geom_linerange(
    aes(xmin = conf.low, xmax = conf.high, color = cohort),
    position = position_dodge(width = 0.6),
    lineend = "round", linewidth = 0.4
  ) +
  geom_point(
    # Colour alone separates the two arms; a second (shape) channel encoding
    # the same variable only halved the apparent size of one arm's dots at
    # this figure's panel width.
    aes(color = cohort),
    position = position_dodge(width = 0.6),
    size = 1.6
  ) +
  geom_text(
    # Colour, not dodge, separates the arms here (see star_x above).
    aes(x = star_x, label = sig_label, colour = cohort),
    # hjust = 1: the labels are right-aligned on their column's anchor, so
    # a one- and a two-asterisk mark end at the same x.
    hjust = 1,
    # The asterisks are the only mark of significance in this panel and at
    # 1.8 they were barely visible against the dots; the x expansion below
    # buys the width the larger glyphs and their wider gap need.
    size = 2.5,
    show.legend = FALSE
  ) +
  scale_colour_manual(name = NULL, values = cohort_colors) +
  # Linear now: the values are already log odds ratios, so a log-transformed
  # axis would be logging them twice. Null is 0, not 1.
  # The star anchors are already inside the scale's range (they are mapped
  # x values), and the labels run leftwards from them, so the right side
  # needs only a hair of expansion -- which is what puts the right-hand
  # column flush against the panel edge.
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.01))) +
  # Shared y now, so no free_y and no scale_y_reordered(): only the leftmost
  # facet draws the cancer-type labels.
  facet_wrap(~outcome_label, nrow = 1) +
  labs(x = "log(OR) of the condition (cancer type vs. all others, adjusted)", y = NULL) +
  theme_std() +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.box.just = "center",
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 6)
  )

ggsave_wrapper(
  plot = forest_combined,
  filename = file.path(out_dir, "myopenia_myosteatosis_liver_steatosis_one_vs_rest_forest.pdf"),
  width = 9,
  height = 3.2
)
