# ============================================================================
# Produces: Figure 5 (A)
# ============================================================================
#######################
# Visualizing the logistic regression results of BMI-pathway associations with a volcano plot,
# highlighting the PI3K/AKT/mTOR pathway
#
#
#
# @author: Ethan Tse
# @date: 2026-04-29
#######################

source("./src/00_shared/load_libraries.R")

bmi_logreg_results_tidy <- read_tsv("./code_outputs/results/tables/bmi/impact/pathway_cancer_bmi_associations/snp_pathway_enrichment_by_bmi_logistic_regression_results.tsv")

# Direction-arrow geometry, same approach as Figure 5D's forest plot
# (src/04_figure4/melanoama_mcc_comparison.R) / Figure 4G's blood-labs
# volcano (pi3k_bloodlabs_volcano.R): positioned just above the highest data
# point via a wider plot margin + coord_cartesian(clip = "off"), rather than
# axis expansion (which just inflates empty white space at the top).
arrow_max_x <- max(abs(bmi_logreg_results_tidy$estimate), na.rm = TRUE) * 1.15
arrow_gap <- arrow_max_x * 0.10 # gap from centre before each arrow starts
arrow_len <- arrow_max_x * 0.40 # equal arrow length on both sides
arrow_y <- max(-log10(bmi_logreg_results_tidy$p.value), na.rm = TRUE) * 1.15
arrow_text_y <- arrow_y * 1.12

# plot a volcano plot of significant BMI-pathways associated with obesity associated cancers vs non-obesity associated cancers
bmi_pathway_volcano_plot <- bmi_logreg_results_tidy %>%
  mutate(
    is_significant = ifelse(padj < 0.05, "Significant (padj < 0.05)", "Not Significant"),
    is_pi3k_pathway = grepl("pi3k", pathway, ignore.case = TRUE),
    cancer_code = cancer_codes[as.character(cancer_type_detailed)],
    # prettify pathway name for labels; specifically format PI3K/AKT/mTOR, Hippo, and RTK MAPK
    pathway_pretty = case_when(
      is_pi3k_pathway ~ "PI3K/AKT/mTOR",
      grepl("hippo", pathway, ignore.case = TRUE) ~ "Hippo",
      grepl("rtk.*mapk|mapk.*rtk", pathway, ignore.case = TRUE) ~ "RTK/RAS/MAPK",
      TRUE ~ gsub("_", " ", pathway)
    ),
    # strwrap(width = 18) breaks "Renal Clear Cell Carcinoma" as "Renal Clear
    # Cell" / "Carcinoma", splitting "Cell" from "Carcinoma"; override with a
    # manual break that keeps "Cell Carcinoma" together on one line instead
    # (matches fig_2b_labels/s1_cancer_labels in figs.qmd /
    # pi3k_bodycomp_forest_plots.R)
    cancer_wrapped = vapply(
      as.character(cancer_type_detailed),
      function(x) {
        if (x == "Renal Clear Cell Carcinoma") {
          "Renal Clear\nCell Carcinoma"
        } else {
          paste(strwrap(x, width = 18), collapse = "\n")
        }
      },
      character(1)
    ),
    labels = ifelse(
      is_significant == "Significant (padj < 0.05)" & is_pi3k_pathway,
      paste0(cancer_wrapped, ":\n", pathway_pretty),
      NA_character_
    )
  ) %>%
  ggplot(aes(x = estimate, y = -log10(p.value))) +
  geom_point(
    aes(color = is_significant),
    size = 2
  ) +
  labs(
    x = "Adj. log(OR) per SD BMI",
    y = expression(-log[10](italic(p) * "-value")),
    color = "Significance"
  ) +
  theme_std() +
  theme(
    legend.position = "top", legend.justification = "center",
    plot.margin = margin(t = 14, r = 6, b = 6, l = 6)
  ) +
  # linewidth/colour matching Figure 5A's dashed line
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.1) +
  # Direction arrows above the points: negative log(OR) = lower BMI, positive
  # = higher BMI. Closed arrowheads, matching Figure 5D / Figure 4G.
  annotate("segment",
    x = -arrow_gap, xend = -(arrow_gap + arrow_len), y = arrow_y, yend = arrow_y,
    arrow = grid::arrow(length = unit(0.12, "cm"), type = "closed"), linewidth = 0.1, colour = "#1f77b4"
  ) +
  annotate("text", x = -arrow_gap, y = arrow_text_y, label = "Lower BMI", size = 2, colour = "#1f77b4", hjust = 1) +
  annotate("segment",
    x = arrow_gap, xend = arrow_gap + arrow_len, y = arrow_y, yend = arrow_y,
    arrow = grid::arrow(length = unit(0.12, "cm"), type = "closed"), linewidth = 0.1, colour = "#d62728"
  ) +
  annotate("text", x = arrow_gap, y = arrow_text_y, label = "Higher BMI", size = 2, colour = "#d62728", hjust = 0) +
  geom_text_repel(
    aes(label = labels),
    box.padding = 1,
    point.padding = 0.3,
    max.overlaps = Inf,
    na.rm = TRUE,
    show.legend = FALSE,
    segment.color = "black",
    force = 5,
    min.segment.length = 0,
    size = 2,
    segment.size = 0.1,
    # ggrepel doesn't know about the annotate() arrows above; keep repelled
    # labels below the arrow row so they can't collide with it.
    ylim = c(NA, arrow_y * 0.92)
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Not Significant" = "#4D4D4D", "Significant (padj < 0.05)" = "#D85A30"),
    breaks = c("Significant (padj < 0.05)", "Not Significant"),
    labels = c("Significant (padj < 0.05)" = "padj < 0.05", "Not Significant" = "NS")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(clip = "off")

ggsave_wrapper(
  plot = bmi_pathway_volcano_plot,
  filename = "./code_outputs/results/figs/bmi/bmi_pathways_analysis/bmi_pathway_volcano_plot_pi3k_highlighted.pdf",
  width = 3,
  height = 3,
  extension = "pdf"
)
