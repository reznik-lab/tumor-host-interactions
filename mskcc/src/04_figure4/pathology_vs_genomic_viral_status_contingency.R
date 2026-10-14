# ============================================================================
# Produces: Figure 4 (C)
# ============================================================================
####################
# 2x2 contingency table: pathology-based vs. genomic (MSK-IMPACT) definition
# of MCC viral status.
#   Pathology: CM2B4 IHC status per path report (Positive = virus-positive,
#     Negative = virus-negative); NOT PERFORMED/AVAILABLE excluded.
#   Genomic (MSK-IMPACT definition): non-viral if (RB1 or TP53 mutant) AND
#     TMB >= 10, else viral -- same definition used for viral_status
#     elsewhere in the MCC codebase (summary_statistics.R).
#
# Styled like Figure 2E/F (fig2ef_contingency_tables.R): 2x2 tile grid, cell
# fill = % of the annotated cohort. No Fisher's exact test annotation -- the
# two methods agree on every case (0 discordant), so OR/p-value are
# degenerate (OR = Inf) and uninformative; the tile grid alone shows the
# perfect concordance.
#
# @author: Ethan Tse
####################

source("./src/00_shared/load_libraries.R")

mcc_df <- read_tsv("./data/private/msk/mcc/nmn_path_annotated_dephi_impact_mcc_table_annotated_08192026.tsv") %>%
  rename(cm2b4_status = `CM2B4 Status Per Path Report`) %>%
  filter(cm2b4_status %in% c("Positive", "Negative")) %>%
  # Drop only cases the genomic rule cannot resolve. TMB is consulted only when
  # RB1 or TP53 is altered, so a missing TMB is irrelevant when both are wild
  # type -- the rule already returns "viral" on the first clause alone. A blanket
  # !is.na(cvr_tmb_score) dropped one such case (RB1/TP53 wild type, TMB missing,
  # CM2B4 Positive) that the rest of the MCC codebase classifies as viral.
  filter(!is.na(rb1), !is.na(tp53)) %>%
  filter(!((rb1 == 1 | tp53 == 1) & is.na(cvr_tmb_score))) %>%
  mutate(
    # y-axis levels reversed vs. x (genomic_status) so the concordant
    # (matching) cells fall on the top-left -> bottom-right diagonal --
    # ggplot places a discrete axis's first level at the bottom, so "MCPyV+"
    # (genomic's first/leftmost level) needs to be pathology's *last* level
    # to land at the top.
    # Labelled by MCPyV status rather than "Viral"/"Non-viral": the axis names
    # the virus directly, and the signed short forms fit the tile labels.
    pathology_status = factor(ifelse(cm2b4_status == "Positive", "MCPyV+", "MCPyV\u2212"), levels = c("MCPyV\u2212", "MCPyV+")),
    genomic_status = factor(ifelse((rb1 == 1 | tp53 == 1) & cvr_tmb_score >= 10, "MCPyV\u2212", "MCPyV+"), levels = c("MCPyV+", "MCPyV\u2212"))
  )

out_dir <- "./code_outputs/results/figs/bmi/merkel_cell_carcinoma"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cell_counts <- mcc_df %>%
  count(pathology_status, genomic_status, .drop = FALSE) %>%
  mutate(pct = n / sum(n))

# No title, no Fisher's exact subtitle (degenerate OR = Inf, uninformative)
# -- but cell n/(%) annotation is kept, matching fig2ef_contingency_tables.R.
# This panel sits in a narrow column (Figure 5, panel E, alongside A/C/D),
# much tighter than fig2ef_contingency_tables.R's standalone 3in-wide
# panels, so the annotation font is shrunk to fit.
pathology_vs_genomic_contingency <- ggplot(cell_counts, aes(x = genomic_status, y = pathology_status)) +
  geom_tile(aes(fill = pct), color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("n = %d\n(%.1f%%)", n, pct * 100)), size = 1.8) +
  scale_fill_gradient(low = "white", high = "#2c7fb8", limits = c(0, 1), guide = "none") +
  labs(
    x = "Genomic status",
    y = "Pathology status"
  ) +
  theme_std(base_size = 6) +
  theme(
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.1),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(2, 2, 2, 2, "mm")
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed()

saveRDS(pathology_vs_genomic_contingency, file.path(out_dir, "pathology_vs_genomic_viral_status_contingency.rds"))

ggsave_wrapper(
  plot = pathology_vs_genomic_contingency,
  filename = file.path(out_dir, "pathology_vs_genomic_viral_status_contingency.pdf"),
  width = 3,
  height = 3.3,
  extension = "pdf"
)
