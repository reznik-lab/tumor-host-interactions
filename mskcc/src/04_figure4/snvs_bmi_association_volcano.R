# ============================================================================
# Produces: Figure 4; Supplementary Table 18
# (a figure listed without a panel letter contributes to a composite
#  panel whose letter could not be attributed unambiguously)
# ============================================================================
###############
# Visualize the results of gene-bmi association using logistic regression
#
# We visualize using a volcano plot, with log(OR) on the x-axis and -log10(p-value) on the y-axis.
# We highlight significant associations and label the most significant ones.
#
# @author: Ethan Tse
# @date: 2026-02-17
###############

source("./src/00_shared/load_libraries.R")
library(purrr)
library(stringr)

gene_bmi_association_results <- read_tsv("./code_outputs/results/tables/bmi/impact/gene_cancer_bmi_associations/somatic_only_ancestry_adjusted_solid/somatic_only_ancestry_adjusted_solid_mv_model_stage_results.tsv")

# PI3K/AKT/MTOR pathway genes (see pathway_constants.R) to always label.
pi3k_akt_mtor_genes <- toupper(c(
  "pik3ca", "pik3cb", "pik3cd", "pik3cg", "pik3r1", "pik3r2", "pik3r3",
  "pten", "akt1", "akt2", "akt3", "mtor", "rptor", "rictor", "tsc1",
  "tsc2", "rheb", "pdpk1", "inpp4b", "inppl1", "ppp2r1a", "stk11"
))

# Build a single-line plotmath label of the form "Cancer Type: GENE",
# with the gene name italicised.
make_label <- function(cancer_name, gene) {
  sprintf("paste('%s: ', italic('%s'))", cancer_name, gene)
}

# Direction-arrow geometry, same approach as Figure 5D's forest plot
# (src/04_figure4/melanoama_mcc_comparison.R) / Figure 4A/4G's volcanoes
# (pathways_bmi_association_volcano.R, pi3k_bloodlabs_volcano.R): positioned
# just above the highest data point via a wider plot margin +
# coord_cartesian(clip = "off"), rather than axis expansion (which just
# inflates empty white space at the top).
full_model_data <- gene_bmi_association_results %>%
  filter(model == "fully_adjusted") %>%
  # Trim extreme outliers so the x-range stays readable. The cut is in the
  # units of the estimate, so it had to move with the switch from per-unit
  # BMI to per-SD BMI (estimates grew by each cancer type's BMI SD, ~5-7x):
  # left at 0.3 it discarded 15 of the 20 padj < 0.05 hits, Merkel Cell
  # Carcinoma RB1 among them. 1.5 drops the same handful of genuine outliers
  # (4 of 2471 rows, |estimate| up to 4.7) and keeps every significant hit.
  filter(abs(estimate) < 1.5)

arrow_max_x <- max(abs(full_model_data$estimate), na.rm = TRUE) * 1.15
arrow_gap <- arrow_max_x * 0.10 # gap from centre before each arrow starts
arrow_len <- arrow_max_x * 0.40 # equal arrow length on both sides
arrow_y <- max(-log10(full_model_data$pvalue), na.rm = TRUE) * 1.15
arrow_text_y <- arrow_y * 1.12

full_model_volcano <- full_model_data %>%
  mutate(cancer_code = cancer_codes[as.character(cancer)]) %>%
  mutate(labels = map2_chr(cancer, gene, make_label)) %>%
  mutate(significance = ifelse(padj < 0.05, "padj < 0.05", "NS")) %>%
  # Baked into a column (rather than referenced as a free variable inside the
  # geom_text_repel data function below) so the saved plot .rds still labels
  # correctly when reloaded in a fresh R session that never sourced this script.
  mutate(is_labelled_gene = (cancer_code == "MCC" & gene == "RB1") | gene %in% pi3k_akt_mtor_genes) %>%
  # Draw NS points first so the significant ones land on top. The rows arrive
  # sorted by padj, and geom_point() paints in row order, so the ~2450 NS
  # points were being drawn over the 20 hits. Pancreatic Adenocarcinoma PIK3CA
  # (padj 0.048) was all but erased by Uterine Endometrioid Carcinoma ESR1
  # (padj 0.052), which sits 0.06 units away -- well inside one point radius --
  # and is the very next row.
  arrange(significance == "padj < 0.05") %>%
  ggplot(aes(x = estimate, y = -log10(pvalue))) +
  geom_point(
    aes(colour = significance),
    size = 2,
  ) +
  labs(
    x = "Adj. log(OR) per SD BMI",
    y = expression(-log[10](italic(p) * "-value")),
    color = "Significance"
  ) +
  theme_std() +
  theme(plot.margin = margin(t = 14, r = 6, b = 6, l = 6)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.1) +
  # Direction arrows above the points: negative log(OR) = enriched at lower
  # BMI (blue), positive = enriched at higher BMI (red). Closed arrowheads,
  # matching Figure 5D / Figure 4A / Figure 4G.
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
    # Only label MCC RB1 and PI3K/AKT/mTOR pathway gene associations; every
    # other significant point is left unlabelled to keep the plot readable.
    data = function(df) {
      df %>% filter(padj < 0.05 & is_labelled_gene)
    },
    aes(label = labels),
    parse = TRUE,
    box.padding = 0.4,
    point.padding = 0.2,
    max.overlaps = Inf,
    show.legend = FALSE,
    segment.color = "black",
    min.segment.length = 0,
    force = 5,
    size = 2,
    segment.size = 0.1,
    # ggrepel doesn't know about the annotate() arrows above; keep repelled
    # labels below the arrow row so they can't collide with it.
    ylim = c(NA, arrow_y * 0.92)
  ) +
  # orange / dark grey, matching Figure 4A/4G's case/control palette (Figure
  # 1C CASE_COLORS: Control = "#4D4D4D", Cancer = "#D85A30")
  scale_color_manual(
    name = NULL,
    values = c("NS" = "#4D4D4D", "padj < 0.05" = "#D85A30"),
    breaks = c("padj < 0.05", "NS"),
  ) +
  scale_x_continuous(n.breaks = 6) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "top",
    legend.justification = "center"
  )

ggsave_wrapper(
  plot = full_model_volcano,
  filename = "./code_outputs/results/figs/bmi/bmi_snvs_analysis/gene_bmi_association_full_model_volcano.pdf",
  width = 4,
  height = 3
)
