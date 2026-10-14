# ============================================================================
# Produces: Figure 5 (C); Supplementary Table 23
# ============================================================================
################
# PI3K/AKT/mTOR Pathway Body Composition — Forest Plot Visualizations
#
# Reads Firth logistic regression results from pi3k_bodycomp_analysis.R and
# renders forest plots for CRC/PDAC/RCC. The LUAD/PIK3CA forest plot lives in
# lung_pik3ca_bodycomp_analysis.R.
# Rows: model type (Marginal, L3 Joint). Columns: cancer type.
#
# @author: Ethan Tse
# @date: May 24th, 2026
################

source("./src/00_shared/load_libraries.R")

pi3k_crc_pdac_rcc <- read_tsv("./code_outputs/results/tables/body_comp/pi3k/pi3k_bodycomp_firth_crc_pdac_rcc.tsv")

all_results <- pi3k_crc_pdac_rcc %>%
  mutate(
    model_label = case_when(
      model_type == "marginal" ~ "Marginal",
      model_type == "joint" ~ "L3 Joint"
    ),
    model_label = factor(model_label, levels = c("Marginal", "L3 Joint")),
    metric_label = case_when(
      grepl("total_muscle_area", var) ~ "SKM",
      grepl("vat_area", var) ~ "VAT",
      grepl("sat_area", var) ~ "SAT",
      var == "bmi_at_scan" ~ "BMI",
      var == "liver_volume" ~ "Liver",
      var == "pancreas_volume" ~ "Pancreas",
      var == "spleen_volume" ~ "Spleen",
      var == "kidney_volume" ~ "Kidney"
    ),
    metric_label = factor(metric_label, levels = c("BMI", "SKM", "VAT", "SAT", "Liver", "Pancreas", "Spleen", "Kidney")),
    level = factor(level, levels = c("L1", "L3", "T10", "T12", "mean", "organ")),
    significance = sig_level(p.value, padj)
  )

forest_theme <- list(
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1),
  geom_linerange(
    aes(
      xmin = conf.low, xmax = conf.high,
      color = significance
    ),
    lineend = "round",
    linewidth = 0.1,
    position = position_dodge(width = 0.3)
  ),
  geom_point(size = 2, shape = 16, aes(colour = significance), position = position_dodge(width = 0.3)),
  scale_colour_manual(
    name = NULL,
    values = SIG_COLOURS,
    drop = FALSE
  ),
  significance_key(size = 2),
  labs(x = "Adj. log(OR) per SD", y = NULL, fill = NULL),
  theme_std(),
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 6),
    panel.spacing.x = unit(2, "pt")
  )
)

# ── Section 1: CRC / PDAC / RCC ───────────────────────────────────────────────

s1_cancers <- c("Colorectal Adenocarcinoma", "Pancreatic Adenocarcinoma", "Renal Clear Cell Carcinoma")

# Combine the L3 body-composition areas with the organ volumes into one panel per
# cancer, with the organ volumes stacked underneath the areas.
combined_crc_pdac_rcc <- all_results %>%
  filter(model_type == "marginal", level %in% c("L3", "organ")) %>%
  filter(cancer_type_detailed %in% s1_cancers) %>%
  mutate(
    metric_label = droplevels(metric_label),
    # reverse so areas sit on top and organ volumes underneath (y position 1 = bottom)
    metric_label = factor(metric_label, levels = rev(levels(metric_label))),
    # split body-composition tissue areas from organ volumes into separate facet rows
    compartment = factor(
      if_else(level == "organ", "Organ volume", "Tissue area"),
      levels = c("Tissue area", "Organ volume")
    )
  )

# label_wrap_gen(14) breaks "Renal Clear Cell Carcinoma" across three lines
# ("Renal Clear" / "Cell" / "Carcinoma"); a manual two-line break keeps "Cell
# Carcinoma" together on one line instead (matches fig_2b_labels, figs.qmd)
s1_cancer_labels <- c(
  "Colorectal Adenocarcinoma" = "Colorectal\nAdenocarcinoma",
  "Pancreatic Adenocarcinoma" = "Pancreatic\nAdenocarcinoma",
  "Renal Clear Cell Carcinoma" = "Renal Clear\nCell Carcinoma"
)

forest_crc_pdac_rcc <- combined_crc_pdac_rcc %>%
  ggplot(aes(x = estimate, y = metric_label)) +
  ggstats::geom_stripped_rows() +
  forest_theme +
  facet_grid(
    compartment ~ cancer_type_detailed,
    scales = "free", space = "free_y",
    labeller = labeller(cancer_type_detailed = as_labeller(s1_cancer_labels))
  )

# ── Save ──────────────────────────────────────────────────────────────────────

create_dir_if_needed("./code_outputs/results/figs/body_comp/pi3k")

ggsave_wrapper(
  plot = forest_crc_pdac_rcc,
  filename = "./code_outputs/results/figs/body_comp/pi3k/pi3k_forest_crc_pdac_rcc.pdf",
  width = 2.5,
  height = 2
)
