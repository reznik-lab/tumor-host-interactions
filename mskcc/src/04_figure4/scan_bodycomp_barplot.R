# ============================================================================
# Produces: Figure 4 (I)
# ============================================================================
#######################
# Figure 3G: body composition (L3 tissue areas) for two example MCC scans
#
# Plots the SAT / VAT / SKM areas measured from the two representative L3
# summary images shown in Figure 3G. The filename prefix of each image is the
# scan `accession`, which keys into the master body composition metadata.
#
#   fig_3g_left  <- "data/private/msk/mcc/mcc_scans/646668_2_L3SummaryImage.jpg"
#   fig_3g_right <- "data/private/msk/mcc/mcc_scans/646736_5_L3SummaryImage.jpg"
#
# @author: Ethan Tse
# @date: 2026-07-03
#######################

source("./src/00_shared/load_libraries.R")

body_comp <- read_tsv("./code_outputs/stable_data/msk/body_comp/metadata/master_body_comp_annotated_metadata.tsv")

# the two representative L3 summary images (filename prefix = accession)
fig_3g_scans <- tibble::tribble(
  ~accession, ~scan_label,
  646668,     "MCPyV\u2212",
  646736,     "MCPyV+"
)

fig_3g_bodycomp <- body_comp %>%
  filter(accession %in% fig_3g_scans$accession) %>%
  transmute(
    accession,
    BMI = bmi_at_scan,
    SAT = l3sat_area,
    VAT = l3vat_area,
    SKM = l3total_muscle_area
  ) %>%
  left_join(fig_3g_scans, by = "accession") %>%
  pivot_longer(c(BMI, SAT, VAT, SKM), names_to = "compartment", values_to = "area") %>%
  mutate(scan_label = factor(scan_label, levels = fig_3g_scans$scan_label))

# fold difference between the two scans per compartment (larger / smaller),
# ordered least -> highest fold change
fig_3g_fold <- fig_3g_bodycomp %>%
  group_by(compartment) %>%
  summarise(fold = max(area) / min(area), y = max(area), .groups = "drop") %>%
  arrange(fold)

# BMI pinned as the first bar; remaining compartments ordered by fold change
compartment_order <- c("BMI", setdiff(as.character(fig_3g_fold$compartment), "BMI"))
fig_3g_bodycomp <- fig_3g_bodycomp %>%
  mutate(compartment = factor(compartment, levels = compartment_order))
fig_3g_fold <- fig_3g_fold %>%
  mutate(
    compartment = factor(compartment, levels = compartment_order),
    label_y = y * 1.15 # per-compartment height for the fold label (free y axis)
  )

scan_cols <- c("MCPyV\u2212" = "#2c7fb8", "MCPyV+" = "#ad3f3f")

# TOP: raw values per compartment (free y axis, since units differ)
fig_3g_barplot <- fig_3g_bodycomp %>%
  ggplot(aes(x = scan_label, y = area, fill = scan_label)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = round(area, 1)), vjust = -0.4, size = 2.4) +
  # fold-difference annotation at the top of each panel, with highlighted background
  geom_label(
    data = fig_3g_fold, inherit.aes = FALSE,
    aes(x = 1.5, y = label_y, label = sprintf("%.1f×", fold)),
    fill = "#F1CE63", label.size = 0, size = 2.8,
    label.padding = unit(0.15, "lines")
  ) +
  facet_wrap(~compartment, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = "Measured value", fill = NULL) +
  scale_fill_manual(values = scan_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  theme_std() +
  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# BOTTOM: standardize each scan against the MCC cohort so compartments are
# comparable on a common scale (SD units / z-score)
qc_removed <- c("P-0036928", "P-0036676", "P-0039150", "P-0054731", "P-0057490", "P-0060039", "P-0058771", "P-0068753")
cohort_stats <- body_comp %>%
  filter(cancer_type_detailed == "Merkel Cell Carcinoma", !dmp_id %in% qc_removed) %>%
  transmute(BMI = bmi_at_scan, SAT = l3sat_area, VAT = l3vat_area, SKM = l3total_muscle_area) %>%
  pivot_longer(everything(), names_to = "compartment", values_to = "area") %>%
  group_by(compartment) %>%
  summarise(mean = mean(area, na.rm = TRUE), sd = sd(area, na.rm = TRUE), .groups = "drop")

fig_3g_z <- fig_3g_bodycomp %>%
  left_join(cohort_stats, by = "compartment") %>%
  mutate(
    z = (area - mean) / sd,
    # re-apply ordering: the join to a character key drops the factor levels
    compartment = factor(compartment, levels = compartment_order)
  )

# viral vs non-viral gap in SD units per compartment (SAT should dwarf BMI):
# this is the comparable "how big is the change" measure across compartments
fig_3g_delta <- fig_3g_z %>%
  select(compartment, scan_label, z) %>%
  pivot_wider(names_from = scan_label, values_from = z) %>%
  mutate(delta = .data[["MCPyV+"]] - .data[["MCPyV\u2212"]])

delta_label_y <- max(fig_3g_z$z) + 0.5

fig_3g_z_plot <- fig_3g_z %>%
  ggplot(aes(x = scan_label, y = z, fill = scan_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey40") +
  geom_col(width = 0.65) +
  # SD gap between the two scans, highlighted (SAT >> BMI)
  geom_label(
    data = fig_3g_delta, inherit.aes = FALSE,
    aes(x = 1.5, y = delta_label_y, label = sprintf("%.1f SD gap", delta)),
    fill = "#F1CE63", label.size = 0, size = 2.8,
    label.padding = unit(0.15, "lines")
  ) +
  facet_wrap(~compartment, nrow = 1) + # shared y: z is comparable across compartments
  labs(x = NULL, y = "SD from MCC\ncohort mean", fill = NULL) +
  scale_fill_manual(values = scan_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  theme_std() +
  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# DELTA: viral vs non-viral gap in SD units, one bar per compartment.
# Directly contrasts the size of the change across compartments (SAT >> BMI).
# Embedded in the manuscript as Figure 5K, a much smaller/squatter cell than
# this script's own 2.6x2.6in save size -- expansion bumped from 0.12 (the
# tallest bar's "1.3 SD" label was clipping against the panel edge) and
# margins/fonts trimmed down so the panel reads clean at that smaller size.
fig_3g_delta_plot <- fig_3g_delta %>%
  mutate(compartment = forcats::fct_reorder(compartment, delta)) %>%
  ggplot(aes(x = compartment, y = delta)) +
  geom_col(width = 0.65, fill = "#4c6b8a") +
  geom_text(aes(label = sprintf("%.1f SD", delta)), vjust = -0.4, size = 2.4) +
  labs(x = NULL, y = "MCPyV+ \u2212 MCPyV\u2212\n(SD units)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
  theme_std() +
  theme(
    plot.margin = margin(t = 9, r = 2, b = 2, l = 2),
    axis.title.y = element_text(size = 6, margin = margin(r = 2)),
    axis.text = element_text(size = 5.5)
  )

# save the figures separately


ggsave_wrapper(
  plot = fig_3g_delta_plot,
  filename = "./code_outputs/results/figs/body_comp/merkel_cell_carcinoma/fig_3g_example_scans_bodycomp_delta_sd.pdf",
  width = 2.6,
  height = 2.6
)
