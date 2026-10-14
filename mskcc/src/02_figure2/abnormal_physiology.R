# ============================================================================
# Produces: Figure 2 (E); Supplementary Table 13
# ============================================================================
####################
# Quantifying the prevalence of extreme physiologies across cancer types
#
# In this, we quantify sarcopenia, myopenia obesity, and myosteatosis across cancer types
# Quantifying
#
# reference: https://pubs.rsna.org/doi/full/10.1148/radiol.222008#sec-2
#
# @author: Ethan Tse
# @date: February 18th, 2026
####################


source("./src/00_shared/load_libraries.R")

body_comp_data <- read_tsv("./code_outputs/stable_data/msk/body_comp/metadata/non_clean_master_body_comp_annotated_metadata.tsv") %>%
  filter(cancer_type_detailed != "Merkel Cell Carcinoma") %>%
  # recalibrate
  mutate(
    l3total_muscle_area_precontrast = 1.0628 * l3total_muscle_area,
    l3total_muscle_median_hu_precontrast = 0.9930 * l3total_muscle_median_hu - 7.5291,
    # the <80 condition is on the pre (converted) value, so branch on the output
    # of the <80 formula and fall back to the >=80 formula when it lands >=80
    liver_median_hu_precontrast = case_when(
      (0.7434 * liver_median_hu - 12.387) < 80 ~ (0.7434 * liver_median_hu) - 12.387,
      TRUE ~ (0.21254 * liver_median_hu) + 31.508,
    ),
  ) %>%
  mutate(cancer_code = cancer_codes[cancer_type_detailed]) %>%
  mutate(SMI = l3total_muscle_area_precontrast / (height_at_scan_cm / 100)^2) %>%
  mutate(myopenia = case_when(
    # Height is missing for ~8% of scans, so SMI is NA for them. Without
    # this guard the final `TRUE ~ 0` silently calls every one of them
    # "not myopenic", inflating the denominator and biasing prevalence down.
    is.na(SMI) ~ NA_real_,
    (gender == "Male") & (SMI < 50) ~ 1,
    (gender == "Female") & (SMI < 39) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(myopenia_obesity = case_when(
    (myopenia == 1) & (bmi_at_scan >= 30) ~ 1,
    # Either component being definitively absent rules the joint phenotype
    # out even when the other is unknown; only genuinely undetermined rows
    # fall through to NA.
    (myopenia == 0) | (bmi_at_scan < 30) ~ 0,
    TRUE ~ NA_real_
  )) %>%
  # myosteatosis is defined relative to the muscle radiodensity distribution,
  # computed within each cancer type (and sex) rather than across the whole cohort
  mutate(myosteatosis = case_when(
    # Same NA guard as myopenia: an unmeasured muscle HU is unknown
    # myosteatosis status, not absence of myosteatosis.
    is.na(l3total_muscle_median_hu_precontrast) ~ NA_real_,
    (gender == "Male") & l3total_muscle_median_hu_precontrast < quantile(l3total_muscle_median_hu_precontrast[gender == "Male"], 0.25, na.rm = TRUE) ~ 1,
    (gender == "Female") & l3total_muscle_median_hu_precontrast < quantile(l3total_muscle_median_hu_precontrast[gender == "Female"], 0.25, na.rm = TRUE) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(liver_steatosis = ifelse(liver_median_hu_precontrast < 57, 1, 0)) %>%
  mutate(obesity = ifelse(bmi_at_scan >= 30, 1, 0))

suff_cancer_types <- body_comp_data %>%
  count(cancer_type_detailed) %>%
  filter(n >= 50) %>%
  pull(cancer_type_detailed)

# plot the prevalence of these extreme physiologies across cancer types
abnormal_physiology_prevalence <- body_comp_data %>%
  mutate(cancer_code = cancer_codes[cancer_type_detailed]) %>%
  filter(cancer_type_detailed %in% suff_cancer_types) %>%
  group_by(cancer_code) %>%
  summarise(
    n = n(),
    myopenia_prevalence = mean(myopenia, na.rm = TRUE),
    obesity_prevalence = mean(obesity, na.rm = TRUE),
    myopenia_obesity_prevalence = mean(myopenia_obesity, na.rm = TRUE),
    myosteatosis_prevalence = mean(myosteatosis, na.rm = TRUE),
    liver_steatosis_prevalence = mean(liver_steatosis, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = -c(cancer_code, n), names_to = "abnormal_physiology", values_to = "prevalence") %>%
  # Wald binomial standard error / 95% CI for each prevalence
  mutate(
    se = sqrt(prevalence * (1 - prevalence) / n),
    ci_lower = pmax(0, prevalence - 1.96 * se),
    ci_upper = pmin(1, prevalence + 1.96 * se)
  ) %>%
  mutate(abnormal_physiology = case_when(
    abnormal_physiology == "myopenia_prevalence" ~ '"Myopenia"',
    abnormal_physiology == "obesity_prevalence" ~ '"Obesity"',
    abnormal_physiology == "myopenia_obesity_prevalence" ~ '"Myopenia/Obesity"',
    abnormal_physiology == "myosteatosis_prevalence" ~ '"Myosteatosis"',
    abnormal_physiology == "liver_steatosis_prevalence" ~ '"Liver Steatosis"',
    TRUE ~ abnormal_physiology
  )) %>%
  mutate(abnormal_physiology = factor(abnormal_physiology, levels = c(
    '"Obesity"',
    '"Myopenia"',
    '"Myopenia/Obesity"',
    '"Myosteatosis"',
    '"Liver Steatosis"'
  ))) %>%
  mutate(cancer_code = gsub(" ", "\n", cancer_code))

# Order cancer_code factor levels by mean BMI value
cancer_bmi_order <- body_comp_data %>%
  filter(cancer_type_detailed %in% suff_cancer_types) %>%
  mutate(cancer_code = cancer_codes[cancer_type_detailed]) %>%
  filter(!is.na(cancer_code)) %>%
  group_by(cancer_code) %>%
  summarize(mean_bmi = mean(bmi_at_scan, na.rm = TRUE)) %>%
  arrange(mean_bmi) %>%
  mutate(cancer_code = gsub(" ", "\n", cancer_code)) %>%
  pull(cancer_code)

# Order cancer_code factor levels by obesity proportion (most obese last, so it
# sits at the top after coord_flip())
obesity_order <- body_comp_data %>%
  mutate(cancer_code = cancer_codes[cancer_type_detailed]) %>%
  filter(cancer_type_detailed %in% suff_cancer_types, !is.na(cancer_code)) %>%
  group_by(cancer_code) %>%
  summarise(obesity_prevalence = mean(obesity, na.rm = TRUE)) %>%
  arrange(obesity_prevalence) %>%
  mutate(cancer_code = gsub(" ", "\n", cancer_code)) %>%
  pull(cancer_code)

# highlight cachexia-prone cancer types in red, the rest in black
cachexia_prone <- c("PDAC", "LUAD", "STAD", "ESCA", "COADREAD")
cancer_category <- setNames(
  ifelse(cancer_bmi_order %in% cachexia_prone, "Cachexia-prone", "Other"),
  cancer_bmi_order
)

category_colors <- c(
  "Cachexia-prone" = "#d62728",
  "Other" = "black"
)

cancer_label_colors <- setNames(
  unname(category_colors[cancer_category[cancer_bmi_order]]),
  cancer_bmi_order
)

# bold the axis labels of cachexia-prone cancer types
cancer_label_faces <- setNames(
  ifelse(cancer_bmi_order %in% cachexia_prone, "bold", "plain"),
  cancer_bmi_order
)

# map cancer codes back to their full names for axis labels;
# for names with a slash, keep only the part before the first slash
cancer_code_to_name <- setNames(gsub("/.*", "", names(cancer_codes)), cancer_codes)

# Precompute axis labels as a static named vector (break -> label) so the saved
# plot object is self-contained; a labels = function(...) closure would lazily
# reference cancer_code_to_name and break when the .rds is reloaded elsewhere.
cancer_x_labels <- setNames(
  stringr::str_wrap(cancer_code_to_name[cancer_bmi_order], width = 100),
  cancer_bmi_order
)

# plot
abnormal_physiology_plot <- abnormal_physiology_prevalence %>%
  filter(abnormal_physiology != '"Myopenia/Obesity"') %>%
  filter(!is.na(cancer_code)) %>%
  mutate(cancer_code = factor(cancer_code, levels = obesity_order)) %>%
  mutate(cancer_category = factor(
    cancer_category[as.character(cancer_code)],
    levels = c("Cachexia-prone", "Other")
  )) %>%
  ggplot(aes(x = cancer_code, y = prevalence, fill = abnormal_physiology)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_linerange(aes(ymin = ci_lower, ymax = ci_upper), linewidth = 0.1) +
  geom_point(aes(color = cancer_category), size = 0, alpha = 0) +
  scale_color_manual(
    values = category_colors,
    guide = guide_legend(override.aes = list(size = 2, alpha = 1, shape = 15))
  ) +
  facet_wrap(abnormal_physiology ~ ., nrow = 1, scales = "free_x", labeller = label_parsed) +
  labs(
    x = NULL,
    y = "Prevalence",
    fill = "Extreme Physiology",
    color = NULL
  ) +
  theme_std() +
  theme(
    # cancer types are on the y-axis after coord_flip(), so colour and bold the
    # cachexia-prone types red
    axis.text.y = element_text(
      color = cancer_label_colors[obesity_order],
      face = cancer_label_faces[obesity_order]
    ),
    panel.spacing = unit(0.75, "lines"),
    legend.position = "bottom",
    legend.justification = "center"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = c(0, 0),
    n.breaks = 3
  ) +
  scale_x_discrete(labels = cancer_x_labels) +
  ggthemes::scale_fill_tableau(palette = "Tableau 20") +
  guides(fill = "none") +
  coord_flip()

ggsave_wrapper(
  plot = abnormal_physiology_plot,
  filename = "./code_outputs/results/figs/body_comp/abnormal_physiology/abnormal_physiology_prevalence_by_cancer_type.pdf",
  width = 6,
  height = 3,
  extension = "pdf"
)

# Scatterplot of myopenia vs myosteatosis prevalence, one point per cancer type.
# figs.qmd loads this object into `fig_2h` but does not draw it into Figure 2 --
# the panel was dropped from the layout and the readRDS() left behind. It is
# still written here because that readRDS() would otherwise fail the chunk.
myopenia_myosteatosis_prevalence <- body_comp_data %>%
  mutate(cancer_code = cancer_codes[cancer_type_detailed]) %>%
  filter(cancer_type_detailed %in% suff_cancer_types) %>%
  filter(!is.na(cancer_code)) %>%
  group_by(cancer_code, cancer_type_detailed) %>%
  summarise(
    n = n(),
    myopenia_prevalence = mean(myopenia, na.rm = TRUE),
    myosteatosis_prevalence = mean(myosteatosis, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  # for names with a slash, keep only the part after the last slash
  mutate(cancer_type_detailed = gsub("/.*", "", cancer_type_detailed)) %>%
  mutate(cachexia_label = ifelse(
    cancer_code %in% c("PDAC", "LUAD", "STAD", "ESCA", "COADREAD"),
    "Cachexia-prone", "Other"
  ))

myopenia_myosteatosis_scatter <- myopenia_myosteatosis_prevalence %>%
  ggplot(aes(x = myopenia_prevalence, y = myosteatosis_prevalence)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkgrey") +
  geom_point(aes(color = cachexia_label), size = 2, alpha = 0.7) +
  ggrepel::geom_text_repel(
    # Pass all points (with blank labels for non-cachexia cancers) so ggrepel
    # repels labels away from every dot, not just the labelled subset.
    aes(
      label = ifelse(cachexia_label == "Cachexia-prone", cancer_type_detailed, ""),
      color = cachexia_label
    ),
    force = 50,
    size = 2,
    box.padding = 0.6,
    point.padding = 0.4,
    min.segment.length = 0,
    max.overlaps = Inf,
    segment.size = 0.2,
    segment.color = "grey50",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("Cachexia-prone" = "#d62728", "Other" = "black"),
    guide = guide_legend(override.aes = list(shape = 15, size = 2))
  ) +
  labs(
    x = "Myopenia Prevalence",
    y = "Myosteatosis Prevalence",
    color = NULL
  ) +
  theme_std() +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 10) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 10) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "left"
  )

ggsave_wrapper(
  plot = myopenia_myosteatosis_scatter,
  filename = "./code_outputs/results/figs/body_comp/abnormal_physiology/myopenia_vs_myosteatosis_prevalence_scatter.pdf",
  width = 4,
  height = 4,
  extension = "pdf"
)
