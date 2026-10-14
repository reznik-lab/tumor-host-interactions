# ============================================================================
# Produces: Figure 4; Supplementary Figure 4 (A)
# (a figure listed without a panel letter contributes to a composite
#  panel whose letter could not be attributed unambiguously)
# ============================================================================
#######################
# Body composition association with tumour mutation burden.
#
# In this file, test if high TMB is associated with differences in body composition
#
#
#
# @author: Ethan Tse
# @date: Jan 13, 2025
#######################

source("./src/00_shared/load_libraries.R")

# Here, we use cutaneous melanoma as a control for whether
# does cutaneous melanoma show the same pattern?
# Load the data
impact_muts <- read_tsv("./code_outputs/stable_data/msk/impact/snvs/dataframes/bmi_oncokb_annotated_first_sample_oncogenic_mutations.tsv")

melanoma_mcc_mutations <- impact_muts %>%
  filter(CANCER_TYPE_DETAILED %in% c("Cutaneous Melanoma", "Merkel Cell Carcinoma"))

# before bodycomp, associate BMI with TMB; pan skin cancers
pan_skin_cancers <- impact_muts %>% filter(grepl("Skin|Melanoma", CANCER_TYPE, ignore.case = TRUE))

bmi_assoc_df <- pan_skin_cancers %>%
  select(
    PATIENT_ID, CANCER_TYPE_DETAILED, CVR_TMB_SCORE,
    BMI_NEAREST_RESULT, CURRENT_AGE_DEID, GENDER,
    ANCESTRY_LABEL, STAGE_CDM_DERIVED
  ) %>%
  distinct() %>%
  group_by(CANCER_TYPE_DETAILED) %>%
  filter(n() >= 100) %>%
  group_modify(~ {
    n_original <- nrow(.x)
    dat <- .x %>% filter(!is.na(CVR_TMB_SCORE), !is.na(BMI_NEAREST_RESULT))
    n_effective <- nrow(dat)
    n_events <- sum(dat$CVR_TMB_SCORE >= 10, na.rm = TRUE)
    # BMI entered per standard deviation, standardized within each cancer
    # type, matching the pan-cancer volcano this panel sits beside. Rescaling a
    # single linear term leaves P-values and significance untouched -- only the
    # estimate and its interval move onto the new scale.
    fit <- glm((CVR_TMB_SCORE >= 10) ~ scale(BMI_NEAREST_RESULT) + CURRENT_AGE_DEID + GENDER + ANCESTRY_LABEL + STAGE_CDM_DERIVED,
      data = dat, family = binomial()
    )
    broom::tidy(fit) %>%
      mutate(
        conf.low = estimate - 1.96 * std.error, conf.high = estimate + 1.96 * std.error,
        n_original = n_original, n_effective = n_effective, n_events = n_events
      )
  }) %>%
  ungroup() %>%
  filter(term == "scale(BMI_NEAREST_RESULT)")

# acral (non uv related) and uveal melanoma (low TMB) cannot be interpreted
skin_bmi_assoc <- bmi_assoc_df %>%
  filter(!CANCER_TYPE_DETAILED %in% c("Acral Melanoma", "Uveal Melanoma"))

skin_bmi_plot_df <- skin_bmi_assoc %>%
  mutate(
    CANCER_TYPE_DETAILED = forcats::fct_reorder(CANCER_TYPE_DETAILED, estimate),
    # BH across the skin cancer types tested in this panel
    padj = p.adjust(p.value, method = "BH"),
    significant = sig_level(p.value, padj)
  )

# directional indicator (coloured, matching Figure 5A's volcano-plot arrows):
# left = lower TMB, right = higher TMB. Rendered in the plot's subtitle slot
# (see element_arrow_subtitle() in styles.R) rather than annotate()d above
# the panel, so it isn't clipped when this plot's top margin gets trimmed
# down for the composed Figure 5 layout.
low_tmb_col <- "#1f77b4"
high_tmb_col <- "#d62728"

# symmetric x-axis around 0 so the estimates sit evenly on both sides
x_lim <- max(abs(c(skin_bmi_plot_df$conf.low, skin_bmi_plot_df$conf.high))) * 1.15

bmi_tmb_forest <- skin_bmi_plot_df %>%
  ggplot(aes(x = estimate, y = CANCER_TYPE_DETAILED, color = significant)) +
  ggstats::geom_stripped_rows(aes(y = CANCER_TYPE_DETAILED), inherit.aes = FALSE) +
  # Dot size and linewidths standardized to the PI3K forest in the next figure
  # (pi3k_forest_crc_pdac_rcc.rds): vline and interval at linewidth 0.1, point
  # at size 2 -- so every forest in the manuscript draws at one weight.
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", linewidth = 0.1) +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high), linewidth = 0.1, lineend = "round") +
  geom_point(size = 2, shape = 16) +
  labs(
    x = "Adj. log(OR) per SD BMI", y = NULL, color = NULL,
    subtitle = "direction"
  ) +
  scale_color_manual(values = SIG_COLOURS, drop = FALSE) +
  significance_key(size = 2) +
  scale_x_continuous(limits = c(-x_lim, x_lim)) +
  theme_std() +
  theme(
    legend.position = "bottom",
    plot.subtitle = element_arrow_subtitle(labels = c("Lower TMB", "Higher TMB"), colors = c(low_tmb_col, high_tmb_col))
  ) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 12))

ggsave_wrapper(
  "code_outputs/figs/impact/body_comp/merkel_cell_carcinoma/bmi_tmb_forest_skin_cancers.pdf",
  bmi_tmb_forest,
  width = 3, height = 2.9
)

# SUPPLEMENTARY: number of patients compared per cancer, split by UV/TMB status
# (n_events = TMB-high >= 10, the UV-high proxy; the rest are TMB-low)
# sorted by total N (largest at top)
tmb_counts <- skin_bmi_assoc %>%
  transmute(
    CANCER_TYPE_DETAILED = forcats::fct_reorder(CANCER_TYPE_DETAILED, n_effective),
    n_effective,
    `TMB-high (≥10)` = n_events,
    `TMB-low (<10)` = n_effective - n_events
  ) %>%
  pivot_longer(
    c(`TMB-high (≥10)`, `TMB-low (<10)`),
    names_to = "tmb_group", values_to = "n"
  ) %>%
  mutate(tmb_group = factor(tmb_group, levels = c("TMB-low (<10)", "TMB-high (≥10)")))

tmb_count_bar <- tmb_counts %>%
  ggplot(aes(x = n, y = CANCER_TYPE_DETAILED, fill = tmb_group)) +
  geom_col(width = 0.7) +
  # per-segment counts (only label non-zero segments)
  geom_text(
    data = ~ dplyr::filter(.x, n > 0),
    aes(label = n),
    position = position_stack(vjust = 0.5), size = 2.2, color = "white"
  ) +
  # total N at the end of each bar
  geom_text(
    data = dplyr::distinct(tmb_counts, CANCER_TYPE_DETAILED, n_effective),
    aes(x = n_effective, y = CANCER_TYPE_DETAILED, label = n_effective),
    inherit.aes = FALSE, hjust = -0.2, size = 2.3
  ) +
  labs(x = "Number of patients", y = NULL, fill = NULL) +
  scale_fill_manual(values = c("TMB-low (<10)" = low_tmb_col, "TMB-high (≥10)" = high_tmb_col)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_std() +
  theme(legend.position = "bottom") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 20))

ggsave_wrapper(
  "./code_outputs/figs/impact/body_comp/merkel_cell_carcinoma/bmi_tmb_counts_skin_cancers_supp.pdf",
  tmb_count_bar,
  width = 3.5, height = 2.6
)


body_comp <- read_tsv("./code_outputs/stable_data/msk/body_comp/metadata/non_clean_master_body_comp_annotated_metadata.tsv")

# patients failing scan QC, excluded from all body composition analyses
qc_removed <- c(
  "P-0036928", "P-0036676", "P-0039150", "P-0054731",
  "P-0057490", "P-0060039", "P-0058771", "P-0068753"
)

bodycomp_melanoma_intersection <- melanoma_mcc_mutations %>%
  inner_join(body_comp, by = c("PATIENT_ID" = "dmp_id")) %>%
  filter(!PATIENT_ID %in% qc_removed) %>%
  select(!contains(".y")) %>%
  # if contains .x remove the .x
  rename_with(~ stringr::str_remove(.x, "\\.x"), contains(".x")) %>%
  mutate(weight_category = case_when(
    BMI_NEAREST_RESULT < 18.5 ~ "Underweight",
    BMI_NEAREST_RESULT >= 18.5 & BMI_NEAREST_RESULT < 25 ~ "Normal weight",
    BMI_NEAREST_RESULT >= 25 & BMI_NEAREST_RESULT < 30 ~ "Overweight",
    BMI_NEAREST_RESULT >= 30 ~ "Obese"
  )) %>%
  mutate(weight_category = factor(weight_category, levels = c("Underweight", "Normal weight", "Overweight", "Obese")))

cutaneous_melanom_tmb_distribution <- bodycomp_melanoma_intersection %>%
  ggplot(aes(x = CVR_TMB_SCORE)) +
  geom_histogram() +
  labs(
    x = "TMB Score",
    y = "Number of Patients",
    title = "TMB Distribution in Cutaneous Melanoma"
  ) +
  theme_bw() +
  facet_wrap(~CANCER_TYPE_DETAILED)

cutaneous_melanom_tmb_distribution

bodycomp_melanoma_intersection %>%
  group_by(CANCER_TYPE_DETAILED) %>%
  summarize(n_distinct(PATIENT_ID))

# verify non-relationship with an adjusted linear model
# repeat with firth's logistic regression (outcome: TMB-high, predictor: body composition)
adjusted_results_firth <- bodycomp_melanoma_intersection %>%
  filter(GENDER == "Male") %>%
  select(
    PATIENT_ID,
    CANCER_TYPE_DETAILED,
    age_at_scan,
    stage_cdm_derived,
    ancestry_label,
    GENDER,
    CVR_TMB_SCORE,
    bmi_at_scan,
    l3sat_area,
    l3vat_area,
    l3total_muscle_area
  ) %>%
  pivot_longer(
    cols = c(bmi_at_scan, l3sat_area, l3vat_area, l3total_muscle_area),
    names_to = "body_compartment",
    values_to = "area"
  ) %>%
  distinct() %>%
  filter(!is.na(CVR_TMB_SCORE)) %>%
  mutate(tmb_group = factor(CVR_TMB_SCORE >= 10, levels = c(FALSE, TRUE))) %>%
  group_by(body_compartment, CANCER_TYPE_DETAILED) %>%
  group_modify(~ {
    ctrl <- logistf::logistf.control(maxit = 1000, maxstep = 0.5)
    model <- logistf::logistf(tmb_group ~ scale(area) + age_at_scan + ancestry_label + stage_cdm_derived,
      data = .x, control = ctrl
    )
    se <- sqrt(diag(model$var))
    # Penalized profile-likelihood limits (logistf's ci.lower/ci.upper), which
    # match the likelihood-ratio p-values in model$prob. The previous
    # pl = FALSE + estimate +/- 1.96 * SE combination gave Wald intervals AND
    # Wald p-values.
    data.frame(
      term      = names(model$coefficients),
      estimate  = model$coefficients,
      std.error = se,
      conf.low  = model$ci.lower,
      conf.high = model$ci.upper,
      p.value   = model$prob
    ) %>%
      filter(term == "scale(area)")
  }) %>%
  ungroup()

adjusted_results_firth_plot <- adjusted_results_firth %>%
  # BH within each cancer type: four compartments tested per disease
  group_by(CANCER_TYPE_DETAILED) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    body_compartment = case_when(
      body_compartment == "bmi_at_scan" ~ "BMI",
      body_compartment == "l3sat_area" ~ "SAT",
      body_compartment == "l3vat_area" ~ "VAT",
      body_compartment == "l3total_muscle_area" ~ "SKM"
    ),
    significance = sig_level(p.value, padj),
    body_compartment = factor(body_compartment, levels = c("BMI", "SAT", "VAT", "SKM"))
  ) %>%
  ggplot(aes(x = estimate, y = body_compartment, color = significance)) +
  ggstats::geom_stripped_rows(aes(y = body_compartment), inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", linewidth = 0.3) +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high), linewidth = 0.3, lineend = "round") +
  geom_point(size = 2, shape = 16) +
  labs(x = "Adjusted log(OR) per SD", y = NULL, color = NULL) +
  scale_color_manual(values = SIG_COLOURS, drop = FALSE) +
  significance_key(size = 2) +
  theme_std() +
  theme(legend.position = "bottom") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) +
  facet_wrap(~CANCER_TYPE_DETAILED, scales = "free_x")

ggsave_wrapper("code_outputs/figs/impact/body_comp/merkel_cell_carcinoma/body_comp/comparison_vs_cutaneous_melanoma/adjusted_results_firth_plot_skin_cancers.pdf",
  adjusted_results_firth_plot,
  width = 3.5,
  height = 3
)
