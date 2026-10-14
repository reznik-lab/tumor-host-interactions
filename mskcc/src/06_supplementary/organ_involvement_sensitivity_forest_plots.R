# ============================================================================
# Produces: Supplementary Figure 2 (B, D); Supplementary Table 12
# ============================================================================
################
# Organ-involvement sensitivity forest plots
#
# Among the three cancers associated with significantly reduced muscle area
# (Pancreatic, Stomach, and Lung Adenocarcinoma; one-vs-rest contrasts from
# bodycomp_vs_background_one_vs_rest.R), compares liver/kidney/spleen volume
# point estimates across the main cohort and the metastatic, treatment-naive
# ("clean scan"), and intersected sensitivity cohorts appended to that script.
#
# @author: Ethan Tse
# @date: July 22nd, 2026
################

library(broom)
source("./src/00_shared/load_libraries.R")

results_dir <- "./code_outputs/results/tables/body_comp/bodycomp_vs_background/lm_w_emmeans/"

analysis_files <- c(
  "Main"                  = paste0(results_dir, "bodycomp_vs_background_population_one_vs_rest_contrasts.tsv"),
  "Metastatic"            = paste0(results_dir, "bodycomp_vs_background_population_one_vs_rest_contrasts_metastatic.tsv"),
  "Treatment-naive"       = paste0(results_dir, "bodycomp_vs_background_population_one_vs_rest_contrasts_clean_scan.tsv"),
  "Metastatic + tx-naive" = paste0(results_dir, "bodycomp_vs_background_population_one_vs_rest_contrasts_metastatic_clean_scan.tsv")
)

target_cancers <- c(
  "Pancreatic Adenocarcinoma effect",
  "Stomach Adenocarcinoma effect",
  "Lung Adenocarcinoma effect"
)
target_organs <- c("liver_volume", "kidney_volume", "spleen_volume")

organ_results <- purrr::imap_dfr(analysis_files, function(path, variant) {
  read_tsv(path, show_col_types = FALSE) %>%
    filter(cancer_type %in% target_cancers, body_compartment %in% target_organs) %>%
    mutate(analysis_variant = variant)
}) %>%
  mutate(
    analysis_variant = factor(analysis_variant, levels = names(analysis_files)),
    cancer_type = factor(
      gsub(" effect", "", cancer_type),
      levels = c("Pancreatic Adenocarcinoma", "Stomach Adenocarcinoma", "Lung Adenocarcinoma")
    ),
    organ = factor(
      case_when(
        body_compartment == "liver_volume" ~ "Liver",
        body_compartment == "kidney_volume" ~ "Kidney",
        body_compartment == "spleen_volume" ~ "Spleen"
      ),
      levels = c("Liver", "Kidney", "Spleen")
    ),
    significance = sig_level(p.value, padj)
  )

organ_forest <- organ_results %>%
  ggplot(aes(x = estimate, y = cancer_type)) +
  ggstats::geom_stripped_rows() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1) +
  geom_linerange(
    aes(xmin = conf.low, xmax = conf.high, color = significance),
    lineend = "round", linewidth = 0.4
  ) +
  geom_point(aes(color = significance), size = 1.8, shape = 16) +
  scale_colour_manual(
    name = NULL,
    values = SIG_COLOURS,
    drop = FALSE
  ) +
  significance_key(size = 1.8) +
  facet_grid(organ ~ analysis_variant, scales = "free_x") +
  labs(x = "Adj. estimate (cm³)", y = NULL) +
  theme_std() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 6),
    panel.spacing.x = unit(2, "pt")
  )

create_dir_if_needed("./code_outputs/results/figs/body_comp/bodycomp_vs_background/")
ggsave_wrapper(
  plot = organ_forest,
  filename = "./code_outputs/results/figs/body_comp/bodycomp_vs_background/organ_involvement_sensitivity_forest.pdf",
  width = 6,
  height = 4
)

# ── Positive control: spleen volume in hematological disease vs. rest ───────
#
# Hematological malignancies (identified by HEMEPACT gene panel assignment)
# are known to be associated with splenomegaly. Testing the same one-vs-rest
# spleen-volume contrast on each individual heme cancer type is a positive
# control: if the approach captures genuine anatomic signal, these patients
# should show a substantially larger adjusted spleen volume than the rest of
# the cohort. Restricted to heme cancer types with >= 30 patients, since the
# main analysis's min_n = 50 threshold excludes all of them.

body_comp_data <- read_tsv("./code_outputs/stable_data/msk/body_comp/metadata/non_clean_master_body_comp_annotated_metadata.tsv") %>%
  filter(cancer_type_detailed != "Merkel Cell Carcinoma")

heme_cancer_types <- body_comp_data %>%
  filter(grepl("HEME", gene_panel)) %>%
  count(cancer_type_detailed) %>%
  filter(n >= 30) %>%
  pull(cancer_type_detailed)

heme_spleen_result <- purrr::map_dfr(heme_cancer_types, function(target) {
  dat <- body_comp_data %>%
    mutate(grp = factor(
      ifelse(cancer_type_detailed == target, "target", "rest"),
      levels = c("rest", "target")
    ))

  # Covariate set matches the main one-vs-rest model in
  # bodycomp_vs_background_one_vs_rest.R.
  mod <- lm(
    spleen_volume ~ grp + age_at_scan + smoking_predictions_3_classes +
      ancestry_label + stage_cdm_derived + gender,
    data = dat
  )

  tidy(mod, conf.int = TRUE) %>%
    filter(term == "grptarget") %>%
    transmute(cancer_type_detailed = target, estimate, conf.low, conf.high, p.value)
}) %>%
  mutate(
    padj = p.adjust(p.value, method = "BH"),
    significance = sig_level(p.value, padj)
  ) %>%
  arrange(estimate)

heme_spleen_forest <- heme_spleen_result %>%
  ggplot(aes(x = estimate, y = reorder(cancer_type_detailed, estimate))) +
  ggstats::geom_stripped_rows() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1) +
  geom_linerange(
    aes(xmin = conf.low, xmax = conf.high, color = significance),
    lineend = "round", linewidth = 0.4
  ) +
  geom_point(aes(color = significance), size = 1.8, shape = 16) +
  scale_colour_manual(
    name = NULL,
    values = SIG_COLOURS,
    drop = FALSE
  ) +
  significance_key(size = 1.8) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 25)) +
  labs(x = "Spleen volume estimate (cm³)", y = NULL) +
  theme_std() +
  theme(
    legend.position = "top",
    axis.ticks.y = element_blank()
  )

ggsave_wrapper(
  plot = heme_spleen_forest,
  filename = "./code_outputs/results/figs/body_comp/bodycomp_vs_background/heme_spleen_positive_control_forest.pdf",
  width = 2.3,
  height = 1.7
)
