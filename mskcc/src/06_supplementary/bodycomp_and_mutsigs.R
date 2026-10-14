# ============================================================================
# Produces: Supplementary Figure 3
# (a figure listed without a panel letter contributes to a composite
#  panel whose letter could not be attributed unambiguously)
# ============================================================================
##################
# Analyze deepsig mutation signatures vs. body composition (L3 areas)
# Mirrors bmi_and_mutsigs.R but uses CT-derived L3 areas as primary exposure.
#
# @author: Ethan Tse
# @date: 2026-06-21
##################

source("./src/00_shared/load_libraries.R")

long_mutsig_df <- read_tsv("./code_outputs/stable_data/msk/impact/mutation_signatures/deepsig/deepsig_exposures_long_format.tsv") %>%
  mutate(CANCER_TYPE_DETAILED = case_when(
    CANCER_TYPE_DETAILED %in% c("Colorectal Adenocarcinoma", "Rectal Adenocarcinoma", "Colon Adenocarcinoma") ~ "Colorectal Adenocarcinoma",
    TRUE ~ CANCER_TYPE_DETAILED
  ))

l3_metrics <- c("l3total_muscle_area", "l3vat_area", "l3sat_area")

body_comp <- read_tsv("./code_outputs/stable_data/msk/body_comp/metadata/non_clean_master_body_comp_annotated_metadata.tsv") %>%
  dplyr::select(sample_id, age_at_scan, all_of(l3_metrics))

bodycomp_mutsigs <- long_mutsig_df %>%
  inner_join(body_comp, by = c(sid = "sample_id"))

# regression for each signature × cancer type × L3 metric
df_res <- bodycomp_mutsigs %>%
  group_by(CANCER_TYPE_DETAILED) %>%
  filter(n_distinct(sid) >= 50) %>%
  ungroup() %>%
  mutate(weighted_exposure = exposure * M) %>%
  mutate(across(all_of(l3_metrics), ~ scale(.x))) %>%
  pivot_longer(cols = all_of(l3_metrics), names_to = "l3_metric", values_to = "l3_value") %>%
  group_by(CANCER_TYPE_DETAILED, signature, l3_metric) %>%
  group_modify(~ {
    cat_vars <- c("GENDER", "STAGE_CDM_DERIVED", "ANCESTRY_LABEL", "SMOKING_PREDICTIONS_3_CLASSES", "PRIOR_MED_TO_MSK")
    varying_cat_vars <- cat_vars[sapply(cat_vars, function(v) length(unique(na.omit(.x[[v]]))) > 1)]
    rhs_vars <- c("l3_value", "age_at_scan", varying_cat_vars)
    formula <- as.formula(paste0("weighted_exposure ~ ", paste(rhs_vars, collapse = " + ")))

    model <- tryCatch(lm(formula = formula, data = .x, na.action = "na.omit"), error = function(e) NULL)
    na_row <- tibble(term = "l3_value", estimate = NA, std.error = NA, statistic = NA, p.value = NA)
    if (is.null(model) || model$df.residual == 0) {
      return(na_row)
    }

    broom::tidy(model, conf.int = TRUE) %>% dplyr::filter(term == "l3_value")
  }) %>%
  ungroup() %>%
  # lenient padj due to smaller sample sizes
  group_by(l3_metric) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  ungroup()


# volcano plot faceted by L3 metric
bodycomp_volcano <- df_res %>%
  mutate(
    signif = ifelse(padj < 0.05, "padj < 0.05", "NS"),
    cancer_code = cancer_codes[CANCER_TYPE_DETAILED]
  ) %>%
  ggplot(aes(x = estimate, y = -log10(p.value))) +
  geom_point(aes(colour = signif), size = 2, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_text_repel(
    data = function(df) df %>% filter(padj < 0.05),
    aes(label = paste0(cancer_code, ":", signature)),
    box.padding = 0.25,
    max.overlaps = 15,
    show.legend = FALSE,
    segment.color = "grey50",
    force = 50,
    size = 2.5
  ) +
  scale_color_manual(
    name = NULL,
    values = c("NS" = "darkgrey", "padj < 0.05" = "#c44c4c"),
    breaks = c("padj < 0.05", "NS"),
    labels = c("padj < 0.05" = "padj < 0.05", "NS" = "NS")
  ) +
  facet_wrap(~l3_metric,
    scales = "free_x", nrow = 1,
    labeller = labeller(l3_metric = c(
      "l3sat_area" = "SAT", "l3vat_area" = "VAT", "l3total_muscle_area" = "SKM"
    ))
  ) +
  labs(
    x = "Adj. estimate of weighted exposure per SD area",
    y = expression(-log[10](italic(p) * "-value")),
    color = "Significance"
  ) +
  theme_std() +
  theme(legend.position = "top", legend.justification = "center")

create_dir_if_needed("./code_outputs/results/figs/bodycomp/mutagenesis/mutation_signatures/deepsig/")
ggsave_wrapper("./code_outputs/results/figs/bodycomp/mutagenesis/mutation_signatures/deepsig/bodycomp_mutsig_volcano.pdf",
  bodycomp_volcano,
  width = 3.7,
  height = 2
)

# significance barplot per metric

# distribution of (unadjusted) p-values per tissue, to check for
# uniformity under the null vs. enrichment near 0 indicating true signal
