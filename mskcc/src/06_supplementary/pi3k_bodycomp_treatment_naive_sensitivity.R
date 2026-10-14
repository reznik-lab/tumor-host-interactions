# ============================================================================
# Produces: Supplementary Figure 5; Supplementary Table 23
# (a figure listed without a panel letter contributes to a composite
#  panel whose letter could not be attributed unambiguously)
# ============================================================================
################
# Sensitivity analysis: does the PI3K/AKT/mTOR ~ body-composition association
# in CRC/PDAC/ccRCC survive restricting the cohort to treatment-naive scans?
#
# pi3k_bodycomp_analysis.R fits these models on ALL scans (its
# `clean_scan == TRUE` filter is commented out) and does not adjust for
# `clean_scan`: the joint SAT+VAT+SKM model cannot carry it, and the marginal
# and joint estimates are reported side by side, so neither carries it. Scan
# curation is handled here instead, by restriction. This script re-fits the
# identical Firth models on the clean_scan == TRUE subset only, which is also
# the stronger claim the manuscript makes.
#
# Mirrors the on-treatment contrast in pi3k_bodycomp_analysis_non_treatment_naive.R
# and the LUAD-side clean-scan restriction already reported from
# lung_pik3ca_bodycomp_analysis.R.
#
# @author: Ethan Tse
# @date: September 4th, 2026
################

source("./src/00_shared/load_libraries.R")
pathway_data <- read_tsv("./code_outputs/stable_data/msk/impact/pathways/sample_level_mutation_w_pathway_annotations.tsv")
bodycomp_data <- read_tsv("./code_outputs/stable_data/msk/body_comp/metadata/non_clean_master_body_comp_annotated_metadata.tsv")

l3_metrics <- c("l3sat_area", "l3vat_area", "l3total_muscle_area")
organ_metrics <- c("liver_volume", "pancreas_volume", "spleen_volume", "kidney_volume")
body_comp_metrics <- c(l3_metrics, organ_metrics)

# Identical to pi3k_bodycomp_analysis.R's body_comp_w_pi3k EXCEPT for the
# clean_scan == TRUE restriction, which is the whole point of this script.
body_comp_w_pi3k <- bodycomp_data %>%
  janitor::clean_names() %>%
  select(
    dmp_id, sample_id, cancer_type_detailed,
    age_at_scan, ancestry_label,
    stage_cdm_derived_granular, stage_cdm_derived, stage_highest_recorded,
    smoking_predictions_3_classes, gender,
    clean_scan, bmi_at_scan,
    gene_panel,
    matches("^l3"),
    any_of(organ_metrics)
  ) %>%
  filter(clean_scan == TRUE) %>%
  left_join(pathway_data %>% select(sample_id, pi3k_akt_mtor_pathway), by = "sample_id") %>%
  filter(!is.na(pi3k_akt_mtor_pathway))

# ── Firth helpers (clean_scan dropped: constant within this cohort) ───────────

firth_tidy <- function(model, terms_of_interest) {
  idx <- which(model$terms %in% terms_of_interest)
  tibble(
    term      = model$terms[idx],
    estimate  = model$coefficients[idx],
    conf.low  = model$ci.lower[idx],
    conf.high = model$ci.upper[idx],
    std.error = sqrt(diag(vcov(model)))[idx],
    p.value   = model$prob[idx]
  )
}

run_firth_marginal <- function(df, metrics = body_comp_metrics) {
  available <- intersect(metrics, colnames(df))
  df %>%
    pivot_longer(cols = all_of(available), names_to = "var", values_to = "value") %>%
    group_by(cancer_type_detailed, var) %>%
    mutate(value = as.numeric(scale(value))) %>%
    group_modify(~ {
      n_samples <- nrow(.x)
      n_altered <- sum(.x$pi3k_akt_mtor_pathway == 1, na.rm = TRUE)
      if (n_altered < 10) {
        return(tibble())
      }
      # n_effective is computed over the variables the formula below actually
      # uses (the main script's version lists stage_cdm_derived_granular here
      # while modelling stage_cdm_derived -- kept consistent here).
      base_vars <- c(
        "value", "pi3k_akt_mtor_pathway", "age_at_scan",
        "ancestry_label", "stage_cdm_derived", "smoking_predictions_3_classes"
      )
      model_vars <- if (n_distinct(.x$gender) > 1) c(base_vars, "gender") else base_vars
      n_effective <- sum(complete.cases(.x[, model_vars]))
      tryCatch(
        {
          if (n_distinct(.x$gender) > 1) {
            model <- logistf::logistf(
              pi3k_akt_mtor_pathway ~ value + age_at_scan + ancestry_label +
                stage_cdm_derived + smoking_predictions_3_classes + gender,
              data = .x
            )
          } else {
            model <- logistf::logistf(
              pi3k_akt_mtor_pathway ~ value + age_at_scan + ancestry_label +
                stage_cdm_derived + smoking_predictions_3_classes,
              data = .x
            )
          }
          firth_tidy(model, "value") %>% mutate(n_samples = n_samples, n_altered = n_altered, n_effective = n_effective)
        },
        error = function(e) {
          message(sprintf("Firth marginal failed: %s / %s: %s", unique(.x$cancer_type_detailed), unique(.x$var), e$message))
          tibble()
        }
      )
    }) %>%
    ungroup()
}

run_firth_joint <- function(df, metrics) {
  available <- intersect(metrics, colnames(df))
  metrics_rhs <- paste(available, collapse = " + ")
  df %>%
    mutate(across(all_of(available), scale)) %>%
    group_by(cancer_type_detailed) %>%
    group_modify(~ {
      n_samples <- nrow(.x)
      n_altered <- sum(.x$pi3k_akt_mtor_pathway == 1, na.rm = TRUE)
      if (n_altered < 10) {
        return(tibble())
      }
      tryCatch(
        {
          covariates <- if (n_distinct(.x$gender) > 1) {
            "age_at_scan + ancestry_label + stage_cdm_derived_granular + smoking_predictions_3_classes + gender"
          } else {
            "age_at_scan + ancestry_label + stage_cdm_derived_granular + smoking_predictions_3_classes"
          }
          model <- logistf::logistf(
            as.formula(paste("pi3k_akt_mtor_pathway ~", metrics_rhs, "+", covariates)),
            data = .x
          )
          firth_tidy(model, available) %>%
            rename(var = term) %>%
            mutate(n_samples = n_samples, n_altered = n_altered)
        },
        error = function(e) {
          message(sprintf("Firth joint failed: %s: %s", unique(.x$cancer_type_detailed), e$message))
          tibble()
        }
      )
    }) %>%
    ungroup()
}

# ── CRC / PDAC / ccRCC, treatment-naive scans only ───────────────────────────

crc_pdac_rcc_cancers <- c("Colorectal Adenocarcinoma", "Pancreatic Adenocarcinoma", "Renal Clear Cell Carcinoma")

crc_pdac_rcc <- body_comp_w_pi3k %>%
  filter(cancer_type_detailed %in% crc_pdac_rcc_cancers)

message("Treatment-naive cohort sizes (scans with pathway calls):")
print(crc_pdac_rcc %>% count(cancer_type_detailed, pi3k_akt_mtor_pathway))

# BH within each metric family, matching pi3k_bodycomp_analysis.R
marginal_l3 <- run_firth_marginal(crc_pdac_rcc, metrics = l3_metrics) %>% mutate(padj = p.adjust(p.value, method = "BH"))
marginal_organ <- run_firth_marginal(crc_pdac_rcc, metrics = organ_metrics) %>% mutate(padj = p.adjust(p.value, method = "BH"))
joint_l3 <- run_firth_joint(crc_pdac_rcc, l3_metrics) %>% mutate(padj = p.adjust(p.value, method = "BH"))

pi3k_results_tx_naive <- bind_rows(
  marginal_l3 %>% mutate(model_type = "marginal", level = "L3"),
  marginal_organ %>% mutate(model_type = "marginal", level = "organ"),
  joint_l3 %>% mutate(model_type = "joint", level = "L3")
) %>%
  mutate(cohort = "Treatment-naive")


# ── Forest plot, same layout as pi3k_forest_crc_pdac_rcc (Figure 4) ──────────

label_results <- function(d) {
  d %>%
    mutate(
      metric_label = case_when(
        grepl("total_muscle_area", var) ~ "SKM",
        grepl("vat_area", var) ~ "VAT",
        grepl("sat_area", var) ~ "SAT",
        var == "liver_volume" ~ "Liver",
        var == "pancreas_volume" ~ "Pancreas",
        var == "spleen_volume" ~ "Spleen",
        var == "kidney_volume" ~ "Kidney"
      ),
      metric_label = factor(metric_label, levels = c("SKM", "VAT", "SAT", "Liver", "Pancreas", "Spleen", "Kidney")),
      significance = sig_level(p.value, padj)
    )
}

s1_cancer_labels <- c(
  "Colorectal Adenocarcinoma" = "Colorectal\nAdenocarcinoma",
  "Pancreatic Adenocarcinoma" = "Pancreatic\nAdenocarcinoma",
  "Renal Clear Cell Carcinoma" = "Renal Clear\nCell Carcinoma"
)

# Full-cohort estimates plotted alongside, so the panel shows the comparison
# the manuscript sentence is actually making rather than the restricted
# estimates in isolation.
full_cohort <- read_tsv("./code_outputs/results/tables/body_comp/pi3k/pi3k_bodycomp_firth_crc_pdac_rcc.tsv") %>%
  mutate(cohort = "All scans")

combined <- bind_rows(pi3k_results_tx_naive, full_cohort) %>%
  filter(model_type == "marginal", level %in% c("L3", "organ")) %>%
  filter(cancer_type_detailed %in% crc_pdac_rcc_cancers) %>%
  label_results() %>%
  mutate(
    metric_label = factor(metric_label, levels = rev(levels(metric_label))),
    compartment = factor(
      if_else(level == "organ", "Organ volume", "Tissue area"),
      levels = c("Tissue area", "Organ volume")
    ),
    cohort = factor(cohort, levels = c("All scans", "Treatment-naive"))
  )

forest_tx_naive <- combined %>%
  ggplot(aes(x = estimate, y = metric_label, shape = cohort)) +
  ggstats::geom_stripped_rows() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1) +
  geom_linerange(
    aes(xmin = conf.low, xmax = conf.high, color = significance),
    lineend = "round", linewidth = 0.1, position = position_dodge(width = 0.6)
  ) +
  geom_point(size = 1.4, aes(colour = significance), position = position_dodge(width = 0.6)) +
  scale_colour_manual(name = NULL, values = SIG_COLOURS, drop = FALSE) +
  significance_key(size = 1.4) +
  scale_shape_manual(name = NULL, values = c("All scans" = 16, "Treatment-naive" = 1)) +
  labs(x = "Adj. log(OR) per SD", y = NULL, fill = NULL) +
  theme_std() +
  theme(
    legend.position = "top",
    legend.box = "vertical",
    legend.spacing.y = unit(1, "pt"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 6),
    panel.spacing.x = unit(2, "pt")
  ) +
  facet_grid(
    compartment ~ cancer_type_detailed,
    scales = "free", space = "free_y",
    labeller = labeller(cancer_type_detailed = as_labeller(s1_cancer_labels))
  )

create_dir_if_needed("./code_outputs/results/figs/body_comp/pi3k")
ggsave_wrapper(
  plot = forest_tx_naive,
  filename = "./code_outputs/results/figs/body_comp/pi3k/pi3k_forest_crc_pdac_rcc_treatment_naive.pdf",
  width = 2.9,
  height = 2.4
)

# ── Console summary for the manuscript sentence ──────────────────────────────

message("\nL3 tissue-area estimates, treatment-naive vs all scans:")
print(
  combined %>%
    filter(level == "L3") %>%
    select(cancer_type_detailed, metric_label, cohort, estimate, conf.low, conf.high, p.value, padj, n_samples, n_altered) %>%
    arrange(cancer_type_detailed, metric_label, cohort) %>%
    as.data.frame()
)
