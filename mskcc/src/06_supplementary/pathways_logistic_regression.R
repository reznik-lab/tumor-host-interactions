# ============================================================================
# Produces: Supplementary Table 22
# ============================================================================
################
# Are there BMI/body composition associated pathways across cancer types?
#
# In this script, we specifically test if oncogenic signaling pathways are
# enriched as a function of BMI and body composition metrics (L3 muscle, VAT, SAT)
# using logistic regression, adjusting for relevant covariates.
#
# We first do a discovery analysis with BMI in the full cohort
# and then do a post-hoc analysis with body composition metrics
#
# @author: Ethan Tse
# @date: Jan 27th, 2026
################

source("./src/00_shared/load_libraries.R")

pathway_data <- read_tsv("./code_outputs/stable_data/msk/impact/pathways/sample_level_mutation_w_pathway_annotations.tsv")
metadata <- read_tsv("./code_outputs/stable_data/msk/metadata/impact_patient_and_sample_metadata_first_surgery.tsv") %>%
  janitor::clean_names()

pathway_data <- pathway_data %>%
  left_join(metadata, by = c("patient_id" = "patient_id")) %>%
  filter(!is.na(bmi_first_result)) %>%
  filter(!is.na(cancer_type_detailed)) %>%
  # look at solid tumours only
  filter(gene_panel %in% c("IMPACT410", "IMPACT468", "IMPACT505", "IMPACT341"))

table(pathway_data$pi3k_akt_mtor_pathway, useNA = "always")

table(pathway_data$pi3k_akt_mtor_pathway, useNA = "always")

# Discovery with full BMI cohort
# do with full cohort (91344 unique patients)
long_pathway_data_full <- pathway_data %>%
  select(
    # cohort,
    patient_id,
    bmi_first_result,
    cancer_type_detailed,
    gender,
    current_age_deid,
    race,
    ancestry_label,
    cvr_tmb_score,
    smoking_predictions_3_classes,
    stage_cdm_derived_granular,
    stage_cdm_derived,
    stage_highest_recorded,
    prior_med_to_msk,
    # height_at_scan,
    # egfr,
    # kras,
    contains("pathway")
  ) %>%
  pivot_longer(cols = contains("pathway"), names_to = "pathway", values_to = "is_pathway_altered") %>%
  # The upstream pathway annotation still carries an immune_evasion_pathway
  # column; the manuscript does not test it, so it is dropped here. This
  # filter is load-bearing -- removing it would re-enter immune evasion into
  # the screen and change every BH-adjusted p-value.
  # Dropped here, before the models, so it also leaves the BH correction denominator.
  filter(pathway != "immune_evasion_pathway")

# enforce minimum number of events
min_event_pathways_full <- long_pathway_data_full %>%
  group_by(cancer_type_detailed, pathway) %>%
  summarize(num_altered = sum(is_pathway_altered == 1, na.rm = TRUE)) %>%
  arrange(desc(num_altered)) %>%
  filter(num_altered > 10)

long_pathway_data_full <- long_pathway_data_full %>%
  semi_join(min_event_pathways_full, by = c("cancer_type_detailed", "pathway"))

table(long_pathway_data_full$is_pathway_altered, useNA = "always")

cancer_types_sufficient_samples <- long_pathway_data_full %>%
  group_by(cancer_type_detailed) %>%
  summarize(num_patients = n_distinct(patient_id)) %>%
  filter(num_patients >= 100) %>%
  pull(cancer_type_detailed)

# bmi-pathway associations across cancer types
bmi_logreg_results <- long_pathway_data_full %>%
  filter(cancer_type_detailed %in% cancer_types_sufficient_samples) %>%
  filter(!is.na(is_pathway_altered)) %>%
  group_by(cancer_type_detailed, pathway) %>%
  group_modify(~ {
    # enforce at least 10 altered samples to run regression
    if (sum(.x$is_pathway_altered == 1, na.rm = TRUE) < 10) {
      return(tibble())
    }

    n_samples <- nrow(.x)
    n_altered <- sum(.x$is_pathway_altered == 1, na.rm = TRUE)
    unique_genders <- unique(.x$gender)
    tryCatch(
      {
        if (length(unique_genders) > 1) {
          logreg_fit <- glm(is_pathway_altered ~ scale(bmi_first_result) + gender + current_age_deid + ancestry_label + smoking_predictions_3_classes + stage_cdm_derived + prior_med_to_msk + cvr_tmb_score, data = .x, family = "binomial")
        } else {
          logreg_fit <- glm(is_pathway_altered ~ scale(bmi_first_result) + current_age_deid + ancestry_label + smoking_predictions_3_classes + stage_cdm_derived + prior_med_to_msk + cvr_tmb_score, data = .x, family = "binomial")
        }
        broom::tidy(logreg_fit) %>%
          mutate(n_samples = n_samples, n_altered = n_altered)
      },
      error = function(e) {
        message(sprintf(
          "GLM failed for cancer type: %s, pathway: %s. Error: %s",
          unique(.x$cancer_type_detailed), unique(.x$pathway), e$message
        ))
        tibble() # return empty tibble on error
        tibble(error_message = e$message)
      }
    )
  }) %>%
  ungroup()

# BH padj correction
bmi_logreg_results_tidy <- bmi_logreg_results %>%
  filter(grepl("scale", term)) %>%
  mutate(
    padj = p.adjust(p.value, method = "BH"),
    OR = exp(estimate),
    OR_ci_lower = exp(estimate - qnorm(0.975) * std.error),
    OR_ci_upper = exp(estimate + qnorm(0.975) * std.error)
  ) %>%
  arrange(p.value)

dir.create("./code_outputs/results/tables/bmi/impact/pathway_cancer_bmi_associations", recursive = TRUE, showWarnings = FALSE)
# NOTE: this unpenalized-glm variant used to write
# snp_pathway_enrichment_by_bmi_logistic_regression_results.tsv -- the same path as the
# Firth version in
# src/05_figure5/pathway_enrichment_by_bodycomp.R
# -- so whichever ran last silently won. The Firth version is the one wired into the
# figure pipelines and reported in the manuscript; this script now writes its own file so
# the two can coexist without clobbering each other.
write_tsv(bmi_logreg_results_tidy, "./code_outputs/results/tables/bmi/impact/pathway_cancer_bmi_associations/snp_pathway_enrichment_by_bmi_logistic_regression_results_glm_unpenalized.tsv")
