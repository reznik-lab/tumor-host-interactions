# ============================================================================
# Produces: Supplementary Table 22
# ============================================================================
################
# Are there BMI/body composition associated pathways across cancer types?
# (SEX-STRATIFIED)
#
# In this script, we specifically test if oncogenic signaling pathways are
# enriched as a function of BMI and body composition metrics (L3 muscle, VAT, SAT)
# using logistic regression, adjusting for relevant covariates.
#
# This is the sex-stratified counterpart to pathways_logistic_regression.R:
# the same discovery model is fit separately within males and females, so
# gender is dropped as a covariate and used as a stratifying variable instead.
#
# @author: Ethan Tse
# @date: Jul 3rd, 2026
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
  filter(pathway != "immune_evasion_pathway") %>%
  filter(pathway == "pi3k_akt_mtor_pathway")

# fit the same model separately within each sex stratum
sexes <- long_pathway_data_full %>%
  filter(!is.na(gender)) %>%
  distinct(gender) %>%
  pull(gender)

bmi_logreg_results_by_sex <- purrr::map_dfr(sexes, function(sex) {
  long_pathway_data_sex <- long_pathway_data_full %>%
    filter(gender == sex)

  # enforce minimum number of events (within stratum)
  min_event_pathways_sex <- long_pathway_data_sex %>%
    group_by(cancer_type_detailed, pathway) %>%
    summarize(num_altered = sum(is_pathway_altered == 1, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(num_altered)) %>%
    filter(num_altered > 50)

  long_pathway_data_sex <- long_pathway_data_sex %>%
    semi_join(min_event_pathways_sex, by = c("cancer_type_detailed", "pathway"))

  cancer_types_sufficient_samples <- long_pathway_data_sex %>%
    group_by(cancer_type_detailed) %>%
    summarize(num_patients = n_distinct(patient_id), .groups = "drop") %>%
    filter(num_patients >= 200) %>%
    pull(cancer_type_detailed)

  # bmi-pathway associations across cancer types, within this sex stratum
  long_pathway_data_sex %>%
    filter(cancer_type_detailed %in% cancer_types_sufficient_samples) %>%
    filter(!is.na(is_pathway_altered)) %>%
    group_by(cancer_type_detailed, pathway) %>%
    group_modify(~ {
      # enforce at least 50 altered samples to run regression
      if (sum(.x$is_pathway_altered == 1, na.rm = TRUE) < 50) {
        return(tibble())
      }

      n_samples <- nrow(.x)
      n_altered <- sum(.x$is_pathway_altered == 1, na.rm = TRUE)
      tryCatch(
        {
          # gender is dropped: it is constant within a sex stratum
          logreg_fit <- glm(is_pathway_altered ~ scale(bmi_first_result) + current_age_deid + ancestry_label + smoking_predictions_3_classes + stage_highest_recorded + prior_med_to_msk + cvr_tmb_score, data = .x, family = "binomial")
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
    ungroup() %>%
    mutate(gender = sex)
})

# BH padj correction (within each sex stratum)
bmi_logreg_results_tidy <- bmi_logreg_results_by_sex %>%
  filter(grepl("scale", term)) %>%
  group_by(gender) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  arrange(padj)

dir.create("./code_outputs/results/tables/bmi/impact/pathway_cancer_bmi_associations", recursive = TRUE, showWarnings = FALSE)
write_tsv(bmi_logreg_results_tidy, "./code_outputs/results/tables/bmi/impact/pathway_cancer_bmi_associations/snp_pathway_enrichment_by_bmi_logistic_regression_results_sex_stratified.tsv")
