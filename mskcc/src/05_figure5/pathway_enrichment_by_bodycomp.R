# ============================================================================
# Produces: Figure 5 (A); Supplementary Figure 5 (A); Supplementary Table 22, 23
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
# @updated: May 24th, 2026
################

source("./src/00_shared/load_libraries.R")
source("./src/00_shared/pathway_constants.R")

pathway_data <- read_tsv("./code_outputs/stable_data/msk/impact/pathways/sample_level_mutation_w_pathway_annotations.tsv")
full_metadata <- read_tsv("./code_outputs/stable_data/msk/metadata/impact_patient_and_sample_metadata_first_surgery.tsv")

pathway_data_bmi <- pathway_data %>%
  left_join(full_metadata, by = "patient_id", suffix = c("", ".meta")) %>%
  # Drop patients with no BMI (the exposure) and no cancer type (the grouping
  # variable). Without these, samples missing a cancer type are pooled into a single
  # NA "cancer type" that is then tested like a real one -- it clears the >= 100
  # patient / > 10 altered thresholds, adds spurious rows to the BH denominator, and
  # breaks downstream plotting code that branches on the cancer-type name.
  filter(!is.na(bmi_nearest_result)) %>%
  filter(!is.na(cancer_type_detailed)) %>%
  filter(gene_panel %in% c("IMPACT341", "IMPACT410", "IMPACT468", "IMPACT505"))


# Discovery with full BMI cohort
# do with full cohort (91344 unique patients)
long_pathway_data_full <- pathway_data_bmi %>%
  select(
    patient_id,
    bmi_nearest_result,
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


cancer_types_sufficient_samples <- long_pathway_data_full %>%
  group_by(cancer_type_detailed) %>%
  summarize(num_patients = n_distinct(patient_id)) %>%
  filter(num_patients >= 100) %>%
  pull(cancer_type_detailed)

# ── Firth helper ──────────────────────────────────────────────────────────────

firth_tidy <- function(model, term) {
  idx <- which(model$terms == term)
  tibble(
    term      = term,
    estimate  = model$coefficients[idx],
    conf.low  = model$ci.lower[idx],
    conf.high = model$ci.upper[idx],
    std.error = sqrt(diag(vcov(model)))[idx],
    p.value   = model$prob[idx]
  )
}

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
    model_df <- .x %>% mutate(bmi_nearest_result = as.numeric(scale(bmi_nearest_result)))

    tryCatch(
      {
        if (length(unique_genders) > 1) {
          logreg_fit <- logistf::logistf(is_pathway_altered ~ bmi_nearest_result + gender + current_age_deid + ancestry_label + smoking_predictions_3_classes + stage_cdm_derived + prior_med_to_msk, data = model_df)
        } else {
          logreg_fit <- logistf::logistf(is_pathway_altered ~ bmi_nearest_result + current_age_deid + ancestry_label + smoking_predictions_3_classes + stage_cdm_derived + prior_med_to_msk, data = model_df)
        }
        firth_tidy(logreg_fit, "bmi_nearest_result") %>%
          mutate(n_samples = n_samples, n_altered = n_altered)
      },
      error = function(e) {
        message(sprintf(
          "Firth logistic regression failed for cancer type: %s, pathway: %s. Error: %s",
          unique(.y$cancer_type_detailed), unique(.y$pathway), e$message
        ))
        tibble() # return empty tibble on error
      }
    )
  }) %>%
  ungroup()

# odds ratio + 95% CI (exponentiated Firth penalized-likelihood log-odds CI), then BH padj correction
bmi_logreg_results_tidy <- bmi_logreg_results %>%
  mutate(
    odds_ratio = exp(estimate),
    or_low = exp(conf.low),
    or_high = exp(conf.high),
    padj = p.adjust(p.value, method = "BH")
  ) %>%
  arrange(p.value)

create_dir_if_needed("./code_outputs/results/tables/bmi/impact/pathway_cancer_bmi_associations")
write_tsv(bmi_logreg_results_tidy, "./code_outputs/results/tables/bmi/impact/pathway_cancer_bmi_associations/snp_pathway_enrichment_by_bmi_logistic_regression_results.tsv")
