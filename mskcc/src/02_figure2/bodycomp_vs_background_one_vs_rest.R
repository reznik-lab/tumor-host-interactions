# ============================================================================
# Produces: Figure 2 (A); Supplementary Figure 2 (B, C, D); Supplementary Table 11, 12
# ============================================================================
# Question: is the average patient of a given cancer type physiologically
# different from the average other cancer patient in the sample, after adjusting
# for covariates?
#
# One-vs-rest variant of bodycomp_vs_background.R. For each tested cancer we fit
# a binary model:
#
#   outcome ~ I(cancer_type == target) + covariates
#
# and report the coefficient on the target indicator. That coefficient is the
# covariate-adjusted difference  mean(target) - mean(all other patients).
#
# The "rest" baseline is PATIENT-weighted (every other patient counts equally),
# so it represents "the average other cancer patient" -- larger cancer types
# contribute more. This is the intended difference from the emmeans "eff"
# (level-equal-weighted grand mean) contrast in bodycomp_vs_background.R, and
# from "del.eff" (which would weight each other cancer TYPE equally).
#
# @author: Ethan Tse
# @date: June 25th, 2026

library(broom)
source("./src/00_shared/load_libraries.R")

# load body composition data
body_comp_data <- read_tsv("./code_outputs/stable_data/msk/body_comp/metadata/non_clean_master_body_comp_annotated_metadata.tsv") %>%
  filter(cancer_type_detailed != "Merkel Cell Carcinoma")

# For each cancer with >= min_n patients, fit outcome ~ (cancer vs rest) + covariates
# and pull the indicator coefficient. Cancers below min_n are not tested as targets
# but still contribute to the "rest" baseline.
run_one_vs_rest <- function(
    outcome_var,
    body_data,
    predictors,
    min_n = 50) {
  covariates <- setdiff(predictors, "cancer_type_detailed")

  # Drop bmi_at_scan from covariates if it's the outcome
  if (outcome_var == "bmi_at_scan") {
    covariates <- setdiff(covariates, "bmi_at_scan")
  }

  # Drop covariates with no variation in this cohort. The stage-restricted
  # sensitivity cohorts below are defined by stage_cdm_derived_granular == "4.0",
  # within which stage_cdm_derived is constant ("Stage 4"), so lm() would fail
  # with "contrasts can be applied only to factors with 2 or more levels".
  # Only constant terms are dropped, so the main cohort keeps the full set.
  constant_covariates <- covariates[vapply(
    covariates,
    function(v) dplyr::n_distinct(body_data[[v]], na.rm = TRUE) < 2,
    logical(1)
  )]
  if (length(constant_covariates) > 0) {
    message(
      "Dropping covariates with no variation in this cohort: ",
      paste(constant_covariates, collapse = ", ")
    )
    covariates <- setdiff(covariates, constant_covariates)
  }

  target_cancers <- body_data %>%
    count(cancer_type_detailed) %>%
    filter(n >= min_n) %>%
    pull(cancer_type_detailed)

  result_list <- lapply(target_cancers, function(target) {
    # 2-level factor: "rest" is the reference, so the "target" coefficient is
    # target - rest.
    dat <- body_data %>%
      mutate(grp = factor(
        ifelse(cancer_type_detailed == target, "target", "rest"),
        levels = c("rest", "target")
      ))

    form <- reformulate(c("grp", covariates), response = outcome_var)
    mod <- lm(form, data = dat)

    tidy(mod, conf.int = TRUE) %>%
      filter(term == "grptarget") %>%
      transmute(
        sex = "All",
        # match the "<Cancer> effect" naming the heatmap script expects
        cancer_type = paste0(target, " effect"),
        estimate,
        std.error,
        conf.low,
        conf.high,
        statistic,
        df = mod$df.residual,
        p.value
      )
  })

  bind_rows(result_list)
}

predictors <- c(
  "cancer_type_detailed",
  "age_at_scan",
  "smoking_predictions_3_classes",
  "ancestry_label",
  "stage_cdm_derived",
  "gender"
)

outcomes <- c(
  "l3vat_area",
  "l3sat_area",
  "l3total_muscle_area",
  "liver_volume",
  "pancreas_volume",
  "spleen_volume",
  "kidney_volume",
  "liver_median_hu",
  "pancreas_median_hu",
  "spleen_median_hu",
  "kidney_median_hu",
  "l3vat_median",
  "l3sat_median",
  "l3total_muscle_median_hu"
)

combined_results <- bind_rows(lapply(outcomes, function(compartment) {
  run_one_vs_rest(compartment, body_comp_data, predictors = predictors) %>%
    mutate(body_compartment = compartment)
})) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  arrange(sex, padj)

create_dir_if_needed("./code_outputs/results/tables/body_comp/bodycomp_vs_background/lm_w_emmeans/")
write_tsv(
  combined_results %>% arrange(padj),
  "./code_outputs/results/tables/body_comp/bodycomp_vs_background/lm_w_emmeans/bodycomp_vs_background_population_one_vs_rest_contrasts.tsv"
)

# ── Sensitivity analyses ──────────────────────────────────────────────────────
# Re-run the same one-vs-rest contrasts restricted to (1) patients with
# metastatic disease at diagnosis, (2) treatment-naive ("clean") scans, and
# (3) the intersection of both, to confirm the main contrasts are not driven
# by stage mix or by scans taken after treatment has already altered body
# composition.

sensitivity_cohorts <- list(
  metastatic = body_comp_data %>% filter(stage_cdm_derived_granular == "4.0"),
  clean_scan = body_comp_data %>% filter(clean_scan == TRUE),
  metastatic_clean_scan = body_comp_data %>% filter(stage_cdm_derived_granular == "4.0", clean_scan == TRUE)
)

for (cohort_name in names(sensitivity_cohorts)) {
  cohort_results <- bind_rows(lapply(outcomes, function(compartment) {
    run_one_vs_rest(compartment, sensitivity_cohorts[[cohort_name]], predictors = predictors) %>%
      mutate(body_compartment = compartment)
  })) %>%
    mutate(padj = p.adjust(p.value, method = "BH")) %>%
    arrange(sex, padj)

  write_tsv(
    cohort_results %>% arrange(padj),
    paste0(
      "./code_outputs/results/tables/body_comp/bodycomp_vs_background/lm_w_emmeans/",
      "bodycomp_vs_background_population_one_vs_rest_contrasts_", cohort_name, ".tsv"
    )
  )
}
