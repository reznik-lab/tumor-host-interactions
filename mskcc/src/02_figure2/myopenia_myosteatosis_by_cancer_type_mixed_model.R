# ============================================================================
# Produces: Figure 2; Supplementary Figure 2 (F); Supplementary Table 15
# (a figure listed without a panel letter contributes to a composite
#  panel whose letter could not be attributed unambiguously)
# ============================================================================
####################
# Is myopenia associated with myosteatosis within each cancer type?
#
# The pooled co-occurrence prevalence in myopenia_cooccurrence.R answers this
# for the cohort as a whole but cannot say whether the association holds (or
# varies in strength) within individual cancer types. Here we fit a separate
# adjusted logistic regression within each major cancer type:
#
#   myosteatosis ~ myopenia + age + sex + stage + ancestry + smoking
#
# and report the myopenia coefficient (log odds ratio, Wald 95% CI, p-value)
# per cancer type, with Benjamini-Hochberg correction across cancer types.
#
# This was previously a single joint model with a myopenia x cancer type
# interaction, with the per-cancer-type contrasts pulled out via emmeans (and,
# before that, a glmer() with a non-identifiable random intercept). Fitting
# one model per cancer type is simpler to state in a paper and gives
# materially the same answer -- every per-cancer-type estimate and
# significance call agrees with the joint model's.
#
# Body composition definitions (myopenia, myosteatosis) are kept consistent
# with abnormal_physiology.R and myopenia_cachexia_mixed_model.R.
#
# @author: Ethan Tse
# @date: August 27th, 2026
####################


source("./src/00_shared/load_libraries.R")

body_comp_data <- read_tsv("./code_outputs/stable_data/msk/body_comp/metadata/non_clean_master_body_comp_annotated_metadata.tsv") %>%
  filter(cancer_type_detailed != "Merkel Cell Carcinoma") %>%
  # recalibrate muscle area and muscle radiodensity to pre-contrast, then
  # compute SMI and myopenia/myosteatosis (consistent with
  # abnormal_physiology.R and myopenia_cachexia_mixed_model.R)
  mutate(
    l3total_muscle_area_precontrast = 1.0628 * l3total_muscle_area,
    l3total_muscle_median_hu_precontrast = 0.9930 * l3total_muscle_median_hu - 7.5291
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
  # myosteatosis is defined relative to the muscle radiodensity distribution,
  # as the sex-specific lowest quartile of pre-contrast muscle median HU
  mutate(myosteatosis = case_when(
    # Same NA guard as myopenia: an unmeasured muscle HU is unknown
    # myosteatosis status, not absence of myosteatosis.
    is.na(l3total_muscle_median_hu_precontrast) ~ NA_real_,
    (gender == "Male") & l3total_muscle_median_hu_precontrast < quantile(l3total_muscle_median_hu_precontrast[gender == "Male"], 0.25, na.rm = TRUE) ~ 1,
    (gender == "Female") & l3total_muscle_median_hu_precontrast < quantile(l3total_muscle_median_hu_precontrast[gender == "Female"], 0.25, na.rm = TRUE) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(myopenia_status = factor(
    ifelse(myopenia == 1, "Myopenia", "No myopenia"),
    levels = c("No myopenia", "Myopenia")
  ))

# Covariates to adjust for (consistent with myopenia_cachexia_mixed_model.R
# and myopenia_one_vs_rest.R): age, sex, stage, genetic ancestry, and smoking.
# BMI is deliberately excluded: it is a body composition measure itself, so
# adjusting for it here would condition out part of the myopenia-myosteatosis
# relationship of interest.
covariates <- c(
  "age_at_scan",
  "gender",
  "stage_cdm_derived",
  "ancestry_label",
  "smoking_predictions_3_classes"
)

# Cancer types need both myopenia levels well represented to support a
# within-cancer-type model of the myopenia effect.
major_cancer_types <- body_comp_data %>%
  filter(!is.na(myopenia), !is.na(myosteatosis)) %>%
  count(cancer_type_detailed, myopenia) %>%
  tidyr::pivot_wider(names_from = myopenia, values_from = n, values_fill = 0) %>%
  filter(`0` >= 10, `1` >= 10) %>%
  pull(cancer_type_detailed)

model_data <- body_comp_data %>%
  filter(cancer_type_detailed %in% major_cancer_types) %>%
  filter(!is.na(myopenia), !is.na(myosteatosis)) %>%
  filter(if_all(all_of(covariates), ~ !is.na(.)))

# Adjusted logistic regression of myosteatosis on myopenia within one cancer
# type. Returns a one-row tidy table for the myopenia term.
fit_one_cancer_type <- function(cancer, data, covariates) {
  d <- data %>% filter(cancer_type_detailed == cancer) %>% droplevels()

  # Drop covariates that are constant within this cancer type (sex in the
  # single-sex cancers), which would otherwise be non-identifiable. The joint
  # interaction model absorbed these automatically; fitting per cancer type
  # makes the drop explicit, so it is recorded in the output.
  usable_covariates <- covariates[vapply(
    covariates,
    function(v) length(unique(d[[v]][!is.na(d[[v]])])) > 1,
    logical(1)
  )]

  model_formula <- as.formula(paste(
    "myosteatosis ~ myopenia_status",
    if (length(usable_covariates)) {
      paste("+", paste(usable_covariates, collapse = " + "))
    } else {
      ""
    }
  ))

  model <- tryCatch(
    glm(model_formula, data = d, family = binomial),
    error = function(e) NULL
  )
  if (is.null(model)) {
    return(NULL)
  }

  cat("\n===== ", cancer, " =====\n", sep = "")
  print(summary(model))

  term <- "myopenia_statusMyopenia"
  coefs <- summary(model)$coefficients
  if (!term %in% rownames(coefs)) {
    return(NULL)
  }

  # Wald CIs, matching what the emmeans contrasts reported before and what
  # the Methods' confidence interval table already states for this analysis.
  ci <- stats::confint.default(model)[term, ]

  tibble(
    cancer_type_detailed = cancer,
    estimate = coefs[term, "Estimate"],
    SE = coefs[term, "Std. Error"],
    conf.low = ci[[1]],
    conf.high = ci[[2]],
    z.ratio = coefs[term, "z value"],
    p.value = coefs[term, "Pr(>|z|)"],
    odds_ratio = exp(coefs[term, "Estimate"]),
    or_conf_low = exp(ci[[1]]),
    or_conf_high = exp(ci[[2]]),
    n = nrow(d),
    n_myopenia = sum(d$myopenia == 1),
    n_myosteatosis = sum(d$myosteatosis == 1),
    dropped_covariates = paste(setdiff(covariates, usable_covariates), collapse = ", ")
  )
}

# Fit every eligible cancer type and BH-correct across them.
fit_by_cancer_type <- function(data, covariates) {
  cancers <- sort(unique(data$cancer_type_detailed))

  results <- purrr::map_dfr(cancers, fit_one_cancer_type,
    data = data, covariates = covariates
  ) %>%
    mutate(p_adj_bh = p.adjust(p.value, method = "BH")) %>%
    arrange(p.value)

  cat("\n===== Per-cancer-type myopenia effect on myosteatosis (BH-adjusted) =====\n")
  print(results)

  results
}

out_dir <- "./code_outputs/results/tables/body_comp/abnormal_physiology"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

results <- fit_by_cancer_type(model_data, covariates)

write_tsv(
  results,
  file.path(out_dir, "myopenia_myosteatosis_by_cancer_type_mixed_model.tsv")
)

# ---------------------------------------------------------------------------
# Sensitivity analysis: restrict to treatment-naive ("clean scan") patients
# ---------------------------------------------------------------------------

clean_scan_major_cancer_types <- model_data %>%
  filter(clean_scan == TRUE) %>%
  count(cancer_type_detailed, myopenia) %>%
  tidyr::pivot_wider(names_from = myopenia, values_from = n, values_fill = 0) %>%
  filter(`0` >= 10, `1` >= 10) %>%
  pull(cancer_type_detailed)

model_data_clean_scan_subset <- model_data %>%
  filter(clean_scan == TRUE, cancer_type_detailed %in% clean_scan_major_cancer_types)

results_clean_scan_subset <- fit_by_cancer_type(model_data_clean_scan_subset, covariates)

write_tsv(
  results_clean_scan_subset,
  file.path(out_dir, "myopenia_myosteatosis_by_cancer_type_mixed_model_clean_scan_subset.tsv")
)
