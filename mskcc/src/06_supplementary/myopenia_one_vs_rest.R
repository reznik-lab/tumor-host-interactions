# ============================================================================
# Produces: Supplementary Figure 2; Supplementary Table 14
# (a figure listed without a panel letter contributes to a composite
#  panel whose letter could not be attributed unambiguously)
# ============================================================================
####################
# One-vs-rest test: in which cancer types is myopenia (and myosteatosis,
# liver steatosis) enriched relative to all other cancer types combined?
#
# The pooled cachexia-prone-vs-other model in myopenia_cachexia_mixed_model.R
# assumes a single, homogeneous effect shared across PDAC/LUAD/STAD/ESCA/
# COADREAD. It cannot say which individual cancer type(s) actually drive that
# effect, or flag enrichment in a cancer type outside that pre-specified list.
#
# Here, for each cancer type with >= min_n patients, we fit:
#
#   outcome ~ I(cancer_type == target) + covariates
#
# on the full cohort (target vs. every other patient pooled as "rest"), and
# report the odds ratio on the target indicator. This mirrors the one-vs-rest
# design in
# src/02_figure2/bodycomp_vs_background_one_vs_rest.R,
# but with glm(family = binomial) in place of lm(), since myopenia/
# myosteatosis/liver_steatosis are binary rather than continuous outcomes.
#
# Body composition definitions (myopenia, myosteatosis, liver_steatosis) are
# kept consistent with abnormal_physiology.R and myopenia_cachexia_mixed_model.R.
#
# @author: Ethan Tse
# @date: August 26th, 2026
####################

library(broom)
source("./src/00_shared/load_libraries.R")

# cancer types considered cachexia-prone (consistent with abnormal_physiology.R
# and myopenia_cachexia_mixed_model.R), used only to flag rows in the output
cachexia_prone <- c("PDAC", "LUAD", "STAD", "ESCA", "COADREAD")

body_comp_data <- read_tsv("./code_outputs/stable_data/msk/body_comp/metadata/non_clean_master_body_comp_annotated_metadata.tsv") %>%
  filter(cancer_type_detailed != "Merkel Cell Carcinoma") %>%
  # recalibrate muscle area, muscle radiodensity, and liver radiodensity to
  # pre-contrast, then compute SMI, myopenia, myosteatosis, and liver steatosis
  # (consistent with abnormal_physiology.R)
  mutate(
    l3total_muscle_area_precontrast = 1.0628 * l3total_muscle_area,
    l3total_muscle_median_hu_precontrast = 0.9930 * l3total_muscle_median_hu - 7.5291,
    liver_median_hu_precontrast = case_when(
      (0.7434 * liver_median_hu - 12.387) < 80 ~ (0.7434 * liver_median_hu) - 12.387,
      TRUE ~ (0.21254 * liver_median_hu) + 31.508,
    )
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
  # liver steatosis: pre-contrast liver median HU below 57 (consistent with
  # abnormal_physiology.R)
  mutate(liver_steatosis = ifelse(liver_median_hu_precontrast < 57, 1, 0))

# Covariates to adjust for (consistent with myopenia_cachexia_mixed_model.R):
# age, sex, stage, genetic ancestry, and smoking.
covariates <- c(
  "age_at_scan",
  "gender",
  "stage_cdm_derived",
  "ancestry_label",
  "smoking_predictions_3_classes"
)

# For each cancer type with >= min_n patients, fit
#   outcome ~ (cancer type vs. every other patient) + covariates
# via logistic regression, and pull the odds ratio on the target indicator.
# Cancer types below min_n are not tested as targets, but still contribute
# to the "rest" baseline for every cancer type that is tested.
run_one_vs_rest_binary <- function(outcome_var, data, covariates, min_n = 50) {
  target_cancers <- data %>%
    count(cancer_type_detailed) %>%
    filter(n >= min_n) %>%
    pull(cancer_type_detailed)

  result_list <- lapply(target_cancers, function(target) {
    dat <- data %>%
      filter(!is.na(.data[[outcome_var]])) %>%
      filter(if_all(all_of(covariates), ~ !is.na(.))) %>%
      # 2-level factor: "rest" is the reference, so the "target" coefficient
      # is the target-vs-rest log-odds difference.
      mutate(grp = factor(
        ifelse(cancer_type_detailed == target, "target", "rest"),
        levels = c("rest", "target")
      ))

    # Guard on the TARGET cancer type's own event counts, not the pooled
    # target+rest cohort (which barely varies across targets and would let
    # this check pass even for a target with almost no events).
    target_rows <- dat[[outcome_var]][dat$grp == "target"]
    n_target_events <- sum(target_rows == 1, na.rm = TRUE)
    n_target_nonevents <- sum(target_rows == 0, na.rm = TRUE)
    if (n_target_events < 10 || n_target_nonevents < 10) {
      return(NULL)
    }

    form <- reformulate(c("grp", covariates), response = outcome_var)
    fit <- tryCatch(glm(form, data = dat, family = binomial), error = function(e) NULL)
    if (is.null(fit)) {
      return(NULL)
    }

    # Wald CIs (estimate +/- 1.96 * SE) throughout, rather than
    # stats::confint()'s profile likelihood with a Wald fallback. The profile
    # can fail inside approx() under near-separation for a small target cancer
    # type, so the fallback meant a single forest plot mixed two interval
    # types across its cancer types. Wald is what the Methods declare for
    # these models and matches
    # myopenia_myosteatosis_by_cancer_type_mixed_model.R.
    ci <- stats::confint.default(fit)

    broom::tidy(fit) %>%
      filter(term == "grptarget") %>%
      mutate(
        conf.low = ci[match(term, rownames(ci)), 1],
        conf.high = ci[match(term, rownames(ci)), 2],
        cancer_type_detailed = target,
        cancer_code = unname(cancer_codes[target]),
        n_total = nrow(dat),
        n_target = length(target_rows),
        n_target_events = n_target_events,
        odds_ratio = exp(estimate),
        or_conf_low = exp(conf.low),
        or_conf_high = exp(conf.high)
      )
  })

  bind_rows(result_list)
}

outcomes <- c("myopenia", "myosteatosis", "liver_steatosis")

# Fit the one-vs-rest test for every outcome x cancer-type combination, then
# BH-adjust p-values pooled across all of them (consistent with
# bodycomp_vs_background_one_vs_rest.R).
run_all_outcomes <- function(data, covariates, min_n = 50) {
  bind_rows(lapply(outcomes, function(outcome) {
    run_one_vs_rest_binary(outcome, data, covariates, min_n = min_n) %>%
      mutate(outcome = outcome)
  })) %>%
    mutate(
      cachexia_prone = cancer_code %in% cachexia_prone,
      padj = p.adjust(p.value, method = "BH")
    ) %>%
    arrange(outcome, padj)
}

out_dir <- "./code_outputs/results/tables/body_comp/abnormal_physiology"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Full cohort
# ---------------------------------------------------------------------------
one_vs_rest_results <- run_all_outcomes(body_comp_data, covariates)

cat("\n===== One-vs-rest enrichment (full cohort): cancer types with padj < 0.05, OR > 1 =====\n")
one_vs_rest_results %>%
  filter(padj < 0.05, odds_ratio > 1) %>%
  select(outcome, cancer_code, odds_ratio, or_conf_low, or_conf_high, p.value, padj, n_target, n_target_events) %>%
  print(n = Inf)

write_tsv(
  one_vs_rest_results,
  file.path(out_dir, "myopenia_myosteatosis_liver_steatosis_one_vs_rest.tsv")
)

# ── Sensitivity analysis: restrict to treatment-naive ("clean scan") patients ──
#
# Consistent with the clean-scan sensitivity analyses elsewhere in
# abnormal_physiology.R.
one_vs_rest_results_clean_scan_subset <- run_all_outcomes(
  body_comp_data %>% filter(clean_scan == TRUE),
  covariates
)

write_tsv(
  one_vs_rest_results_clean_scan_subset,
  file.path(out_dir, "myopenia_myosteatosis_liver_steatosis_one_vs_rest_clean_scan_subset.tsv")
)
