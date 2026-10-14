# ============================================================================
# Produces: Supplementary Figure 5 (G, H); Supplementary Table 24
# ============================================================================
################
# Sensitivity analysis: does the PI3K/AKT/mTOR-VAT association in EGFR-mutant
# LUAD survive excluding STK11 from the pathway definition?
#
# lung_pik3ca_bodycomp_analysis.R currently defines "pathway altered" using
# the FULL PI3K/AKT/mTOR gene list, i.e. INCLUDING STK11 (as of the Aug 17
# 2026 "update figures" commit, which deliberately disabled the prior
# `setdiff(pi3k_genes, "stk11")` exclusion — see that script's header). STK11
# loss is independently associated with cachexia after tumor diagnosis, so a
# VAT association driven by the PI3K/AKT/mTOR pathway as a whole could in
# principle just be STK11 riding along. This script re-fits the same Firth
# models, restricted to EGFR-mutant LUAD, with STK11 removed from the
# pathway-altered gene set (the pre-Aug-17 definition), to test whether the
# association is robust to that exclusion.
#
# Four arms are fit, crossing the two axes the manuscript claims separately,
# and mirroring lung_pik3ca_bodycomp_analysis.R's Plots 2-4 on the
# STK11-inclusive definition:
#   unadjusted             -- matches that script's results_l3_unadj
#   bmi_adjusted           -- matches its results_l3_bmiadj (adds bmi_at_scan)
#   clean_scan             -- matches its results_l3_clean (treatment-proximal only)
#   clean_scan_bmi_adjusted-- treatment-proximal AND adjusted for bmi_at_scan
# Without these, "independent of STK11", "independent of overall body size" and
# "holds in treatment-proximal scans" are each only established one at a time
# against the STK11-inclusive definition, and whether they hold jointly has no
# answer in the outputs.
#
# READ THE clean_scan ARMS WITH CARE. Restricting to treatment-proximal scans
# leaves n = 133 with only 14 pathway-altered patients, against ~12 model
# parameters (ancestry alone contributes 5 dummies) -- roughly one event per
# parameter. Firth penalization keeps the estimates finite, but at that ratio
# the covariate-adjusted estimates are close to unidentified and the intervals
# are correspondingly wide. n_params and events_per_param are written to the
# output so this is visible in the table rather than only in this comment. The
# same caveat applies to the STK11-inclusive clean_scan arm that the manuscript
# already reports (n = 133, 15 altered).
#
# @author: Ethan Tse
# @date: August 31st, 2026
################

source("./src/00_shared/load_libraries.R")
source("./src/00_shared/obesity_constants.R")
pi3k_genes <- pi3k_akt_mtor_pathway
pi3k_genes_no_stk11 <- setdiff(pi3k_genes, "stk11")

pathway_data <- read_tsv("./code_outputs/stable_data/msk/impact/pathways/sample_level_mutation_w_pathway_annotations.tsv")
bodycomp_data <- read_tsv("./code_outputs/stable_data/msk/body_comp/metadata/non_clean_master_body_comp_annotated_metadata.tsv")
mutation_mtx <- read_tsv("./code_outputs/stable_data/msk/impact/snvs/matrix/panel_corrected/bmi_oncokb_annotated_first_sample_oncogenic_mutation_mtx.tsv") %>%
  select(SAMPLE_ID, KRAS, EGFR) %>%
  janitor::clean_names()
mutation_data <- read_tsv("./code_outputs/stable_data/msk/impact/snvs/dataframes/bmi_oncokb_annotated_first_sample_oncogenic_mutations.tsv")

l3_metrics <- c("l3sat_area", "l3vat_area", "l3total_muscle_area")
organ_metrics <- c("liver_volume", "pancreas_volume", "spleen_volume", "kidney_volume")

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
  left_join(pathway_data %>% select(sample_id, pi3k_akt_mtor_pathway), by = "sample_id") %>%
  left_join(mutation_mtx, by = "sample_id") %>%
  filter(!is.na(pi3k_akt_mtor_pathway))

pi3k_altered_samples_no_stk11 <- mutation_data %>%
  filter(tolower(Hugo_Symbol) %in% pi3k_genes_no_stk11) %>%
  distinct(SAMPLE_ID) %>%
  pull(SAMPLE_ID)

# EGFR-mutant LUAD only (the subgroup this sensitivity analysis is for),
# pathway status recomputed excluding STK11.
luad_egfr_no_stk11 <- body_comp_w_pi3k %>%
  filter(cancer_type_detailed == "Lung Adenocarcinoma", egfr == 1) %>%
  mutate(
    pi3k_akt_mtor_pathway = if_else(sample_id %in% pi3k_altered_samples_no_stk11, 1, 0),
    cancer_type_detailed = "EGFR-mut"
  )

# ── Firth helper (matches lung_pik3ca_bodycomp_analysis.R) ────────────────────

firth_tidy <- function(model, terms_of_interest) {
  idx <- which(model$terms %in% terms_of_interest)
  tibble(
    term = model$terms[idx],
    estimate = model$coefficients[idx],
    conf.low = model$ci.lower[idx],
    conf.high = model$ci.upper[idx],
    std.error = sqrt(diag(vcov(model)))[idx],
    p.value = model$prob[idx],
    or = exp(estimate),
    or.conf.low = exp(conf.low),
    or.conf.high = exp(conf.high)
  )
}

# adjust_bmi = TRUE adds bmi_at_scan as a covariate, except when bmi_at_scan is
# itself the metric being tested (which would be collinear) -- the same rule as
# lung_pik3ca_bodycomp_analysis.R, so the STK11-excluded estimates are directly
# comparable to the STK11-inclusive ones arm for arm.
run_firth_marginal <- function(df, metrics, adjust_bmi = FALSE) {
  df %>%
    pivot_longer(cols = all_of(metrics), names_to = "var", values_to = "value") %>%
    group_by(cancer_type_detailed, var) %>%
    mutate(value = as.numeric(scale(value))) %>%
    group_modify(~ {
      this_var <- .y$var
      covars <- c(
        "age_at_scan", "ancestry_label",
        "stage_cdm_derived", "smoking_predictions_3_classes"
      )
      if (n_distinct(.x$gender) > 1) covars <- c(covars, "gender")
      if (n_distinct(.x$clean_scan) > 1) covars <- c(covars, "clean_scan")
      if (adjust_bmi && this_var != "bmi_at_scan") covars <- c(covars, "bmi_at_scan")
      model_vars <- c("value", "pi3k_akt_mtor_pathway", covars)
      n_effective <- sum(complete.cases(.x[, model_vars]))
      model <- logistf::logistf(
        as.formula(paste("pi3k_akt_mtor_pathway ~ value +", paste(covars, collapse = " + "))),
        data = .x
      )
      n_altered <- sum(.x$pi3k_akt_mtor_pathway == 1, na.rm = TRUE)
      n_params <- length(model$coefficients)
      firth_tidy(model, "value") %>%
        mutate(
          n_samples = nrow(.x),
          n_altered = n_altered,
          n_effective = n_effective,
          n_params = n_params,
          # Events per parameter. Below ~5 the adjusted estimate is carrying
          # little information regardless of its p-value; see header.
          events_per_param = round(n_altered / n_params, 2)
        )
    }) %>%
    ungroup()
}

# No multiple-testing correction is applied in this script, deliberately. These
# are post-hoc, confirmatory re-tests of one association that was already
# identified and corrected for in the discovery analysis
# (lung_pik3ca_bodycomp_analysis.R, where BH spans all compartments and driver
# subtypes). Re-testing that single hit under alternative pathway definitions,
# covariate sets and cohort restrictions is not a new discovery family, and
# BH-correcting a confirmatory re-test penalises the estimate for the number of
# ways it was probed rather than for the number of hypotheses screened. Raw
# p-values are reported; the discovery-stage correction is what governs the
# claim, and the Results cite raw p throughout this paragraph.
results_egfr_no_stk11_unadj <- run_firth_marginal(
  luad_egfr_no_stk11,
  metrics = l3_metrics, adjust_bmi = FALSE
)

results_egfr_no_stk11_bmiadj <- run_firth_marginal(
  luad_egfr_no_stk11,
  metrics = l3_metrics, adjust_bmi = TRUE
)

# Treatment-proximal ("clean scan") restriction, matching Plot 4 of the main
# script. clean_scan is constant after this filter, so the helper's
# n_distinct(clean_scan) > 1 test drops it from the covariates automatically.
luad_egfr_no_stk11_clean <- luad_egfr_no_stk11 %>% filter(clean_scan == TRUE)

results_egfr_no_stk11_clean <- run_firth_marginal(
  luad_egfr_no_stk11_clean,
  metrics = l3_metrics, adjust_bmi = FALSE
)

results_egfr_no_stk11_clean_bmiadj <- run_firth_marginal(
  luad_egfr_no_stk11_clean,
  metrics = l3_metrics, adjust_bmi = TRUE
)

results_egfr_no_stk11 <- bind_rows(
  results_egfr_no_stk11_unadj %>% mutate(model = "unadjusted"),
  results_egfr_no_stk11_bmiadj %>% mutate(model = "bmi_adjusted"),
  results_egfr_no_stk11_clean %>% mutate(model = "clean_scan"),
  results_egfr_no_stk11_clean_bmiadj %>% mutate(model = "clean_scan_bmi_adjusted")
) %>%
  # BH within each model arm, matching run_firth_marginal()'s convention in
  # lung_pik3ca_bodycomp_analysis.R so the two STK11 definitions are comparable.
  group_by(model) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  arrange(model, cancer_type_detailed, p.value)

cat("\n===== EGFR-mutant LUAD, PI3K/AKT/mTOR pathway (STK11 excluded) ~ L3 body composition =====\n")
cat("-- unadjusted --\n")
print(results_egfr_no_stk11_unadj)
cat("-- additionally adjusted for BMI at scan --\n")
print(results_egfr_no_stk11_bmiadj)
cat("-- treatment-proximal scans only (see header: ~1 event per parameter) --\n")
print(results_egfr_no_stk11_clean)
cat("-- treatment-proximal scans only, adjusted for BMI at scan --\n")
print(results_egfr_no_stk11_clean_bmiadj)

out_dir <- "./code_outputs/results/tables/body_comp/pi3k"
create_dir_if_needed(out_dir)
write_tsv(
  results_egfr_no_stk11,
  file.path(out_dir, "pi3k_bodycomp_firth_luad_egfr_stk11_excluded.tsv")
)
