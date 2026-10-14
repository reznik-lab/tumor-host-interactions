# ============================================================================
# Produces: Supplementary Figure 5 (A); Supplementary Table 22
# ============================================================================
################
# Are BMI/pathway associations robust to using a treatment-naive BMI label?
#
# Sensitivity analysis for pathway_enrichment_by_bodycomp.R: restricts the
# cohort to patients with a treatment-naive BMI (first BMI recorded strictly
# before their first treatment start date; patients without one are dropped)
# and to patients with no prior medications to MSK, then re-tests ONLY the
# cancer/pathway associations that were significant (padj < 0.05) in the
# full-cohort discovery analysis. Also re-tests the same hits with an
# additional, stricter restriction: on top of the tx-naive BMI + no-prior-meds
# cohort, samples must have a surgery/specimen-acquisition date that precedes
# the patient's first treatment date. Full-cohort vs tx-naive estimates are
# then plotted side-by-side in one forest plot.
#
# @author: Ethan Tse
# @date: 2026-07-06
################

source("./src/00_shared/load_libraries.R")
source("./src/00_shared/pathway_constants.R")

# ── Discovery hits to re-test ─────────────────────────────────────────────────

discovery_results <- read_tsv(
  "./code_outputs/results/tables/bmi/impact/pathway_cancer_bmi_associations/snp_pathway_enrichment_by_bmi_logistic_regression_results.tsv"
)

discovery_hits <- discovery_results %>%
  filter(padj < 0.05) %>%
  distinct(cancer_type_detailed, pathway)

pathway_data <- read_tsv("./code_outputs/stable_data/msk/impact/pathways/sample_level_mutation_w_pathway_annotations.tsv")
full_metadata <- read_tsv("./code_outputs/stable_data/msk/metadata/impact_patient_and_sample_metadata_first_surgery.tsv")

pathway_data_bmi_full <- pathway_data %>%
  left_join(full_metadata, by = "patient_id", suffix = c("", ".meta")) %>%
  filter(gene_panel %in% c("IMPACT341", "IMPACT410", "IMPACT468", "IMPACT505")) %>%
  filter(cancer_type_detailed %in% discovery_hits$cancer_type_detailed)

# ── Compute first treatment-naive BMI per patient ────────────────────────────

longitudinal_bmi <- read_tsv("./data/private/impact/msk_solid_heme/data_timeline_bmi.txt") %>%
  janitor::clean_names() %>%
  dplyr::rename("bmi" = "result") %>%
  filter(bmi > 10) # filter out implausible BMI values

tx_data <- read_tsv("./data/private/impact/msk_solid_heme/data_timeline_treatment.txt") %>%
  janitor::clean_names()

surgery_data <- read_tsv("data/private/impact/msk_solid_heme/data_timeline_specimen_surgery.txt") %>%
  janitor::clean_names()

first_tx_date <- tx_data %>%
  group_by(patient_id) %>%
  summarise(first_tx_date = min(start_date, na.rm = TRUE), .groups = "drop")

tx_naive_bmi <- longitudinal_bmi %>%
  left_join(first_tx_date, by = "patient_id") %>%
  filter(is.na(first_tx_date) | start_date < first_tx_date) %>%
  group_by(patient_id) %>%
  slice_min(start_date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(patient_id, tx_naive_bmi = bmi)

# ── Restrict to tx-naive BMI + no prior medications to MSK ───────────────────

pathway_data_bmi_tx_naive <- pathway_data_bmi_full %>%
  inner_join(tx_naive_bmi, by = "patient_id") %>%
  filter(prior_med_to_msk == "No prior medications")

message(
  "Patients retained after tx-naive BMI + no-prior-medications filter: ",
  n_distinct(pathway_data_bmi_tx_naive$patient_id)
)

long_pathway_data_sensitivity <- pathway_data_bmi_tx_naive %>%
  select(
    patient_id,
    tx_naive_bmi,
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
    contains("pathway")
  ) %>%
  pivot_longer(cols = contains("pathway"), names_to = "pathway", values_to = "is_pathway_altered") %>%
  # The upstream pathway annotation still carries an immune_evasion_pathway
  # column; the manuscript does not test it, so it is dropped here. This
  # filter is load-bearing -- removing it would re-enter immune evasion into
  # the screen and change every BH-adjusted p-value.
  # Dropped here, before the models, so it also leaves the BH correction denominator.
  filter(pathway != "immune_evasion_pathway") %>%
  # only re-test the specific cancer/pathway pairs significant in the discovery analysis
  semi_join(discovery_hits, by = c("cancer_type_detailed", "pathway"))

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

# re-test only the discovery hits, in the tx-naive/no-prior-meds cohort
# (prior_med_to_msk is dropped from the model since it's constant by construction)
sensitivity_results <- long_pathway_data_sensitivity %>%
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
    model_df <- .x %>% mutate(tx_naive_bmi = as.numeric(scale(tx_naive_bmi)))

    tryCatch(
      {
        if (length(unique_genders) > 1) {
          logreg_fit <- logistf::logistf(is_pathway_altered ~ tx_naive_bmi + gender + current_age_deid + ancestry_label + smoking_predictions_3_classes + stage_cdm_derived, data = model_df)
        } else {
          logreg_fit <- logistf::logistf(is_pathway_altered ~ tx_naive_bmi + current_age_deid + ancestry_label + smoking_predictions_3_classes + stage_cdm_derived, data = model_df)
        }
        firth_tidy(logreg_fit, "tx_naive_bmi") %>%
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
sensitivity_results_tidy <- sensitivity_results %>%
  mutate(
    odds_ratio = exp(estimate),
    or_low = exp(conf.low),
    or_high = exp(conf.high),
    padj = p.adjust(p.value, method = "BH")
  ) %>%
  arrange(p.value)


# ── Re-test the same discovery hits in the full cohort, additionally adjusting
# for TMB (checks whether the BMI-pathway association is confounded by TMB) ──

long_pathway_data_tmb_adj <- pathway_data_bmi_full %>%
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
    contains("pathway")
  ) %>%
  pivot_longer(cols = contains("pathway"), names_to = "pathway", values_to = "is_pathway_altered") %>%
  # The upstream pathway annotation still carries an immune_evasion_pathway
  # column; the manuscript does not test it, so it is dropped here. This
  # filter is load-bearing -- removing it would re-enter immune evasion into
  # the screen and change every BH-adjusted p-value.
  # Dropped here, before the models, so it also leaves the BH correction denominator.
  filter(pathway != "immune_evasion_pathway") %>%
  semi_join(discovery_hits, by = c("cancer_type_detailed", "pathway"))

full_cohort_tmb_results <- long_pathway_data_tmb_adj %>%
  filter(!is.na(is_pathway_altered)) %>%
  group_by(cancer_type_detailed, pathway) %>%
  group_modify(~ {
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
          logreg_fit <- logistf::logistf(is_pathway_altered ~ bmi_nearest_result + gender + current_age_deid + ancestry_label + smoking_predictions_3_classes + stage_cdm_derived + prior_med_to_msk + cvr_tmb_score, data = model_df)
        } else {
          logreg_fit <- logistf::logistf(is_pathway_altered ~ bmi_nearest_result + current_age_deid + ancestry_label + smoking_predictions_3_classes + stage_cdm_derived + prior_med_to_msk + cvr_tmb_score, data = model_df)
        }
        firth_tidy(logreg_fit, "bmi_nearest_result") %>%
          mutate(n_samples = n_samples, n_altered = n_altered)
      },
      error = function(e) {
        message(sprintf(
          "Firth logistic regression (full cohort + TMB adj.) failed for cancer type: %s, pathway: %s. Error: %s",
          unique(.y$cancer_type_detailed), unique(.y$pathway), e$message
        ))
        tibble()
      }
    )
  }) %>%
  ungroup()

full_cohort_tmb_results_tidy <- full_cohort_tmb_results %>%
  mutate(
    odds_ratio = exp(estimate),
    or_low = exp(conf.low),
    or_high = exp(conf.high),
    padj = p.adjust(p.value, method = "BH")
  ) %>%
  arrange(p.value)


# ── Re-test the same discovery hits, further restricted to samples whose
# surgery/specimen-acquisition date precedes the patient's first treatment date
# (on top of the tx-naive BMI + no-prior-medications restriction above) ──────

tx_naive_sample_ids <- surgery_data %>%
  select(patient_id, sample_id, surgery_start_date = start_date) %>%
  left_join(first_tx_date, by = "patient_id") %>%
  filter(is.na(first_tx_date) | surgery_start_date < first_tx_date) %>%
  pull(sample_id)

pathway_data_bmi_tx_naive_sample <- pathway_data_bmi_tx_naive %>%
  filter(sample_id %in% tx_naive_sample_ids)

message(
  "Patients retained after tx-naive BMI + tx-naive sample + no-prior-medications filter: ",
  n_distinct(pathway_data_bmi_tx_naive_sample$patient_id)
)

long_pathway_data_tx_naive_sample <- pathway_data_bmi_tx_naive_sample %>%
  select(
    patient_id,
    tx_naive_bmi,
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
    contains("pathway")
  ) %>%
  pivot_longer(cols = contains("pathway"), names_to = "pathway", values_to = "is_pathway_altered") %>%
  # The upstream pathway annotation still carries an immune_evasion_pathway
  # column; the manuscript does not test it, so it is dropped here. This
  # filter is load-bearing -- removing it would re-enter immune evasion into
  # the screen and change every BH-adjusted p-value.
  # Dropped here, before the models, so it also leaves the BH correction denominator.
  filter(pathway != "immune_evasion_pathway") %>%
  semi_join(discovery_hits, by = c("cancer_type_detailed", "pathway"))

# (prior_med_to_msk is dropped from the model since it's constant by construction)
tx_naive_sample_results <- long_pathway_data_tx_naive_sample %>%
  filter(!is.na(is_pathway_altered)) %>%
  group_by(cancer_type_detailed, pathway) %>%
  group_modify(~ {
    if (sum(.x$is_pathway_altered == 1, na.rm = TRUE) < 10) {
      return(tibble())
    }

    n_samples <- nrow(.x)
    n_altered <- sum(.x$is_pathway_altered == 1, na.rm = TRUE)
    unique_genders <- unique(.x$gender)
    model_df <- .x %>% mutate(tx_naive_bmi = as.numeric(scale(tx_naive_bmi)))

    tryCatch(
      {
        if (length(unique_genders) > 1) {
          logreg_fit <- logistf::logistf(is_pathway_altered ~ tx_naive_bmi + gender + current_age_deid + ancestry_label + smoking_predictions_3_classes + stage_cdm_derived, data = model_df)
        } else {
          logreg_fit <- logistf::logistf(is_pathway_altered ~ tx_naive_bmi + current_age_deid + ancestry_label + smoking_predictions_3_classes + stage_cdm_derived, data = model_df)
        }
        firth_tidy(logreg_fit, "tx_naive_bmi") %>%
          mutate(n_samples = n_samples, n_altered = n_altered)
      },
      error = function(e) {
        message(sprintf(
          "Firth logistic regression (tx-naive BMI + tx-naive sample) failed for cancer type: %s, pathway: %s. Error: %s",
          unique(.y$cancer_type_detailed), unique(.y$pathway), e$message
        ))
        tibble()
      }
    )
  }) %>%
  ungroup()

tx_naive_sample_results_tidy <- tx_naive_sample_results %>%
  mutate(
    odds_ratio = exp(estimate),
    or_low = exp(conf.low),
    or_high = exp(conf.high),
    padj = p.adjust(p.value, method = "BH")
  ) %>%
  arrange(p.value)


# ── Forest plot: full cohort vs full cohort + TMB adj. vs tx-naive/no-prior-meds ─

comparison_df <- bind_rows(
  discovery_results %>%
    semi_join(discovery_hits, by = c("cancer_type_detailed", "pathway")) %>%
    mutate(model = "Full cohort"),
  full_cohort_tmb_results_tidy %>%
    mutate(model = "Full cohort + TMB adj."),
  sensitivity_results_tidy %>%
    mutate(model = "Tx-naive BMI, no prior meds"),
  tx_naive_sample_results_tidy %>%
    mutate(model = "Tx-naive BMI + sample, no prior meds")
) %>%
  mutate(
    # prettify pathway name, consistent with pathways_bmi_association_volcano.R
    pathway_pretty = case_when(
      grepl("pi3k", pathway, ignore.case = TRUE) ~ "PI3K/AKT/mTOR",
      grepl("hippo", pathway, ignore.case = TRUE) ~ "Hippo",
      grepl("rtk.*mapk|mapk.*rtk", pathway, ignore.case = TRUE) ~ "RTK/RAS/MAPK",
      grepl("cell_cycle", pathway, ignore.case = TRUE) ~ "Cell Cycle",
      grepl("notch", pathway, ignore.case = TRUE) ~ "Notch",
      grepl("wnt", pathway, ignore.case = TRUE) ~ "Wnt",
      grepl("tgfb", pathway, ignore.case = TRUE) ~ "TGF-β",
      grepl("p53", pathway, ignore.case = TRUE) ~ "p53",
      grepl("myc", pathway, ignore.case = TRUE) ~ "Myc",
      grepl("nrf2", pathway, ignore.case = TRUE) ~ "NRF2",
      TRUE ~ gsub("_", " ", pathway)
    ),
    cancer_wrapped = vapply(
      as.character(cancer_type_detailed),
      function(x) paste(strwrap(x, width = 18), collapse = "\n"),
      character(1)
    ),
    row_label = paste0(cancer_wrapped, ":\n", pathway_pretty),
    model = factor(model, levels = c("Tx-naive BMI + sample, no prior meds", "Tx-naive BMI, no prior meds", "Full cohort + TMB adj.", "Full cohort")),
    # figs.qmd re-maps this panel's colour from `model` onto `significance`
    # (each model already has its own facet there), so the three-level scale
    # has to exist in the saved object.
    significance = sig_level(p.value, padj),
    # Stars are dropped by figs.qmd for the manuscript panel but kept here so
    # the standalone PDF is readable: one star = nominal, two = survives BH.
    sig_label = case_when(
      padj < 0.05 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# order rows by the full-cohort effect size
row_order <- comparison_df %>%
  filter(model == "Full cohort") %>%
  arrange(estimate) %>%
  pull(row_label)
comparison_df <- comparison_df %>%
  mutate(row_label = factor(row_label, levels = row_order))

tx_naive_bmi_forest <- comparison_df %>%
  ggplot(aes(x = estimate, y = row_label, group = model)) +
  ggstats::geom_stripped_rows() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1) +
  geom_linerange(
    aes(xmin = conf.low, xmax = conf.high, color = model),
    lineend = "round",
    linewidth = 0.1,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    size = 1.5, aes(colour = model, shape = model),
    position = position_dodge(width = 0.6)
  ) +
  geom_text(
    aes(x = conf.high, label = sig_label, group = model),
    position = position_dodge(width = 0.6),
    hjust = -0.3,
    size = 2
  ) +
  scale_colour_manual(
    name = NULL,
    values = c(
      "Full cohort" = "black",
      "Full cohort + TMB adj." = "#1B7837",
      "Tx-naive BMI, no prior meds" = "darkred",
      "Tx-naive BMI + sample, no prior meds" = "#2166AC"
    )
  ) +
  scale_shape_manual(
    name = NULL,
    values = c(
      "Full cohort" = 16,
      "Full cohort + TMB adj." = 15,
      "Tx-naive BMI, no prior meds" = 17,
      "Tx-naive BMI + sample, no prior meds" = 18
    )
  ) +
  labs(
    x = paste(strwrap("Adj. log(OR) per SD BMI", width = 30), collapse = "\n"),
    y = NULL
  ) +
  guides(
    colour = guide_legend(nrow = 2),
    shape = guide_legend(nrow = 2)
  ) +
  theme_std() +
  theme(
    legend.position = "top",
    axis.ticks.y = element_blank()
  )

create_dir_if_needed("./code_outputs/results/figs/bmi/impact/pathway_cancer_bmi_associations")
ggsave_wrapper(
  "./code_outputs/results/figs/bmi/impact/pathway_cancer_bmi_associations/tx_naive_bmi_pathway_forest_plot.pdf",
  tx_naive_bmi_forest,
  width = 3.5,
  height = 3
)
