# ============================================================================
# Produces: Figure 5 (C); Supplementary Figure 5; Supplementary Table 23
# (a figure listed without a panel letter contributes to a composite
#  panel whose letter could not be attributed unambiguously)
# ============================================================================
################
# PI3K/AKT/MTOR Pathway and Body Composition — Dedicated Analysis
#
# In a pan-cancer screen of oncogenic signaling pathways, we found that
# alterations in the PI3K/AKT/MTOR pathway were significantly negatively associated with
# BMI and muscle mass in pancreatic ductal adenocarcinoma (PDAC) patients.
#
# Firth's penalized logistic regression of PI3K/AKT/MTOR pathway alteration
# against L1 and L3 body composition metrics (muscle, VAT, SAT).
# Both marginal (each metric separately) and joint (all metrics together) models.
#
# Covers CRC, PDAC, RCC. The lung adenocarcinoma analysis (PIK3CA-specific,
# with a driver-subtype interaction) lives in lung_pik3ca_bodycomp_analysis.R.
#
# @author: Ethan Tse
# @date: May 24th, 2026
################

source("./src/00_shared/load_libraries.R")

pathway_data <- read_tsv("./code_outputs/stable_data/msk/impact/pathways/sample_level_mutation_w_pathway_annotations.tsv")
bodycomp_data <- read_tsv("./code_outputs/stable_data/msk/body_comp/metadata/non_clean_master_body_comp_annotated_metadata.tsv")

l3_metrics <- c("l3sat_area", "l3vat_area", "l3total_muscle_area")
organ_metrics <- c("liver_volume", "pancreas_volume", "spleen_volume", "kidney_volume")
body_comp_metrics <- c(l3_metrics, organ_metrics)

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
  filter(!is.na(pi3k_akt_mtor_pathway))

# ── Firth helpers ─────────────────────────────────────────────────────────────

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
      # clean_scan is deliberately NOT a covariate here. run_firth_joint() below
      # cannot carry it, and the marginal and joint estimates are reported side by
      # side in the same figure and table -- keeping it in only one of them meant
      # the two were not nested and the single Methods formula that covers both
      # was wrong for one of them. Treatment-naive status is instead handled by
      # the clean-scan-restricted sensitivity refit.
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
          message(sprintf("Firth marginal failed: %s / %s: %s", unique(.y$cancer_type_detailed), unique(.y$var), e$message))
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
    group_by(cancer_type_detailed) %>%
    # Standardize WITHIN cancer type, matching run_firth_marginal above. Scaling
    # before the group_by() pooled the SD across CRC/PDAC/RCC, so it absorbed
    # between-cancer variation in body composition and the resulting "per SD"
    # estimates were not on the same footing as the marginal estimates they are
    # reported alongside in the same output table.
    mutate(across(all_of(available), ~ as.numeric(scale(.x)))) %>%
    group_modify(~ {
      n_samples <- nrow(.x)
      n_altered <- sum(.x$pi3k_akt_mtor_pathway == 1, na.rm = TRUE)
      if (n_altered < 10) {
        return(tibble())
      }
      tryCatch(
        {
          covariates <- if (n_distinct(.x$gender) > 1) {
            "age_at_scan + ancestry_label + stage_cdm_derived + smoking_predictions_3_classes + gender"
          } else {
            "age_at_scan + ancestry_label + stage_cdm_derived + smoking_predictions_3_classes"
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
          message(sprintf("Firth joint failed: %s: %s", unique(.y$cancer_type_detailed), e$message))
          tibble()
        }
      )
    }) %>%
    ungroup()
}

# ── Section 1: CRC / PDAC / RCC ──────────────────────────────────────────────

crc_pdac_rcc_cancers <- c("Colorectal Adenocarcinoma", "Pancreatic Adenocarcinoma", "Renal Clear Cell Carcinoma")

crc_pdac_rcc <- body_comp_w_pi3k %>%
  filter(cancer_type_detailed %in% crc_pdac_rcc_cancers)

marginal_l3_crc_pdac_rcc <- run_firth_marginal(crc_pdac_rcc, metrics = l3_metrics) %>%
  mutate(padj = p.adjust(p.value, method = "BH"))
marginal_organ_crc_pdac_rcc <- run_firth_marginal(crc_pdac_rcc, metrics = organ_metrics) %>%
  mutate(padj = p.adjust(p.value, method = "BH"))

joint_l3_crc_pdac_rcc <- run_firth_joint(crc_pdac_rcc, l3_metrics) %>%
  mutate(padj = p.adjust(p.value, method = "BH"))

pi3k_results_crc_pdac_rcc <- bind_rows(
  marginal_l3_crc_pdac_rcc %>% mutate(model_type = "marginal", level = "L3"),
  marginal_organ_crc_pdac_rcc %>% mutate(model_type = "marginal", level = "organ"),
  joint_l3_crc_pdac_rcc %>% mutate(model_type = "joint", level = "L3")
)

# ── Write outputs ─────────────────────────────────────────────────────────────

create_dir_if_needed("./code_outputs/results/tables/body_comp/pi3k")
write_tsv(pi3k_results_crc_pdac_rcc %>% arrange(p.value), "./code_outputs/results/tables/body_comp/pi3k/pi3k_bodycomp_firth_crc_pdac_rcc.tsv")
