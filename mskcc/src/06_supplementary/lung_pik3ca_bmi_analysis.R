# ============================================================================
# Produces: Supplementary Figure 5; Supplementary Table 24
# (a figure listed without a panel letter contributes to a composite
#  panel whose letter could not be attributed unambiguously)
# ============================================================================
################
# PI3K/AKT/mTOR Pathway and BMI in Lung Adenocarcinoma (full IMPACT cohort)
#
# BMI is recorded for the entire IMPACT LUAD cohort (~8-9k patients), not just the
# small CT-scan subset with body-composition metrics, so all BMI analyses run here
# off the full mutation matrix at full power. The body-composition (scan cohort)
# analyses live in lung_pik3ca_bodycomp_analysis.R.
#
# As in the body-comp script, the PI3K/AKT/mTOR pathway is defined INCLUDING STK11,
# so that the BMI and body-composition halves of the manuscript paragraph share one
# consistent pathway definition. Robustness to STK11 — a common, lineage-specific
# lung event independently linked to cachexia — is assessed by a dedicated
# STK11-excluded sensitivity analysis (lung_pik3ca_bmi_stk11_exclusion_sensitivity.R).
# A sample is "pathway altered" if it carries an oncogenic mutation in any
# PI3K/AKT/mTOR gene (STK11 included).
#
# Produces an adjusted multivariable Firth forest of pathway alteration ~ BMI per SD.
#
# @author: Ethan Tse
# @date: July 1st, 2026
################

source("./src/00_shared/load_libraries.R")
library(ggtext)
source("./src/00_shared/obesity_constants.R")
# rename to avoid collision with the same-named binary column in the data
pi3k_genes <- pi3k_akt_mtor_pathway
# Main pathway definition includes STK11 (see header)
pi3k_pathway_genes <- pi3k_genes

# Full IMPACT cohort matrix (one row per first sample). Pathway alteration
# (STK11 included) is recomputed from the per-gene oncogenic-mutation columns.
mutation_mtx_full <- read_tsv("./code_outputs/stable_data/msk/impact/snvs/matrix/panel_corrected/bmi_oncokb_annotated_first_sample_oncogenic_mutation_mtx.tsv") %>%
  select(
    SAMPLE_ID, KRAS, EGFR, CANCER_TYPE_DETAILED, BMI_NEAREST_RESULT,
    ANCESTRY_LABEL, STAGE_CDM_DERIVED, GENDER, CURRENT_AGE_DEID,
    SMOKING_PREDICTIONS_3_CLASSES, GENE_PANEL,
    any_of(toupper(pi3k_pathway_genes))
  ) %>%
  janitor::clean_names() %>%
  mutate(pi3k_akt_mtor_pathway = as.integer(
    rowSums(across(any_of(pi3k_pathway_genes)), na.rm = TRUE) > 0
  ))

# ── Firth helper ──────────────────────────────────────────────────────────────

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

# ── LUAD cohort with mutually-exclusive driver subtypes ───────────────────────
# EGFR+ | KRAS+ (not EGFR) | WT (neither). BMI must be recorded.

luad_bmi <- mutation_mtx_full %>%
  filter(cancer_type_detailed == "Lung Adenocarcinoma") %>%
  mutate(driver_subtype = case_when(
    egfr == 1 ~ "EGFR+",
    kras == 1 ~ "KRAS+",
    egfr == 1 & kras == 1 ~ "Double",
    egfr == 0 & kras == 0 ~ "WT"
  )) %>%
  filter(!is.na(driver_subtype), !is.na(bmi_nearest_result))

# ── Adjusted multivariable Firth forest, pathway ~ BMI per SD ─────────────────
# Firth logistic (pathway alteration as the outcome, BMI per SD as the predictor),
# adjusted for age, ancestry, stage, smoking (+ gender where it varies), run for
# all LUAD then within each driver subset.

luad_bmi_grouped <- bind_rows(
  luad_bmi %>% mutate(group = paste("LUAD", driver_subtype))
)

run_bmi_firth <- function(df) {
  df %>%
    group_by(group) %>%
    group_modify(~ {
      n_samples <- nrow(.x)
      n_altered <- sum(.x$pi3k_akt_mtor_pathway == 1, na.rm = TRUE)
      if (n_altered < 10) {
        return(tibble())
      }
      .x <- .x %>% mutate(bmi_z = as.numeric(scale(bmi_nearest_result)))
      # stage_cdm_derived (stage at diagnosis) is the study-wide stage covariate;
      # it matches the companion body-composition analysis in
      # lung_pik3ca_bodycomp_analysis.R, which tests the same driver subtypes.
      covariates <- if (n_distinct(.x$gender) > 1) {
        "current_age_deid + ancestry_label + stage_cdm_derived + smoking_predictions_3_classes + gender"
      } else {
        "current_age_deid + ancestry_label + stage_cdm_derived + smoking_predictions_3_classes"
      }
      tryCatch(
        {
          model <- logistf::logistf(
            as.formula(paste("pi3k_akt_mtor_pathway ~ bmi_z +", covariates)),
            data = .x
          )
          firth_tidy(model, "bmi_z") %>%
            mutate(n_samples = n_samples, n_altered = n_altered)
        },
        error = function(e) {
          message(sprintf("BMI Firth failed: %s: %s", unique(.x$group), e$message))
          tibble()
        }
      )
    }) %>%
    ungroup()
}

pi3k_bmi_luad <- run_bmi_firth(luad_bmi_grouped) %>%
  mutate(padj = p.adjust(p.value, method = "BH"))

# Laid out as a single "BMI" row faceted by driver subtype (KRAS/EGFR WT,
# EGFR-mut, KRAS-mut) so it aligns with the scan-cohort body-composition forests
# (SKM/VAT/SAT) in lung_pik3ca_bodycomp_analysis.R.
pi3k_bmi_forest_luad <- pi3k_bmi_luad %>%
  mutate(
    group = factor(case_when(
      group == "LUAD WT" ~ "*KRAS*<sup>WT</sup>/*EGFR*<sup>WT</sup>",
      group == "LUAD EGFR+" ~ "*EGFR*<sup>mut</sup>",
      group == "LUAD KRAS+" ~ "*KRAS*<sup>mut</sup>"
    ), levels = c("*KRAS*<sup>WT</sup>/*EGFR*<sup>WT</sup>", "*EGFR*<sup>mut</sup>", "*KRAS*<sup>mut</sup>")),
    metric_label = factor("BMI", levels = "BMI"),
    significance = sig_level(p.value, padj)
  ) %>%
  ggplot(aes(x = estimate, y = metric_label)) +
  ggstats::geom_stripped_rows() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1) +
  geom_linerange(
    aes(
      xmin = conf.low, xmax = conf.high,
      color = significance
    ),
    lineend = "round", linewidth = 0.1, position = position_dodge(width = 0.3)
  ) +
  geom_point(size = 2, shape = 16, aes(colour = significance), position = position_dodge(width = 0.3)) +
  # Dark grey / red -- red matches Figure 5D's significance colour
  # ("#ad3f3f"), applied here to highlight (red) vs non-significant (grey) points.
  scale_colour_manual(
    name = NULL,
    values = SIG_COLOURS,
    drop = FALSE
  ) +
  significance_key(size = 2) +
  facet_wrap(~group, nrow = 1, drop = FALSE) +
  labs(x = "Adj. log(OR) per SD BMI", y = NULL, fill = NULL) +
  theme_std() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.y = element_blank(),
    strip.text = element_markdown(size = 6),
    panel.spacing.x = unit(1, "pt")
  )

# ── Write outputs ─────────────────────────────────────────────────────────────

create_dir_if_needed("./code_outputs/results/tables/body_comp/pi3k")
write_tsv(
  pi3k_bmi_luad %>% arrange(p.value),
  "./code_outputs/results/tables/body_comp/pi3k/pi3k_bmi_firth_luad_full_cohort.tsv"
)

create_dir_if_needed("./code_outputs/results/figs/body_comp/pi3k")
# Adjusted BMI forest, full LUAD cohort, single BMI row faceted by subtype
ggsave_wrapper(
  plot = pi3k_bmi_forest_luad,
  filename = "./code_outputs/results/figs/body_comp/pi3k/pi3k_bmi_forest_luad_full_cohort.pdf",
  width = 3,
  height = 1.3
)
