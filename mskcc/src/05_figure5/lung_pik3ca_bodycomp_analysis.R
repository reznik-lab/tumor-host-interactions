# ============================================================================
# Produces: Figure 5; Supplementary Figure 5 (G, H); Supplementary Table 24
# (a figure listed without a panel letter contributes to a composite
#  panel whose letter could not be attributed unambiguously)
# ============================================================================
################
# PI3K/AKT/mTOR Pathway and Body Composition in Lung Adenocarcinoma
#
# Split out from pi3k_bodycomp_analysis.R (which covers CRC/PDAC/RCC).
#
# The PI3K/AKT/mTOR pathway is defined INCLUDING STK11 here, matching the full
# oncogenic-pathway definition used in the pan-cancer / BMI analyses
# (lung_pik3ca_bmi_analysis.R) so that the BMI and body-composition halves of the
# manuscript paragraph share one consistent pathway definition. Because STK11 loss
# is a common, lineage-specific lung event independently linked to cachexia, the
# robustness of the main association to STK11 is assessed in a dedicated
# STK11-excluded sensitivity analysis, lung_pik3ca_bodycomp_stk11_exclusion_sensitivity.R.
# A sample is "pathway altered" if it carries an oncogenic mutation in any
# PI3K/AKT/mTOR gene (STK11 included).
#
# @author: Ethan Tse
# @date: June 26th, 2026
################

source("./src/00_shared/load_libraries.R")
source("./src/00_shared/obesity_constants.R")
# rename to avoid collision with the same-named binary column in the data
pi3k_genes <- pi3k_akt_mtor_pathway
# Main pathway definition includes STK11 (see header)
pi3k_pathway_genes <- pi3k_genes
library(ggtext)

# Markdown facet labels: italic gene names with superscript mutation status,
# rendered via element_markdown(). The underlying group strings are kept plain
# (they also key the grouping / results tables); this only relabels the strips.
pi3k_facet_labeller <- as_labeller(c(
  "KRAS/EGFR WT" = "*KRAS*<sup>WT</sup>/*EGFR*<sup>WT</sup>",
  "EGFR-mut"     = "*EGFR*<sup>mut</sup>",
  "KRAS-mut"     = "*KRAS*<sup>mut</sup>"
))

pathway_data <- read_tsv("./code_outputs/stable_data/msk/impact/pathways/sample_level_mutation_w_pathway_annotations.tsv")
bodycomp_data <- read_tsv("./code_outputs/stable_data/msk/body_comp/metadata/non_clean_master_body_comp_annotated_metadata.tsv")
mutation_mtx <- read_tsv("./code_outputs/stable_data/msk/impact/snvs/matrix/panel_corrected/bmi_oncokb_annotated_first_sample_oncogenic_mutation_mtx.tsv") %>%
  select(SAMPLE_ID, KRAS, EGFR) %>%
  janitor::clean_names()
mutation_data <- read_tsv("./code_outputs/stable_data/msk/impact/snvs/dataframes/bmi_oncokb_annotated_first_sample_oncogenic_mutations.tsv")

# NOTE: the full-cohort BMI analyses (preliminary BMI boxplot + adjusted BMI
# forest) live in lung_pik3ca_bmi_analysis.R — BMI is recorded cohort-wide, so
# those run on the full mutation matrix rather than the CT-scan subset used here.

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
    matches("^l[13]"),
    matches("^t1[02]"),
    any_of(organ_metrics)
  ) %>%
  # join keeps only IMPACT-sequenced samples (pathway annotation present); the
  # pathway alteration value itself is recomputed below from pi3k_pathway_genes
  # (full pathway, STK11 included).
  left_join(pathway_data %>% select(sample_id, pi3k_akt_mtor_pathway), by = "sample_id") %>%
  left_join(mutation_mtx, by = "sample_id") %>%
  filter(!is.na(pi3k_akt_mtor_pathway))

# ── Firth helpers ─────────────────────────────────────────────────────────────

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

# Firth logistic of pathway alteration ~ each (scaled) metric, one model per
# metric. Covariates are built per group, dropping any that collapse to a single
# level (e.g. clean_scan in the clean-only cohort, gender in single-sex subsets).
# adjust_bmi = TRUE adds bmi_at_scan as a covariate, except when bmi_at_scan is
# itself the metric being tested (which would be collinear).
run_firth_marginal <- function(df, metrics = body_comp_metrics, adjust_bmi = FALSE) {
  available <- intersect(metrics, colnames(df))
  # typed empty result so columns (incl. p.value) survive even when every group is
  # skipped — otherwise downstream p.adjust(p.value) errors with "object not found"
  empty_result <- tibble(
    term = character(), estimate = numeric(), conf.low = numeric(),
    conf.high = numeric(), std.error = numeric(), p.value = numeric(),
    or = numeric(), or.conf.low = numeric(), or.conf.high = numeric(),
    n_samples = integer(), n_altered = integer(), n_effective = integer()
  )
  df %>%
    pivot_longer(cols = all_of(available), names_to = "var", values_to = "value") %>%
    group_by(cancer_type_detailed, var) %>%
    mutate(value = as.numeric(scale(value))) %>%
    group_modify(~ {
      n_samples <- nrow(.x)
      n_altered <- sum(.x$pi3k_akt_mtor_pathway == 1, na.rm = TRUE)
      if (n_altered < 10) {
        return(empty_result)
      }
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
      tryCatch(
        {
          model <- logistf::logistf(
            as.formula(paste("pi3k_akt_mtor_pathway ~ value +", paste(covars, collapse = " + "))),
            data = .x
          )
          firth_tidy(model, "value") %>% mutate(n_samples = n_samples, n_altered = n_altered, n_effective = n_effective)
        },
        error = function(e) {
          message(sprintf("Firth marginal failed: %s / %s: %s", unique(.x$cancer_type_detailed), unique(.x$var), e$message))
          empty_result
        }
      )
    }) %>%
    ungroup()
}

# ── Pathway alteration status (STK11 included) & LUAD driver subtypes ──────────
# Mutually exclusive: EGFR+ | KRAS+ (not EGFR) | WT (neither). A sample is
# pathway-altered if it has an oncogenic mutation in any PI3K/AKT/mTOR gene
# (STK11 included).

pi3k_altered_samples <- mutation_data %>%
  filter(tolower(Hugo_Symbol) %in% pi3k_pathway_genes) %>%
  distinct(SAMPLE_ID) %>%
  pull(SAMPLE_ID)

luad <- body_comp_w_pi3k %>%
  filter(cancer_type_detailed == "Lung Adenocarcinoma") %>%
  mutate(
    pi3k_akt_mtor_pathway = if_else(sample_id %in% pi3k_altered_samples, 1, 0),
    driver_subtype = case_when(
      egfr == 1 ~ "EGFR+",
      kras == 1 ~ "KRAS+",
      egfr == 0 & kras == 0 ~ "WT"
    )
  ) %>%
  filter(!is.na(driver_subtype))

# ── Label the mutually-exclusive driver subtypes ──────────────────────────────
# Each plot below contrasts the LUAD driver subtypes (KRAS/EGFR WT, EGFR-mut,
# KRAS-mut) — the pooled "All" cohort is not shown. Modeled per subtype so
# coefficients are directly comparable within a forest plot.

luad_grouped <- luad %>%
  mutate(cancer_type_detailed = case_when(
    driver_subtype == "WT" ~ "KRAS/EGFR WT",
    driver_subtype == "EGFR+" ~ "EGFR-mut",
    driver_subtype == "KRAS+" ~ "KRAS-mut"
  ))
# Arm 4 cohort: the high-quality ("clean") scan subset only.
luad_grouped_clean <- luad_grouped %>% filter(clean_scan == TRUE)

# Shared forest-plot builder: facets the driver subtypes (KRAS/EGFR WT, EGFR-mut,
# KRAS-mut), one row per metric (SKM / VAT / SAT), points colored by nominal
# significance. BMI is not shown here (see lung_pik3ca_bmi_analysis.R).
make_pi3k_forest <- function(res, xlab = "Adj. log(OR) per SD") {
  res %>%
    # organ volumes are still fit in the association (kept in the results tables),
    # but are not shown in the forest — only BMI + tissue compartments are plotted.
    filter(!var %in% organ_metrics) %>%
    mutate(
      metric_label = case_when(
        grepl("total_muscle_area", var) ~ "SKM",
        grepl("vat_area", var) ~ "VAT",
        grepl("sat_area", var) ~ "SAT"
      ),
      # reverse so the rows read top-to-bottom SKM -> VAT -> SAT (y pos 1 = bottom)
      metric_label = factor(metric_label, levels = rev(c("SKM", "VAT", "SAT"))),
      group = factor(cancer_type_detailed, levels = c("KRAS/EGFR WT", "EGFR-mut", "KRAS-mut")),
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
    scale_colour_manual(
      name = NULL,
      # red matches Figure 5D's significance colour ("#ad3f3f")
      values = SIG_COLOURS,
      drop = FALSE
    ) +
    significance_key(size = 2) +
    # drop = FALSE keeps all four panels (incl. empty ones). A driver subtype with
    # < 10 pathway-altered samples is skipped in run_firth_marginal (notably LUAD
    # KRAS+ in the clean-scan cohort), so its panel renders empty rather than
    # disappearing — keeping panel layout consistent across Plots 2-3.
    facet_wrap(~group, nrow = 1, drop = FALSE, labeller = pi3k_facet_labeller) +
    labs(x = xlab, y = NULL, fill = NULL) +
    theme_std() +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks.y = element_blank(),
      strip.text = element_markdown(size = 6),
      panel.spacing.x = unit(1, "pt")
    )
}

# ── Plot 2: L3 body composition + organ volumes ~ pathway (scan) ───────────────
# Marginal Firth models, NOT adjusted for BMI. L3 SKM/VAT/SAT and the organ
# volumes each tested separately so all rows are directly comparable. BMI is not
# shown here — it is plotted at full-cohort power in lung_pik3ca_bmi_analysis.R.
results_l3_unadj <- run_firth_marginal(
  luad_grouped,
  metrics = c(l3_metrics, organ_metrics), adjust_bmi = FALSE
) %>%
  mutate(padj = p.adjust(p.value, method = "BH"))

pi3k_forest_l3_unadj <- make_pi3k_forest(
  results_l3_unadj,
  xlab = "Adj. log(OR) per SD"
)

# ── Plot 3: L3 body composition ~ pathway, adjusted for bmi_at_scan ────────────
# Same scan cohort, but each L3 metric is additionally adjusted for BMI at scan
# (isolating body-composition signal independent of overall obesity). BMI itself
# is dropped here since it cannot be adjusted for itself.
results_l3_bmiadj <- run_firth_marginal(
  luad_grouped,
  metrics = c(l3_metrics, organ_metrics), adjust_bmi = TRUE
) %>%
  mutate(padj = p.adjust(p.value, method = "BH"))

pi3k_forest_l3_bmiadj <- make_pi3k_forest(
  results_l3_bmiadj,
  xlab = "Adj. (incl. BMI) log(OR) per SD"
)

# ── Arm 4: L3 body composition ~ pathway (clean scans only) ───────────────────
# Identical to Plot 2 (unadjusted) but restricted to clean_scan == TRUE scans.
# Table-only: this arm has no panel of its own, it is one row block of
# pi3k_bodycomp_firth_luad_subtypes.tsv.
results_l3_clean <- run_firth_marginal(
  luad_grouped_clean,
  metrics = c(l3_metrics, organ_metrics), adjust_bmi = FALSE
) %>%
  mutate(padj = p.adjust(p.value, method = "BH"))

# ── Arm 4b: clean scans, additionally adjusted for BMI at scan ───────────────
# Table-only, as Arm 4. The fourth cell of the {BMI adjustment} x
# {treatment-naive} grid. Plots 2-3 and Arm 4 fit three of the four -- all
# scans unadjusted, all scans + BMI, clean scans unadjusted -- so "independent
# of overall body size" and "holds in treatment-naive scans" were each only
# ever established one at a time. The
# STK11-excluded sensitivity script
# (lung_pik3ca_bodycomp_stk11_exclusion_sensitivity.R) already fits all four on
# its own pathway definition; this closes the same grid on the main,
# STK11-inclusive definition.
#
# READ WITH THE SAME CARE AS THE OTHER clean_scan ARM. The treatment-naive LUAD
# scan cohort is small (see n_samples/n_altered in the output table) and this
# arm spends one more parameter than Arm 4 on bmi_at_scan, so its intervals
# are the widest of the four. It is here to show the point estimate does not
# move, not to add power.
results_l3_clean_bmiadj <- run_firth_marginal(
  luad_grouped_clean,
  metrics = c(l3_metrics, organ_metrics), adjust_bmi = TRUE
) %>%
  mutate(padj = p.adjust(p.value, method = "BH"))

# Combined results table: all four scan-cohort model arms.
pi3k_results_luad <- bind_rows(
  results_l3_unadj %>% mutate(model = "unadjusted"),
  results_l3_bmiadj %>% mutate(model = "bmi_adjusted"),
  results_l3_clean %>% mutate(model = "clean_scan"),
  results_l3_clean_bmiadj %>% mutate(model = "clean_scan_bmi_adjusted")
) %>% arrange(model, cancer_type_detailed, p.value)

# ── Write outputs ─────────────────────────────────────────────────────────────
# (Full-cohort BMI outputs are written by lung_pik3ca_bmi_analysis.R.)

create_dir_if_needed("./code_outputs/results/tables/body_comp/pi3k")
write_tsv(pi3k_results_luad, "./code_outputs/results/tables/body_comp/pi3k/pi3k_bodycomp_firth_luad_subtypes.tsv")

create_dir_if_needed("./code_outputs/results/figs/body_comp/pi3k")
# Plot 2: L3 body comp + organ volumes ~ pathway, scan cohort, unadjusted, by subtype
ggsave_wrapper(
  plot = pi3k_forest_l3_unadj,
  filename = "./code_outputs/results/figs/body_comp/pi3k/pi3k_forest_luad_l3_unadjusted.pdf",
  width = 3,
  height = 2
)
# Plot 3: L3 body comp + organ volumes ~ pathway, scan cohort, adjusted for BMI, by subtype
ggsave_wrapper(
  plot = pi3k_forest_l3_bmiadj,
  filename = "./code_outputs/results/figs/body_comp/pi3k/pi3k_forest_luad_l3_bmi_adjusted.pdf",
  width = 3,
  height = 2
)

