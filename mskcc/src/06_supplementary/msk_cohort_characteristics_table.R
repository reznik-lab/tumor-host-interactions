# ============================================================================
# Produces: Supplementary Table 8
# ============================================================================
#####################
# Supplementary Table 8: patient characteristics (age, sex, ancestry, smoking,
# stage, prior medications at MSK) of the whole solid-tumor MSK-IMPACT cohort.
#
# Cohort: solid-tumor IMPACT panels only, no further restriction -- N = 78,213,
# every one of whom has a BMI. This is the cohort the manuscript opens with.
# The >= 400-patients-per-cancer-type restriction used elsewhere belongs to the
# Figure 1g-h forest (oac_vs_noac_lm_msk.R), not to a description of the cohort,
# so it is deliberately not applied here.
#
# @author: Ethan Tse
# @date: 2026-08-05
#####################

source("./src/00_shared/load_libraries.R")
library(gtsummary)

full_impact_cohort <- read_tsv("./code_outputs/stable_data/msk/metadata/impact_patient_and_sample_metadata_first_surgery.tsv") %>%
  janitor::clean_names() %>%
  select(
    patient_id,
    tumor_code = oncotree_code,
    prior_med_to_msk,
    current_age_deid,
    gender,
    ancestry_label,
    stage_highest_recorded,
    smoking_predictions_3_classes,
    gene_panel
  ) %>%
  mutate(tumor_code = if_else(tumor_code %in% c("COAD", "READ", "COADREAD"), "COADREAD", tumor_code)) %>%
  filter(gene_panel %in% c("IMPACT341", "IMPACT410", "IMPACT468", "IMPACT505"))

stopifnot(nrow(full_impact_cohort) == 78213)
cat("N patients (full solid-tumor IMPACT cohort):", nrow(full_impact_cohort), "\n")

full_summary_tbl <- full_impact_cohort %>%
  select(
    current_age_deid,
    gender,
    ancestry_label,
    smoking_predictions_3_classes,
    stage_highest_recorded,
    prior_med_to_msk
  ) %>%
  tbl_summary(
    label = list(
      current_age_deid ~ "Age at MSK-IMPACT sequencing",
      gender ~ "Sex",
      ancestry_label ~ "Ancestry",
      smoking_predictions_3_classes ~ "Smoking status",
      stage_highest_recorded ~ "Stage",
      prior_med_to_msk ~ "Prior medications at MSK"
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing_text = "Unknown"
  ) %>%
  bold_labels()

out_dir <- "./code_outputs/results/tables/bmi/impact/msk_cohort_characteristics"
create_dir_if_needed(out_dir)

# Raw label + stat_* columns, matching the AoU table1_summary.tsv convention
# (rows with all-NA stat_* are variable-label header rows; blank/indented
# rows below them are the category levels).
write_tsv(
  full_summary_tbl$table_body %>% select(label, starts_with("stat_")),
  file.path(out_dir, "supp_table_full_impact_cohort_summary.tsv")
)
write_tsv(
  tibble(cohort = "Full solid-tumor MSK-IMPACT cohort", n = nrow(full_impact_cohort)),
  file.path(out_dir, "supp_table_full_impact_cohort_n.tsv")
)
