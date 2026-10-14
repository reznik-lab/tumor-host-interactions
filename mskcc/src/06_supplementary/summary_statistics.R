# ============================================================================
# Produces: Supplementary Table 19
# ============================================================================
#######################
# Summary statistics (age, sex, stage, smoking status, ancestry) of the cohort, stratified by viral status. This is to be used in the supplement.
# What fraction are obese in each group? What fraction are smokers in each group? What is the age distribution in each group?
#
#
#
# @author: Ethan Tse
# @date: May 6th, 2026
#######################

source("./src/00_shared/load_libraries.R")

# Load the data
mutation_metadata <- read_tsv("./code_outputs/stable_data/msk/impact/snvs/dataframes/bmi_oncokb_annotated_first_sample_oncogenic_mutations.tsv") %>%
  filter(CANCER_TYPE_DETAILED == "Merkel Cell Carcinoma") %>%
  janitor::clean_names() %>%
  select(sample_id, bmi_nearest_result, current_age_deid, gender, stage_highest_recorded, smoking_predictions_3_classes, prior_med_to_msk, so_comments, cvr_tmb_score, path_slide_exists, msk_slide_id) %>%
  distinct()

mutation_mtx <- read_tsv("./code_outputs/stable_data/msk/impact/snvs/matrix/panel_corrected/bmi_oncokb_annotated_first_sample_oncogenic_mutation_mtx.tsv") %>%
  filter(CANCER_TYPE_DETAILED == "Merkel Cell Carcinoma") %>%
  select(SAMPLE_ID, RB1, TP53) %>%
  janitor::clean_names()

# Body composition subcohort, as defined in src/04_figure4/bodycomp_associations.R:
# MCC patients with body comp scan data, excluding QC-failed scans
bodycomp_qc_removed <- c("P-0036928", "P-0036676", "P-0039150", "P-0054731", "P-0057490", "P-0060039", "P-0058771", "P-0068753")
bodycomp_subcohort <- read_tsv("./code_outputs/stable_data/msk/body_comp/metadata/master_body_comp_annotated_metadata.tsv") %>%
  filter(cancer_type_detailed == "Merkel Cell Carcinoma") %>%
  filter(!dmp_id %in% bodycomp_qc_removed) %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  select(sample_id, age_at_scan)

bodycomp_subcohort_sample_ids <- bodycomp_subcohort$sample_id

# joint data
mcc_df <- mutation_mtx %>%
  left_join(mutation_metadata, by = c("sample_id" = "sample_id")) %>%
  mutate(viral_status = ifelse((rb1 == 1 | tp53 == 1) & cvr_tmb_score >= 10, "non_viral", "viral")) %>%
  mutate(in_bodycomp_subcohort = sample_id %in% bodycomp_subcohort_sample_ids) %>%
  # age_at_scan is defined only for the body-composition subcohort; the
  # cohort-characteristics table (Supplementary Table 6) reports it there,
  # since that table describes the scan subcohort.
  left_join(bodycomp_subcohort, by = "sample_id")

mcc_df <- mcc_df %>%
  arrange(desc(path_slide_exists))

create_dir_if_needed("./code_outputs/results/tables/bmi/impact/merkel_cell_carcinoma")
write_tsv(mcc_df, "./code_outputs/results/tables/bmi/impact/merkel_cell_carcinoma/msk_impact_mcc_df.tsv")
