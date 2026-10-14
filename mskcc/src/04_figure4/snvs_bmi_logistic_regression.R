# ============================================================================
# Produces: Figure 4; Supplementary Table 18
# (a figure listed without a panel letter contributes to a composite
#  panel whose letter could not be attributed unambiguously)
# ============================================================================
#####################
# Logistic regression to detect gene-cancer associations with BMI
#
# This is a claude refactor of a prior script.
#
# In this, we fit a base unadjusted model, a partial (age, sex), and a
# fully adjusted model (age, sex, ancestry, stage, smoking, prior treatment),
# using glm(mut ~ var + covariates, data = adf, family = "binomial").
#
# Self-reported race was dropped from the fully adjusted model (see
# logistic_regression_helpers.R): genetic ancestry already covers it, and no
# other adjusted model in the manuscript carries it.
#
# @author: Ethan Tse
# @date: 2026-04-29
#####################

source("./src/00_shared/load_libraries.R")
source("./src/00_shared/utils.R")
source("./src/00_shared/logistic_regression_helpers.R")

start_time <- Sys.time()

# ── 1. Load data ───────────────────────────────────────────────────────────────

seq_bmi_primary_mut <- read_tsv(
  "./code_outputs/stable_data/msk/impact/snvs/dataframes/bmi_oncokb_annotated_first_sample_oncogenic_mutations.tsv"
) %>%
  filter(!is.na(ANCESTRY_LABEL)) %>%
  mutate(Tumor_Sample_Barcode = SAMPLE_ID) %>%
  filter(GENE_PANEL %in% c("IMPACT410", "IMPACT468", "IMPACT505", "IMPACT341")) %>%
  mutate(CANCER_TYPE_DETAILED = ifelse(CANCER_TYPE_DETAILED %in% c("Colorectal Adenocarcinoma", "Colon Adenocarcinoma", "Rectal Adenocarcinoma"), "Colorectal Adenocarcinoma", CANCER_TYPE_DETAILED))

# filter(grepl("Heme", GENE_PANEL, ignore.case = TRUE)) # for heme-only analysis
length(unique(seq_bmi_primary_mut$Tumor_Sample_Barcode))
mut_mtx <- read_tsv(
  "./code_outputs/stable_data/msk/impact/snvs/matrix/panel_corrected/bmi_oncokb_annotated_first_sample_oncogenic_mutation_mtx.tsv"
) %>%
  dplyr::rename(Tumor_Sample_Barcode = SAMPLE_ID) %>%
  filter(Tumor_Sample_Barcode %in% seq_bmi_primary_mut$Tumor_Sample_Barcode)

# what's the cancer sample size
cancers_n_geq_50 <- mut_mtx %>%
  group_by(CANCER_TYPE_DETAILED) %>%
  summarise(n = n()) %>%
  arrange(n) %>%
  filter(n >= 100)

seq_bmi_primary_mut <- seq_bmi_primary_mut %>%
  filter(CANCER_TYPE_DETAILED %in% cancers_n_geq_50$CANCER_TYPE_DETAILED)

# BMI is entered per standard deviation rather than per unit, standardized
# WITHIN each cancer type -- the same specification the pathway-level volcano
# uses (scale(bmi_nearest_result) inside each cancer type's own fit), so the two
# volcanoes are on one scale. Standardizing here rather than inside
# impact_logistic_regression() leaves every other analysis sharing that helper
# on its existing scale. Rescaling a single linear term does not move P-values
# or significance -- only the estimate and its interval change units.
seq_bmi_primary_mut <- seq_bmi_primary_mut %>%
  group_by(CANCER_TYPE_DETAILED) %>%
  mutate(BMI_NEAREST_RESULT_Z = as.numeric(scale(BMI_NEAREST_RESULT))) %>%
  ungroup()

# ── 2. Gene-panel lookup ───────────────────────────────────────────────────────

genes_in_panels <- read_tsv(
  "./code_outputs/stable_data/msk/impact/panels/synonymous_genes_merged/impact_panels_new.tsv"
) %>%
  dplyr::rename(
    "impact-heme-400" = "hemepact_400",
    "impact-heme-468" = "impact_heme_468",
    "access-heme-115" = "access_heme_115"
  ) %>%
  mutate(panels = apply(
    ., 1,
    function(row) list(names(row)[!is.na(row) & !names(row) %in% c("gene_name", "in_panels")])
  ))

# ── 3. Derived metadata ────────────────────────────────────────────────────────

cancer_types <- unique(seq_bmi_primary_mut$CANCER_TYPE_DETAILED)

single_sex_cancers <- seq_bmi_primary_mut %>%
  group_by(CANCER_TYPE_DETAILED) %>%
  summarise(n_sexes = n_distinct(GENDER), .groups = "drop") %>%
  filter(CANCER_TYPE_DETAILED %in% cancer_types, n_sexes == 1)

# ── 4. Build joined matrix (mutations + clinical covariates) ───────────────────

regression_meta <- seq_bmi_primary_mut %>%
  select(
    Tumor_Sample_Barcode,
    age = CURRENT_AGE_DEID,
    sex = GENDER,
    ancestry = ANCESTRY_LABEL,
    # BMI_NEAREST_RESULT,
    BMI_NEAREST_RESULT_Z,
    cancer = CANCER_TYPE_DETAILED,
    stage = STAGE_CDM_DERIVED,
    smoking = SMOKING_PREDICTIONS_3_CLASSES,
    prior_treatment = PRIOR_MED_TO_MSK,
    sample_type = SAMPLE_TYPE
  ) %>%
  distinct()

cancer_data_matrix <- mut_mtx %>%
  filter(Tumor_Sample_Barcode %in% seq_bmi_primary_mut$Tumor_Sample_Barcode) %>%
  inner_join(regression_meta, by = "Tumor_Sample_Barcode") %>%
  tibble::column_to_rownames("Tumor_Sample_Barcode")

# ── 5. Run logistic regression ─────────────────────────────────────────────────

logistic_reg_results <- impact_logistic_regression(
  variable_of_interest = "BMI_NEAREST_RESULT_Z",
  mutation_mtx = cancer_data_matrix,
  seq_bmi_primary_mut = seq_bmi_primary_mut,
  single_sex_cancers = single_sex_cancers,
  gene_in_panel = genes_in_panels,
  min_mutations = 10
)

# ── 6. Save results ────────────────────────────────────────────────────────────

save_impact_logistic_regression(
  impact_logistic_regression_results_obj = logistic_reg_results,
  analysis_name = "somatic_only_ancestry_adjusted_solid"
  # analysis_name = "impact_tmb_adjusted"
)

message(
  "Done | elapsed: ",
  round(difftime(Sys.time(), start_time, units = "mins"), 2), " mins"
)
