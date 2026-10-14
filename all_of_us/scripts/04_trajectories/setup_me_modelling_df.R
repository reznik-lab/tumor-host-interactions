#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 04.2 -- build the input table for the mixed-effects landmark models.
##
## Produces no manuscript figure or table: assembles the input frame for the
## landmark mixed-effects models.
##
## Reads:
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm_matched_usable_data/solid_heme/full_cohort/all_cancers_geq1_before_diagnosis_after_psm.tsv
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm_matched_usable_data/solid_heme/full_cohort/all_controls_geq1_before_diagnosis_after_psm.tsv
## Writes:
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/temp/mixed_effects_analysis_input.tsv
## ---------------------------------------------------------------------------

## --- locations -------------------------------------------------------------
## Found from this script's own path, so it runs from any directory.
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (!startsWith(.f, "/")) .f <- file.path(Sys.getenv("PWD", getwd()), .f)
REPO <- normalizePath(file.path(dirname(.f), "..", ".."), mustWork = TRUE)

## Data root: two levels above the repo. Override with $AOU_WORKSPACE.
WORKSPACE <- normalizePath(Sys.getenv(
  "AOU_WORKSPACE", unset = file.path(REPO, "..", "..", "workspace")),
  mustWork = TRUE)

source(file.path(REPO, "src", "setup.R"))


## --- psm data --------------------------------------------------------------

cancers_psm_cohort <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm_matched_usable_data/solid_heme/full_cohort/all_cancers_geq1_before_diagnosis_after_psm.tsv"))

controls_psm_cohort <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm_matched_usable_data/solid_heme/full_cohort/all_controls_geq1_before_diagnosis_after_psm.tsv"))

combined_psm_chort <- bind_rows(cancers_psm_cohort, controls_psm_cohort) %>%
  # for each case, assign their cancer_type to be the same as their matched case
  # as defined by subclass
  group_by(subclass) %>%
  mutate(
    cancer_type = if_else(
      case == 0,                        # or case == "control", however it's coded
      first(cancer_type[case == 1]),    # pull the case's cancer_type
      cancer_type
    )
  ) %>%
  ungroup()

print(head(combined_psm_chort))


## ---------------------------------------------------------------------------
## 1. Prep: compute days_from_dx, filter to BMI, clean range
## ---------------------------------------------------------------------------
bmi_long <- combined_psm_chort %>%
  filter(t_disp_standard_concept_name == "Body mass index (BMI) [Ratio]") %>%
  mutate(
    measurement_datetime = ymd_hms(measurement_datetime, quiet = TRUE) %>%
      coalesce(ymd(measurement_datetime, quiet = TRUE)),
    pseudo_diagnosis_datetime = ymd_hms(pseudo_diagnosis_datetime, quiet = TRUE) %>%
      coalesce(ymd(pseudo_diagnosis_datetime, quiet = TRUE)),
    days_from_dx = as.numeric(difftime(measurement_datetime, pseudo_diagnosis_datetime, units = "days")),
    bmi = value_as_number
  ) %>%
  filter_usable_bmi("bmi") %>%   # value column was renamed above
  # restrict to a reasonable pre-dx window so distant outlier measurements don't distort the smooth
  filter(days_from_dx >= -1825 & days_from_dx <= 30) %>%
  mutate(
    case      = factor(case, levels = c(0, 1), labels = c("control", "case")),
    person_id = factor(person_id)
  ) %>%
  select(person_id, days_from_dx, bmi, case, subclass, match_cancer_type) %>%
  filter(!is.na(days_from_dx), !is.na(bmi))


## ---------------------------------------------------------------------------
## 2. Only fit cancer types with enough patients AND enough repeated measures
## ---------------------------------------------------------------------------
## Counted on CASES ONLY. This gate used to be n_distinct(person_id) >= 50
## across both arms, which at 10:1 matching let a type through on ~5 cases
## (CODE_REVIEW_FINDINGS #17); it now uses the pipeline-wide floor from
## setup.R, so this frame holds the same cancer types every other stage does.
cancer_types_to_test <- bmi_long %>%
  group_by(match_cancer_type) %>%
  filter(n_distinct(person_id[case == "case"]) >= MIN_CASES_PER_TYPE) %>%
  pull(match_cancer_type) %>%
  unique()

bmi_long_filtered <- bmi_long %>%
  filter(match_cancer_type %in% cancer_types_to_test)


## --- write -----------------------------------------------------------------

out_path <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/temp/mixed_effects_analysis_input.tsv")

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write_tsv(bmi_long_filtered, out_path)

message("wrote ", nrow(bmi_long_filtered), " rows -> ", out_path)
