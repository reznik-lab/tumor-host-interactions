#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 00.3 -- first primary malignancy per person.
##
## Produces no manuscript figure or table: supplies each case's first primary
## malignancy, which defines their cancer type and diagnosis date.
##
## Reads:
##   <workspace>/cdrv9/cancer_patients/aou_cancer_patients_data/bq_exports/conditionOccurrence/20260712/conditionOccurrence/conditionOccurrence_000000000000.csv
## Writes:
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/first_cancer_diagnosis_dates_grouped.tsv
## ---------------------------------------------------------------------------

## --- locations -------------------------------------------------------------
## Found from this script's own path, so it runs from any directory.
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (!startsWith(.f, "/")) .f <- file.path(Sys.getenv("PWD", getwd()), .f)
REPO <- normalizePath(file.path(dirname(.f), "..", ".."), mustWork = TRUE)

## Data root: two levels above the repo. Override with $AOU_WORKSPACE.
WORKSPACE <- normalizePath(
  Sys.getenv(
    "AOU_WORKSPACE",
    unset = file.path(REPO, "..", "..", "workspace")
  ),
  mustWork = TRUE
)

source(file.path(REPO, "src", "setup.R"))


## --- read ------------------------------------------------------------------

cancer_diagnosis_meta <- read_csv(file.path(
  WORKSPACE,
  "cdrv9/cancer_patients/aou_cancer_patients_data/bq_exports/conditionOccurrence/20260712/conditionOccurrence/conditionOccurrence_000000000000.csv"
)) %>%
  janitor::clean_names()


## --- first diagnosis per person --------------------------------------------

first_cancer <- cancer_diagnosis_meta %>%
  group_by(person_id) %>%
  arrange(condition_start_datetime) %>%
  slice(1) %>%
  ungroup()

cancer_diagnosis_meta_grouped <- first_cancer %>%
  mutate(cancer_type_grouped = assign_tissue_group(t_disp_standard_concept_name))

testing <- cancer_diagnosis_meta_grouped %>%
  select(person_id,
    cancer_type = t_disp_standard_concept_name,
    cancer_type_grouped,
    first_dx_datetime = condition_start_datetime
  )


## --- write -----------------------------------------------------------------

out_path <- file.path(
  WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/first_cancer_diagnosis_dates_grouped.tsv"
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write_tsv(testing, out_path)

message("wrote ", nrow(testing), " rows -> ", out_path)
