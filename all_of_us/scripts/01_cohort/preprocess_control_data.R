#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 01.2 -- build the PSM matching table for control individuals.
##
## Produces no manuscript figure or table: builds the control pool's matching
## table.
##
## Reads:
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/control_individuals/raw_concatenated/measurementOccurrence.tsv
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/control_individuals/raw_concatenated/person.tsv
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/control_individuals/raw_concatenated/zipcodeSocioeconomic.tsv
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/aou_enrollment_date.tsv
## Writes:
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/control_individuals/pre_psm_matched_metadata.tsv
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


## --- BMI measurements ------------------------------------------------------

controls_bmi_measurements <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/control_individuals/raw_concatenated/measurementOccurrence.tsv"),
  progress = TRUE) %>%
  janitor::clean_names()

print(head(controls_bmi_measurements))

## filter_usable_bmi() (src/setup.R) replaces the provenance exclusion that used
## to sit here -- see the note there, and the matching change in
## preprocess_cancer_data.R.
controls_bmi_measurements <- controls_bmi_measurements %>%
  filter_usable_bmi() %>%
  group_by(person_id) %>%
  summarise(
    n_bmi_measurements = n_distinct(measurement_datetime),
    earliest_bmi_year  = year(min(measurement_datetime, na.rm = TRUE))
  ) %>%
  ungroup()


## --- demographics, SES, enrollment -----------------------------------------

aou_enrollment_date <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/aou_enrollment_date.tsv"))

print(dim(aou_enrollment_date))
print(head(aou_enrollment_date))

controls_demographics <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/control_individuals/raw_concatenated/person.tsv"))

# Socioeconomic status
# for SES, use a composite metric, deprivation_index
control_ses <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/control_individuals/raw_concatenated/zipcodeSocioeconomic.tsv")) %>%
  janitor::clean_names() %>%
  select(person_id, deprivation_index)


## --- matching table --------------------------------------------------------
## Demographics, socioeconomic status, and the year of BMI measurements.
## Matching by the number of measurements is not feasible, because the matching
## of interest is on pre-diagnosis BMI measurements, which can only be computed
## with a matching cancer individual.

control_bmi_summary <- controls_bmi_measurements %>%
  left_join(controls_demographics, by = "person_id") %>%
  left_join(control_ses,           by = "person_id") %>%
  left_join(aou_enrollment_date,   by = "person_id") %>%
  mutate(
    age_at_enrollment = as.numeric(
      difftime(as.Date(primary_consent_date), as.Date(date_of_birth), units = "days")
    ) / 365.25,
    case = 0
  )


## --- write -----------------------------------------------------------------

out_path <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/control_individuals/pre_psm_matched_metadata.tsv")

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write_tsv(control_bmi_summary %>% janitor::clean_names(), out_path)

message("wrote ", nrow(control_bmi_summary), " rows -> ", out_path)
