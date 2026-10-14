#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 01.1 -- build the PSM matching table(s) for cancer individuals.
##
## Produces no manuscript figure or table: builds the case arm's matching table.
##
## Reads:
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/aou_enrollment_date.tsv
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/raw_concatenated/person.tsv
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/raw_concatenated/zipcodeSocioeconomic.tsv
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/raw_concatenated/measurementOccurrence.tsv
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/first_cancer_diagnosis_dates_grouped.tsv
## Writes:
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/pre_psm_matched_cancer_bmi_summary_prediagnosis.tsv
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/pre_psm_matched_cancer_bmi_summary_anytime.tsv
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


## ---------------------------------------------------------------------------
## 1. Load core reference tables
## ---------------------------------------------------------------------------
aou_enrollment_date <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/aou_enrollment_date.tsv"))

cancer_demographics <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/raw_concatenated/person.tsv")) %>%
  janitor::clean_names()

# Socioeconomic status: composite deprivation_index
cancer_ses <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/raw_concatenated/zipcodeSocioeconomic.tsv")) %>%
  janitor::clean_names() %>%
  select(person_id, deprivation_index)


## ---------------------------------------------------------------------------
## 2. First cancer diagnosis (any type) per person
## ---------------------------------------------------------------------------
# this df has been filtered to the first diagnosis, as well as grouped to tissue groups
first_dx_date <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/first_cancer_diagnosis_dates_grouped.tsv")) %>%
  janitor::clean_names()


## ---------------------------------------------------------------------------
## 3. Load and clean raw BMI measurements
## ---------------------------------------------------------------------------
cancer_bmi_measurements_raw <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/raw_concatenated/measurementOccurrence.tsv"),
  progress = TRUE) %>%
  janitor::clean_names()

## filter_usable_bmi() (src/setup.R) replaces the provenance exclusion that used
## to sit here. It dropped self-report / unmapped / no-information rows -- 8.5%
## of the table -- that every later stage kept, so this covariate counted a
## different set of measurements than the analyses ran on. See the note there.
cancer_bmi_measurements <- cancer_bmi_measurements_raw %>%
  filter_usable_bmi() %>%
  distinct(person_id, measurement_datetime, .keep_all = TRUE) %>%   # one row per person-datetime
  mutate(measurement_date = as.Date(measurement_datetime)) %>%
  left_join(first_dx_date %>% select(person_id, first_dx_datetime), by = "person_id")


## ---------------------------------------------------------------------------
## 4. Shared helper: count / gap / span summary from a (pre-filtered) long
##    BMI table
## ---------------------------------------------------------------------------
summarise_bmi_history <- function(bmi_df) {
  bmi_df %>%
    arrange(person_id, measurement_date) %>%
    group_by(person_id) %>%
    mutate(
      gap_days = as.numeric(measurement_date - lag(measurement_date))  # NA for each person's first measurement
    ) %>%
    summarise(
      n_bmi_measurements = n_distinct(measurement_date),
      earliest_bmi_year  = year(min(measurement_date, na.rm = TRUE)),
      mean_gap_days      = mean(gap_days, na.rm = TRUE),   # NaN if only 1 measurement
      overall_span_days  = as.numeric(max(measurement_date) - min(measurement_date)),
      .groups = "drop"
    ) %>%
    mutate(
      mean_gap_years     = mean_gap_days / 365.25,
      overall_span_years = overall_span_days / 365.25
    )
}


## ---------------------------------------------------------------------------
## 5a. PRE-DIAGNOSIS table: only measurements strictly before first_dx_datetime
## ---------------------------------------------------------------------------
bmi_history_prediagnosis <- cancer_bmi_measurements %>%
  filter(!is.na(first_dx_datetime),
         measurement_datetime < first_dx_datetime) %>%
  summarise_bmi_history()

## ---------------------------------------------------------------------------
## 5b. ANY-TIME table: all measurements, no diagnosis filter
## ---------------------------------------------------------------------------
bmi_history_anytime <- cancer_bmi_measurements %>%
  summarise_bmi_history()


## ---------------------------------------------------------------------------
## 7. Shared demographic/SES/enrollment base (identical for both output
##    tables -- only the BMI-history columns differ)
## ---------------------------------------------------------------------------
cancer_base <- cancer_demographics %>%
  left_join(cancer_ses,          by = "person_id") %>%
  left_join(first_dx_date,       by = "person_id") %>%
  left_join(aou_enrollment_date, by = "person_id") %>%
  mutate(
    age_at_enrollment = as.numeric(
      difftime(as.Date(primary_consent_date), as.Date(date_of_birth), units = "days")
    ) / 365.25,
    case = 1
  )


## ---------------------------------------------------------------------------
## 8. Final output tables
##    (join BMI-history metrics; each is later filtered to
##    n_bmi_measurements >= 4 at whichever scope is appropriate -- pre-diagnosis
##    for table 1, anytime for table 2 -- as a separate explicit step, not baked
##    in here, so you can see drop counts before committing to the filter)
## ---------------------------------------------------------------------------
cancer_bmi_summary_prediagnosis <- cancer_base %>%
  inner_join(bmi_history_prediagnosis, by = "person_id")   # inner: drops anyone with 0 pre-diagnosis measurements

cancer_bmi_summary_anytime <- cancer_base %>%
  inner_join(bmi_history_anytime, by = "person_id")        # inner: drops anyone with 0 measurements at all


## ---------------------------------------------------------------------------
## 9. Sanity checks
## ---------------------------------------------------------------------------
message(sprintf("\nPre-diagnosis table: %d individuals (of %d in cancer_base).",
                nrow(cancer_bmi_summary_prediagnosis), nrow(cancer_base)))
message(sprintf("Any-time table:       %d individuals (of %d in cancer_base).",
                nrow(cancer_bmi_summary_anytime), nrow(cancer_base)))

message(sprintf("\nPre-diagnosis: %d individuals (%.1f%%) have >= 4 pre-diagnosis BMI measurements.",
                sum(cancer_bmi_summary_prediagnosis$n_bmi_measurements >= 4),
                100 * mean(cancer_bmi_summary_prediagnosis$n_bmi_measurements >= 4)))
message(sprintf("Any-time:      %d individuals (%.1f%%) have >= 4 BMI measurements at any time.",
                sum(cancer_bmi_summary_anytime$n_bmi_measurements >= 4),
                100 * mean(cancer_bmi_summary_anytime$n_bmi_measurements >= 4)))


## --- write -----------------------------------------------------------------

prediagnosis_path <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/pre_psm_matched_cancer_bmi_summary_prediagnosis.tsv")

anytime_path <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/pre_psm_matched_cancer_bmi_summary_anytime.tsv")

dir.create(dirname(prediagnosis_path), recursive = TRUE, showWarnings = FALSE)

write_tsv(cancer_bmi_summary_prediagnosis %>% janitor::clean_names(), prediagnosis_path)
write_tsv(cancer_bmi_summary_anytime      %>% janitor::clean_names(), anytime_path)
