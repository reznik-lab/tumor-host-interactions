#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 06.2 -- cohort demographics, cases vs. matched controls.
##
## Produces: Supplementary Table 1 -- cohort demographics, cases vs. matched
##           controls.
##
## Reads:
##   <workspace>/.../stable_data/psm/solid_heme/full_cohort/aggregate/psm_matched_cohort_all_cancer_types.tsv
##   <workspace>/.../stable_data/{cancer_patients,control_individuals}/ cohort tables
## Writes:
##   .../results/tables/data_overview/table1_summary.tsv (+ .xlsx)
##       [Supplementary Table 1]
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
library(gtsummary)
library(gt)


## --- load ------------------------------------------------------------------

# load PSM-matched samples; full cohort PSM was matched using all samples, and
# then post-filtered to a subset with at least one pre-diagnosis BMI measurement
matched_controls <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm_matched_usable_data/solid_heme/full_cohort/all_controls_geq1_before_diagnosis_after_psm.tsv"))

cancers <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm_matched_usable_data/solid_heme/full_cohort/all_cancers_geq1_before_diagnosis_after_psm.tsv"))


## --- restrict to the reported cancer types ---------------------------------
## Table 1 describes the cohort the manuscript analyses, so it has to be the
## same cohort: matched, post-geq1, and restricted to the cancer types
## clearing MIN_CASES_PER_TYPE. It previously covered every matched type, the
## sub-floor ones included, so its totals were larger than any figure's.
##
## ONE DEFINITION OF "the 33", read from the matched cohort rather than
## recomputed from the geq1 frames: those have smaller case counts, so gating
## on them would admit a different set. Same source supplementary_psm_quantification.R
## and Stage 05 use.
reportable_types <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm/solid_heme/full_cohort/aggregate/psm_matched_cohort_all_cancer_types.tsv")) %>%
  filter(case == 1) %>%
  distinct(person_id, match_cancer_type) %>%
  count(match_cancer_type, name = "n_matched_cases") %>%
  filter(n_matched_cases >= MIN_CASES_PER_TYPE) %>%
  pull(match_cancer_type)

stopifnot(length(reportable_types) > 0)

## A CONTROL IS KEPT IF ANY OF ITS MATCHES IS TO A REPORTED TYPE. Filtering
## rows, not people: a control matched to both a reported and a sub-floor type
## keeps its reported match and stays in the table once.
##
## IT CAN STILL LOSE MEASUREMENTS. The geq1 frames hold a match's PRE-DIAGNOSIS
## measurements, and the cutoff is that match's own pseudo-diagnosis date, so
## the rows differ between a control's matches rather than repeating. Dropping
## its sub-floor match therefore drops any measurement that was pre-diagnosis
## only under that match. That is the intended reading: count_measurements()
## below is a union over the person's retained matches, so after this filter it
## describes the person as they appear in the REPORTED cohort, which is what
## this table is for. It does mean the "BMI Measurements per Person" row can
## fall relative to the unrestricted version for reasons other than people
## leaving.
restrict_to_reported <- function(df, arm) {
  before <- n_distinct(df$person_id)
  out    <- df %>% filter(match_cancer_type %in% reportable_types)
  after  <- n_distinct(out$person_id)
  message(sprintf("%s: %s of %s kept after restricting to the %d reported cancer types (%s dropped).",
                  arm, scales::comma(after), scales::comma(before),
                  length(reportable_types), scales::comma(before - after)))
  out
}

matched_controls <- restrict_to_reported(matched_controls, "Controls")
cancers          <- restrict_to_reported(cancers,          "Cases")

stopifnot(nrow(matched_controls) > 0, nrow(cancers) > 0)

# load person level metadata
matched_controls_person_metadata <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/control_individuals/pre_psm_matched_metadata.tsv"))

## THE PREDIAGNOSIS FILE, NOT THE ANYTIME ONE. preprocess_cancer_data.R writes
## both, and they differ in exactly the BMI-history columns -- earliest_bmi_year
## and n_bmi_measurements -- because one summarises a case's pre-diagnosis
## measurements and the other all of them. psm.R matched on the PREDIAGNOSIS
## table, so the anytime table's earliest_bmi_year is not the value the match
## saw, and a Table 1 built on it would report a covariate the match never used
## while claiming balance on it. Controls have one metadata table and are
## unaffected.
cancers_person_metadata <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/pre_psm_matched_cancer_bmi_summary_prediagnosis.tsv"))

metadata <- bind_rows(matched_controls_person_metadata, cancers_person_metadata) %>%
  ## psm.R derives this from primary_consent_date and calipers on it; it is not
  ## a stored column, so it has to be rebuilt the same way here.
  mutate(enrollment_year = year(as.Date(primary_consent_date)))

## Every covariate the match used has to be present, or the table silently
## drops a row of itself. Named explicitly rather than discovered at select().
PSM_MATCH_COLS <- c("t_disp_sex_at_birth", "t_disp_self_reported_category",
                    "earliest_bmi_year", "enrollment_year", "age_at_enrollment")
missing_cols <- setdiff(PSM_MATCH_COLS, names(metadata))
if (length(missing_cols) > 0) {
  stop("metadata is missing PSM match variable(s): ", paste(missing_cols, collapse = ", "),
       ". Table 1 must report the variables psm.R matched on -- see its header.")
}

print(head(metadata))


## --- per-person measurement counts -----------------------------------------
##
## Count DISTINCT (person_id, measurement_datetime), not rows: a control
## matched to k cancer types carries k copies of each measurement, so counting
## rows would multiply their measurement count by roughly k while leaving the
## cases, which have one cancer type each, unaffected.

count_measurements <- function(df) {
  df %>%
    distinct(person_id, measurement_datetime) %>%
    count(person_id, name = "n_measurements")
}

ctrl_measurement_counts    <- count_measurements(matched_controls)
cancers_measurement_counts <- count_measurements(cancers)

# Combine the two count tables into one, one row per person
all_counts <- bind_rows(ctrl_measurement_counts, cancers_measurement_counts)

## The cohorts are disjoint, so no person should be counted twice here. If one were, the
## right_join below would fan their metadata row out and double-count them in
## every row of the table, which is exactly the kind of thing this table is
## supposed to be trusted not to do -- so stop rather than warn.
if (any(duplicated(all_counts$person_id))) {
  stop(sum(duplicated(all_counts$person_id)),
       " person_id(s) appear in BOTH the cancer and control cohort tables. ",
       "The cancer and control cohorts must be disjoint.")
}

# Merge onto metadata by person_id
metadata_merged <- metadata %>%
  right_join(all_counts, by = "person_id")


## --- Table 1 ---------------------------------------------------------------

# Helper function to collapse non-response variants into one category
collapse_unknown <- function(x) {
  x <- as.character(x)
  x <- case_when(
    str_detect(str_to_lower(x), "prefer not to answer|prefer not|pna") ~ "Unknown/Not Reported",
    str_detect(str_to_lower(x), "unknown") ~ "Unknown/Not Reported",
    str_detect(str_to_lower(x), "skip") ~ "Unknown/Not Reported",
    is.na(x) ~ "Unknown/Not Reported",
    TRUE ~ x
  )
  x
}

## ONE ROW FOR THREE CATEGORIES: American Indian or Alaska Native + Middle
## Eastern or North African + Native Hawaiian or Other Pacific Islander. THOSE
## THREE ONLY.
##
## THE LABEL NAMES ITS MEMBERS. This row was "Other", a label that did not say
## what was in it, so the membership lived only in this comment and the row
## could not be read from the table alone. Naming the three categories in the
## label costs a long row header and removes the need for a Methods sentence to
## carry the definition.
##
## IT IS STILL NOT A CATCH-ALL, and should not be read as one: "None of these",
## "More than one population" and "Unknown/Not Reported" remain their own rows.
## Nothing else joins this one.
##
## PRESENTATION ONLY. psm.R matched EXACTLY on self_reported_category with its
## full set of levels, so the arms are balanced within each member category
## separately -- collapsing here does not pool anything the match pooled, it
## merges three display rows, and the merged counts remain the sum of three
## separately balanced strata.
##
## Matched on substrings rather than the full concept name so a relabelling
## upstream (the AoU concept names carry prefixes in some releases) does not
## silently drop the rule and leave the categories reported apart.
AIAN_MENA_NHPI <- paste(
  "American Indian or Alaska Native",
  "Middle Eastern or North African",
  "Native Hawaiian or Other Pacific Islander",
  sep = "/")

collapse_small_race <- function(x) {
  x <- as.character(x)
  case_when(
    str_detect(str_to_lower(x), "american indian|alaska native") ~ AIAN_MENA_NHPI,
    str_detect(str_to_lower(x), "middle eastern|north african") ~ AIAN_MENA_NHPI,
    str_detect(str_to_lower(x), "hawaiian|pacific islander") ~ AIAN_MENA_NHPI,
    TRUE ~ x
  )
}

## A YEAR IS A LABEL, NOT A QUANTITY. gtsummary runs continuous statistics
## through style_number(), which carries big.mark = "," and picks a decimal
## place, so the year 2017 rendered as "2,017.0" -- a thousands separator and a
## tenth of a year, on a calendar year. Both year rows are formatted as plain
## integers instead.
##
## The variables STAY NUMERIC. This is a formatting change, not a change of
## summary: the median/IQR/min/max is still what these rows report, for the
## reason in the note below (an exact factor with ~45 levels would otherwise
## be ~45 table rows). Only the printing of the numbers changes.
##
## A MEDIAN CAN LAND ON A HALF-YEAR when the arm has an even number of people,
## and round() is banker's rounding, so 2014.5 prints as 2014 while 2015.5
## prints as 2016. The alternative is a decimal place on a year, which is what
## this is fixing. Quote the IQR alongside it if the exact midpoint matters.
style_year <- function(x) {
  ifelse(is.na(x), NA_character_,
         format(round(as.numeric(x)), trim = TRUE, scientific = FALSE,
                big.mark = ""))
}

## THE ROWS ARE THE PSM COVARIATES, IN THE ORDER psm.R USES THEM.
##
## All five matched covariates are reported: sex at birth, self_reported_category
## (the field psm.R matches on -- not ethnicity, which is a different AoU field
## and was never matched on), age at enrollment, earliest_bmi_year and
## enrollment_year. The last two are summarised as continuous for presentation;
## psm.R matches earliest_bmi_year as an exact factor with ~45 levels.
##
## The BMI measurement count is NOT matched on and is labelled as such, so the
## table does not imply a balance that was never enforced.

table1 <- metadata_merged %>%
  mutate(
    case = factor(case, levels = c(0, 1), labels = c("Control", "Case")),
    t_disp_sex_at_birth           = collapse_unknown(t_disp_sex_at_birth),
    t_disp_self_reported_category = collapse_small_race(
                                      collapse_unknown(t_disp_self_reported_category))
  ) %>%
  select(case,
         ## matched on
         t_disp_sex_at_birth, t_disp_self_reported_category,
         earliest_bmi_year, enrollment_year, age_at_enrollment,
         ## not matched on
         n_measurements) %>%
  tbl_summary(
    by = case,
    ## STATED, NOT GUESSED. tbl_summary() infers a variable's type, and its rule
    ## makes a numeric with fewer than 10 distinct values CATEGORICAL.
    ## enrollment_year spans 2017-2024 -- 8 values -- so it was being split into
    ## one row per year, while earliest_bmi_year (~45 values) stayed continuous:
    ## two year covariates summarised two different ways, neither of them chosen.
    ##
    ## It also broke the percentages. A categorical row's statistics are {n} and
    ## {p}, and {p} is held as a FRACTION -- gtsummary's default style_percent()
    ## is what multiplies it by 100. The digits = style_year below applies to
    ## every statistic of the variable it names, so {p} went through round()
    ## instead and every enrollment-year row printed "(0%)", the large ones
    ## included.
    ##
    ## Declaring the type fixes both and keeps the summary the comment above
    ## describes. Do not drop this line because the data "looks" continuous: the
    ## cohort's enrollment span is what decides, and it is under the threshold.
    type = list(
      earliest_bmi_year ~ "continuous",
      enrollment_year   ~ "continuous"
    ),
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75}), min: {min}, max: {max}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    ## Applies to every statistic of these two variables -- median, p25, p75,
    ## min and max all print as bare years. age_at_enrollment and
    ## n_measurements keep the default formatting: those ARE quantities.
    digits = list(
      earliest_bmi_year ~ style_year,
      enrollment_year   ~ style_year
    ),
    label = list(
      t_disp_sex_at_birth           ~ "Sex at Birth (exact match)",
      t_disp_self_reported_category ~ "Self-Reported Race/Ethnicity (exact match)",
      earliest_bmi_year             ~ "Year of Earliest BMI (exact match)",
      enrollment_year               ~ "Year of Enrollment (caliper +/- 2 years)",
      age_at_enrollment             ~ "Age at Enrollment (caliper +/- 2 years)",
      n_measurements                ~ "BMI Measurements per Person (not matched on)"
    ),
    missing = "no"
  )


## --- export ----------------------------------------------------------------

# Clean column names: remove embedded newlines and markdown bold
table1_flat <- table1
names(table1_flat) <- names(table1_flat) %>%
  str_replace_all("\n", " ") %>%
  str_remove_all("\\*\\*")

out_dir <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/results/tables/data_overview")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

table1_flat %>%
  gtsummary::as_tibble(col_labels = FALSE) %>%
  write_table(file.path(out_dir, "table1_summary.tsv"))
