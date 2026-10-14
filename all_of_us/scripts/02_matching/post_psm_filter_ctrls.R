#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 02.3 -- build the analytic cohort tables every figure script reads.
##
## Produces no manuscript figure or table: applies the pre-diagnosis BMI
## coverage gates to the matched cohort.
##
## Reads:
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm/solid_heme/full_cohort/aggregate/psm_matched_cohort_all_cancer_types.tsv
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/control_individuals/raw_concatenated/measurementOccurrence.tsv
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/raw_concatenated/measurementOccurrence.tsv
## Writes (all into .../stable_data/psm_matched_usable_data/solid_heme/full_cohort/):
##   all_{controls,cancers}_geq4_before_diagnosis_after_psm.tsv
##   all_{controls,cancers}_geq1_before_diagnosis_after_psm.tsv
##   all_{controls,cancers}_all_bmi_measurements_geq4_eligible.tsv
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


## --- load ------------------------------------------------------------------

merged_df <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm/solid_heme/full_cohort/aggregate/psm_matched_cohort_all_cancer_types.tsv"))

# we have an issue where control BMI also has the person ids of the cancer group
control_bmi <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/control_individuals/raw_concatenated/measurementOccurrence.tsv")) %>%
  janitor::clean_names()

cancer_bmi <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/raw_concatenated/measurementOccurrence.tsv")) %>%
  janitor::clean_names()

controls_pseudo_diagnosis_df <- merged_df %>%
  filter(case == 0) %>%
  select(person_id, case, subclass, cancer_type, match_cancer_type, pseudo_diagnosis_datetime)

cancer_pseudo_diagnosis_df <- merged_df %>%
  filter(case == 1) %>%
  select(person_id, case, subclass, cancer_type, match_cancer_type, pseudo_diagnosis_datetime)

## The eligibility gates below key on (person_id, subclass) as the identity of a
## single match. That only means anything if it is unique -- if a person could
## appear twice in one subclass, the gates would count their measurements twice.
## psm.R guarantees it (replace = FALSE within a cancer type, and subclass is
## prefixed with the cancer type so it is globally unique), but assert it here,
## since every downstream cohort is built on it.
stopifnot(
  !any(duplicated(controls_pseudo_diagnosis_df[c("person_id", "subclass")])),
  !any(duplicated(cancer_pseudo_diagnosis_df[c("person_id", "subclass")]))
)


## --- prefilter -------------------------------------------------------------

# prefilter bmi values (physiologically plausible range), dedupe
# The range gate is filter_usable_bmi() in src/setup.R, shared with Stage 01 and
# Stage 05 so every stage selects the same rows.
prefilter_bmi <- function(bmi_df) {
  bmi_df %>%
    select(person_id, measurement_datetime, t_disp_standard_concept_name, value_as_number) %>%
    filter_usable_bmi() %>%
    distinct()
}

ctrl_bmi   <- prefilter_bmi(control_bmi)
cancer_bmi <- prefilter_bmi(cancer_bmi)


#' Join pre-filtered BMI measurements to pseudo-diagnosis dates, restrict to
#' BMI measurements taken before the pseudo-diagnosis date, and require at
#' least `min_n` such measurements per MATCH (see MATCH, NOT PERSON below).
#'
#' @param bmi_df prefiltered BMI df (person_id, measurement_datetime, t_disp_standard_concept_name, value_as_number)
#' @param diagnosis_df df with person_id, subclass, pseudo_diagnosis_datetime (+ metadata cols)
#' @param min_n minimum number of pre-diagnosis BMI measurements required per match.
#'   Use 4 for the "geq4" analysis subset, 1 for the "full cohort" (any measurement) subset.
#' @param join_by column(s) to join bmi_df and diagnosis_df on
get_pre_diagnosis_bmi_samples <- function(bmi_df, diagnosis_df, min_n = 1, join_by = "person_id") {
  bmi_df %>%
    inner_join(diagnosis_df, by = join_by) %>%
    mutate(
      measurement_datetime      = as_datetime(measurement_datetime),
      pseudo_diagnosis_datetime = as_datetime(pseudo_diagnosis_datetime)
    ) %>%
    filter(
      str_detect(t_disp_standard_concept_name, regex("BMI", ignore_case = TRUE)),
      measurement_datetime < pseudo_diagnosis_datetime
    ) %>%
    ## group_by(person_id, subclass) -- the MATCH, not the person. See the
    ## MATCH, NOT PERSON note on get_all_bmi_samples below; the same pooling
    ## defect was here, and the fix is the same.
    group_by(person_id, subclass) %>%
    ## n_distinct(measurement_datetime), NOT n(): prefilter_bmi() keeps a row per
    ## distinct (datetime, value), so a person with two different BMI values
    ## recorded at the same instant would otherwise count twice toward min_n.
    filter(n_distinct(measurement_datetime) >= min_n) %>%
    ungroup()
}


#' Join pre-filtered BMI measurements to pseudo-diagnosis dates, restrict to
#' MATCHES with at least `min_n` BMI measurements before that match's
#' pseudo-diagnosis date, but then return ALL of that person's BMI measurements
#' (not just the pre-diagnosis ones). Used for cachexia episode detection, where
#' an episode can start pre-diagnosis (e.g. up to 6 months prior) and continue
#' into/past diagnosis -- the min_n requirement establishes eligibility, but
#' episode detection needs the full trajectory.
#'
#' @param bmi_df prefiltered BMI df (person_id, measurement_datetime, t_disp_standard_concept_name, value_as_number)
#' @param diagnosis_df df with person_id, subclass, pseudo_diagnosis_datetime (+ metadata cols)
#' @param min_n minimum number of pre-diagnosis BMI measurements a match must have
#'   for that match to be included at all (default 4, matching the geq4 analysis)
#' @param join_by column(s) to join bmi_df and diagnosis_df on
get_all_bmi_samples <- function(bmi_df, diagnosis_df, min_n = 4, join_by = "person_id") {
  joined <- bmi_df %>%
    inner_join(diagnosis_df, by = join_by) %>%
    mutate(
      measurement_datetime      = as_datetime(measurement_datetime),
      pseudo_diagnosis_datetime = as_datetime(pseudo_diagnosis_datetime)
    ) %>%
    filter(
      str_detect(t_disp_standard_concept_name, regex("BMI", ignore_case = TRUE))
    )

  ## -------------------------------------------------------------------------
  ## ELIGIBILITY IS PER MATCH, NOT PER PERSON.
  ##
  ## psm.R matches each cancer type against the full control pool (replace =
  ## FALSE only *within* a type), so one control can sit in one subclass per
  ## cancer type, each with a DIFFERENT pseudo-diagnosis date. "Before
  ## diagnosis" therefore means a different set of measurements for each of a
  ## control's matches, and min_n has to be tested against each of them
  ## separately. (person_id, subclass) is that unit: a control belongs to at
  ## most one subclass per cancer type, and a case to exactly one overall.
  ##
  ## !! FIXED: THIS USED TO COUNT PER PERSON, POOLED OVER ALL THEIR MATCHES !!
  ## It was distinct(person_id, measurement_datetime) |> count(person_id), i.e.
  ## the UNION of the pre-diagnosis measurements across every match, and the
  ## resulting verdict was applied to all of that person's matches via
  ## `person_id %in% eligible_ids`. Because "before diagnosis" is a prefix of a
  ## person's time-ordered measurements, that union is exactly the prefix of the
  ## match with the LATEST pseudo-diagnosis date -- so the gate was decided by
  ## that one match and then granted to all the others. A control with
  ##
  ##     measurements 2015, 2016, 2017, 2018
  ##     match A pseudo-dx 2015-06  -> 1 pre-dx measurement
  ##     match B pseudo-dx 2016-06  -> 2
  ##     match D pseudo-dx 2019-01  -> 4
  ##
  ## pooled to 4 and entered the geq4 cohort under ALL THREE matches, carrying
  ## one measurement under A and two under B. Cases have one match each and were
  ## never affected, so the admitted rows were all controls -- and those rows
  ## have sparse pre-diagnosis coverage, so a qualifying weight loss is less
  ## likely to be observed in them. They landed in the denominator as non-calls,
  ## deflating the control rate and inflating the case-vs-control contrast in
  ## Figure 1B and Supp Fig 1D.
  ##
  ## (A sensitivity analysis that is not reported here keyed its spacing gate
  ## on (person_id, pseudo_diagnosis_datetime) and so compensated for this
  ## downstream; the main path did not. The pairs that only existed because of
  ## this defect no longer do.)
  ## -------------------------------------------------------------------------
  eligible_matches <- joined %>%
    filter(measurement_datetime < pseudo_diagnosis_datetime) %>%
    ## distinct before counting: prefilter_bmi() keeps a row per distinct
    ## (datetime, value), so two values recorded at the same instant are one
    ## measurement, not two.
    distinct(person_id, subclass, measurement_datetime) %>%
    count(person_id, subclass) %>%
    filter(n >= min_n) %>%
    select(person_id, subclass)

  # return ALL measurements (pre- and post-diagnosis) for those eligible matches
  joined %>%
    semi_join(eligible_matches, by = c("person_id", "subclass")) %>%
    mutate(
      months_from_diagnosis = as.numeric(difftime(measurement_datetime, pseudo_diagnosis_datetime, units = "days")) / 30.44
    )
}


## --- run and write ---------------------------------------------------------

out_dir <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm_matched_usable_data/solid_heme/full_cohort")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#' Run the pipeline for one (bmi_df, diagnosis_df) pair at a given threshold
#' and write the result out with a threshold-specific filename.
run_and_write <- function(bmi_df, diagnosis_df, min_n, label, out_dir) {
  result <- get_pre_diagnosis_bmi_samples(bmi_df, diagnosis_df, min_n = min_n)
  fname <- glue::glue("{label}_geq{min_n}_before_diagnosis_after_psm.tsv")
  write_tsv(result, file.path(out_dir, fname))
  result
}

# --- geq4 subset (existing analysis) ---
controls_geq4 <- run_and_write(ctrl_bmi,   controls_pseudo_diagnosis_df, min_n = 4, label = "all_controls", out_dir = out_dir)
cancer_geq4   <- run_and_write(cancer_bmi, cancer_pseudo_diagnosis_df,   min_n = 4, label = "all_cancers",  out_dir = out_dir)

# --- full cohort subset (any measurement pre-diagnosis, no >=4 requirement) ---
controls_full <- run_and_write(ctrl_bmi,   controls_pseudo_diagnosis_df, min_n = 1, label = "all_controls", out_dir = out_dir)
cancer_full   <- run_and_write(cancer_bmi, cancer_pseudo_diagnosis_df,   min_n = 1, label = "all_cancers",  out_dir = out_dir)

# --- all BMI measurements (pre- AND post-diagnosis) for people meeting the
#     geq4-before-diagnosis eligibility criterion. Used for cachexia episode
#     detection, since an episode can start pre-diagnosis and continue past it. ---
all_bmi_controls <- get_all_bmi_samples(ctrl_bmi,   controls_pseudo_diagnosis_df, min_n = 4)
all_bmi_cancer   <- get_all_bmi_samples(cancer_bmi, cancer_pseudo_diagnosis_df,   min_n = 4)

write_tsv(all_bmi_controls, file.path(out_dir, "all_controls_all_bmi_measurements_geq4_eligible.tsv"))
write_tsv(all_bmi_cancer,   file.path(out_dir, "all_cancers_all_bmi_measurements_geq4_eligible.tsv"))


## --- what the per-match gate cost -------------------------------------------
## The gates above used to pool a person's pre-diagnosis measurements across all
## of their matches (see the MATCH, NOT PERSON note on get_all_bmi_samples).
## Report the difference so the size of the correction is visible in the log
## rather than having to be re-derived: how many (person, match) pairs the
## pooled gate would have admitted that the per-match gate does not. Cases have
## one match each, so their two numbers are identical by construction -- if they
## are not, something is wrong with the subclass keys.

pooled_gate_pairs <- function(bmi_df, diagnosis_df, min_n = 4) {
  joined <- bmi_df %>%
    inner_join(diagnosis_df, by = "person_id") %>%
    mutate(
      measurement_datetime      = as_datetime(measurement_datetime),
      pseudo_diagnosis_datetime = as_datetime(pseudo_diagnosis_datetime)
    ) %>%
    filter(measurement_datetime < pseudo_diagnosis_datetime)

  ids <- joined %>%
    distinct(person_id, measurement_datetime) %>%
    count(person_id) %>%
    filter(n >= min_n) %>%
    pull(person_id)

  joined %>% filter(person_id %in% ids) %>% distinct(person_id, subclass) %>% nrow()
}

report_gate <- function(arm, bmi_df, diagnosis_df, kept) {
  n_now    <- nrow(distinct(kept, person_id, subclass))
  n_pooled <- pooled_gate_pairs(bmi_df, diagnosis_df, min_n = 4)
  message(sprintf(
    "geq4 gate [%s]: %d (person, match) pairs kept; the old pooled gate kept %d (%+d, %.1f%%). %d distinct people.",
    arm, n_now, n_pooled, n_now - n_pooled,
    if (n_pooled > 0) 100 * (n_now - n_pooled) / n_pooled else NA_real_,
    n_distinct(kept$person_id)))
}

report_gate("controls", ctrl_bmi,   controls_pseudo_diagnosis_df, all_bmi_controls)
report_gate("cancers",  cancer_bmi, cancer_pseudo_diagnosis_df,   all_bmi_cancer)

# sanity check for cancer subset -- should be one-to-one, so row count into the
# join should not expand relative to cancer_bmi
print(dim(cancer_pseudo_diagnosis_df))
print(head(cancer_pseudo_diagnosis_df))
