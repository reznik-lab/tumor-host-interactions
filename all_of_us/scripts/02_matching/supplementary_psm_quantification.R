#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 02.4 -- cohort size and match quality per cancer type.
##
## Produces: Supplementary Table 2 -- cohort size and match quality per cancer
##           type.
##
## Reads:
##   <workspace>/.../stable_data/psm_matched_usable_data/solid_heme/full_cohort/ geq1 tables
## Writes:
##   .../results/tables/psm_quantification/supp_fig1c_by_cancer_type_summary.tsv (+ .xlsx)
##       [Supplementary Table 2]
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


## --- 1. Data ---------------------------------------------------------------

psm_cancer_matches <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm_matched_usable_data/solid_heme/full_cohort/all_cancers_geq1_before_diagnosis_after_psm.tsv"))

psm_control_matches <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm_matched_usable_data/solid_heme/full_cohort/all_controls_geq1_before_diagnosis_after_psm.tsv"))

print(head(psm_cancer_matches))

psm_cancer_matches_subset <- psm_cancer_matches %>%
  select(person_id, case, subclass, match_cancer_type)

psm_control_matches_subset <- psm_control_matches %>%
  select(person_id, case, subclass, match_cancer_type)

psm_stats_df <- bind_rows(psm_cancer_matches_subset, psm_control_matches_subset) %>%
  distinct()

print(head(psm_stats_df))


## --- the 33: types clearing the pipeline-wide floor ------------------------
## ONE definition of "the 33", shared with every other stage rather than
## recomputed here. The floor is >= MIN_CASES_PER_TYPE MATCHED cases, counted
## on the raw PSM output; this script's own cohort is the geq1 one, whose case
## counts are smaller, so gating on those would admit a different set and
## quietly stop being the 33 the rest of the manuscript means.
##
## This is a list of type NAMES only. Every count below is still computed over
## this script's own geq1 cohort.
reportable_types <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm/solid_heme/full_cohort/aggregate/psm_matched_cohort_all_cancer_types.tsv")) %>%
  filter(case == 1) %>%
  distinct(person_id, match_cancer_type) %>%
  count(match_cancer_type, name = "n_matched_cases") %>%
  filter(n_matched_cases >= MIN_CASES_PER_TYPE) %>%
  pull(match_cancer_type)

message(sprintf("%d cancer types clear MIN_CASES_PER_TYPE = %d.",
                length(reportable_types), MIN_CASES_PER_TYPE))


## --- 0. Setup --------------------------------------------------------------

## This used to be code_outputs_scripts/ itself, which dropped six TSVs and a
## PDF at the top level beside stable_data/ and temp/. Tables and figures now
## go to the same results/ tree as every other script.
tab_dir <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/results/tables/psm_quantification")
fig_dir <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/results/figs/psm_quantification")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_std())

# psm_stats_df expected columns: person_id, case (1 = case, 0 = control), subclass, match_cancer_type
psm_stats_df <- psm_stats_df %>%
  mutate(case = as.integer(case))


## --- 2. Overall supplementary values ---------------------------------------

n_controls_total <- psm_stats_df %>%
  filter(case == 0) %>%
  distinct(person_id) %>%
  nrow()

n_cases_total <- psm_stats_df %>%
  filter(case == 1) %>%
  distinct(person_id) %>%
  nrow()

# Matches per case (controls per subclass)
matches_per_case <- psm_stats_df %>%
  filter(case == 0) %>%
  count(subclass, name = "n_controls")

avg_matches_per_case    <- mean(matches_per_case$n_controls)
median_matches_per_case <- median(matches_per_case$n_controls)
range_matches_per_case  <- range(matches_per_case$n_controls)

## The quartiles are here because the manuscript reports the match ratio as a
## median with an IQR, not as a mean: the match is variable-ratio 1:k with a
## caliper, so the distribution is capped at k and left-skewed, and a mean
## sitting between two quartiles describes no actual case.
q_matches_per_case <- quantile(matches_per_case$n_controls, c(0.25, 0.75), type = 7)

summary_overall <- tibble(
  n_cases                 = n_cases_total,
  n_controls              = n_controls_total,
  n_cancer_types          = n_distinct(psm_stats_df$match_cancer_type),
  avg_matches_per_case    = round(avg_matches_per_case, 2),
  median_matches_per_case = median_matches_per_case,
  q1_matches_per_case     = unname(q_matches_per_case[1]),
  q3_matches_per_case     = unname(q_matches_per_case[2]),
  min_matches_per_case    = range_matches_per_case[1],
  max_matches_per_case    = range_matches_per_case[2]
)

print(summary_overall)

print(glue::glue(
  "yielding a total of {n_controls_total} controls and an average of ",
  "{round(avg_matches_per_case, 1)} matches per case"
))

## The sentence the manuscript actually quotes, printed in its own shape so a
## re-run can be checked against the text without re-deriving it by hand.
print(glue::glue(
  "{format(n_cases_total, big.mark = ',')} cases across ",
  "{summary_overall$n_cancer_types} cancer types, matched to ",
  "{format(n_controls_total, big.mark = ',')} controls, with a median of ",
  "{median_matches_per_case} matches per case ",
  "(IQR {summary_overall$q1_matches_per_case}-{summary_overall$q3_matches_per_case})"
))


## --- 2b. The same overall values, restricted to the reportable types -------
## What the figures and analyses actually cover. summary_overall above is the
## match as a whole, including types too small for any stage to report; this is
## the subset the manuscript draws, and it is the pair of numbers to quote
## beside a per-type figure.
##
## CONTROLS ARE RECOUNTED, NOT SUMMED. A control matched to several cancer
## types is one person; adding the per-type control counts over the reportable
## types would count them once per type. n_distinct over the subset's own rows
## is the only correct way to do this.

psm_stats_reportable <- psm_stats_df %>%
  filter(match_cancer_type %in% reportable_types)

matches_per_case_reportable <- psm_stats_reportable %>%
  filter(case == 0) %>%
  count(subclass, name = "n_controls")

q_reportable <- quantile(matches_per_case_reportable$n_controls, c(0.25, 0.75), type = 7)

summary_overall_reportable <- tibble(
  n_cases                 = n_distinct(psm_stats_reportable$person_id[psm_stats_reportable$case == 1]),
  n_controls              = n_distinct(psm_stats_reportable$person_id[psm_stats_reportable$case == 0]),
  n_cancer_types          = n_distinct(psm_stats_reportable$match_cancer_type),
  avg_matches_per_case    = round(mean(matches_per_case_reportable$n_controls), 2),
  median_matches_per_case = median(matches_per_case_reportable$n_controls),
  q1_matches_per_case     = unname(q_reportable[1]),
  q3_matches_per_case     = unname(q_reportable[2]),
  min_matches_per_case    = min(matches_per_case_reportable$n_controls),
  max_matches_per_case    = max(matches_per_case_reportable$n_controls)
)

print(summary_overall_reportable)

## Said as the sentence it would appear as, and beside the unrestricted one so
## the cost of the restriction is visible rather than inferred.
print(glue::glue(
  "restricted to the {summary_overall_reportable$n_cancer_types} cancer types with >= ",
  "{MIN_CASES_PER_TYPE} matched cases: ",
  "{format(summary_overall_reportable$n_cases, big.mark = ',')} cases and ",
  "{format(summary_overall_reportable$n_controls, big.mark = ',')} controls ",
  "(of {format(n_cases_total, big.mark = ',')} and ",
  "{format(n_controls_total, big.mark = ',')} across all ",
  "{summary_overall$n_cancer_types} types)"
))


## --- 3. Breakdown by cancer type -------------------------------------------
## Controls per case, one row per cancer type, with a pan-cancer row on top.
##
## matches_per_case_by_type is built FIRST: the per-type distribution stats
## below are computed from it, so it has to exist before the summary does.

matches_per_case_by_type <- psm_stats_df %>%
  filter(case == 0) %>%
  count(subclass, match_cancer_type, name = "n_controls") %>%
  mutate(reportable = match_cancer_type %in% reportable_types)

## The headline row, over the cancer types clearing MIN_CASES_PER_TYPE.
PAN_CANCER_LABEL <- sprintf("All cancers (%d types)", length(reportable_types))

#' Cases, distinct controls and the controls-per-case distribution, for one
#' set of rows of psm_stats_df.
#'
#' POOLED, NOT AVERAGED. The quantiles come from every subclass in `stats` at
#' once, so the pan-cancer row is the distribution over all cases of the
#' reportable types -- not the mean of the per-type medians, which would weight
#' a 107-case type the same as a 3,000-case one.
summarise_match_ratio <- function(stats, label) {
  per_case <- stats %>% filter(case == 0) %>% count(subclass, name = "k")
  tibble(
    cancer_type    = label,
    n_cases        = n_distinct(stats$person_id[stats$case == 1]),
    ## DISTINCT people. For a single type this equals the number of pairs
    ## (replace = FALSE within a type); pooled over types it does NOT, because
    ## one control serves several. Recounted here either way, never summed.
    n_controls     = n_distinct(stats$person_id[stats$case == 0]),
    mean_matches   = round(mean(per_case$k), 2),
    median_matches = median(per_case$k),
    q1_matches     = quantile(per_case$k, 0.25, type = 7),
    q3_matches     = quantile(per_case$k, 0.75, type = 7),
    iqr_matches    = IQR(per_case$k, type = 7)
  )
}

## EVERY ROW IS REPORTABLE, so there is no `reportable` column: this table is
## the reported set and nothing else. The types below the floor are still
## visible elsewhere in this script -- summary_by_type_detail,
## matches_per_case_by_type and the boxplot carry every type with the flag --
## so the whole-match view is not lost by restricting here.
summary_by_type <- bind_rows(
  ## the pan-cancer row, pooled over the same types the rows below list
  summarise_match_ratio(psm_stats_reportable, PAN_CANCER_LABEL) %>%
    mutate(is_pan_cancer = TRUE, match_cancer_type = NA_character_),
  psm_stats_reportable %>%
    group_split(match_cancer_type) %>%
    map_dfr(~ summarise_match_ratio(.x, display_type_label(.x$match_cancer_type[1])) %>%
              mutate(is_pan_cancer = FALSE,
                     match_cancer_type = .x$match_cancer_type[1])) %>%
    arrange(desc(n_cases))
) %>%
  relocate(cancer_type, n_cases, n_controls,
           mean_matches, median_matches, q1_matches, q3_matches, iqr_matches,
           is_pan_cancer, match_cancer_type)

## The pan-cancer row covers exactly the rows beneath it.
stopifnot(summary_by_type$is_pan_cancer[1],
          sum(summary_by_type$is_pan_cancer) == 1,
          nrow(summary_by_type) - 1 == n_distinct(psm_stats_reportable$match_cancer_type))

## Summing n_controls down the table is wrong and the pan-cancer row shows by
## how much: it recounts distinct people, so it is smaller than the column sum
## over the same types. Said in the log rather than left to be rediscovered.
message(sprintf(
  "Pan-cancer row: %s cases, %s distinct controls. Summing the per-type control column over those %d types instead gives %s -- %.2fx, because one control serves several types.",
  scales::comma(summary_by_type$n_cases[1]),
  scales::comma(summary_by_type$n_controls[1]),
  length(reportable_types),
  scales::comma(sum(summary_by_type$n_controls[!summary_by_type$is_pan_cancer])),
  sum(summary_by_type$n_controls[!summary_by_type$is_pan_cancer]) /
    summary_by_type$n_controls[1]))

## The floor is defined on MATCHED cases while this table is built from the
## geq1 cohort, so a type could in principle clear the floor and have no row
## here. Say so rather than silently writing a short table.
missing_here <- setdiff(reportable_types, summary_by_type$match_cancer_type)
if (length(missing_here) > 0) {
  warning("clears the floor but has no row in summary_by_type: ",
          paste(missing_here, collapse = ", "), call. = FALSE)
}

print(as.data.frame(summary_by_type %>% select(-match_cancer_type)), row.names = FALSE)

summary_by_type_detail <- matches_per_case_by_type %>%
  group_by(match_cancer_type) %>%
  summarise(
    avg_matches_per_case    = round(mean(n_controls), 2),
    median_matches_per_case = median(n_controls),
    min_matches_per_case    = min(n_controls),
    max_matches_per_case    = max(n_controls),
    n_cases_with_matches    = n(),
    .groups = "drop"
  ) %>%
  mutate(reportable = match_cancer_type %in% reportable_types)

print(summary_by_type_detail)


## --- 4. Write tables (TSV) -------------------------------------------------

write_table(summary_by_type, file.path(tab_dir, "supp_fig1c_by_cancer_type_summary.tsv"))
