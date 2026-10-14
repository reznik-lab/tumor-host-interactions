#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 03.2 -- pre-wrangle cachexia calls: attach cancer type + diagnosis date.
##
## Produces no manuscript figure or table: attaches cancer type and diagnosis
## date to the called cachexia episodes.
##
## Loops over all three weight-loss thresholds, so each output file holds the
## threshold its name claims.
##
## Writes (one per threshold):
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/temp/cachexia_calls/wl_w_meta_{5,10,15}_pct.tsv
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


## --- which episode files to postprocess -------------------------------------
## calling_cachexia_episodes.py stamps its outputs with the date it ran:
##   df_episodes_all_precomp_<WL label>_<YYYYMMDD>.csv
## Set EPISODE_DATE to that run's stamp. Pinned explicitly rather than globbed
## for the newest file, so a stale episode file can never be picked up silently.
EPISODE_DATE <- "20260929"

## label in the caller's filename -> suffix in this script's output filename
WL_RUNS <- c("WL5" = "5", "WL10" = "10", "WL15" = "15")


## --- read ------------------------------------------------------------------

# get cancer type and diagnosis date labels:
# NOTE: these data files have been prefiltered. they are all BMI measurements of
# individuals who have at least 4 pre-diagnosis BMI measurements.
cancer_meta <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm_matched_usable_data/solid_heme/full_cohort/all_cancers_all_bmi_measurements_geq4_eligible.tsv")) %>%
  select(person_id, case, subclass, cancer_type, match_cancer_type, pseudo_diagnosis_datetime) %>%
  distinct()

control_meta <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm_matched_usable_data/solid_heme/full_cohort/all_controls_all_bmi_measurements_geq4_eligible.tsv")) %>%
  select(person_id, case, subclass, cancer_type, match_cancer_type, pseudo_diagnosis_datetime) %>%
  distinct()

merge_df <- bind_rows(cancer_meta, control_meta)


## --- join + write, one file per weight-loss threshold -----------------------
## Each output keeps the threshold it was actually computed at, so the WL label
## in the filename matches the contents.

out_dir <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/temp/cachexia_calls")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (wl_label in names(WL_RUNS)) {
  pct <- WL_RUNS[[wl_label]]

  in_path <- file.path(WORKSPACE, sprintf(
    "pre_diagnosis_bmi_analysis/code_outputs_scripts/cachexia_episodes/processed_data/df_episodes_all_precomp_%s_%s.csv",
    wl_label, EPISODE_DATE))

  if (!file.exists(in_path)) {
    stop("No episode file for ", wl_label, " at ", in_path,
         "\n  Check that calling_cachexia_episodes.py has run and that",
         " EPISODE_DATE (", EPISODE_DATE, ") matches its output stamp.")
  }

  calls_w_meta <- read_csv(in_path) %>%
    right_join(merge_df, by = "person_id")

  out_path <- file.path(out_dir, sprintf("wl_w_meta_%s_pct.tsv", pct))
  write_tsv(calls_w_meta, out_path)

  message("wrote ", nrow(calls_w_meta), " rows -> ", out_path)
}
