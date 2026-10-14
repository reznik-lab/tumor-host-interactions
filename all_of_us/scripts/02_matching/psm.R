#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 02.1 -- propensity score matching: cancer cases vs. controls.
##
## Produces no manuscript figure or table: the propensity score match every
## reported analysis is built on.
##
## Outputs per cancer type (dir_cancer_specific):
##   - match_object_<slug>.rds               : full matchit object (RDS only)
##   - matched_cohort_<slug>.rds / .tsv      : matched case+control rows
##   - balance_table_prematch_<slug>.tsv     : pre-match covariate balance (TSV only)
##   - balance_table_matched_<slug>.tsv      : post-match covariate balance (TSV only)
##
## Aggregate outputs (dir_aggregate):
##   - matched_cohort_all_cancer_types.rds / .tsv : all matched cohorts stacked
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
library(MatchIt)


## ---------------------------------------------------------------------------
## 0. Set up output directories
## ---------------------------------------------------------------------------
dir_cancer_specific <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm/solid_heme/full_cohort/cancer_specific")
dir_aggregate <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm/solid_heme/full_cohort/aggregate")

dir.create(dir_cancer_specific, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_aggregate,       showWarnings = FALSE, recursive = TRUE)

## helper to make a filesystem-safe version of a cancer_type_grouped value
sanitize_name <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}


## ---------------------------------------------------------------------------
## 1. Load / assemble data
## ---------------------------------------------------------------------------
cancer_df <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/pre_psm_matched_cancer_bmi_summary_prediagnosis.tsv"))
    # filter(n_bmi_measurements >= 4)

control_df <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/control_individuals/pre_psm_matched_metadata.tsv"))

DIAGNOSIS_DATE_COL <- "first_dx_datetime"
CANCER_TYPE_COL    <- "cancer_type_grouped"

stopifnot(DIAGNOSIS_DATE_COL %in% names(cancer_df))
stopifnot(CANCER_TYPE_COL %in% names(cancer_df))
stopifnot(all(cancer_df$case == 1))
stopifnot(all(control_df$case == 0))


## ---------------------------------------------------------------------------
## 2. Derive enrollment year and enforce types needed for exact matching /
##    caliper matching -- do this once, up front, on cases and controls
##    separately, since controls don't have cancer_type_grouped
## ---------------------------------------------------------------------------
prep <- function(df) {
  df %>%
    mutate(
      primary_consent_date = as.Date(primary_consent_date),
      enrollment_year      = year(primary_consent_date),

      sex_at_birth           = factor(sex_at_birth),
      self_reported_category = factor(self_reported_category),
      earliest_bmi_year      = factor(earliest_bmi_year),

      age_at_enrollment = as.numeric(age_at_enrollment),
      deprivation_index = as.numeric(deprivation_index)
    )
}

cancer_df  <- prep(cancer_df)
control_df <- prep(control_df)

# match_vars <- c("case", "sex_at_birth", "self_reported_category",
#                  "earliest_bmi_year", "enrollment_year",
#                  "age_at_enrollment", "deprivation_index", "n_bmi_measurements")

match_vars <- c("case", "sex_at_birth", "self_reported_category",
                "earliest_bmi_year", "enrollment_year",
                "age_at_enrollment")

n_cancer_before <- nrow(cancer_df)
cancer_df <- cancer_df %>% filter(if_all(all_of(match_vars), ~ !is.na(.)))
message(sprintf("Cases: dropped %d rows with missing match variables (of %d).",
                n_cancer_before - nrow(cancer_df), n_cancer_before))

n_control_before <- nrow(control_df)
control_df <- control_df %>% filter(if_all(all_of(match_vars), ~ !is.na(.)))
message(sprintf("Controls: dropped %d rows with missing match variables (of %d).",
                n_control_before - nrow(control_df), n_control_before))

cancer_types <- cancer_df %>%
  pull(.data[[CANCER_TYPE_COL]]) %>%
  unique() %>%
  sort()

message(sprintf("Found %d cancer_type_grouped values: %s",
                length(cancer_types), paste(cancer_types, collapse = ", ")))


## ---------------------------------------------------------------------------
## 3. Run one match per cancer type against the full control pool, then stitch
##    the results back together.
##
## The types are independent -- each matches against the same control pool and
## writes its own files -- so they run in parallel rather than one after the
## other. This is a scheduling change only: the matchit() call, its arguments
## and everything derived from it are unchanged, and nearest-neighbour matching
## on a fitted propensity score draws no random numbers, so the per-type
## outputs do not depend on which worker runs them or in what order. The
## combined table is rebuilt in cancer_types order below, so its rows are in
## the same order the serial loop produced.
##
## PSM_CORES sets the worker count (default: one fewer than the machine
## reports, leaving a core for everything else). PSM_CORES=1 restores the
## original serial run -- worth using if memory is tight, since each worker
## holds its own combined_ct and matchit object.
## ---------------------------------------------------------------------------
library(parallel)

PSM_CORES <- local({
  env <- Sys.getenv("PSM_CORES")
  n <- if (nzchar(env)) suppressWarnings(as.integer(env)) else max(1L, detectCores() - 1L)
  if (is.na(n) || n < 1L) n <- 1L
  ## mclapply forks, which is Unix-only; elsewhere fall back to serial.
  if (.Platform$OS.type != "unix") n <- 1L
  min(n, length(cancer_types))
})

## A threaded BLAS inside each forked worker would oversubscribe the cores the
## workers are already using. Numerically irrelevant, but it costs wall time.
if (PSM_CORES > 1L) Sys.setenv(OMP_NUM_THREADS = "1")

message(sprintf("Matching %d cancer types across %d worker(s).",
                length(cancer_types), PSM_CORES))

## The body below is the serial loop's body verbatim, returning matched_ct
## instead of assigning into matched_list.
match_one_cancer_type <- function(ct) {

  ct_slug <- sanitize_name(ct)   # computed up front; needed for match + balance table filenames

  cases_ct <- cancer_df %>% filter(.data[[CANCER_TYPE_COL]] == ct)

  combined_ct <- bind_rows(cases_ct, control_df)

  message(sprintf("\n--- Matching cancer_type_grouped = '%s' (n cases = %d) ---",
                  ct, nrow(cases_ct)))

  m.out <- tryCatch(
    matchit(
      case ~ age_at_enrollment +
             sex_at_birth + self_reported_category + earliest_bmi_year,
      data        = combined_ct,
      method      = "nearest",
      distance    = "glm",
      exact       = ~ sex_at_birth + self_reported_category + earliest_bmi_year,
      caliper     = c(enrollment_year = 2, age_at_enrollment = 2),
      std.caliper = FALSE,
      ratio       = 10,
      replace     = FALSE
    ),
    error = function(e) {
      warning(sprintf("Match failed for cancer_type_grouped = '%s': %s", ct, e$message))
      NULL
    }
  )

  if (is.null(m.out)) return(NULL)

  m.sum <- summary(m.out)

  ## -------------------------------------------------------------------------
  ## Save balance tables (pre-match and post-match) as TSV only
  ## (not printed to console, not saved as RDS)
  ## -------------------------------------------------------------------------
  if (is.null(m.sum$sum.matched)) {
    warning(sprintf("No post-match balance table available for '%s' (no matches?)", ct))
  } else {
    balance_matched <- as.data.frame(m.sum$sum.matched) %>%
      rownames_to_column("covariate")
    write_tsv(
      balance_matched,
      file.path(dir_cancer_specific, sprintf("balance_table_matched_%s.tsv", ct_slug))
    )
    rm(balance_matched)
  }

  if (!is.null(m.sum$sum.all)) {
    balance_all <- as.data.frame(m.sum$sum.all) %>%
      rownames_to_column("covariate")
    write_tsv(
      balance_all,
      file.path(dir_cancer_specific, sprintf("balance_table_prematch_%s.tsv", ct_slug))
    )
    rm(balance_all)
  }

  rm(m.sum)

  matched_ct <- match.data(m.out) %>%
    mutate(
      !!DIAGNOSIS_DATE_COL := as_datetime(.data[[DIAGNOSIS_DATE_COL]]),
      ## make subclass globally unique across cancer-type loops
      subclass = paste(ct, subclass, sep = "__"),
      ## tag every row (case AND its matched controls) with the
      ## cancer type this match was run for
      match_cancer_type = ct
    )

  ## -------------------------------------------------------------------------
  ## Assign pseudo tumor diagnosis date to matched controls
  ## = the diagnosis datetime of the case in the same subclass
  ## -------------------------------------------------------------------------
  diag_lookup <- matched_ct %>%
    filter(case == 1) %>%
    select(subclass, pseudo_diagnosis_datetime = all_of(DIAGNOSIS_DATE_COL))

  matched_ct <- matched_ct %>%
    left_join(diag_lookup, by = "subclass") %>%
    mutate(
      pseudo_diagnosis_datetime = if_else(
        case == 1,
        .data[[DIAGNOSIS_DATE_COL]],
        pseudo_diagnosis_datetime
      )
    )

  ## -------------------------------------------------------------------------
  ## Save this cancer type's match object + matched data
  ## (written immediately so a later failure doesn't lose
  ## already-completed types)
  ## -------------------------------------------------------------------------
  saveRDS(
    m.out,
    file.path(dir_cancer_specific, sprintf("match_object_%s.rds", ct_slug))
  )
  saveRDS(
    matched_ct,
    file.path(dir_cancer_specific, sprintf("matched_cohort_%s.rds", ct_slug))
  )
  write_tsv(
    matched_ct,
    file.path(dir_cancer_specific, sprintf("psm_matched_cohort_%s.tsv", ct_slug))
  )

  message(sprintf("Saved cancer-specific outputs for '%s' (slug: %s)", ct, ct_slug))

  # remove the cancer type match object to clear memory; can combine later
  rm(m.out)
  gc()

  matched_ct
}

## mc.preschedule = FALSE hands out one type at a time: the types differ hugely
## in size, and fixed up-front blocks would leave workers idle behind whichever
## block drew the big ones. Worker messages interleave in the log as a result.
matched_list <- if (PSM_CORES > 1L) {
  mclapply(cancer_types, match_one_cancer_type,
           mc.cores = PSM_CORES, mc.preschedule = FALSE)
} else {
  lapply(cancer_types, match_one_cancer_type)
}
names(matched_list) <- cancer_types

## A worker that raises comes back as a try-error rather than stopping the run,
## so it has to be checked for explicitly.
crashed <- vapply(matched_list, function(x) inherits(x, "try-error"), logical(1))
if (any(crashed)) {
  for (ct in cancer_types[crashed]) {
    warning(sprintf("Worker for cancer_type_grouped = '%s' failed: %s",
                    ct, conditionMessage(attr(matched_list[[ct]], "condition"))))
  }
}

## Drop crashes and the NULLs returned by a failed matchit(), exactly as the
## serial loop's `next` did.
matched_list <- matched_list[!crashed]
matched_list <- matched_list[!vapply(matched_list, is.null, logical(1))]

## Serially, a worker killed outright (the OOM killer, say) took the whole run
## down with it. A forked one does not: mclapply just returns NULL for it, which
## is indistinguishable here from a matchit() that declined to match. Either way
## the type is missing from the combined table, so say which ones and how many,
## rather than leaving it to be inferred from the counts printed below.
dropped <- setdiff(cancer_types, names(matched_list))
if (length(dropped) > 0) {
  warning(sprintf(
    "%d of %d cancer types produced no matched cohort and are absent from the combined table: %s",
    length(dropped), length(cancer_types), paste(dropped, collapse = ", ")))
}


## ---------------------------------------------------------------------------
## 4. Combine all per-cancer-type matched datasets
##    Note: the same control can legitimately appear in more than one cancer
##    type's matched set (it was matched independently each time). If you need
##    controls to be exclusive to a single cancer type, that requires a
##    different design (e.g. splitting the control pool up front, or
##    de-duplicating post hoc).
## ---------------------------------------------------------------------------
matched <- bind_rows(matched_list)


## ---------------------------------------------------------------------------
## 5. Sanity checks
## ---------------------------------------------------------------------------
cat("\n=== Overall ===\n")
cat("N cancer types successfully matched:", length(matched_list), "of", length(cancer_types), "\n")
cat("N case-rows in combined matched data:   ", sum(matched$case == 1), "\n")
cat("N control-rows in combined matched data:", sum(matched$case == 0), "\n")
cat("N distinct controls used (may be reused across cancer types):",
    n_distinct(matched$person_id[matched$case == 0]), "\n")

cat("\n=== Per cancer type ===\n")
matched %>%
  count(match_cancer_type, case) %>%
  pivot_wider(names_from = case, values_from = n, names_prefix = "case_") %>%
  print(n = Inf)

cat("\nControls per case (summary), pooled across all types:\n")
print(matched %>% filter(case == 0) %>% count(subclass) %>% summary())

cat("\nControls missing a pseudo diagnosis datetime:",
    sum(matched$case == 0 & is.na(matched$pseudo_diagnosis_datetime)), "\n")


## ---------------------------------------------------------------------------
## 6. Save aggregate output (combined across all cancer types)
## ---------------------------------------------------------------------------
saveRDS(matched,   file.path(dir_aggregate, "matched_cohort_all_cancer_types.rds"))
write_tsv(matched, file.path(dir_aggregate, "psm_matched_cohort_all_cancer_types.tsv"))

message(sprintf("\nSaved aggregate outputs to %s", dir_aggregate))
message(sprintf("Saved %d sets of cancer-specific outputs to %s",
                length(matched_list), dir_cancer_specific))
