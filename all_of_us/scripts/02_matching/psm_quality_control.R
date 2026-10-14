#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 02.2 -- covariate balance of the propensity score match.
##
## Produces: Supplementary Figure 1b -- each cancer type's worst covariate
##           imbalance, before and after matching.
##           Supplementary Table 4 -- the same, per covariate.
##
## Reads one matchit object per cancer type and reports, for each covariate,
## the standardized mean difference before and after matching.
##
## Reads:
##   <workspace>/.../stable_data/psm/solid_heme/full_cohort/cancer_specific/match_object_<type>.rds
## Writes:
##   <workspace>/.../results/tables/psm/cancer_specific/supp_fig1a_smd_across_cancer_types.tsv (+ .xlsx)
##       [Supplementary Table 4]
##   <workspace>/.../results/figs/psm/cancer_specific/supp_fig1a_max_abs_smd_by_cancer_type.pdf (+ .rds)
##       [Supplementary Figure 1b]
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
library(cobalt)      # balance tables + love plots
library(ggthemes)    # tableau_color_pal("Tableau20")
library(openxlsx)

## theme_std() comes from setup.R. tableau20 is local to this script.
tableau20 <- ggthemes::tableau_color_pal("Tableau 20")(20)

## helper -- same slugification used when the PSM outputs were saved, so
## filenames stay consistent across the two scripts
sanitize_name <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}


## ---------------------------------------------------------------------------
## 1. Input locations (written by psm.R) and output locations (this script
##    writes QC results, so these go under results/, not stable_data/)
## ---------------------------------------------------------------------------

## These must stay in step with dir_cancer_specific / dir_aggregate in psm.R.
psm_cancer_specific_dir <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm/solid_heme/full_cohort/cancer_specific")
psm_aggregate_dir <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm/solid_heme/full_cohort/aggregate")

## results/figs/<analysis>/ and results/tables/<analysis>/ -- the layout every
## other script uses, so every figure in the pipeline sits under one parent.
results_root <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/results")
figs_dir     <- file.path(results_root, "figs",   "psm", "cancer_specific")
tabs_dir     <- file.path(results_root, "tables", "psm", "cancer_specific")
figs_dir_agg <- file.path(results_root, "figs",   "psm", "aggregate")
tabs_dir_agg <- file.path(results_root, "tables", "psm", "aggregate")

dir.create(figs_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(tabs_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(figs_dir_agg, showWarnings = FALSE, recursive = TRUE)
dir.create(tabs_dir_agg, showWarnings = FALSE, recursive = TRUE)


## ---------------------------------------------------------------------------
## 2. Find the saved matchit objects
## ---------------------------------------------------------------------------

## Find the saved matchit objects. The same directory also holds
## matched_cohort_*.rds (plain data frames), so match the prefix rather than
## every .rds -- bal.tab() has nothing to do with those.
matchit_files <- list.files(
  psm_cancer_specific_dir,
  pattern = "^match_object_.*\\.rds$",
  full.names = TRUE,
  recursive = TRUE
)
if (length(matchit_files) == 0) {
  stop("No matchit .rds files found in ", psm_cancer_specific_dir)
}

## Derive a cancer-type label from each filename
cancer_labels <- matchit_files %>%
  basename() %>%
  str_remove("^match_object_") %>%
  str_remove("\\.rds$")

## Output dir for per-cancer-type balance workbooks

## The per-type loop that used to live here wrote a workbook for EVERY cancer
## type, including ones too small to report. It now runs in section 4, gated on
## MIN_CASES_PER_TYPE and alongside that type's love plot, so the workbooks and
## the plots always describe the same set of cancer types.


## ---------------------------------------------------------------------------
## 3. Readable covariate labels, and the love plot itself
## ---------------------------------------------------------------------------
## cobalt expands each factor covariate to one row per level, named
## "<var>_<level>", where <level> is the raw value stored in the data: an
## OMOP/AoU concept id for sex_at_birth and self_reported_category, a calendar
## year for earliest_bmi_year. Printed as-is those codes are unreadable, and
## earliest_bmi_year's ~45 levels crowd out every other covariate. This section
## fixes both, then draws one plot per cancer type.

## The factors psm.R matches on, longest name first so that prefix matching
## below cannot mis-split a covariate whose name starts with another's.
PSM_FACTOR_VARS <- c("self_reported_category", "earliest_bmi_year", "sex_at_birth")

## Factors collapsed to a single row. earliest_bmi_year is an EXACT matching
## variable (psm.R's `exact = ~ ...`), so every one of its adjusted SMDs is 0 by
## construction -- the 45 rows cost two thirds of the y axis and carry no
## post-match information. The collapsed row keeps the single worst level, so a
## real pre-match imbalance is still visible rather than averaged away.
PSM_COLLAPSE_VARS <- c("earliest_bmi_year")

## Labels for covariates that are not factor expansions. Anything absent here
## is title-cased from its own name rather than left in snake_case.
PSM_COVARIATE_LABELS <- c(
  distance           = "Propensity score",
  age_at_enrollment  = "Age at enrollment",
  enrollment_year    = "Enrollment year",
  deprivation_index  = "Deprivation index",
  n_bmi_measurements = "BMI measurements per person"
)

## Readable name for each factor, used as the prefix on its level labels.
PSM_FACTOR_LABELS <- c(
  self_reported_category = "Race/ethnicity",
  sex_at_birth           = "Sex at birth",
  earliest_bmi_year      = "Earliest BMI year"
)


## --- concept id -> display name --------------------------------------------
## Built from the same files psm.R matched on, which carry AoU's own display
## column beside each coded column, rather than from a hard-coded crosswalk: a
## guessed concept id would mislabel a race or sex category silently, and a
## mislabelled category on a balance figure is worse than a visible code.
##
## Candidate display-column names, tried in order; the first present wins. If
## none is present the code is kept and the message below says which variables
## fell back, so a missing crosswalk is always visible.
display_col_candidates <- function(var) {
  c(paste0("t_disp_", var),
    paste0(var, "_name"),
    paste0(var, "_concept_name"),
    paste0("t_disp_", var, "_concept_id"))
}

## psm.R's inputs. Kept in step with the read_tsv() calls at the top of psm.R.
psm_input_files <- c(
  file.path(WORKSPACE,
    "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/pre_psm_matched_cancer_bmi_summary_prediagnosis.tsv"),
  file.path(WORKSPACE,
    "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/control_individuals/pre_psm_matched_metadata.tsv")
)

#' Crosswalk from stored code to display name, for one covariate.
#'
#' Reads only the two columns it needs -- the headers are inspected first, so a
#' table lacking a display column costs one header read rather than a full load.
#'
#' @return a named character vector (names = codes, values = labels), or NULL
#'   when no source table carries a display column for `var`.
build_level_lookup <- function(var, paths) {
  pairs <- NULL

  for (p in paths) {
    if (!file.exists(p)) {
      warning("PSM input not found, skipping for label lookup: ", p)
      next
    }
    hdr <- names(read_tsv(p, n_max = 0, show_col_types = FALSE, progress = FALSE))
    if (!(var %in% hdr)) next

    disp <- intersect(display_col_candidates(var), hdr)
    if (length(disp) == 0) next

    d <- read_tsv(p, col_select = all_of(c(var, disp[1])),
                  show_col_types = FALSE, progress = FALSE)
    names(d) <- c("code", "label")
    pairs <- bind_rows(pairs, d)
  }

  if (is.null(pairs)) return(NULL)

  pairs <- pairs %>%
    mutate(code = as.character(code), label = as.character(label)) %>%
    filter(!is.na(code), !is.na(label), label != "") %>%
    distinct()

  ## One code must mean one thing. If a code carries two display names the
  ## crosswalk is not trustworthy for that code, so keep the code itself.
  ambiguous <- pairs %>% count(code) %>% filter(n > 1) %>% pull(code)
  if (length(ambiguous) > 0) {
    warning(var, ": ", length(ambiguous),
            " code(s) map to more than one display name; left as codes.")
    pairs <- pairs %>% filter(!code %in% ambiguous)
  }
  if (nrow(pairs) == 0) return(NULL)

  setNames(pairs$label, pairs$code)
}

level_lookups <- set_names(PSM_FACTOR_VARS) %>%
  map(build_level_lookup, paths = psm_input_files)

resolved <- names(level_lookups)[!map_lgl(level_lookups, is.null)]
unresolved <- setdiff(PSM_FACTOR_VARS, resolved)
message("Label crosswalk resolved for: ",
        if (length(resolved)) paste(resolved, collapse = ", ") else "(none)")
if (length(unresolved) > 0) {
  message("  NO display column found for: ", paste(unresolved, collapse = ", "),
          " -- these keep their raw codes on the y axis. Add the display ",
          "column name to display_col_candidates() if the source tables have one.")
}


## --- covariate name -> (variable, level) -----------------------------------

#' Split a cobalt covariate name into its variable and, for a factor
#' expansion, its level. Non-factor covariates come back with level = NA.
split_covariate <- function(cov) {
  var   <- rep(NA_character_, length(cov))
  level <- rep(NA_character_, length(cov))

  for (v in PSM_FACTOR_VARS) {
    hit <- is.na(var) & startsWith(cov, paste0(v, "_"))
    var[hit]   <- v
    level[hit] <- substring(cov[hit], nchar(v) + 2)
  }
  var[is.na(var)] <- cov[is.na(var)]
  tibble(covariate = cov, var = var, level = level)
}

#' Human-readable label for one covariate row.
pretty_label <- function(var, level) {
  ifelse(
    is.na(level),
    ## not a factor expansion
    ifelse(var %in% names(PSM_COVARIATE_LABELS),
           PSM_COVARIATE_LABELS[var],
           str_to_sentence(str_replace_all(var, "_", " "))),
    ## factor expansion: "<Factor>: <level name>"
    paste0(
      ifelse(var %in% names(PSM_FACTOR_LABELS), PSM_FACTOR_LABELS[var],
             str_to_sentence(str_replace_all(var, "_", " "))),
      ": ",
      map2_chr(var, level, function(v, l) {
        lk <- level_lookups[[v]]
        if (!is.null(lk) && !is.na(l) && l %in% names(lk)) lk[[l]] else l
      })
    )
  )
}


## --- build the plot frame --------------------------------------------------

## ---------------------------------------------------------------------------
## 3b. The numbers behind the balance table
## ---------------------------------------------------------------------------
## The love plot shows the SMD before and after matching. The table that
## accompanies it has to show what those SMDs are a difference OF: the mean
## (SD) of each continuous covariate and the percentage in each factor level,
## for cases and for controls, in the unmatched pool and in the matched sample.
##
## All of it comes from the one bal.tab() call that already feeds the plot,
## with disp = c("means", "sds") added, so the table and the figure cannot
## disagree -- same matchit object, same weights, same covariate expansion.
## The post-match means are weighted by the matching weights, which is what
## makes them the right comparison for psm.R's variable-ratio 1:k match.

#' Which treatment level is which, read off cobalt's own column names.
#'
#' cobalt names its mean/SD columns after the treatment values found in the
#' data -- "M.0.Un"/"M.1.Un" for psm.R's 0/1 `case`, but "M.Control.Un"/
#' "M.Case.Un" had the treatment been a factor -- and emits them control-first,
#' the same orientation as its Diff = treated - control. Reading the levels
#' back off the names rather than hard-coding "0"/"1" keeps a relabelled
#' treatment from silently swapping the case and control columns.
balance_group_levels <- function(bal) {
  hits <- grep("^M\\..+\\.Un$", colnames(bal), value = TRUE)
  if (length(hits) != 2) {
    stop("expected 2 unadjusted mean columns from cobalt, found ", length(hits),
         " -- was bal.tab() called with disp = c(\"means\", \"sds\")?")
  }
  lv <- sub("^M\\.(.+)\\.Un$", "\\1", hits)
  ## psm.R codes case = 1, control = 0. When that is what we are looking at,
  ## pin the orientation outright rather than trusting column order alone.
  if (setequal(lv, c("0", "1"))) lv <- c("0", "1")
  list(control = lv[1], treated = lv[2])
}

#' bal.tab() asked for the arm means and SDs, across cobalt's two spellings.
#'
#' cobalt 4 takes disp = c("means", "sds"); cobalt 3 spelled the same request
#' disp.means / disp.sds, and absorbs an unknown `disp` into its `...` without
#' complaint -- which would leave every mean column out of the table and every
#' cancer type failing in balance_group_levels(). Ask the modern way, check the
#' columns actually arrived, and fall back rather than depend on the version.
bal_tab_with_means <- function(m) {
  bt <- cobalt::bal.tab(m, un = TRUE, disp = c("means", "sds"))
  has_means <- any(grepl("^M\\..+\\.Un$",
                         colnames(as.data.frame(bt[["Balance"]]))))
  if (has_means) return(bt)
  cobalt::bal.tab(m, un = TRUE, disp.means = TRUE, disp.sds = TRUE)
}

#' One statistic column out of a cobalt Balance table; NA if cobalt did not
#' emit it (an older version, or a `disp` that left the SDs out).
bal_stat <- function(bal, stat, lvl, sample) {
  nm <- paste(stat, lvl, sample, sep = ".")
  if (nm %in% colnames(bal)) as.numeric(bal[[nm]]) else rep(NA_real_, nrow(bal))
}

#' Display form of one cell. A binary covariate's "mean" is the proportion of
#' the arm in that level, so it reads as a percentage; a continuous one reads
#' as mean (SD).
fmt_balance_value <- function(type, m, s) {
  ifelse(is.na(m), NA_character_,
         ifelse(type == "Binary",
                sprintf("%.1f%%", 100 * m),
                ifelse(is.na(s), sprintf("%.2f", m),
                       sprintf("%.2f (%.2f)", m, s))))
}

#' Flag a percentage whose implied cell count -- or that count's complement --
#' could fall in 1-20, the range AoU's dissemination policy forbids exporting.
#' A percentage sitting beside its denominator IS such a count, so the table
#' carries the flag rather than leaving it to be worked out by hand.
#'
#' The post-match percentages are weighted, so the implied count is not exact;
#' floor/ceiling bracket it and the smaller side is taken, which errs towards
#' flagging. Nothing is blanked here -- the table is written for use inside the
#' workbench, and which rows need suppressing before it leaves is a judgement
#' for whoever exports it.
small_cell <- function(type, prop, n) {
  k_lo <- floor(prop * n)
  k_hi <- ceiling(prop * n)
  cnt  <- pmin(k_lo, n - k_hi)
  !is.na(prop) & !is.na(n) & type == "Binary" & cnt >= 1 & cnt <= 20
}

#' Per-covariate balance for one cancer type: means/percentages for each arm
#' and the SMD, before and after matching, one row per covariate (factors
#' expanded to one row per level, as cobalt expands them).
#'
#' @param bal_tab a cobalt bal.tab() result built with
#'   `un = TRUE, disp = c("means", "sds")`.
#' @param stats that type's row from match_stats(), for the denominators.
balance_detail <- function(bal_tab, cancer_type, stats) {
  bal <- as.data.frame(bal_tab[["Balance"]])
  lv  <- balance_group_levels(bal)

  cov_type <- if ("Type" %in% colnames(bal)) as.character(bal$Type) else NA_character_
  parts    <- split_covariate(rownames(bal))

  tibble(
    cancer_type = cancer_type,
    covariate   = rownames(bal),
    label       = pretty_label(parts$var, parts$level),
    var         = parts$var,
    level       = parts$level,
    type        = cov_type,

    cases_pre_mean     = bal_stat(bal, "M", lv$treated, "Un"),
    cases_pre_sd       = bal_stat(bal, "S", lv$treated, "Un"),
    controls_pre_mean  = bal_stat(bal, "M", lv$control, "Un"),
    controls_pre_sd    = bal_stat(bal, "S", lv$control, "Un"),
    smd_pre_match      = as.numeric(bal$Diff.Un),

    cases_post_mean    = bal_stat(bal, "M", lv$treated, "Adj"),
    cases_post_sd      = bal_stat(bal, "S", lv$treated, "Adj"),
    controls_post_mean = bal_stat(bal, "M", lv$control, "Adj"),
    controls_post_sd   = bal_stat(bal, "S", lv$control, "Adj"),
    smd_post_match     = as.numeric(bal$Diff.Adj),

    n_cases_pre        = stats$n_cases,
    n_controls_pre     = stats$n_controls,
    n_cases_post       = stats$n_cases_matched,
    n_controls_post    = stats$n_controls_matched
  ) %>%
    mutate(
      cases_pre     = fmt_balance_value(type, cases_pre_mean,     cases_pre_sd),
      controls_pre  = fmt_balance_value(type, controls_pre_mean,  controls_pre_sd),
      cases_post    = fmt_balance_value(type, cases_post_mean,    cases_post_sd),
      controls_post = fmt_balance_value(type, controls_post_mean, controls_post_sd),
      low_cell_count =
        small_cell(type, cases_pre_mean,      n_cases_pre)     |
        small_cell(type, controls_pre_mean,   n_controls_pre)  |
        small_cell(type, cases_post_mean,     n_cases_post)    |
        small_cell(type, controls_post_mean,  n_controls_post),
      across(c(ends_with("_mean"), ends_with("_sd"),
               smd_pre_match, smd_post_match), ~ round(.x, 4))
    ) %>%
    ## The love plot's own row order, so figure and table read side by side:
    ## propensity score first (FALSE sorts before TRUE), then by how imbalanced
    ## the covariate started out.
    arrange(covariate != "distance", desc(abs(smd_pre_match))) %>%
    select(cancer_type, covariate, label, var, level, type,
           cases_pre, controls_pre, smd_pre_match,
           cases_post, controls_post, smd_post_match,
           everything())
}

## ---------------------------------------------------------------------------
## 4. Per-cancer-type balance, for every type with at least
##    MIN_CASES_PER_TYPE cases.
## ---------------------------------------------------------------------------


#' Case/control counts and matching ratios, read from the matchit object
#' itself rather than from cobalt's printed Observations table.
#'
#' `treat` (0/1 per unit), `weights` (0 for unmatched) and `subclass` are
#' documented elements of a matchit object, so they do not move between cobalt
#' versions the way the Observations table's row and column labels can.
match_stats <- function(m, cancer_type) {
  treat <- m$treat
  if (is.null(treat)) {
    stop("matchit object for ", cancer_type, " has no $treat element")
  }
  w <- m$weights
  if (is.null(w)) w <- rep(1, length(treat))   # no weights => everything kept
  matched <- w > 0

  sub <- m$subclass
  ## controls per matched case, one row per subclass
  per_case <- if (!is.null(sub)) {
    tibble(subclass = as.character(sub), treat = treat, matched = matched) %>%
      filter(matched, !is.na(subclass)) %>%
      count(subclass, treat) %>%
      tidyr::pivot_wider(names_from = treat, values_from = n,
                         names_prefix = "t", values_fill = 0) %>%
      transmute(n_controls = if ("t0" %in% names(.)) t0 else 0L)
  } else {
    tibble(n_controls = integer())
  }

  tibble(
    cancer_type        = cancer_type,
    n_cases            = sum(treat == 1),
    n_cases_matched    = sum(treat == 1 & matched),
    n_controls         = sum(treat == 0),
    n_controls_matched = sum(treat == 0 & matched),
    matched_sets       = nrow(per_case),
    controls_per_case_mean   = if (nrow(per_case)) round(mean(per_case$n_controls), 2) else NA_real_,
    controls_per_case_median = if (nrow(per_case)) median(per_case$n_controls) else NA_real_,
    controls_per_case_min    = if (nrow(per_case)) min(per_case$n_controls) else NA_real_,
    controls_per_case_max    = if (nrow(per_case)) max(per_case$n_controls) else NA_real_
  )
}

#' Balance summary: how many covariates clear the |SMD| < 0.1 rule of thumb,
#' before and after matching. `distance` is the propensity score, not a
#' covariate, so it is excluded from the counts.
balance_summary <- function(bal_df) {
  cov <- bal_df %>% filter(covariate != "distance")
  tibble(
    n_covariates        = nrow(cov),
    n_imbalanced_pre    = sum(abs(cov$smd_pre_match)  >= 0.1, na.rm = TRUE),
    n_imbalanced_post   = sum(abs(cov$smd_post_match) >= 0.1, na.rm = TRUE),
    max_abs_smd_pre     = suppressWarnings(max(abs(cov$smd_pre_match),  na.rm = TRUE)),
    max_abs_smd_post    = suppressWarnings(max(abs(cov$smd_post_match), na.rm = TRUE))
  ) %>%
    mutate(across(starts_with("max_abs"), ~ ifelse(is.finite(.x), round(.x, 4), NA_real_)))
}

## Each element is list(stats =, summary =, detail =): the counts for every
## type (gate passed or not), and -- for the types that passed -- the balance
## summary and the per-covariate detail section 5 stacks.
per_type <- map(seq_along(matchit_files), function(i) {
  cancer_type <- cancer_labels[i]

  out <- tryCatch({
    matchit_obj <- readRDS(matchit_files[i])

    stats <- match_stats(matchit_obj, cancer_type)

    ## THE GATE. Counted on cases, like MIN_CASES_PER_TYPE everywhere else.
    if (stats$n_cases < MIN_CASES_PER_TYPE) {
      message(sprintf("  %-28s %5d cases -- below %d, skipped",
                      cancer_type, stats$n_cases, MIN_CASES_PER_TYPE))
      return(list(stats = stats %>% mutate(drawn = FALSE)))
    }

    ## The means and SDs are what turn the balance table from SMDs alone into
    ## the mean/percentage comparison the supplementary table reports. The love
    ## plot reads the same object, so figure and table cannot drift apart.
    bal_tab <- bal_tab_with_means(matchit_obj)

    detail <- balance_detail(bal_tab, cancer_type, stats)

    ## the two SMD columns on their own, for the balance summary below
    bal_df <- detail %>%
      select(covariate, label, smd_pre_match, smd_post_match)

    type_summary <- bind_cols(stats, balance_summary(bal_df))

    list(stats   = stats %>% mutate(drawn = TRUE),
         summary = type_summary,
         detail  = detail)
  },
  error = function(e) {
    ## One bad cancer type must not cost the whole run.
    warning(cancer_type, ": ", conditionMessage(e))
    message(sprintf("  %-28s FAILED: %s", cancer_type, conditionMessage(e)))
    NULL
  })

  out
})

## tryCatch() returns NULL for a type that failed outright; compact() drops
## those, and pulling by name drops the types the gate skipped from `details`
## and `summaries` while keeping them in `results`.
per_type  <- purrr::compact(per_type)
results   <- map_dfr(per_type, "stats")
summaries <- map_dfr(per_type, "summary")
details   <- map_dfr(per_type, "detail")

if (nrow(results) == 0) {
  stop("No cancer type produced any output -- see the warnings above.")
}

## ---------------------------------------------------------------------------
## 6. Supp Fig 1A -- the balance supplementary table, drawn
## ---------------------------------------------------------------------------
## Section 4 draws one love plot per cancer type: 40-odd separate figures,
## which is the right form for checking a single match and the wrong one for
## the claim the manuscript makes, that balance holds ACROSS types. This
## section draws the stacked table from section 5 instead, so the figure and
## the supplementary table are the same numbers twice -- both built from
## `details` / `summaries` in memory, not re-read from the .tsv.
##
## Two files rather than one grid, because they are read at different widths:
##
##   a. covariate-level. One row per covariate, one point per cancer type,
##      before and after matching. The spread of the "after" cloud inside the
##      +/- 0.1 lines IS the balance claim; the "before" cloud says what the
##      match had to remove.
##   b. type-level. One row per cancer type, its single worst |SMD| before and
##      after. The panel that answers "did any type fail?" -- which the
##      covariate panel, pooling types, cannot.

if (nrow(details) == 0) {
  warning("No cancer type cleared the gate, so Supp Fig 1A was not drawn.")
} else {

  ## --- readable cancer-type labels -----------------------------------------
  ## `cancer_type` here is psm.R's filename slug ("lymphoma_non_hodgkin_b_cell"),
  ## not the raw cancer_type_grouped value display_type_label() keys on. The
  ## slug is lossy -- "&", "/" and "-" all became "_" -- so it is un-slugged
  ## against the real values, read from psm.R's own input, rather than
  ## title-cased back, which would print "Liver biliary tract".
  slug_to_label <- local({
    p <- psm_input_files[1]
    raw <- character()
    if (file.exists(p)) {
      hdr <- names(read_tsv(p, n_max = 0, show_col_types = FALSE, progress = FALSE))
      if ("cancer_type_grouped" %in% hdr) {
        raw <- read_tsv(p, col_select = all_of("cancer_type_grouped"),
                        show_col_types = FALSE, progress = FALSE) %>%
          pull(1) %>% as.character() %>% unique()
        raw <- raw[!is.na(raw) & nzchar(raw)]
      }
    }
    slug <- sanitize_name(raw)
    ## sanitize_name() is injective in practice but not by construction. If two
    ## raw groups ever slug the same, psm.R has already overwritten one's files
    ## with the other's, so drop the slug rather than label it with whichever
    ## happened to win.
    keep <- !(slug %in% slug[duplicated(slug)])
    setNames(display_type_label(raw[keep]), slug[keep])
  })

  #' Slug -> figure label. An unmatched slug keeps its title-cased form, so a
  #' renamed cancer group shows up under a readable name rather than as NA.
  unslug_type <- function(slug) {
    hit <- match(slug, names(slug_to_label))
    ifelse(is.na(hit), str_to_sentence(str_replace_all(slug, "_", " ")),
           unname(slug_to_label[hit]))
  }

  n_unresolved <- sum(!unique(details$cancer_type) %in% names(slug_to_label))
  if (n_unresolved > 0) {
    message(sprintf(
      "Supp Fig 1A: %d of %d cancer-type slugs had no cancer_type_grouped match; title-cased instead.",
      n_unresolved, n_distinct(details$cancer_type)))
  }

  ## Sample labels and colours.
  as_sample <- function(x) factor(x, levels = c("smd_pre_match", "smd_post_match"),
                                  labels = c("Unadjusted", "Adjusted"))

  ## --- a. covariate-level panel --------------------------------------------
  ## The same collapse the love plots use: earliest_bmi_year is an exact
  ## matching variable, so every one of its ~45 levels is 0 after matching and
  ## they would be 45 rows of zeros crowding out every covariate that moved.
  ## One row per type, holding that type's worst level before matching.
  collapsed <- details %>%
    filter(var %in% PSM_COLLAPSE_VARS, !is.na(level)) %>%
    group_by(cancer_type, var) %>%
    slice_max(abs(smd_pre_match), n = 1, with_ties = FALSE, na_rm = TRUE) %>%
    ungroup() %>%
    mutate(label = paste0(ifelse(var %in% names(PSM_FACTOR_LABELS),
                                 PSM_FACTOR_LABELS[var], var),
                          " (worst level)"))

  smd_by_cov <- details %>%
    filter(!(var %in% PSM_COLLAPSE_VARS & !is.na(level))) %>%
    bind_rows(collapsed) %>%
    select(cancer_type, covariate, var, label, smd_pre_match, smd_post_match)

  ## Row order: propensity score at the top, then the covariate that started
  ## out worst. Ordered on the MEDIAN across types
  ## rather than the max, so one outlying cancer type cannot promote a
  ## covariate that is fine everywhere else.
  cov_rank <- smd_by_cov %>%
    group_by(label) %>%
    summarise(is_ps          = any(covariate == "distance"),
              median_abs_pre = median(abs(smd_pre_match), na.rm = TRUE),
              .groups = "drop") %>%
    arrange(desc(is_ps), desc(median_abs_pre))
  cov_order <- cov_rank$label

  ## The layer the points are drawn around, and a table in its own right: what
  ## the stacked balance table says once it is read across types rather than
  ## down one. n_types_ge_0.1 is the column to look at -- the count of cancer
  ## types where this covariate misses the 0.1 rule of thumb.
  smd_cov_summary <- smd_by_cov %>%
    pivot_longer(c(smd_pre_match, smd_post_match),
                 names_to = "sample", values_to = "smd") %>%
    group_by(label, sample) %>%
    summarise(
      n_cancer_types = sum(!is.na(smd)),
      median_smd     = median(smd, na.rm = TRUE),
      median_abs_smd = median(abs(smd), na.rm = TRUE),
      min_smd        = suppressWarnings(min(smd, na.rm = TRUE)),
      max_smd        = suppressWarnings(max(smd, na.rm = TRUE)),
      max_abs_smd    = suppressWarnings(max(abs(smd), na.rm = TRUE)),
      n_types_ge_0.1 = sum(abs(smd) >= 0.1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric), ~ ifelse(is.finite(.x), round(.x, 4), NA_real_))) %>%
    arrange(match(label, cov_order), sample != "smd_pre_match")

  ## The figure's own numbers, written out: the supplementary table is per
  ## covariate PER TYPE, and this is the across-type summary the figure draws.
  ## as.character(), not the factor: openxlsx writes a factor's level labels
  ## but readr and openxlsx disagree often enough about factor columns that the
  ## two copies of a table are not worth risking on it.
  write_table(smd_cov_summary %>% mutate(sample = as.character(as_sample(sample))),
              file.path(tabs_dir, "supp_fig1a_smd_across_cancer_types.tsv"),
              sheet = "smd_across_types")

  ## --- b. type-level panel --------------------------------------------------
  ## The worst covariate each type has, before and after matching. The
  ## propensity score is excluded already: balance_summary() drops `distance`,
  ## which is what the match balanced ON and would flatter every row.
  max_smd_by_type <- summaries %>%
    transmute(type_label     = unslug_type(cancer_type),
              smd_pre_match  = max_abs_smd_pre,
              smd_post_match = max_abs_smd_post) %>%
    pivot_longer(c(smd_pre_match, smd_post_match),
                 names_to = "sample", values_to = "max_abs_smd") %>%
    mutate(sample = as_sample(sample))

  type_order_fig <- max_smd_by_type %>%
    filter(sample == "Unadjusted") %>%
    arrange(max_abs_smd) %>%
    pull(type_label)

  p_max_smd_by_type <- max_smd_by_type %>%
    mutate(type_label = factor(type_label, levels = type_order_fig)) %>%
    ggplot(aes(x = max_abs_smd, y = type_label, colour = sample)) +
    geom_vline(xintercept = 0.1, linetype = "dashed", colour = "black",
               linewidth = 0.25) +
    ## the connector makes the before/after pair one row to the eye, so the
    ## distance travelled is what the row shows
    geom_line(aes(group = type_label), colour = "grey70", linewidth = 0.2,
              orientation = "y") +
    geom_point(size = 1.4) +
    scale_colour_manual(values = tableau20[c(2, 1)], name = "Sample") +
    labs(
      x = "Largest |SMD| over the type's covariates", y = NULL,
      title = "Worst-balanced covariate, by cancer type",
      subtitle = sprintf("%d of %d types have every covariate under 0.1 after matching",
                         sum(summaries$n_imbalanced_post == 0, na.rm = TRUE),
                         nrow(summaries))
    ) +
    theme_std() +
    theme(legend.position = "bottom")

  ggsave_wrapper(file.path(figs_dir, "supp_fig1a_max_abs_smd_by_cancer_type.pdf"),
                 p_max_smd_by_type, width = 5,
                 height = max(3, 1.3 + 0.13 * nrow(summaries)))

  worst_post <- suppressWarnings(max(summaries$max_abs_smd_post, na.rm = TRUE))
  message("\nSupplementary Figure 1b:")
  message("  ", file.path(figs_dir, "supp_fig1a_max_abs_smd_by_cancer_type.pdf"))
  message(sprintf("  %d covariate rows x %d cancer types; worst post-match |SMD| anywhere = %s",
                  n_distinct(smd_by_cov$label), n_distinct(smd_by_cov$cancer_type),
                  if (is.finite(worst_post)) sprintf("%.3f", worst_post) else "NA"))
}
