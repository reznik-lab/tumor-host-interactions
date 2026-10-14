#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 04.4 -- the landmark case-control BMI difference for every cancer
## type above the floor, at five pre-diagnosis landmarks.
##
## Produces: Supplementary Table 6 -- case-control BMI difference at five
##           pre-diagnosis landmarks, every cancer type above the floor.
##
## The p-values here are BH-adjusted over this script's own family of cells,
## which is larger than the main figure's, so they are not the main figure's.
## Reads:
##   <workspace>/.../temp/mixed_effects_analysis_input.tsv
## Writes:
##   .../results/tables/bmi_landmarks/supp_landmarks_all_types.tsv (+ .xlsx)
##       [Supplementary Table 6]
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
library(lmerTest)
library(broom.mixed)

fig_dir <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/results/figs/bmi_landmarks/supplementary")
table_dir <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/results/tables/bmi_landmarks")
dir.create(fig_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

## Per-cell floor. NOT set from MIN_CASES_PER_TYPE: that constant is the
## cohort-level case floor for reporting a cancer type at all, while this gates
## how many cases must sit inside one landmark window for that single cell to
## be fitted. They answer different questions and should move independently --
## the same reasoning the main script's per-timebin gate is written under.
MIN_CASES_PER_CELL <- 20


## --- load ------------------------------------------------------------------

bmi_long <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/temp/mixed_effects_analysis_input.tsv"))

## control = reference, so the fixed effect is named "casecase" and its sign
## reads as case minus control. The main script does the same; done here too
## because read_tsv() returns a character column with no level order at all,
## and alphabetical order would silently make "case" the reference and flip
## every estimate on both figures.
bmi_long <- bmi_long %>%
  filter(case %in% c("case", "control")) %>%
  mutate(case = factor(case, levels = c("control", "case")))

stopifnot(identical(levels(bmi_long$case), c("control", "case")))


## --- which cancer types clear the cohort-level floor -----------------------
## Counted on CASES only, never cases + controls pooled: at 10:1 matching a
## pooled gate passes a type on a tenth of the cases it looks like.

types_ok <- bmi_long %>%
  group_by(match_cancer_type) %>%
  summarise(n_cases = n_distinct(person_id[case == "case"]), .groups = "drop") %>%
  filter(n_cases >= MIN_CASES_PER_TYPE) %>%
  arrange(desc(n_cases))

message(sprintf("%d cancer types with >= %d cases:", nrow(types_ok), MIN_CASES_PER_TYPE))
print(types_ok, n = Inf)

stopifnot(nrow(types_ok) > 0)


## --- one row per (cancer type x landmark) cell, with its n -----------------

cells <- map_dfr(seq_len(nrow(LANDMARK_TIMEPOINTS)), function(i) {
  tp <- LANDMARK_TIMEPOINTS$label[i]
  bmi_long %>%
    filter(match_cancer_type %in% types_ok$match_cancer_type,
           days_from_dx >= LANDMARK_TIMEPOINTS$window_lo[i],
           days_from_dx <= LANDMARK_TIMEPOINTS$window_hi[i]) %>%
    mutate(timepoint = tp)
})

cell_n <- cells %>%
  group_by(match_cancer_type, timepoint) %>%
  summarise(
    n_cases      = n_distinct(person_id[case == "case"]),
    n_controls   = n_distinct(person_id[case == "control"]),
    n_obs        = n(),
    n_subclass   = n_distinct(subclass),
    .groups = "drop"
  ) %>%
  ## every type x landmark combination exists as a row, so a cell that is
  ## EMPTY is visibly empty rather than quietly absent from the grid
  complete(match_cancer_type = types_ok$match_cancer_type,
           timepoint = LANDMARK_TIMEPOINTS$label,
           fill = list(n_cases = 0, n_controls = 0, n_obs = 0, n_subclass = 0)) %>%
  mutate(fitted = n_cases >= MIN_CASES_PER_CELL & n_controls >= 2)

message(sprintf("%d of %d cells clear >= %d cases and will be fitted.",
                sum(cell_n$fitted), nrow(cell_n), MIN_CASES_PER_CELL))


## --- fit each cell ---------------------------------------------------------
## One model per cell. A cell that fails to converge is reported as not fitted
## rather than dropped silently: `note` carries the reason, so a blank cell on
## the figure can always be traced to either too few cases or a failed fit.

fit_cell <- function(ct, tp) {
  d <- cells %>% filter(match_cancer_type == ct, timepoint == tp)
  note <- NA_character_

  fit <- withCallingHandlers(
    tryCatch(
      lmerTest::lmer(bmi ~ case + (1 | subclass) + (1 | person_id),
                     data = d, REML = TRUE),
      error = function(e) { note <<- paste("error:", conditionMessage(e)); NULL }),
    warning = function(w) {
      note <<- paste(na.omit(c(note, paste("warning:", conditionMessage(w)))),
                     collapse = "; ")
      invokeRestart("muffleWarning")
    })
  if (is.null(fit)) return(tibble(note = note))

  out <- broom.mixed::tidy(fit, effects = "fixed") %>%
    filter(term == "casecase") %>%
    transmute(estimate, std.error, df, statistic, p.value)

  if (nrow(out) != 1) {
    return(tibble(note = paste(na.omit(c(note, "no casecase term")), collapse = "; ")))
  }
  out %>% mutate(note = note)
}

stopifnot(any(cell_n$fitted))

res <- cell_n %>%
  filter(fitted) %>%
  mutate(fit = map2(match_cancer_type, timepoint, fit_cell)) %>%
  unnest(fit)

## unnest() only creates the columns the fits actually returned: if every fit
## failed, fit_cell() returned `note` alone and `estimate` does not exist, so
## the select() below would fail with a missing-column error rather than an
## empty figure. Make the shape of this frame independent of how the fits went.
for (col in c("estimate", "std.error", "df", "statistic", "p.value")) {
  if (!col %in% names(res)) res[[col]] <- NA_real_
}
if (!"note" %in% names(res)) res$note <- NA_character_

## Cells that were gated out, or that failed, carried through with NA estimates
## so the grid stays complete.
res <- cell_n %>%
  left_join(res %>% select(match_cancer_type, timepoint,
                           estimate, std.error, df, statistic, p.value, note),
            by = c("match_cancer_type", "timepoint"))

n_failed <- sum(res$fitted & is.na(res$estimate))
if (n_failed > 0) {
  warning(sprintf("%d cell(s) cleared the case floor but produced no estimate; see `note` in the TSV.",
                  n_failed))
}


## --- BH, over every cell THIS script fitted --------------------------------
## Its own family -- see the header. Recorded in the table so the number behind
## the correction travels with it instead of living only in this comment.

bh_family_size <- sum(!is.na(res$p.value))

res <- res %>%
  mutate(
    padj = p.adjust(p.value, method = "BH"),
    bh_family_size = bh_family_size,
    sig = !is.na(padj) & padj < 0.05
  )

message(sprintf("BH family: %d fitted cells across %d cancer types x %d landmarks.",
                bh_family_size, nrow(types_ok), nrow(LANDMARK_TIMEPOINTS)))


## --- the 95% interval ------------------------------------------------------
## Computed once, here, so the supplementary table and the forest panel cannot
## drift apart: the figure below reads these columns rather than recomputing
## them from std.error a second time.
##
## Wald, matching me_timepoint_landmarks_diff.pdf. The p-value is Satterthwaite
## and the interval is not, so the supplementary legend must not call them the
## same approximation.

res <- res %>%
  mutate(
    conf_low  = estimate - 1.96 * std.error,
    conf_high = estimate + 1.96 * std.error
  )


## --- order and labels ------------------------------------------------------
## Rows ordered by the 48-60 month estimate, so the figure opens on the types
## that are most elevated furthest from diagnosis -- the claim's starting
## point. A type not fitted there sorts last rather than being dropped.

type_order_raw <- res %>%
  filter(timepoint == "5_years") %>%
  mutate(sort_key = coalesce(estimate, -Inf)) %>%
  arrange(desc(sort_key)) %>%
  pull(match_cancer_type)
type_order_raw <- c(type_order_raw,
                    setdiff(types_ok$match_cancer_type, type_order_raw))

res <- res %>%
  mutate(
    type_label = factor(display_type_label(match_cancer_type),
                        levels = display_type_label(type_order_raw)),
    landmark   = factor(LANDMARK_LABELS[timepoint],
                        levels = unname(LANDMARK_LABELS))
  )

## display_type_label() must not collapse two raw types onto one label, or rows
## would merge here. It is asserted injective in setup.R; asserted again
## against THIS cohort's types, since only a subset is ever present.
stopifnot(!any(is.na(res$type_label)), !any(is.na(res$landmark)),
          n_distinct(res$type_label) == n_distinct(res$match_cancer_type))

## --- the supplementary table -----------------------------------------------
## write_table() writes both .tsv and .xlsx; the .xlsx is the supplementary
## file, the .tsv the one anything downstream would read.
##
## The columns are the ones the table is meant to show, in reading order:
## which cell, how much data stood behind it, the effect with its interval,
## and the evidence. What is NOT here is recoverable or intermediate --
## std.error from the interval, `sig` from p_adj, `fitted` from `note` --
## with two exceptions kept deliberately:
##
##   note            why a cell that cleared the case floor still has no
##                   estimate. Nothing else records it, and without it a blank
##                   cell on the figure cannot be told from a failed fit.
##   bh_family_size  the family the correction was taken over. It is constant
##                   down the column, but these p_adj values mean nothing
##                   without it -- and mixed_effects_modelling.R reports a
##                   DIFFERENT family for some of the same cells (see header).
##
## cancer_type is the display label, which display_type_label() is asserted
## injective on above, so no two raw types collapse onto one row.

write_table(res %>%
            transmute(
              cancer_type = type_label,
              timepoint,
              landmark,
              n_cases, n_controls, n_obs, n_subclass,
              estimate,
              conf_low, conf_high,
              p_value = p.value,
              p_adj   = padj,
              bh_family_size,
              note) %>%
            arrange(cancer_type, landmark),
          file.path(table_dir, "supp_landmarks_all_types.tsv"))


## --- what the run did ------------------------------------------------------

message("\n-- cells by landmark --")
print(res %>%
        group_by(landmark) %>%
        summarise(types_fitted = sum(!is.na(estimate)),
                  types_gated  = sum(!fitted),
                  sig_padj     = sum(sig, na.rm = TRUE),
                  .groups = "drop"))

message("\n-- significant cells, strongest first --")
print(res %>%
        filter(sig) %>%
        arrange(padj) %>%
        select(type_label, landmark, n_cases, estimate, p.value, padj),
      n = Inf)
