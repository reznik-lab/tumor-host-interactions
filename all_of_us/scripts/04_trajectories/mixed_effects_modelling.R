#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 04.3 -- landmark mixed-effects models, drawn as case/control boxplots.
##
## Produces: Figure 1c-f, the landmark boxplot beside each trajectory --
##           1c stomach, 1d pancreas, 1e kidney, 1f endometrial.
##           The trajectory and sample-size row come from analysis.R.
##
## BMI ~ case + (1 | person_id) + (1 | subclass), fitted per cancer type and
## landmark window; p-values BH-adjusted over every cell fitted.
## Reads:
##   <workspace>/.../temp/mixed_effects_analysis_input.tsv
## Writes:
##   .../results/figs/bmi_landmarks/bmi_diff_boxplot_<type>.pdf (+ .rds)
##       [Figure 1c-f, the boxplot beside each trajectory]
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
library(ggthemes)

landmark_dir <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/results/figs/bmi_landmarks")
## This script used to write figures only, so the landmark estimates quoted in
## the Results existed nowhere as a table. They do now.
table_dir <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/results/tables/bmi_landmarks")
dir.create(landmark_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir,    recursive = TRUE, showWarnings = FALSE)


## --- load ------------------------------------------------------------------

bmi_long_filtered <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/temp/mixed_effects_analysis_input.tsv"))

print(head(bmi_long_filtered))
print(unique(bmi_long_filtered$match_cancer_type))


## --- landmark windows ------------------------------------------------------

## MOVED to src/setup.R as LANDMARK_TIMEPOINTS: supp_landmarks_all_types.R
## windows the same measurements the same way, and a second literal copy of the
## table is exactly the kind of drift this repo keeps having to fix. The month
## ranges each window really covers, and why "5_years" is the 48-60 month
## window, are documented there.
timepoints <- LANDMARK_TIMEPOINTS


## --- cohort restriction ----------------------------------------------------

# Restrict to cancer types with enough individuals overall
suff_cancer_types <- bmi_long_filtered %>%
  group_by(match_cancer_type) %>%
  summarize(n_cases = n_distinct(person_id[case == "case"]), .groups = "drop") %>%
  filter(n_cases >= MIN_CASES_PER_TYPE)

# Count cases per (cancer type, timepoint) using the same windowing logic as `res`
timebin_counts <- map_dfr(1:nrow(timepoints), function(i) {
  bmi_long_filtered %>%
    filter(days_from_dx >= timepoints$window_lo[i],
           days_from_dx <= timepoints$window_hi[i],
           case == "case") %>%
    mutate(timepoint = timepoints$label[i])
}) %>%
  group_by(match_cancer_type, timepoint) %>%
  summarize(n_cases = n_distinct(person_id), .groups = "drop")

# A cancer type qualifies only if it has >= 20 cases in EVERY timepoint
suff_by_timebin <- timebin_counts %>%
  complete(match_cancer_type, timepoint = timepoints$label, fill = list(n_cases = 0)) %>%
  group_by(match_cancer_type) %>%
  summarize(min_n = min(n_cases), .groups = "drop") %>%
  filter(min_n >= 20)

# Apply both the overall MIN_CASES_PER_TYPE filter and the per-timebin n>=20 filter
bmi_long_filtered <- bmi_long_filtered %>%
  filter(match_cancer_type %in% suff_cancer_types$match_cancer_type) %>%
  filter(match_cancer_type %in% suff_by_timebin$match_cancer_type) %>%
  filter(case %in% c("case", "control")) %>%
  mutate(case = factor(case, levels = c("control", "case")))  # control = reference

print(timebin_counts)
print(levels(bmi_long_filtered$case))


## --- landmark models -------------------------------------------------------

res <- map_dfr(1:nrow(timepoints), function(i) {
  bmi_long_filtered %>%
    filter(days_from_dx >= timepoints$window_lo[i],
           days_from_dx <= timepoints$window_hi[i]) %>%
    mutate(timepoint = timepoints$label[i])
}) %>%
  group_by(match_cancer_type, timepoint) %>%
  group_modify(~ {
    fit <- lmerTest::lmer(bmi ~ case + (1 | subclass) + (1 | person_id),
                          data = .x, REML = TRUE)
    broom.mixed::tidy(fit, effects = "fixed")
  })

print(res %>% filter(term == "casecase", p.value < 0.05))


## --- per-person averaged BMI within each landmark window -------------------

bmi_fitered_subset <- bmi_long_filtered %>%
  filter(grepl("pancrea|stomach|kidney|uter", match_cancer_type, ignore.case = TRUE))

# 1. Assign each raw measurement to a timepoint window, then average per person
#    within that window (so a person with several ~3yr measurements gets
#    collapsed to one averaged BMI for "3_years")
bmi_avg_by_tp <- map_dfr(timepoints$label, function(tp) {
  lo <- timepoints$window_lo[timepoints$label == tp]
  hi <- timepoints$window_hi[timepoints$label == tp]

  bmi_fitered_subset %>%
    filter(days_from_dx >= lo, days_from_dx <= hi) %>%
    group_by(person_id, case, subclass, match_cancer_type) %>%
    summarize(bmi = mean(bmi, na.rm = TRUE), .groups = "drop") %>%
    mutate(timepoint = tp)
}) %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("5_years", "3_years", "1_year", "6_months", "at_diagnosis")))

print(colnames(res))


## ===========================================================================
## BOXPLOTS -- Tukey five-number summary, two timepoints.  [Figure 1c-f]
## ===========================================================================

# ---- helper: format BH-adjusted p-values ----
## Formatting comes from fmt_p() in src/setup.R -- the paper-wide convention
## (italic capital P, "adj" subscripted, exact value, m x 10^e below 0.001).
## It returns plotmath source, so the geom_text() drawing it passes parse = TRUE.

# facet label lookup: raw timepoint value -> display label
timepoint_labels <- c(
  "5_years"      = "48-60 months",
  "at_diagnosis" = "At Diagnosis"
)

## Pull the case-vs-control fixed-effect row per cancer type x timepoint.
##
## Labels come from res_plot, which is BH-adjusted over every (cancer type x
## timepoint) case effect fitted -- not just the cells drawn below, so the
## adjustment does not depend on the figure. "LME estimate" names the `case`
## fixed effect; the boxes are a Tukey summary of per-person window means.
model_labels <- res_plot %>%
  filter(as.character(timepoint) %in% names(timepoint_labels)) %>%
  transmute(
    match_cancer_type,
    timepoint = as.character(timepoint),
    label = fmt_p(
      padj, adjusted = TRUE,
      prefix = paste0("LME estimate = ", formatC(estimate, digits = 2, format = "f"), ",")
    )
  )

# ---- summarize raw BMI into standard Tukey five-number summary per group ----
# replicates exactly what geom_boxplot() computes internally, so the
# rendered box + whiskers are identical to the original plot
bmi_summary_by_tp <- bmi_avg_by_tp %>%
  filter(timepoint %in% c("5_years", "at_diagnosis")) %>%
  group_by(match_cancer_type, timepoint, case) %>%
  summarise(
    q1     = quantile(bmi, 0.25, na.rm = TRUE),
    q3     = quantile(bmi, 0.75, na.rm = TRUE),
    iqr    = q3 - q1,
    middle = median(bmi, na.rm = TRUE),
    lower  = q1,
    upper  = q3,
    ymin   = suppressWarnings(min(bmi[bmi >= q1 - 1.5 * iqr], na.rm = TRUE)),
    ymax   = suppressWarnings(max(bmi[bmi <= q3 + 1.5 * iqr], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  select(-q1, -q3, -iqr)

# ---- function to build one cancer type's plot ----
make_case_control_boxplot <- function(ct) {
  df  <- bmi_summary_by_tp %>% filter(match_cancer_type == ct)
  lbl <- model_labels %>% filter(match_cancer_type == ct)
  lbl <- df %>%
    group_by(timepoint) %>%
    summarise(y_pos = max(ymax, na.rm = TRUE) * 1.08, .groups = "drop") %>%
    right_join(lbl, by = "timepoint")

  ggplot(df, aes(x = case, fill = case)) +
    geom_boxplot(
      aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
      stat = "identity",
      width = 0.5, color = "black", alpha = 0.5
    ) +
    geom_text(
      data = lbl,
      aes(x = 1.5, y = y_pos, label = label),
      inherit.aes = FALSE,
      parse = TRUE,          # label is plotmath from fmt_p(); see src/setup.R
      size = 2,
      lineheight = 0.9
    ) +
    facet_wrap(~ timepoint, scales = "free_y", nrow = 1,
               labeller = as_labeller(timepoint_labels)) +
    scale_fill_tableau() +
    scale_color_tableau() +
    labs(
      title = paste0("Averaged BMI by timepoint — ", ct),
      x = NULL, y = "BMI"
    ) +
    theme_std() +
    theme(legend.position = "none")
}

# ---- generate + save the 4 plots ----
cancer_types <- c("Kidney", "Pancreas", "Stomach", "Uterine & endometrial")
plots <- setNames(lapply(cancer_types, make_case_control_boxplot), cancer_types)

walk2(
  plots, names(plots),
  ~ ggsave_wrapper(
      filename = file.path(landmark_dir,
                           paste0("bmi_diff_boxplot_", janitor::make_clean_names(.y), ".pdf")),
      plot = .x,
      width = 4, height = 3
    )
)
