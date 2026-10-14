#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 05 -- downward weight-class shift between one year before diagnosis
## and diagnosis, Obese and Overweight pooled.
##
## Produces: Supplementary Figure 1d -- downward weight-class shift between one
##           year before diagnosis and diagnosis, per cancer type.
##           Supplementary Table 7 is built from that panel's .rds.
##
## Reads:
##   <workspace>/.../stable_data/cancer_patients/first_cancer_diagnosis_dates_grouped.tsv
##   <workspace>/.../stable_data/cancer_patients/raw_concatenated/measurementOccurrence.tsv
##   <workspace>/.../stable_data/psm/solid_heme/full_cohort/aggregate/psm_matched_cohort_all_cancer_types.tsv
##     -- for the cancer-type list only; no matched patient enters any number here
## Writes:
##   .../results/figs/weight_class/prop_downward_class_change_barplot_pooled.pdf (+ .rds)
##       [Supplementary Figure 1d; Supplementary Table 7 is built from the .rds]
##   .../results/tables/weight_class/weight_class_shift_pooled.tsv (+ .xlsx)
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
library(ggalluvial)

landmark_dir <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/results/figs/weight_class")

## Where the plots are written.

weight_colors <- c(
  "Underweight" = "#4477AA",
  "Normal"      = "#66CCAA",
  "Overweight"  = "#EEBB44",
  "Obese"       = "#CC6677"
)


## --- load data -------------------------------------------------------------

# cancer BMI data
cancer_dx_data <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/first_cancer_diagnosis_dates_grouped.tsv"))

cancer_bmi_data <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/raw_concatenated/measurementOccurrence.tsv"))

cancer_bmi_df <- cancer_bmi_data %>%
  inner_join(., cancer_dx_data, by = "person_id")

print(head(cancer_bmi_df))


## ---------------------------------------------------------------------------
## ASSUMPTIONS TO VERIFY:
## 1. BMI concept name match below - update BMI_CONCEPT_NAME if needed
## 2. Grouping variable = cancer_type_grouped (change to cancer_type if preferred)
## 3. Timepoint windows come from LANDMARK_TIMEPOINTS (src/setup.R), shared
##    with the Stage 04 landmark models. See the note at `targets` below for
##    how the bins are defined.
## ---------------------------------------------------------------------------

GROUP_VAR <- "cancer_type_grouped"   # change to "cancer_type" if preferred

# ---- target timepoints: the SHARED landmark windows ------------------------
## From LANDMARK_TIMEPOINTS in src/setup.R, the same windows the Stage 04
## landmark models use. This script used to carry its own: point targets at
## 5 * 365.25 / 365.25 / 182.625 / 0 days before diagnosis, with a flat
## MAX_DIST_DAYS = 90 tolerance either side. Two of those differed from the
## shared definition and nothing said so --
##
##        this script was        LANDMARK_TIMEPOINTS is
##   5y   -1826 +/- 90           -1642.5 +/- 180
##   1y   -365  +/- 90           -365    +/- 120
##   6mo  -182.6 +/- 90          -182.5  +/- 120
##   dx    0    +/- 90            0      +/- 30
##
## -- so "at diagnosis" here admitted a measurement up to 89 days AFTER
## diagnosis, which the shared window does not, and the 5-year bin was centred
## half a year earlier than every other script's. Supp Fig 1E's numbers move
## as a result; that is the point of the change, not a side effect.
##
## WHAT DOES NOT CHANGE: one value per (person, timepoint), because the
## alluvial and the class-change comparison need a single weight class per
## landmark. HOW that value is formed did change -- it is now the mean of
## every measurement in the window, as in Stage 04, rather than the single
## measurement nearest the centre. See the note at landmark_bmi below.

LANDMARK_LABELS <- c(
  "5_years"      = "5 years",
  "1_year"       = "1 year",
  "6_months"     = "6 months",
  "at_diagnosis" = "At Diagnosis"
)

## 3_years is in the shared table but has never been one of this script's bins;
## the alluvial is drawn over four. Dropped by the filter below rather than by
## editing the shared table, which Stage 04 also reads.
targets <- LANDMARK_TIMEPOINTS %>%
  filter(label %in% names(LANDMARK_LABELS)) %>%
  mutate(time_bin = factor(unname(LANDMARK_LABELS[label]),
                           levels = unname(LANDMARK_LABELS))) %>%
  select(time_bin, center, window_lo, window_hi) %>%
  arrange(time_bin)

## Every label this script expects has to be in the shared table -- if one is
## renamed there, fail here rather than silently drawing three bins.
stopifnot(nrow(targets) == length(LANDMARK_LABELS),
          !any(is.na(targets$time_bin)))

print(targets)

# ---- usable BMI measurements ----
## This stage previously had NO value gate, unlike every other stage. That
## matters here more than anywhere else: the weight_category case_when below
## ends in an open-ended `TRUE ~ "Obese"`, so an implausibly high BMI was
## classified as Obese rather than dropped, and an implausibly low one as
## Underweight. filter_usable_bmi() (src/setup.R) is the shared gate.
##
## The table is already BMI-only -- t_disp_standard_concept_name has the single
## value "Body mass index (BMI) [Ratio]" -- so no concept filter is needed here.
## (The ASSUMPTIONS block above refers to a BMI_CONCEPT_NAME that was never
## defined; there is nothing to define, and the note is left as-is only so the
## numbering is unchanged.)
bmi_measurements <- cancer_bmi_df %>%
  mutate(
    first_dx_datetime    = as.POSIXct(first_dx_datetime, tz = "UTC"),
    measurement_datetime = as.POSIXct(measurement_datetime, tz = "UTC")
  ) %>%
  filter_usable_bmi() %>%
  filter(!is.na(measurement_datetime), !is.na(first_dx_datetime))

# ---- for each person x target timepoint, the mean BMI in that window ----
## ---- restrict to the reportable cancer types ------------------------------
## THE SAME SET OF TYPES EVERY OTHER PANEL DRAWS, not merely the same rule.
##
## The floor is >= MIN_CASES_PER_TYPE (src/setup.R) counted on CASES, and every
## other stage counts it on cases in the PSM-MATCHED cohort. This script runs
## on the raw, unmatched cohort, which is a superset: a type can have 120 raw
## cases and 85 matched ones. Applying the floor to this script's own counts
## would therefore be the same rule over a different denominator, and could put
## a cancer type in Supp Fig 1E that appears in no other figure. So the type
## list is read from the matched cohort rather than recomputed here.
##
## This is the one input this script takes from the PSM pipeline; everything
## else below is still the raw cohort, and the DENOMINATORS are still raw-cohort
## patients. Only WHICH TYPES are drawn comes from the match.
reportable_types <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm/solid_heme/full_cohort/aggregate/psm_matched_cohort_all_cancer_types.tsv"),
  show_col_types = FALSE) %>%
  filter(case == 1) %>%
  distinct(person_id, match_cancer_type) %>%
  count(match_cancer_type, name = "n_matched_cases") %>%
  filter(n_matched_cases >= MIN_CASES_PER_TYPE) %>%
  pull(match_cancer_type)

stopifnot(length(reportable_types) > 0)

message(sprintf("Type gate: %d cancer types with >= %d matched cases.",
                length(reportable_types), MIN_CASES_PER_TYPE))

## What the gate costs HERE, in this script's own cohort -- the types dropped
## and the patients with them. Reported rather than silent, since this cohort
## is not the one the gate was computed on.
dropped_types <- setdiff(unique(bmi_measurements$cancer_type_grouped), reportable_types)
if (length(dropped_types) > 0) {
  dropped_n <- bmi_measurements %>%
    filter(cancer_type_grouped %in% dropped_types) %>%
    distinct(person_id, cancer_type_grouped) %>%
    count(cancer_type_grouped, name = "n") %>%
    arrange(desc(n))
  message(sprintf("  dropped %d types / %s patients from this cohort: %s",
                  nrow(dropped_n), scales::comma(sum(dropped_n$n)),
                  paste(sprintf("%s (%d)", dropped_n$cancer_type_grouped, dropped_n$n),
                        collapse = ", ")))
}

## A type that clears the floor in the matched cohort but is absent here means
## the two cohorts have diverged -- worth seeing, not worth stopping for.
missing_here <- setdiff(reportable_types, unique(bmi_measurements$cancer_type_grouped))
if (length(missing_here) > 0) {
  warning("reportable in the matched cohort but absent from this raw cohort: ",
          paste(missing_here, collapse = ", "), call. = FALSE)
}

bmi_measurements <- bmi_measurements %>%
  filter(cancer_type_grouped %in% reportable_types)


## MEAN BMI OVER THE WINDOW, one value per (person, timepoint).
##
## This took the single measurement nearest the window's centre. Every other
## within-person aggregation in the pipeline is a mean over the window --
## mixed_effects_modelling.R's per-person landmark BMI, analysis.R's
## per-person bin BMI -- and nothing said why this one differed.
##
## It also discarded most of the data: the windows are up to +/- 180 days
## wide, so a person with several measurements in one contributed exactly one,
## chosen by how close their visit happened to fall to the centre. That made
## the class assignment turn on visit timing, which is not a property of the
## patient. It matters more here than elsewhere because the downstream step is
## a hard threshold -- a single measurement at BMI 30.2 rather than 29.8 flips
## a patient between Obese and Overweight, and the whole analysis is about
## crossing those boundaries.
##
## MEAN BMI, THEN CLASSIFY -- not the modal class of individually classified
## measurements, which would throw away exactly the information that matters
## near a cut point.
##
## n_measurements_in_window is kept so a value can be read against how much it
## rests on: a mean over one measurement is the old behaviour for that person.
landmark_bmi <- bmi_measurements %>%
  select(person_id, cancer_type, cancer_type_grouped, first_dx_datetime,
         measurement_datetime, value_as_number) %>%
  tidyr::crossing(targets) %>%
  mutate(
    ## signed, diagnosis-relative, as LANDMARK_TIMEPOINTS is defined: negative
    ## is before diagnosis. The old code worked in a positive days_before_dx
    ## and an absolute tolerance, which cannot express an asymmetric window.
    days_from_dx = as.numeric(difftime(measurement_datetime, first_dx_datetime,
                                       units = "days"))
  ) %>%
  filter(days_from_dx >= window_lo, days_from_dx <= window_hi) %>%
  group_by(person_id, cancer_type, cancer_type_grouped, first_dx_datetime, time_bin) %>%
  summarise(
    n_measurements_in_window = n_distinct(measurement_datetime),
    value_as_number          = mean(value_as_number, na.rm = TRUE),
    .groups = "drop"
  )

## One row per (person, timepoint), or every count downstream is a measurement
## count wearing a person count's name.
stopifnot(!any(duplicated(landmark_bmi[c("person_id", "time_bin")])))

message(sprintf("Landmark BMI: %s (person, timepoint) values from %s measurements; median %g per window, %.1f%% rest on a single measurement.",
                scales::comma(nrow(landmark_bmi)),
                scales::comma(sum(landmark_bmi$n_measurements_in_window)),
                median(landmark_bmi$n_measurements_in_window),
                100 * mean(landmark_bmi$n_measurements_in_window == 1)))

# ---- categorize BMI into weight category ----
landmark_bmi <- landmark_bmi %>%
  mutate(
    weight_category = case_when(
      value_as_number < 18.5 ~ "Underweight",
      value_as_number < 25   ~ "Normal",
      value_as_number < 30   ~ "Overweight",
      TRUE                   ~ "Obese"
    ),
    weight_category = factor(weight_category,
                             levels = c("Underweight", "Normal", "Overweight", "Obese"))
  )

# ---- denominator: people with a matched BMI measurement at that timepoint ----
n_totals_bmi <- landmark_bmi %>%
  distinct(person_id, .data[[GROUP_VAR]], time_bin) %>%
  count(.data[[GROUP_VAR]], time_bin, name = "n_total")

# ---- proportion in each weight category, per time bin ----
weight_category_by_bin <- landmark_bmi %>%
  count(.data[[GROUP_VAR]], time_bin, weight_category, name = "n") %>%
  complete(
    weight_category = c("Underweight", "Normal", "Overweight", "Obese"),
    nesting(!!sym(GROUP_VAR), time_bin),
    fill = list(n = 0)
  ) %>%
  left_join(n_totals_bmi, by = c(GROUP_VAR, "time_bin")) %>%
  mutate(
    weight_category = factor(weight_category,
                             levels = c("Underweight", "Normal", "Overweight", "Obese")),
    proportion = n / n_total
  ) %>%
  arrange(.data[[GROUP_VAR]], time_bin, weight_category)

print(weight_category_by_bin)


## --- two-bin cohort: 1 year prior -> at diagnosis --------------------------

BIN_LEVELS_2 <- c("1 year", "At Diagnosis")

# ---- restrict to individuals with a matched measurement in BOTH bins ----
complete_persons_2bin <- landmark_bmi %>%
  filter(time_bin %in% BIN_LEVELS_2) %>%
  distinct(person_id, time_bin) %>%
  count(person_id, name = "n_bins_present") %>%
  filter(n_bins_present == length(BIN_LEVELS_2)) %>%
  pull(person_id)

# carry cancer_type_grouped through
alluvial_data_2bin <- landmark_bmi %>%
  filter(person_id %in% complete_persons_2bin, time_bin %in% BIN_LEVELS_2) %>%
  select(person_id, cancer_type_grouped, time_bin, weight_category) %>%
  distinct(person_id, cancer_type_grouped, time_bin, .keep_all = TRUE)

# ---- pivot wide: one row per person, one column per time bin ----
alluvial_wide_2bin <- alluvial_data_2bin %>%
  pivot_wider(names_from = time_bin, values_from = weight_category) %>%
  select(person_id, cancer_type_grouped, all_of(BIN_LEVELS_2))

# ---- collapse into counts of each flow pattern, WITHIN each cancer_type_grouped ----
alluvial_counts_2bin <- alluvial_wide_2bin %>%
  count(cancer_type_grouped, `1 year`, `At Diagnosis`, name = "freq")




## ---------------------------------------------------------------------------
## Downward weight-class change, one year before diagnosis -> at diagnosis.
##
## The cohort is already restricted to types clearing MIN_CASES_PER_TYPE. On
## top of that, MIN_COHORT_PER_STRATUM is a per-cell sample-size floor: it
## decides which bars a covered type gets, over a much smaller denominator
## (patients with a usable BMI in BOTH landmark windows), so the two floors
## move independently.
## ---------------------------------------------------------------------------

weight_levels <- c("Underweight", "Normal", "Overweight", "Obese")  # low -> high

## Minimum patients in one (cancer type x starting class) cell for that bar to
## be drawn. At 25 an exact 95% CI is still ~ +/- 20 points wide, which is why
## the figure carries the interval and the table carries n.
MIN_COHORT_PER_STRATUM <- 25

ALL_TYPES_LABEL <- "All cancers (pan-cancer)"

## ---- the two-bin cohort, keyed by where the patient started ---------------
class_change_data_all <- alluvial_wide_2bin %>%
  filter(!is.na(`1 year`), !is.na(`At Diagnosis`)) %>%
  filter(`1 year` %in% c("Obese", "Overweight")) %>%
  mutate(
    baseline_class = factor(as.character(`1 year`), levels = c("Obese", "Overweight")),
    rank_prior     = match(`1 year`, weight_levels),
    rank_diagnosis = match(`At Diagnosis`, weight_levels),
    class_down     = rank_diagnosis < rank_prior
  )

## ---- which (type x starting class) cells clear the floor ------------------
cell_n <- class_change_data_all %>%
  count(cancer_type_grouped, baseline_class, name = "n_cohort")

cells_ok <- cell_n %>% filter(n_cohort >= MIN_COHORT_PER_STRATUM)

message(sprintf(
  "Supp Fig 1E: %d of %d (cancer type x starting class) cells clear n >= %d, spanning %d cancer types.",
  nrow(cells_ok), nrow(cell_n), MIN_COHORT_PER_STRATUM,
  n_distinct(cells_ok$cancer_type_grouped)))

stopifnot(nrow(cells_ok) > 0)

## The per-type rows and the pan-cancer row are the SAME computation over two
## groupings, so they are stacked here and summarised once below rather than
## written twice. The pan-cancer row is ungated on purpose: it is the reference
## line, and should be the whole cohort rather than the plotted subset of it.
stratum_input <- bind_rows(
  class_change_data_all %>%
    semi_join(cells_ok, by = c("cancer_type_grouped", "baseline_class")),
  class_change_data_all %>%
    mutate(cancer_type_grouped = ALL_TYPES_LABEL)
)

## ---- rate of any downward class change, within each cell ------------------
## Exact (Clopper-Pearson) interval, so a bar drawn off 25 patients cannot be
## read as a bar drawn off 500.
exact_ci <- function(x, n, bound) {
  vapply(seq_along(x), function(i) {
    if (is.na(n[i]) || n[i] == 0) return(NA_real_)
    stats::binom.test(x[i], n[i])$conf.int[bound]
  }, numeric(1))
}

class_down_by_stratum <- stratum_input %>%
  group_by(cancer_type_grouped, baseline_class) %>%
  summarise(n_cohort = n(), n_class_down = sum(class_down), .groups = "drop") %>%
  mutate(
    proportion = n_class_down / n_cohort,
    ci_low     = exact_ci(n_class_down, n_cohort, 1),
    ci_high    = exact_ci(n_class_down, n_cohort, 2)
  )

## ---- where they landed: destination class at diagnosis, within each cell --
## `At Diagnosis` is carried as CHARACTER through the complete(): the column is
## a factor on the way in, and complete() replaces it with the character vector
## it is given, so the rejoin would be factor-against-character. Re-levelled
## once below, after the fill.
dest_by_stratum <- stratum_input %>%
  mutate(dx_class = as.character(`At Diagnosis`)) %>%
  count(cancer_type_grouped, baseline_class, dx_class, name = "n") %>%
  complete(
    dx_class = weight_levels,
    nesting(cancer_type_grouped, baseline_class),
    fill = list(n = 0)
  ) %>%
  left_join(
    class_down_by_stratum %>% select(cancer_type_grouped, baseline_class, n_cohort),
    by = c("cancer_type_grouped", "baseline_class")
  ) %>%
  mutate(
    dx_class   = factor(dx_class, levels = weight_levels),
    proportion = n / n_cohort,
    rank_dest  = match(as.character(dx_class), weight_levels),
    rank_base  = match(as.character(baseline_class), weight_levels),
    moved_down = rank_dest < rank_base
  )

## the stacked segments: destinations BELOW where the patient started, so the
## segments of one bar sum to that cell's downward rate above
dest_down <- dest_by_stratum %>% filter(moved_down, n_cohort > 0)

overall_by_stratum <- class_down_by_stratum %>%
  filter(cancer_type_grouped == ALL_TYPES_LABEL)

print(overall_by_stratum)
for (i in seq_len(nrow(overall_by_stratum))) {
  message(sprintf(
    "Supp Fig 1E pan-cancer: %.1f%% of patients %s 1 year before diagnosis (%d of %d) had dropped at least one weight class by diagnosis.",
    100 * overall_by_stratum$proportion[i],
    tolower(as.character(overall_by_stratum$baseline_class[i])),
    overall_by_stratum$n_class_down[i], overall_by_stratum$n_cohort[i]))
}

## the destination split, pan-cancer -- the sentence this panel is for
pan_dest <- dest_down %>% filter(cancer_type_grouped == ALL_TYPES_LABEL)
for (i in seq_len(nrow(pan_dest))) {
  message(sprintf(
    "  ... of those %s at 1 year, %.1f%% were %s at diagnosis (%d of %d).",
    tolower(as.character(pan_dest$baseline_class[i])),
    100 * pan_dest$proportion[i],
    tolower(as.character(pan_dest$dx_class[i])),
    pan_dest$n[i], pan_dest$n_cohort[i]))
}


## --- continuity panel: the pooled version, same cohort and same gate -------
## This is what Supp Fig 1E was before the stratification above: one bar per
## cancer type, Obese and Overweight patients pooled. Kept so the published
## panel still regenerates, but note it is NOT identical to the published one:
## it is now drawn over every type clearing both the MIN_CASES_PER_TYPE cohort
## gate and MIN_COHORT_PER_STRATUM, rather than the six matched by the old
## hard-coded grepl(), and its dashed line is the computed pan-cancer rate
## rather than the former literal 0.14.

## ---- which types get a bar -----------------------------------------------
## A type is drawn only where all four of these clear MIN_POOLED_BAR_N:
##   n_cohort           the denominator
##   n_obese_1yr        the Obese half of it
##   n_overweight_1yr   the Overweight half
##   n_dropped_a_class  the numerator
##
## The first and last are the two numbers the bar publishes. The two baseline
## counts are here because this panel POOLS them, and a "pooled" bar built
## from one baseline is not the same quantity as one built from both: an Obese
## patient has two classes to fall through and an Overweight patient one, so a
## type contributing only Overweight patients sits on a mechanically lower
## ceiling and is not comparable with the bars beside it. The per-cell gate
## above drops such a cell's patients but leaves the type, which is exactly
## the case this catches.
##
## 20, as asked. Note the AoU band is 1-20 INCLUSIVE, so a cell of exactly 20
## still discloses a count; 21 would be the strictly compliant floor. The
## complement is not gated here either -- see the check after the plot.
MIN_POOLED_BAR_N <- 20

pooled_all_types <- stratum_input %>%
  filter(cancer_type_grouped != ALL_TYPES_LABEL) %>%
  group_by(cancer_type_grouped) %>%
  summarise(
    n_cohort         = n(),
    n_obese_1yr      = sum(baseline_class == "Obese"),
    n_overweight_1yr = sum(baseline_class == "Overweight"),
    n_class_down     = sum(class_down),
    .groups = "drop"
  ) %>%
  mutate(proportion = n_class_down / n_cohort,
         plotted = n_cohort         >= MIN_POOLED_BAR_N &
                   n_obese_1yr      >= MIN_POOLED_BAR_N &
                   n_overweight_1yr >= MIN_POOLED_BAR_N &
                   n_class_down     >= MIN_POOLED_BAR_N)

class_down_pooled <- pooled_all_types %>%
  filter(plotted) %>%
  arrange(desc(proportion))

dropped_bars <- pooled_all_types %>% filter(!plotted)
if (nrow(dropped_bars) > 0) {
  message(sprintf("Pooled panel: %d of %d types drawn; %d below MIN_POOLED_BAR_N = %d.",
                  nrow(class_down_pooled), nrow(pooled_all_types),
                  nrow(dropped_bars), MIN_POOLED_BAR_N))
  message("  excluded: ", paste(sprintf(
    "%s (cohort %d, obese %d, overweight %d, dropped %d)",
    display_type_label(dropped_bars$cancer_type_grouped),
    dropped_bars$n_cohort, dropped_bars$n_obese_1yr,
    dropped_bars$n_overweight_1yr, dropped_bars$n_class_down), collapse = "; "))
}

stopifnot(nrow(class_down_pooled) > 0)

## The bar publishes n_cohort and n_class_down, so n_cohort - n_class_down is
## recoverable by subtraction and is disclosed whether or not it is printed.
## Not part of the gate above, which is the four counts as specified -- flagged
## rather than filtered, since dropping a bar is an editorial call.
complement_small <- class_down_pooled %>%
  filter(dplyr::between(n_cohort - n_class_down, 1, MIN_POOLED_BAR_N))
if (nrow(complement_small) > 0) {
  warning("drawn bars whose NON-dropped count is 1-", MIN_POOLED_BAR_N,
          " and so disclosed by subtraction: ",
          paste(sprintf("%s (%d)", display_type_label(complement_small$cancer_type_grouped),
                        complement_small$n_cohort - complement_small$n_class_down),
                collapse = ", "), call. = FALSE)
}

overall_class_down <- class_change_data_all %>%
  summarise(
    n_cohort     = n(),
    n_class_down = sum(class_down),
    proportion   = n_class_down / n_cohort
  )

overall_rate <- overall_class_down$proportion  # pan-cancer, computed not pinned

message(sprintf(
  "Supp Fig 1E pooled overall rate: %.4f (%d of %d in the pan-cancer two-bin cohort). Was hard-coded at 0.14.",
  overall_rate, overall_class_down$n_class_down, overall_class_down$n_cohort))

print(class_down_pooled)

prop_downward_class_change_pooled <- ggplot(
  class_down_pooled,
  aes(x = reorder(display_type_label(cancer_type_grouped), proportion),
      y = proportion)) +
  geom_col(fill = "#4477AA") +
  geom_hline(yintercept = overall_rate, linetype = "dashed", color = "black", linewidth = 0.1) +
  annotate(
    "text",
    x = 0.65, y = overall_rate,
    label = paste0("Overall: ", scales::percent(overall_rate, accuracy = 1)),
    vjust = -0.6, hjust = 0, size = 2, color = "black"
  ) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA), expand = c(0, 0)) +
  labs(x = NULL, y = "% with lower weight class at diagnosis") +
  theme_std()

ggsave_wrapper(file.path(landmark_dir, "prop_downward_class_change_barplot_pooled.pdf"),
               prop_downward_class_change_pooled,
               width = 2.5, height = max(1.5, 0.12 * nrow(class_down_pooled) + 0.6))


table_dir <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/results/tables/weight_class")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


## --- supplementary table, POOLED -------------------------------------------
## The same question with the starting classes collapsed: of the patients who
## were Obese OR Overweight a year before diagnosis, how many had fallen to a
## lower class by diagnosis, per cancer type. One row per type instead of two.
## This is the table for the pooled panel
## (prop_downward_class_change_barplot_pooled.pdf).
##
## BUILT FROM stratum_input, THE SAME FRAME THE POOLED FIGURE USES, so the two
## cannot disagree. That frame is gated per (type x starting class), which has
## a consequence worth stating: a type whose Overweight cell cleared
## MIN_COHORT_PER_STRATUM but whose Obese cell did not contributes only its
## Overweight patients here, so n_cohort is NOT that type's full pooled cohort.
## n_cohort_ungated carries the full figure beside it, and n_excluded_by_cell_gate
## the difference, so a row where the pooling is partial is visible rather than
## silently understated. The pan-cancer row is ungated by construction (see the
## stratum_input note above), so for it the two agree.
##
## "class_down" IS PER PATIENT, AGAINST THEIR OWN BASELINE. Obese -> Overweight
## is a drop; Overweight -> Overweight is not. That is why the destination
## columns below are not interchangeable with the drop count: ending Overweight
## means something different depending on where the patient started.
##
## baseline composition is reported for the same reason. The pooled rate
## depends on the mix: an Obese patient has two classes to fall through and an
## Overweight patient one, so a type with proportionally more Obese patients
## will show a higher pooled rate for that reason alone. Two types' pooled
## rates are not comparable without n_obese_1yr / n_overweight_1yr beside them.

## Full pooled cohort per type, before the per-cell gate -- the denominator the
## gated one is compared against above. class_change_data_all carries real
## cancer types only; the pan-cancer row is added here the same way
## stratum_input adds it, or the join below would leave it NA.
pooled_ungated <- bind_rows(
  class_change_data_all %>% count(cancer_type_grouped, name = "n_cohort_ungated"),
  tibble(cancer_type_grouped = ALL_TYPES_LABEL,
         n_cohort_ungated    = nrow(class_change_data_all))
)

pooled_core <- stratum_input %>%
  group_by(cancer_type_grouped) %>%
  summarise(
    n_cohort        = n(),
    n_obese_1yr     = sum(baseline_class == "Obese"),
    n_overweight_1yr = sum(baseline_class == "Overweight"),
    n_dropped_a_class = sum(class_down),
    .groups = "drop"
  ) %>%
  mutate(
    prop_dropped_a_class = n_dropped_a_class / n_cohort,
    ci_low  = exact_ci(n_dropped_a_class, n_cohort, 1),
    ci_high = exact_ci(n_dropped_a_class, n_cohort, 2)
  )

## Where they ended up, pooled: one column pair per class at diagnosis. These
## are destinations, NOT drops -- see the note above.
pooled_dest <- stratum_input %>%
  mutate(dx_class = as.character(`At Diagnosis`)) %>%
  count(cancer_type_grouped, dx_class, name = "n") %>%
  complete(dx_class = weight_levels, nesting(cancer_type_grouped),
           fill = list(n = 0)) %>%
  left_join(pooled_core %>% select(cancer_type_grouped, n_cohort),
            by = "cancer_type_grouped") %>%
  mutate(dest_key = str_to_lower(dx_class), proportion = n / n_cohort) %>%
  select(cancer_type_grouped, dest_key, n, proportion) %>%
  pivot_wider(names_from = dest_key, values_from = c(n, proportion),
              names_glue = "{dest_key}_at_dx_{.value}")

supp_weight_class_shift_pooled <- pooled_core %>%
  left_join(pooled_dest,     by = "cancer_type_grouped") %>%
  left_join(pooled_ungated,  by = "cancer_type_grouped") %>%
  mutate(
    cancer_type_label = display_type_label(cancer_type_grouped),
    is_pan_cancer     = cancer_type_grouped == ALL_TYPES_LABEL,
    n_excluded_by_cell_gate = n_cohort_ungated - n_cohort
  ) %>%
  select(
    cancer_type = cancer_type_grouped, cancer_type_label, is_pan_cancer,
    n_cohort, n_obese_1yr, n_overweight_1yr,
    n_dropped_a_class, prop_dropped_a_class, ci_low, ci_high,
    starts_with("obese_at_dx_"), starts_with("overweight_at_dx_"),
    starts_with("normal_at_dx_"), starts_with("underweight_at_dx_"),
    n_cohort_ungated, n_excluded_by_cell_gate
  ) %>%
  ## which rows the pooled panel actually draws, so the table and the figure
  ## beside it cannot disagree about that. The pan-cancer row is the reference
  ## line rather than a bar, so it is not gated.
  left_join(pooled_all_types %>% select(cancer_type = cancer_type_grouped, plotted),
            by = "cancer_type") %>%
  mutate(plotted = if_else(is_pan_cancer, NA, plotted)) %>%
  arrange(desc(is_pan_cancer), desc(prop_dropped_a_class))

## The pan-cancer row is ungated, so its two cohort figures must agree. If this
## fires, stratum_input's pan-cancer arm has started being filtered.
pan <- supp_weight_class_shift_pooled %>% filter(is_pan_cancer)
stopifnot(nrow(pan) == 1, pan$n_excluded_by_cell_gate == 0)

## The drop count cannot exceed the cohort, and the two baselines must account
## for all of it.
stopifnot(
  all(supp_weight_class_shift_pooled$n_dropped_a_class <=
        supp_weight_class_shift_pooled$n_cohort),
  all(supp_weight_class_shift_pooled$n_obese_1yr +
        supp_weight_class_shift_pooled$n_overweight_1yr ==
        supp_weight_class_shift_pooled$n_cohort)
)

write_table(supp_weight_class_shift_pooled,
            file.path(table_dir, "weight_class_shift_pooled.tsv"),
            sheet = "pooled")

print(supp_weight_class_shift_pooled, n = Inf, width = Inf)

partial <- supp_weight_class_shift_pooled %>%
  filter(!is_pan_cancer, n_excluded_by_cell_gate > 0)
if (nrow(partial) > 0) {
  message(sprintf(
    "%d type(s) pool only one starting class, because the other cell is below MIN_COHORT_PER_STRATUM = %d: %s. Their n_cohort is not the full pooled cohort -- see n_cohort_ungated.",
    nrow(partial), MIN_COHORT_PER_STRATUM,
    paste(sprintf("%s (%d of %d)", partial$cancer_type_label,
                  partial$n_cohort, partial$n_cohort_ungated), collapse = ", ")))
}

message(sprintf("Wrote %d rows to %s",
                nrow(supp_weight_class_shift_pooled),
                file.path(table_dir, "weight_class_shift_pooled.tsv")))
