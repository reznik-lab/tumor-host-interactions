#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 03.3 -- prevalence of pre-diagnosis cachexia per cancer type.
##
## Produces: Figure 1b -- prevalence of pre-diagnosis cachexia per cancer
##           type, episodes starting >= 180 days before diagnosis.
##
## Reads:
##   <workspace>/.../temp/cachexia_calls/wl_w_meta_<WL_PCT>_pct.tsv
## Writes:
##   <workspace>/.../results/figs/cachexia/cachexia_prevalence_180d.pdf (+ .rds)
##       [Figure 1b]
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


## --- read postprocessed calls ----------------------------------------------

## "5", "10" or "15".
WL_PCT <- "10"

## An episode must start at least this many days before diagnosis to count.
MAIN_LEAD_DAYS <- 180

wl_10_calls_w_meta <- read_tsv(file.path(WORKSPACE, sprintf(
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/temp/cachexia_calls/wl_w_meta_%s_pct.tsv",
  WL_PCT)))

## No pre-diagnosis date filter here, deliberately: people whose episodes are
## all post-diagnosis belong in the denominator as FALSE. The pre-diagnosis
## restriction is applied in before_dx below.



## --- one row per person ----------------------------------------------------

# collapse to one row per person (earliest call, if multiple calls exist per person)
# if a person has no call, start_date/end_date are NA and stay NA here
cach_per_person <- wl_10_calls_w_meta %>%
  group_by(person_id, case, match_cancer_type, cancer_type, pseudo_diagnosis_datetime) %>%
  slice_min(start_date, n = 1, with_ties = FALSE, na_rm = TRUE) %>%
  ungroup() %>%
  select(person_id, case, match_cancer_type, cancer_type, pseudo_diagnosis_datetime,
         start_date, end_date)

# slice_min(na_rm = TRUE) drops people whose start_date is always NA (no call at all)
# add them back explicitly so they remain in the denominator, with NA dates
no_call_people <- wl_10_calls_w_meta %>%
  distinct(person_id, case, match_cancer_type, cancer_type, pseudo_diagnosis_datetime) %>%
  anti_join(cach_per_person, by = "person_id") %>%
  mutate(start_date = as.Date(NA), end_date = as.Date(NA))

cach_per_person <- bind_rows(cach_per_person, no_call_people)

# flag before_dx: call must start before diagnosis AND the episode must end
# within 6 months (180 days) of diagnosis; no-call people are FALSE, not missing
#' Flag people whose episode started at least `lead_days` before diagnosis.
#'
#' FIXED: the original was `- 180`, which subtracts 180 SECONDS because
#' pseudo_diagnosis_datetime is POSIXct -- so the 6-month lead the Figure 1B
#' y-axis label claims was never applied, and before_dx meant simply "started
#' before diagnosis". lubridate::days() makes the unit explicit. Bars can only
#' fall relative to the published figure, since anything qualifying at a real
#' 180-day lead also qualified at the old 180-second one.
#'
#' Longer leads are strictly more restrictive: the flagged set at 180 days is a
#' subset of the set at 90, which is a subset of 60, and so on.
flag_before_dx <- function(df, lead_days) {
  df %>%
    mutate(
      before_dx = if_else(
        !is.na(start_date) & !is.na(end_date) &
          start_date < pseudo_diagnosis_datetime - lubridate::days(lead_days),
          # end_date < pseudo_diagnosis_datetime + 30,
        TRUE, FALSE
      )
    )
}

cach_per_person <- flag_before_dx(cach_per_person, MAIN_LEAD_DAYS)

## --- the gate --------------------------------------------------------------
## A cancer type is drawn only where it has >= MIN_CASES_PER_TYPE cases,
## >= MIN_CACHEXIA_CALLS of them with an episode, and >= MIN_CACHEXIA_CALLS
## without one -- the complement too, since the denominator is published
## beside the bar and the unreported cell is otherwise recoverable.
MIN_CACHEXIA_CALLS <- 20

summarise_at_lead <- function(df, lead_days) {
  flag_before_dx(df, lead_days) %>%
    group_by(match_cancer_type, case) %>%
    summarise(
      n_total              = n(),
      n_before_dx          = sum(before_dx),
      proportion_before_dx = mean(before_dx),
      .groups = "drop"
    ) %>%
    mutate(
      lead_days   = lead_days,
      publishable = case == 1 &
        n_total                 >= MIN_CASES_PER_TYPE &
        n_before_dx             >= MIN_CACHEXIA_CALLS &
        (n_total - n_before_dx) >= MIN_CACHEXIA_CALLS
    )
}

prevalence_main    <- summarise_at_lead(cach_per_person, MAIN_LEAD_DAYS)


## --- local theme override --------------------------------------------------
## Not the theme_std() in setup.R: theme_classic/ArialMT, not theme_cowplot.
theme_std <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "ArialMT") %+replace%
    theme(
      line = element_line(colour = "black", linewidth = 0.1),
      axis.line.x = element_line(colour = "black", linewidth = 0.1),
      axis.line.y = element_line(colour = "black", linewidth = 0.1),
      axis.ticks = element_line(linewidth = 0.1),
      axis.ticks.length = unit(1, "pt"))
}

## --- shared row order ------------------------------------------------------
## Ascending prevalence over every cancer type, before the gate.
type_order_all <- prevalence_main %>%
  filter(case == 1) %>%
  mutate(match_cancer_type = display_type_label(match_cancer_type)) %>%
  arrange(proportion_before_dx) %>%
  pull(match_cancer_type)


## --- Figure 1B ------------------------------------------------------------

plot_prevalence <- function(df, y_label) {
  df %>%
    filter(publishable) %>%
    mutate(match_cancer_type = factor(display_type_label(match_cancer_type),
                                      levels = type_order_all)) %>%
    droplevels() %>%
    ggplot(aes(x = match_cancer_type, y = proportion_before_dx)) +
    geom_col(width = 0.8, fill = "#B07AA1", linewidth = 0.3) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = y_label) +
    theme_std() +
    theme(plot.margin = margin(t = 5, r = 15, b = 5, l = 5)) +
    coord_flip(clip = "off")
}

fig_dir <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/results/figs/cachexia")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

ggsave_wrapper(file.path(fig_dir, sprintf("cachexia_prevalence_%dd.pdf", MAIN_LEAD_DAYS)),
  plot_prevalence(prevalence_main,
                  paste0("% with start of ", WL_PCT, "% WL >= ", MAIN_LEAD_DAYS,
                         " days before diagnosis")),
  width = 5, height = 3.5)
