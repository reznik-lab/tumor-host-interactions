#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 03.3d -- cachexia prevalence per cancer type at 5, 10 and 15 % weight
## loss, three facets sharing an x scale and a row order.
##
## Produces: Supplementary Figure 1c -- cachexia prevalence per cancer type at
##           5, 10 and 15 % weight loss.
##           Supplementary Table 5 -- the counts behind those bars.
##
## Reads:
##   <workspace>/.../temp/cachexia_calls/wl_w_meta_{5,10,15}_pct.tsv
## Writes:
##   <workspace>/.../results/figs/cachexia/supp_cachexia_prevalence_by_threshold_180d.pdf (+ .rds)
##       [Supplementary Figure 1c]
##   <workspace>/.../results/tables/cachexia/supp_cachexia_prevalence_by_threshold_180d.tsv (+ .xlsx)
##       [Supplementary Table 5]
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


## --- constants -------------------------------------------------------------

## Thresholds, in the order the facets are drawn; the shared row order comes
## from WL_PCT_MAIN, whose facet reproduces Figure 1B.
WL_PCTS     <- c("5", "10", "15")
WL_PCT_MAIN <- "10"
stopifnot(WL_PCT_MAIN %in% WL_PCTS)

MAIN_LEAD_DAYS <- 180

## The AoU 1-20 disclosure band: 20 is inside it, 21 is strictly compliant.
MIN_CACHEXIA_CALLS <- 20


## --- one row per person, per threshold -------------------------------------
## cachexia_frequency.R's construction, parameterised by threshold. See the
## notes there: no pre-diagnosis date filter on the read (it would silently
## delete people whose episodes are all post-diagnosis, who belong in the
## denominator as FALSE), and no_call_people added back because slice_min(na_rm)
## drops anyone with no call at all.

cach_per_person_at <- function(pct) {
  calls <- read_tsv(file.path(WORKSPACE, sprintf(
    "pre_diagnosis_bmi_analysis/code_outputs_scripts/temp/cachexia_calls/wl_w_meta_%s_pct.tsv",
    pct)), show_col_types = FALSE)

  per_person <- calls %>%
    group_by(person_id, case, match_cancer_type, cancer_type, pseudo_diagnosis_datetime) %>%
    slice_min(start_date, n = 1, with_ties = FALSE, na_rm = TRUE) %>%
    ungroup() %>%
    select(person_id, case, match_cancer_type, cancer_type, pseudo_diagnosis_datetime,
           start_date, end_date)

  no_call_people <- calls %>%
    distinct(person_id, case, match_cancer_type, cancer_type, pseudo_diagnosis_datetime) %>%
    anti_join(per_person, by = "person_id") %>%
    mutate(start_date = as.Date(NA), end_date = as.Date(NA))

  bind_rows(per_person, no_call_people) %>%
    mutate(threshold_pct = as.integer(pct))
}

cach_per_person <- map_dfr(WL_PCTS, cach_per_person_at)


## --- flag and gate ---------------------------------------------------------
## Both copied from cachexia_frequency.R. flag_before_dx uses lubridate::days(),
## so the lead is in days; summarise_at_lead gains threshold_pct in its
## grouping and nothing else.

flag_before_dx <- function(df, lead_days) {
  df %>%
    mutate(
      before_dx = if_else(
        !is.na(start_date) & !is.na(end_date) &
          start_date < pseudo_diagnosis_datetime - lubridate::days(lead_days),
        TRUE, FALSE
      )
    )
}

summarise_at_lead <- function(df, lead_days) {
  flag_before_dx(df, lead_days) %>%
    group_by(threshold_pct, match_cancer_type, case) %>%
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


## --- the nesting assertion -------------------------------------------------
## Counts must fall as the threshold tightens. A violation is a defect
## upstream, not a result, so it stops the run.

check_monotone <- function(df, what) {
  bad <- df %>%
    arrange(match_cancer_type, case, threshold_pct) %>%
    group_by(match_cancer_type, case) %>%
    filter(any(diff(n_before_dx) > 0)) %>%
    ungroup()
  if (nrow(bad) > 0) {
    print(bad %>% select(match_cancer_type, case, threshold_pct, n_total, n_before_dx))
    stop(sprintf(
      "%s: call counts RISE as the weight-loss threshold tightens for %d (type, arm) pair(s). The thresholds nest, so this cannot happen.",
      what, n_distinct(bad[c("match_cancer_type", "case")])), call. = FALSE)
  }
}

check_monotone(prevalence_main, sprintf("%dd", MAIN_LEAD_DAYS))


## --- local theme override --------------------------------------------------
## NOT the theme_std() in setup.R -- theme_classic/ArialMT rather than
## theme_cowplot, as in cachexia_frequency.R, so these panels match Figure 1B.
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
## Ascending prevalence at WL_PCT_MAIN and MAIN_LEAD_DAYS over EVERY cancer
## type, before any gate -- the identical construction cachexia_frequency.R
## uses for type_order_all. Taking subsets of one ordering is what makes the
## three facets, the two figures and Figure 1B itself line up row for row.
type_order_all <- prevalence_main %>%
  filter(case == 1, threshold_pct == as.integer(WL_PCT_MAIN)) %>%
  mutate(match_cancer_type = display_type_label(match_cancer_type)) %>%
  arrange(proportion_before_dx) %>%
  pull(match_cancer_type)


## --- the two figures -------------------------------------------------------
## A type is drawn in a facet only where it clears the gate AT THAT THRESHOLD.
## It keeps its row in the figure as long as it clears the gate somewhere, so a
## facet with no bar for a row means "suppressed here", not "zero" -- said in
## the caption note printed below. Dropping the type from the whole figure
## instead would discard the thresholds at which it IS publishable.

## Light to dark as the threshold tightens, from the same two anchors
## cachexia_frequency.R ramps its lead-time colours between.
threshold_colours <- setNames(
  colorRampPalette(c("#E4CFE0", "#7F4F72"))(length(WL_PCTS)),
  paste0(WL_PCTS, "% WL"))

plot_by_threshold <- function(df, y_label) {
  drawn <- df %>% filter(publishable)

  ## Rows: every type publishable in at least one facet, in the shared order.
  rows <- intersect(type_order_all,
                    display_type_label(unique(drawn$match_cancer_type)))

  drawn %>%
    mutate(
      match_cancer_type = factor(display_type_label(match_cancer_type),
                                 levels = rows),
      threshold_label   = factor(paste0(threshold_pct, "% WL"),
                                 levels = paste0(WL_PCTS, "% WL"))
    ) %>%
    ggplot(aes(x = match_cancer_type, y = proportion_before_dx,
               fill = threshold_label)) +
    geom_col(width = 0.8) +
    ## drop = FALSE: a threshold whose facet is entirely suppressed still gets
    ## an (empty) panel, so the reader sees it was asked and not reported.
    facet_wrap(~ threshold_label, nrow = 1, drop = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.15))) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = threshold_colours, guide = "none") +
    labs(x = NULL, y = y_label) +
    theme_std() +
    theme(plot.margin = margin(t = 5, r = 15, b = 5, l = 5),
          panel.spacing.x = unit(8, "pt"),
          strip.background = element_blank(),
          strip.text = element_text(size = 7)) +
    coord_flip(clip = "off")
}

fig_dir <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/results/figs/cachexia")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

ggsave_wrapper(
  file.path(fig_dir, sprintf("supp_cachexia_prevalence_by_threshold_%dd.pdf",
                             MAIN_LEAD_DAYS)),
  plot_by_threshold(prevalence_main,
                    sprintf("%% with start of weight loss >= %d days before diagnosis",
                            MAIN_LEAD_DAYS)),
  width = 9, height = 3.5)


## --- the numbers behind the bars -------------------------------------------
## Gated rows only, so every count printed clears MIN_CACHEXIA_CALLS and so
## does its complement.

supp_table_for <- function(df) {
  df %>%
    filter(publishable) %>%
    transmute(
      cancer_type   = display_type_label(match_cancer_type),
      threshold_pct,
      n_cases       = n_total,
      n_cachexia    = n_before_dx,
      pct_cachexia  = round(100 * n_before_dx / n_total, 1)
    ) %>%
    arrange(cancer_type, threshold_pct)
}

supp_main    <- supp_table_for(prevalence_main)

tab_dir <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/results/tables/cachexia")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

write_table(supp_main,
            file.path(tab_dir, sprintf("supp_cachexia_prevalence_by_threshold_%dd.tsv",
                                       MAIN_LEAD_DAYS)),
            sheet = sprintf("%dd", MAIN_LEAD_DAYS))

message(sprintf("Types drawn: %d of %d (summed over facets).",
                nrow(supp_main), 3 * n_distinct(prevalence_main$match_cancer_type)))
