#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 04.1 -- BMI trajectories before diagnosis, one panel per cancer type.
##
## Produces: Figure 1c-f, the BMI trajectory and the sample-size row of each
##           panel -- 1c stomach, 1d pancreas, 1e kidney, 1f endometrial.
##           The boxplot beside each comes from mixed_effects_modelling.R.
##
## Reads:
##   <workspace>/.../stable_data/psm_matched_usable_data/solid_heme/full_cohort/all_cancers_all_bmi_measurements_geq4_eligible.tsv
##   <workspace>/.../stable_data/psm_matched_usable_data/solid_heme/full_cohort/all_controls_all_bmi_measurements_geq4_eligible.tsv
## Writes (under .../results/figs/bmi_trajectories/):
##   trajectory/<type>_bmi_trajectory.pdf (+ .rds)             [Figure 1c-f]
##   sample_size_table/<type>_bmi_sample_size_table.pdf (+ .rds)
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
library(broom)


## --- load ------------------------------------------------------------------

cancers <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm_matched_usable_data/solid_heme/full_cohort/all_cancers_geq1_before_diagnosis_after_psm.tsv"))

controls <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm_matched_usable_data/solid_heme/full_cohort/all_controls_geq1_before_diagnosis_after_psm.tsv"))

print(unique(cancers$match_cancer_type))

# for each case, assign their cancer_type to be the same as their matched case
# as defined by subclass
matched_df <- bind_rows(cancers, controls) %>%
  group_by(subclass) %>%
  mutate(
    cancer_type = if_else(
      case == 0,                        # or case == "control", however it's coded
      first(cancer_type[case == 1]),    # pull the case's cancer_type
      cancer_type
    )
  ) %>%
  ungroup()


## --- shared bin definitions ------------------------------------------------
## Bin edges in months before diagnosis: 0-1, 1-3, 3-6, 6-9, 9-12, 12-18,
## 18-24, 24-36, 36-48, 48-60, then an open-ended "60+ months" bin pooling
## everything beyond 5 years.
BIN_LABELS <- c("0-1", "1-3", "3-6", "6-9",
                "9-12", "12-18", "18-24",
                "24-36", "36-48", "48-60",
                "60+")
BIN_BREAKS <- c(0, 1, 3, 6, 9, 12, 18, 24, 36, 48, 60, Inf)

# x-axis order: furthest from diagnosis on the left -> closest on the right
BIN_ORDER_LEFT_TO_RIGHT <- c("60+", "48-60", "36-48",
                             "24-36", "18-24", "12-18",
                             "9-12", "6-9", "3-6",
                             "1-3", "0-1")

MIN_N_PER_GROUP <- 2

base_dir <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/results/figs/bmi_trajectories")


## ===========================================================================
## One trajectory panel per cancer type.
## ===========================================================================

## Cancer types to run, with optional per-type exclusion regex applied against
## `match_cancer_type`.
cancer_types <- tibble::tribble(
  ~ct_of_interest, ~exclude_pattern,
  "kidney",         NA_character_,
  "uterine",        NA_character_,
  "pancrea",        "neuroendocrine",   # excludes PanNETs, keeps exocrine/adenocarcinoma
  "stomach",        NA_character_,
  "lung",           NA_character_
)

## Builds bin_summary, star_positions, p, p_table for one cancer type and
## returns them + saves outputs.
run_bmi_trajectory <- function(ct_of_interest, exclude_pattern = NA_character_,
                               matched_df, base_dir) {

  message("Running: ", ct_of_interest,
          if (!is.na(exclude_pattern)) paste0(" (excluding '", exclude_pattern, "')") else "")

  ## -- y-axis limits: uterine gets a higher range, everything else shares 25-32 --
  y_limits <- if (grepl("uterine", ct_of_interest, ignore.case = TRUE)) {
    c(28, 35)
  } else {
    c(25, 32)
  }

  ## -- 0. Filter to cancer type, optionally excluding a histology pattern --
  df <- matched_df %>%
    mutate(
      measurement_datetime      = as_datetime(measurement_datetime),
      pseudo_diagnosis_datetime = as_datetime(pseudo_diagnosis_datetime),
      case = factor(case, levels = c(0, 1), labels = c("Control", "Cancer"))
    ) %>%
    filter(grepl(ct_of_interest, match_cancer_type, ignore.case = TRUE))

  if (!is.na(exclude_pattern)) {
    df <- df %>%
      filter(!grepl(exclude_pattern, match_cancer_type, ignore.case = TRUE))
  }

  ## -- 0b. Pipeline-wide floor (MIN_CASES_PER_TYPE, src/setup.R) --
  ## This type list is hand-written, so a type below the floor is skipped
  ## LOUDLY rather than quietly drawn: the message names it and the count.
  n_cases <- n_distinct(df$person_id[df$case == "Cancer"])
  if (n_cases < MIN_CASES_PER_TYPE) {
    message(sprintf("  SKIPPED -- %d cases, below MIN_CASES_PER_TYPE = %d.",
                    n_cases, MIN_CASES_PER_TYPE))
    return(NULL)
  }

  ## -- 1. BMI measurements before diagnosis, binned --
  bmi_pre <- df %>%
    filter(
      str_detect(t_disp_standard_concept_name, regex("BMI", ignore_case = TRUE)),
      measurement_datetime <= pseudo_diagnosis_datetime,
      !is.na(value_as_number)
    ) %>%
    mutate(
      months_before_dx = as.numeric(
        difftime(pseudo_diagnosis_datetime, measurement_datetime, units = "days")
      ) / (365.25 / 12),
      time_bin = cut(
        months_before_dx,
        breaks = BIN_BREAKS,
        labels = BIN_LABELS,
        right = FALSE,
        include.lowest = TRUE
      )
    ) %>%
    filter(!is.na(time_bin))

  ## -- 2. Per-person mean BMI within each bin --
  person_bin_means <- bmi_pre %>%
    group_by(person_id, case, time_bin) %>%
    summarise(person_mean_bmi = mean(value_as_number, na.rm = TRUE), .groups = "drop")

  ## -- 3. Across-person mean + SEM per bin, per group --
  bin_summary <- person_bin_means %>%
    group_by(case, time_bin) %>%
    summarise(
      n_people = n(),
      mean_bmi = mean(person_mean_bmi, na.rm = TRUE),
      sd_bmi   = sd(person_mean_bmi, na.rm = TRUE),
      sem_bmi  = sd_bmi / sqrt(n_people),
      .groups = "drop"
    )

  ## -- 4. Significance per bin (Wilcoxon) --
  sig_by_bin <- person_bin_means %>%
    group_by(time_bin) %>%
    filter(
      n_distinct(person_id[case == "Cancer"])  >= MIN_N_PER_GROUP,
      n_distinct(person_id[case == "Control"]) >= MIN_N_PER_GROUP
    ) %>%
    group_modify(~ {
      .x <- droplevels(.x)
      if (nlevels(.x$case) != 2) {
        return(tibble(p.value = NA_real_))
      }
      tidy(wilcox.test(person_mean_bmi ~ case, data = .x, exact = FALSE)) %>%
        select(p.value)
    }) %>%
    ungroup()

  if (nrow(sig_by_bin) == 0) {
    sig_by_bin <- tibble(time_bin  = factor(character(), levels = BIN_LABELS),
                         p.value   = numeric(0),
                         sig_label = character(0))
    message("  No time_bin had >=", MIN_N_PER_GROUP,
            " people in both groups -- no significance tests run.")
  } else {
    sig_by_bin <- sig_by_bin %>%
      mutate(
        sig_label = case_when(
          p.value < 0.001 ~ "***",
          p.value < 0.01  ~ "**",
          p.value < 0.05  ~ "*",
          TRUE            ~ ""
        )
      ) %>%
      select(time_bin, p.value, sig_label)
  }

  star_positions <- bin_summary %>%
    group_by(time_bin) %>%
    summarise(y_star = max(mean_bmi + sem_bmi, na.rm = TRUE) + 0.4, .groups = "drop") %>%
    left_join(sig_by_bin, by = "time_bin") %>%
    filter(sig_label != "")

  ## -- 5. x-axis setup --
  bin_summary <- bin_summary %>%
    mutate(time_bin = factor(time_bin, levels = BIN_ORDER_LEFT_TO_RIGHT),
           x_num = as.numeric(time_bin))

  star_positions <- star_positions %>%
    mutate(time_bin = factor(time_bin, levels = BIN_ORDER_LEFT_TO_RIGHT),
           x_num = as.numeric(time_bin))

  ## -- Plot p --
  p <- ggplot(bin_summary, aes(x = x_num, y = mean_bmi, color = case, group = case)) +
    geom_hline(yintercept = 30, linetype = "dashed", color = "red", linewidth = 0.1) +
    geom_smooth(method = "loess", span = 0.2, linewidth = 0.1, se = FALSE) +
    geom_linerange(aes(ymin = mean_bmi - sem_bmi, ymax = mean_bmi + sem_bmi), linewidth = 0.1) +
    geom_point(size = 2, shape = 19) +
    scale_color_manual(values = c("Control" = "#4477AA", "Cancer" = "#D85A30")) +
    scale_x_continuous(
      breaks = seq_along(BIN_ORDER_LEFT_TO_RIGHT),
      labels = BIN_ORDER_LEFT_TO_RIGHT,
      expand = expansion(mult = 0.03)
    ) +
    labs(x = "Time to cancer diagnosis (Months)", y = "Mean BMI", color = NULL, fill = NULL) +
    theme_std() +
    theme(
      axis.ticks      = element_blank(),
      legend.position = "none",
      axis.title.x    = element_text(margin = margin(t = 6)),
      plot.margin     = margin(t = 2, r = 2, b = 0, l = 2, unit = "pt")
    ) +
    scale_y_continuous(expand = c(0.02, 0), n.breaks = 5, limits = y_limits)

  ## -- Table p_table --
  p_table <- bin_summary %>%
    select(case, x_num, n_people) %>%
    complete(case, x_num, fill = list(n_people = 0)) %>%
    ggplot(aes(x = x_num, y = case, label = n_people)) +
    geom_text(size = 2, color = "black") +
    scale_y_discrete(
      limits = c("Control", "Cancer"),
      labels = c("Control" = "N controls", "Cancer" = "N cancer")
    ) +
    scale_x_continuous(
      breaks = seq_along(BIN_ORDER_LEFT_TO_RIGHT),
      labels = BIN_ORDER_LEFT_TO_RIGHT,
      expand = expansion(mult = 0.03)
    ) +
    labs(x = NULL, y = NULL) +
    theme_std() +
    theme(
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks  = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(color = c("#4477AA", "#D85A30")),
      plot.margin = margin(t = 0, r = 2, b = 2, l = 2, unit = "pt")
    )

  ## -- Save --
  trajectory_dir <- file.path(base_dir, "trajectory")
  table_dir      <- file.path(base_dir, "sample_size_table")
  for (d in c(trajectory_dir, table_dir)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }

  ## Use ct_of_interest as filename stem (sanitized); pancreas gets a distinct
  ## stem so the exocrine-only run doesn't overwrite anything unexpectedly
  stem <- gsub("[^a-zA-Z0-9]", "_", ct_of_interest)
  if (!is.na(exclude_pattern)) {
    stem <- paste0(stem, "_excl_", gsub("[^a-zA-Z0-9]", "_", exclude_pattern))
  }

  ggsave_wrapper(file.path(trajectory_dir, paste0(stem, "_bmi_trajectory.pdf")),
                 p, width = 4, height = 2.25)
  saveRDS(p, file.path(trajectory_dir, paste0(stem, "_bmi_trajectory.rds")))
  ggsave_wrapper(file.path(table_dir, paste0(stem, "_bmi_sample_size_table.pdf")),
                 p_table, width = 4, height = 0.75)
  saveRDS(p_table, file.path(table_dir, paste0(stem, "_bmi_sample_size_table.rds")))

  list(bin_summary = bin_summary, p = p, p_table = p_table)
}


## --- run across all cancer types -------------------------------------------

results <- pmap(
  cancer_types,
  ~ run_bmi_trajectory(
      ct_of_interest  = ..1,
      exclude_pattern = ..2,
      matched_df      = matched_df,
      base_dir        = base_dir
    )
)
names(results) <- cancer_types$ct_of_interest
