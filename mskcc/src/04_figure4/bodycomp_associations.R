# ============================================================================
# Produces: Figure 4; Supplementary Figure 4; Supplementary Table 20
# (a figure listed without a panel letter contributes to a composite
#  panel whose letter could not be attributed unambiguously)
# ============================================================================
#######################
# Association between body composition and RB1 mutants in Merkel Cell Carcinom (MCC)
#
# In this file, we explore the relationship between RB1 mutants and SAT/VAT/SKM area
#
#
# @author: Ethan Tse
# @date: Jan 13, 2025
#######################

source("./src/00_shared/load_libraries.R")

ct_scan_features <- c("l3sat_area", "l3vat_area", "l3total_muscle_area")

# Load the data
mutation_metadata <- read_tsv("./code_outputs/stable_data/msk/impact/snvs/dataframes/bmi_oncokb_annotated_first_sample_oncogenic_mutations.tsv") %>%
  filter(CANCER_TYPE_DETAILED == "Merkel Cell Carcinoma") %>%
  janitor::clean_names() %>%
  select(sample_id, so_comments, cvr_tmb_score) %>%
  distinct()

mutation_mtx <- read_tsv("./code_outputs/stable_data/msk/impact/snvs/matrix/panel_corrected/bmi_oncokb_annotated_first_sample_oncogenic_mutation_mtx.tsv") %>%
  filter(CANCER_TYPE_DETAILED == "Merkel Cell Carcinoma") %>%
  select(SAMPLE_ID, RB1, TP53) %>%
  janitor::clean_names()

qc_removed <- c("P-0036928", "P-0036676", "P-0039150", "P-0054731", "P-0057490", "P-0060039", "P-0058771", "P-0068753")
body_comp <- read_tsv("./code_outputs/stable_data/msk/body_comp/metadata/master_body_comp_annotated_metadata.tsv") %>%
  filter(cancer_type_detailed == "Merkel Cell Carcinoma") %>%
  filter(!dmp_id %in% qc_removed)

body_comp_mutation_mtx <- body_comp %>%
  left_join(mutation_mtx, by = c("sample_id" = "sample_id")) %>%
  select(sample_id, dmp_id, cancer_type_detailed, age_at_scan, gender, stage_cdm_derived, ancestry_label, surgery_to_diagnosis_days, scan_to_diagnosis_days, smoking_predictions_3_classes, all_of(ct_scan_features), bmi_at_scan, rb1, tp53, clean_scan) %>%
  # add notes
  left_join(mutation_metadata, by = c("sample_id" = "sample_id")) %>%
  mutate(viral_status = ifelse((rb1 == 1 | tp53 == 1) & cvr_tmb_score >= 10, "non_viral", "viral")) %>%
  mutate(l3vatsat_ratio = l3vat_area / l3sat_area)

# Tableau 20 palette
tableau20 <- c(
  "#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F",
  "#8CD17D", "#B6992D", "#F1CE63", "#499894", "#86BCB6",
  "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", "#D37295",
  "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6"
)

viral_status_colours <- c(viral = tableau20[1], non_viral = tableau20[3])

# Barplot: proportion with RB1 mutation, TP53 mutation, TMB >= 10 by BMI category
# Use full mutation matrix (not body_comp subset) so all MCC patients with genomic data are included
mutation_metadata_bmi <- read_tsv("./code_outputs/stable_data/msk/impact/snvs/dataframes/bmi_oncokb_annotated_first_sample_oncogenic_mutations.tsv") %>%
  filter(CANCER_TYPE_DETAILED == "Merkel Cell Carcinoma") %>%
  janitor::clean_names() %>%
  # Covariates come along for the adjusted model below. They are complete on
  # this cohort apart from ancestry_label (6 of 176 missing).
  select(sample_id, cvr_tmb_score, bmi_nearest_result, current_age_deid, gender,
         stage_cdm_derived, smoking_predictions_3_classes, ancestry_label) %>%
  distinct()

# BMI by viral status, split on the genomically-derived viral status this
# analysis uses everywhere else --
# (oncogenic RB1 or TP53) AND TMB >= 10 => non-viral/UV -- rather than on RB1
# alone, which is only one of the two mutations and ignores the TMB criterion.
# Same full-genomics cohort (not the CT subcohort), so it is the whole-cohort
# counterpart of the viral-status body-composition panels below. The RB1
# version above is kept: it is the narrower result this one supersedes in the
# main figure, and is still worth having on disk.
viral_bmi_boxplot_df <- mutation_mtx %>%
  left_join(mutation_metadata_bmi, by = "sample_id") %>%
  filter(!is.na(bmi_nearest_result), !is.na(rb1), !is.na(tp53), !is.na(cvr_tmb_score)) %>%
  mutate(viral_status = factor(
    ifelse((rb1 == 1 | tp53 == 1) & cvr_tmb_score >= 10, "non_viral", "viral"),
    levels = c("non_viral", "viral")
  ))

# The P reported on this panel is the BMI term of a multivariable logistic
# regression of viral status, not the unadjusted two-group test this panel
# used to carry: the comparison is observational, and age and sex in
# particular differ between the virus-positive and virus-negative groups.
#
# Covariates are the project-standard set (age, sex, stage_cdm_derived,
# smoking, ancestry). The Merkel-cell exception that drops smoking and
# ancestry applies to the BODY-COMPOSITION models further down, which run on
# 40 men with 10 in the smaller class; it does not apply here, where the
# cohort is the full genomics set (n = 176, 49 virus-negative) and carries
# both covariates.
#
# logistf (Firth penalized), not glm: with ancestry_label in the model the
# unpenalized fit separates -- ADM/AFR/EAS have 3/2/3 patients and glm
# returns |coef| ~ 15 with SE ~ 2255 and "fitted probabilities numerically 0
# or 1". Firth is also what every other regression in this file uses. The
# BMI estimate itself is stable either way (glm 0.179, P = 1.1e-4; Firth
# 0.164, P = 2.7e-5); the penalization is for the ancestry terms.
viral_bmi_adjusted_fit <- logistf::logistf(
  viral_status ~ bmi_nearest_result + current_age_deid + gender +
    stage_cdm_derived + smoking_predictions_3_classes + ancestry_label,
  data = viral_bmi_boxplot_df,
  family = binomial()
)
viral_bmi_adjusted_res <- tibble(
  term = names(viral_bmi_adjusted_fit$coefficients),
  estimate = unname(viral_bmi_adjusted_fit$coefficients),
  ci_lower = viral_bmi_adjusted_fit$ci.lower,
  ci_upper = viral_bmi_adjusted_fit$ci.upper,
  p_value = viral_bmi_adjusted_fit$prob
) %>%
  mutate(odds_ratio = exp(estimate), n = viral_bmi_adjusted_fit$n)
# First write of the run into this directory. readr does not create missing
# directories, and the figure wrapper -- which does, via create.dir -- may not
# have run yet, so create it here rather than relying on a pre-existing tree.
create_dir_if_needed("./code_outputs/results/tables/body_comp/merkel_cell_carcinoma/")
write_tsv(
  viral_bmi_adjusted_res,
  "./code_outputs/results/tables/body_comp/merkel_cell_carcinoma/viral_status_bmi_adjusted_logistic_regression.tsv"
)
viral_bmi_adjusted_bmi <- viral_bmi_adjusted_res[viral_bmi_adjusted_res$term == "bmi_nearest_result", ]
viral_bmi_adjusted_p <- viral_bmi_adjusted_bmi$p_value
print(viral_bmi_adjusted_res)

# The jitter is seeded so re-running this script does not reshuffle the
# points under a figure that is already in the manuscript.
set.seed(20260101)
viral_bmi_boxplot <- viral_bmi_boxplot_df %>%
  ggplot(aes(x = viral_status, y = bmi_nearest_result, fill = viral_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.65, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.4) +
  # Written as plain "p = ..." text; figs.qmd's relabel_p() rewrites in-panel
  # text labels matching that pattern into the manuscript's italic capital P.
  # Effect size above the P: odds ratio per unit of BMI with its Firth
  # penalized profile-likelihood interval. Two SEPARATE single-line
  # annotations rather than one label with a "\n" in it -- figs.qmd's
  # relabel_p() routes a multi-line in-panel label through a data-column path,
  # which an annotate() layer (whose label is a fixed aes_param, not a mapped
  # column) does not have, and errors with "replacement has 0 rows".
  # size 1.8, not 2: at size 2 the OR line measures 0.834in against this
  # panel's 0.902in of panel width, which the hjust inset then overruns.
  annotate(
    "text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.4, size = 1.8,
    label = sprintf(
      "OR = %.2f (%.2f-%.2f)",
      viral_bmi_adjusted_bmi$odds_ratio,
      exp(viral_bmi_adjusted_bmi$ci_lower),
      exp(viral_bmi_adjusted_bmi$ci_upper)
    )
  ) +
  annotate(
    "text", x = -Inf, y = Inf, hjust = -0.08, vjust = 2.9, size = 1.8,
    label = sprintf("p = %.1e", viral_bmi_adjusted_p)
  ) +
  scale_fill_manual(values = viral_status_colours) +
  # Axis labelled by MCPyV status rather than "UV"/"Viral": the grouping is
  # the virus, and the UV wording named the presumed cause of the
  # virus-negative group rather than the group itself. Signed short forms
  # (not "MCPyV positive"/"negative") because these are tick labels on a
  # panel under an inch wide.
  scale_x_discrete(labels = c(non_viral = "MCPyV\u2212", viral = "MCPyV+")) +
  labs(x = "MCPyV status", y = "BMI", title = "MCC: MCPyV status") +
  theme_std() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

print(table(viral_bmi_boxplot_df$viral_status))

ggsave_wrapper(
  plot = viral_bmi_boxplot,
  filename = "./code_outputs/results/figs/body_comp/merkel_cell_carcinoma/viral_status_bmi_boxplot.pdf",
  width = 2,
  height = 2.5
)

# show the dichotomy between sexes — body composition at L3 by viral status.
#
# Helper to build a faceted viral-status boxplot for an arbitrary set of
# measures. `measures` maps raw column names to display labels (in facet order).
# When `adjust_bmi = TRUE`, each measure is replaced by the residuals of a linear
# fit on bmi_at_scan, i.e. the compartment variation not explained by BMI.
# `show_wilcox = FALSE` suppresses the in-panel two-group test, for panels that
# carry a covariate-adjusted estimate annotated by the caller instead (the BMI
# panel below).
plot_bodycomp_by_viral_status <- function(measures, title, adjust_bmi = FALSE,
                                          show_wilcox = TRUE) {
  df <- body_comp_mutation_mtx %>%
    filter(gender == "Male")

  if (adjust_bmi) {
    df <- df %>%
      mutate(dplyr::across(
        dplyr::all_of(names(measures)),
        ~ stats::resid(stats::lm(.x ~ bmi_at_scan, na.action = stats::na.exclude))
      ))
  }

  df %>%
    select(viral_status, gender, dplyr::all_of(names(measures))) %>%
    pivot_longer(
      cols = dplyr::all_of(names(measures)),
      names_to = "measure",
      values_to = "value"
    ) %>%
    mutate(measure = factor(measures[measure], levels = unname(measures))) %>%
    ggplot(aes(x = viral_status, y = value, fill = viral_status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65, width = 0.6) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.3) +
    (if (show_wilcox) {
      ggpubr::stat_compare_means(
        method = "wilcox.test",
        label = "p.format",
        size = 2,
        label.y.npc = 0.97
      )
    } else {
      NULL
    }) +
    scale_fill_manual(values = viral_status_colours) +
    scale_x_discrete(labels = c(non_viral = "MCPyV\u2212", viral = "MCPyV+")) +
    facet_wrap(~measure, scales = "free_y", nrow = 1) +
    labs(x = NULL, y = "Area", fill = NULL, title = title) +
    theme_std() +
    theme(
      legend.position = "none"
    )
}

# Main figure: SAT, VAT, and SKM
# Supplementary figure: VAT/SAT ratio only
# Supplementary figure: BMI within the CT scan cohort.
#
# This panel is the CT-subcohort check on the whole-cohort BMI association of
# Figure 4e, so it carries the same KIND of statistic: the BMI term of a
# multivariable Firth logistic regression of viral status, not the unadjusted
# two-group test it used to show. The covariate set is this file's MCC
# body-composition set (age + stage), NOT the project-standard set used for
# Figure 4e -- smoking, ancestry and clean_scan are excluded here for exactly
# the reason documented at the adjusted body-composition models below: this is
# the 40-male subcohort with 10 non-viral cases. Sex is not a covariate because
# the subcohort is male by construction.
viral_bmi_scan_fit <- logistf::logistf(
  viral_status ~ bmi_at_scan + age_at_scan + stage_cdm_derived,
  data = body_comp_mutation_mtx %>%
    filter(gender == "Male") %>%
    mutate(
      viral_status = ifelse(viral_status == "viral", 1, 0),
      age_at_scan = as.numeric(age_at_scan)
    ),
  control = logistf::logistf.control(maxit = 1000, maxstep = 0.5),
  plcontrol = logistf::logistpl.control(maxit = 1000)
)
viral_bmi_scan_res <- tibble(
  term = names(viral_bmi_scan_fit$coefficients),
  estimate = unname(viral_bmi_scan_fit$coefficients),
  ci_lower = viral_bmi_scan_fit$ci.lower,
  ci_upper = viral_bmi_scan_fit$ci.upper,
  p_value = viral_bmi_scan_fit$prob
) %>%
  mutate(odds_ratio = exp(estimate), n = viral_bmi_scan_fit$n)
write_tsv(
  viral_bmi_scan_res,
  "./code_outputs/results/tables/body_comp/merkel_cell_carcinoma/viral_status_bmi_at_scan_logistic_regression.tsv"
)
viral_bmi_scan_bmi <- viral_bmi_scan_res[viral_bmi_scan_res$term == "bmi_at_scan", ]
print(viral_bmi_scan_res)

bmi_scan_boxplot <- plot_bodycomp_by_viral_status(
  measures = c(bmi_at_scan = "BMI (kg/m²)"),
  title = "BMI at scan of male MCC patients",
  show_wilcox = FALSE
) +
  # Same two-annotation form, inset and sizes as Figure 4e's panel: one
  # single-line annotate() each, because figs.qmd's relabel_p() cannot rewrite
  # a multi-line annotate() label (see the note at that panel).
  annotate(
    "text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.4, size = 1.8,
    label = sprintf(
      "OR = %.2f (%.2f-%.2f)",
      viral_bmi_scan_bmi$odds_ratio,
      exp(viral_bmi_scan_bmi$ci_lower),
      exp(viral_bmi_scan_bmi$ci_upper)
    )
  ) +
  annotate(
    "text", x = -Inf, y = Inf, hjust = -0.08, vjust = 2.9, size = 1.8,
    label = sprintf("p = %.1e", viral_bmi_scan_bmi$p_value)
  )

ggsave_wrapper(
  plot = bmi_scan_boxplot,
  filename = "./code_outputs/results/figs/body_comp/merkel_cell_carcinoma/bmi_at_scan_viral_status_boxplot_supplementary.pdf",
  width = 1.3,
  height = 1.5
)

# BMI-adjusted figure: SAT, VAT, SKM residuals after regressing out BMI
# we focus our analysis on males only
# fit one logistic regression model per tissue, then a final joint model
male_data <- body_comp_mutation_mtx %>%
  filter(gender == "Male") %>%
  select(viral_status, all_of(ct_scan_features), age_at_scan, stage_cdm_derived, clean_scan) %>%
  mutate(age_at_scan = as.numeric(age_at_scan)) %>%
  mutate(viral_status = ifelse(viral_status == "viral", 1, 0))

# try firth's regression
male_data_firths <- male_data %>%
  mutate(
    # Standardized (SD) rather than per 10 cm2. A fixed 10 cm2 step is a very
    # different fraction of each compartment's spread in this cohort (SAT SD
    # 85.0, VAT SD 112.6, SKM SD 26.9 cm2), so per-10cm2 log-ORs are not
    # comparable across compartments -- it compresses VAT toward zero and
    # inflates SKM ~4x relative to VAT purely through the unit choice. Per SD
    # also matches every other body-composition model in this study.
    across(
      all_of(ct_scan_features),
      ~ as.numeric(scale(.x)),
      .names = "{.col}_z"
    )
  )

fit_firth <- function(vars) {
  formula <- as.formula(
    paste(
      "viral_status ~",
      paste(c(vars, "age_at_scan", "stage_cdm_derived"), collapse = " + ")
    )
  )

  logistf::logistf(
    formula = formula,
    data = male_data_firths,
    control = logistf::logistf.control(maxit = 1000, maxstep = 0.5),
    plcontrol = logistf::logistpl.control(maxit = 1000)
  )
}

models <- list(
  l3sat = fit_firth("l3sat_area_z"),
  l3vat = fit_firth("l3vat_area_z"),
  l3muscle = fit_firth("l3total_muscle_area_z"),
  joint_l3 = fit_firth(c("l3sat_area_z", "l3vat_area_z", "l3total_muscle_area_z"))
)

tidy_logistf <- function(model) {
  data.frame(
    term = names(model$coefficients),
    estimate = model$coefficients,
    std.error = sqrt(diag(model$var)),
    # Penalized profile-likelihood limits on the log-odds scale, taken from
    # logistf's own ci.lower/ci.upper. These are the intervals that correspond
    # to the penalized likelihood-ratio p-values in model$prob; Wald limits
    # (estimate +/- 1.96 * SE) are symmetric by construction and understate the
    # upper bound, because at this event count (viral = 30, non-viral = 10) the
    # likelihood is strongly right-skewed. Where the profile does not close on
    # the upper side (near-separation), ci.upper is very large or infinite --
    # that is the intended behaviour, not a failure to converge.
    ci.lower = model$ci.lower,
    ci.upper = model$ci.upper,
    p.value = model$prob
  )
}

# Shared forest-plot style, matching the manuscript body-composition forest
# plots: coefficients shown on the log-odds scale, coloured by nominal
# significance, with marginal and joint models separated into their own panels.
# Expects columns: estimate, ci.lower, ci.upper, metric_label, model_type,
# significance (ci.lower/ci.upper on the log-odds scale, from tidy_logistf).
forest_sig_colours <- SIG_COLOURS

# Penalized profile limits are not guaranteed to close on both sides: under
# near-separation the profile never falls by qchisq(0.95, 1) / 2 in one
# direction and logistf returns a bound of +/-Inf (or an astronomically large
# finite value). Rescaling the panel to such a bound would flatten every other
# interval into a dot, so intervals are drawn clipped at +/-log(clip_or) and the
# clipped end is capped with an arrowhead, which is the usual forest-plot
# convention for "this bound is open". The unclipped values are what get written
# to the results TSVs.
plot_firth_forest <- function(df, facet = facet_wrap(~model_type, ncol = 1),
                              x_lab = "log(OR) per SD", arrow_labels = NULL,
                              clip_or = 5, center_zero = TRUE) {
  # clip_or is a floor, not a fixed cut. The clip must never fall inside a point
  # estimate: the open-bound arrow is drawn as a segment from estimate to the
  # clipped bound, so if the estimate lies beyond the clip that segment reverses
  # and renders as an arrowhead pointing *inward*, which reads as the interval
  # running the wrong way. (This is exactly what happened when these models moved
  # from per-10cm2 to per-SD units: log-ORs grew from ~0.1-0.5 to ~0.9-3.8, while
  # the default clip stayed at log(5) = 1.61 and truncated every row.) Expand the
  # clip past the largest finite estimate (with headroom, so an arrowhead on the
  # largest row still has room to render) and let the panel scale with the
  # units rather than silently mis-drawing.
  # The headroom is a fraction of the panel's own extent, not of the estimate,
  # so the gap between a clipped row's point and its arrowhead stays visible
  # regardless of scale. A multiplicative margin on the estimate alone is not
  # enough: when the largest estimate is also the clipped row (as here for SAT),
  # the marker and the arrowhead end up a few points apart and read as one blob.
  # Clip ONLY bounds that are not finite. Clipping a large-but-finite bound is
  # actively misleading here: these profile-likelihood intervals are strongly
  # right-skewed under near-separation (SAT: estimate 3.81, CI [1.32, 10.37],
  # i.e. 2.5 left and 6.6 right of the point), and truncating the long side
  # makes the drawn interval look *left*-skewed -- the opposite of the truth.
  # A finite bound, however wide, is real information about how uncertain the
  # upper end is, and at these magnitudes the other rows still occupy ~17-20%
  # of the panel rather than collapsing to dots.
  finite_bounds <- c(df$ci.lower, df$ci.upper)
  finite_bounds <- finite_bounds[is.finite(finite_bounds)]
  clip_at <- max(log(clip_or), if (length(finite_bounds)) max(abs(finite_bounds)) else 0)

  df <- df %>%
    mutate(
      model_type = factor(
        model_type,
        levels = c("Marginal", "Joint", "Marginal (BMI adj)", "Joint (BMI adj)")
      ),
      ci_lower_draw = pmax(ci.lower, -clip_at),
      ci_upper_draw = pmin(ci.upper, clip_at),
      open_upper = ci.upper > clip_at,
      open_lower = ci.lower < -clip_at
    )

  p <- df %>%
    ggplot(aes(x = estimate, y = metric_label)) +
    ggstats::geom_stripped_rows() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1) +
    geom_linerange(
      aes(
        xmin = ci_lower_draw, xmax = ci_upper_draw,
        color = significance
      ),
      lineend = "round", linewidth = 0.1
    ) +
    geom_segment(
      data = ~ dplyr::filter(.x, open_upper),
      aes(x = estimate, xend = ci_upper_draw, y = metric_label, yend = metric_label,
          colour = significance),
      linewidth = 0.1, arrow = arrow(length = unit(0.04, "in"), type = "closed")
    ) +
    geom_segment(
      data = ~ dplyr::filter(.x, open_lower),
      aes(x = estimate, xend = ci_lower_draw, y = metric_label, yend = metric_label,
          colour = significance),
      linewidth = 0.1, arrow = arrow(length = unit(0.04, "in"), type = "closed")
    ) +
    # size 2, matching the PI3K forest in the next figure
    # (pi3k_forest_crc_pdac_rcc.rds), so every forest uses one dot size.
    geom_point(aes(colour = significance), size = 2) +
    scale_colour_manual(name = NULL, values = forest_sig_colours, drop = FALSE) +
    significance_key(size = 2) +
    facet +
    labs(x = x_lab, y = NULL) +
    theme_std() +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks.y = element_blank(),
      strip.text = element_text(size = 6),
      panel.spacing.x = unit(1, "pt")
    )

  if (!is.null(arrow_labels)) {
    # Directional indicator (styled like Figure 5A's volcano-plot arrows) in
    # the plot's subtitle slot -- see element_arrow_subtitle() in styles.R --
    # rather than forest_arrows()'s in-panel annotations, so it isn't clipped
    # when H/I's top margin gets trimmed down for the composed Figure 5
    # layout, and can't crowd the cowplot panel letters on either side.
    # Drawn (clipped) bounds set the panel width, so one open-ended interval
    # cannot squash the rest of the forest.
    # x range, either centred on zero or data-driven.
    # center_zero = TRUE forces +/- max|bound|, putting the reference line in
    # the middle of the panel so the Lower|Higher indicator reads as a balanced
    # two-sided scale and the eye can judge distance from the null the same way
    # on both sides. The cost is empty panel: these effects are almost entirely
    # one-sided (l3 panel spans -0.97 to +15.69), so centring nearly doubles
    # the drawn width and pushes every interval into the right half.
    # center_zero = FALSE lets the range follow the data and hands the true
    # zero position to element_arrow_subtitle()'s zero_frac, so the Lower|Higher
    # split still lands on the reference line without reserving space for a
    # side that has no data.
    est_bounds <- c(df$ci_lower_draw, df$ci_upper_draw, df$estimate)
    est_bounds <- est_bounds[is.finite(est_bounds)]
    lo <- min(c(est_bounds, 0), na.rm = TRUE)
    hi <- max(c(est_bounds, 0), na.rm = TRUE)
    pad <- 0.08 * (hi - lo)
    if (center_zero) {
      x_lim <- max(abs(c(lo, hi))) + pad
      x_lo <- -x_lim
      x_hi <- x_lim
    } else {
      x_lo <- lo - pad
      x_hi <- hi + pad
    }
    zero_frac <- (0 - x_lo) / (x_hi - x_lo)

    p <- p +
      labs(subtitle = "direction") +
      scale_x_continuous(limits = c(x_lo, x_hi)) +
      theme(
        # H/I sit in a narrow third-width column (see figs.qmd), so the
        # subtitle is sized down from the default to fit without the two
        # halves overlapping.
        plot.subtitle = element_arrow_subtitle(labels = arrow_labels, text_size = 5.5, height = unit(0.28, "in"), zero_frac = zero_frac)
      )
  }

  p
}


res_df <- purrr::imap_dfr(models, function(model, name) {
  tidy_logistf(model) %>%
    dplyr::filter(grepl("l3|l1|t10|t12", term, ignore.case = TRUE)) %>%
    dplyr::mutate(
      model = name,
      odds_ratio = exp(estimate),
      conf.low = exp(ci.lower),
      conf.high = exp(ci.upper)
    )
}) %>%
  dplyr::mutate(
    level = dplyr::case_when(
      grepl("l1", model, ignore.case = TRUE) ~ "L1",
      grepl("l3", model, ignore.case = TRUE) ~ "L3",
      grepl("t10", model, ignore.case = TRUE) ~ "T10",
      grepl("t12", model, ignore.case = TRUE) ~ "T12",
      TRUE ~ "Other"
    )
  ) %>%
  # BH within each anatomical level: Figure 4f is the L3 panel, so its family is
  # exactly the six tests it shows -- the three marginal compartment models plus
  # the three terms of the joint model. Figure 4g's BMI-adjusted models are a
  # separate family, corrected on their own below.
  dplyr::group_by(level) %>%
  dplyr::mutate(padj = p.adjust(p.value, method = "BH")) %>%
  dplyr::ungroup()

# save results of firth's regression
create_dir_if_needed("./code_outputs/results/tables/body_comp/merkel_cell_carcinoma/")
write_tsv(res_df, "./code_outputs/results/tables/body_comp/merkel_cell_carcinoma/firths_logistic_regression_results.tsv")

# despite estimation instability due to low N = 40, viral = 30, non-viral = 10

# visualize firth's regression results
plot_df <- res_df %>%
  mutate(
    model_type = ifelse(grepl("^joint_", model), "Joint", "Marginal"),
    significance = sig_level(p.value, padj),
    metric_label = dplyr::case_when(
      grepl("sat", term) ~ "SAT",
      grepl("vat", term) ~ "VAT",
      grepl("muscle", term) ~ "SKM"
    )
  )

# primary plot visualizing L3 features (marginal vs joint separated by panel)
firth_logistic_regression_forest_plot_l3 <- plot_df %>%
  filter(level == "L3") %>%
  plot_firth_forest(x_lab = "Adj. log(OR) per SD", arrow_labels = c("Lower", "Higher"))

ggsave_wrapper(
  plot = firth_logistic_regression_forest_plot_l3,
  filename = "./code_outputs/results/figs/body_comp/merkel_cell_carcinoma/firths_logistic_regression_forest_plot_l3.pdf",
  width = 2.0,
  height = 3.3
)

# Fully adjusted L3 model: the same covariates as fit_firth above (age, stage;
# sex is fixed as this cohort is male-only) plus BMI at scan, so this is exactly
# the marginal model with one term added -- which is what the figure caption
# claims ("as in h, additionally adjusted for BMI at scan").
#
# clean_scan is deliberately NOT a covariate in any MCC model. It is missing for
# 2 of the 40 male MCC patients, so listwise deletion would drop them -- one of
# only 10 non-viral cases -- costing 5% of an analysis that is already
# separation-limited. Smoking and ancestry are excluded for the same reason of
# scale: with viral = 30 and non-viral = 10 the model cannot support them, and
# neither the methods nor the caption claims them. Firth penalization is
# retained given that small, imbalanced N.
male_data_adjusted <- body_comp_mutation_mtx %>%
  filter(gender == "Male") %>%
  select(
    viral_status, l3sat_area, l3vat_area, l3total_muscle_area,
    age_at_scan, stage_cdm_derived, smoking_predictions_3_classes,
    ancestry_label, bmi_at_scan, clean_scan
  ) %>%
  mutate(
    age_at_scan = as.numeric(age_at_scan),
    viral_status = ifelse(viral_status == "viral", 1, 0),
    l3sat_area_z = as.numeric(scale(l3sat_area)),
    l3vat_area_z = as.numeric(scale(l3vat_area)),
    l3total_muscle_area_z = as.numeric(scale(l3total_muscle_area))
  )

adjust_covariates <- c(
  "age_at_scan", "stage_cdm_derived", "bmi_at_scan"
)

fit_firth_adjusted <- function(vars) {
  formula <- as.formula(
    paste(
      "viral_status ~",
      paste(c(vars, adjust_covariates), collapse = " + ")
    )
  )
  logistf::logistf(
    formula = formula,
    data = male_data_adjusted,
    control = logistf::logistf.control(maxit = 1000, maxstep = 0.5),
    plcontrol = logistf::logistpl.control(maxit = 1000)
  )
}

models_adjusted <- list(
  l3sat = fit_firth_adjusted("l3sat_area_z"),
  l3vat = fit_firth_adjusted("l3vat_area_z"),
  l3muscle = fit_firth_adjusted("l3total_muscle_area_z"),
  joint_l3 = fit_firth_adjusted(c("l3sat_area_z", "l3vat_area_z", "l3total_muscle_area_z"))
)

res_df_adjusted <- purrr::imap_dfr(models_adjusted, function(model, name) {
  tidy_logistf(model) %>%
    dplyr::filter(grepl("l3", term)) %>%
    dplyr::mutate(
      model = name,
      model_type = ifelse(grepl("^joint_", name), "Joint (BMI adj)", "Marginal (BMI adj)"),
      odds_ratio = exp(estimate),
      conf.low = exp(ci.lower),
      conf.high = exp(ci.upper),
      metric_label = dplyr::case_when(
        grepl("sat", term) ~ "SAT",
        grepl("vat", term) ~ "VAT",
        grepl("muscle", term) ~ "SKM"
      )
    )
}) %>%
  dplyr::mutate(padj = p.adjust(p.value, method = "BH"),
                significance = sig_level(p.value, padj))

write_tsv(res_df_adjusted, "./code_outputs/results/tables/body_comp/merkel_cell_carcinoma/firths_logistic_regression_l3_fully_adjusted_results.tsv")

firth_adjusted_forest_plot_l3 <- res_df_adjusted %>%
  plot_firth_forest(x_lab = "Adj. log(OR) per SD", arrow_labels = c("Lower", "Higher"))

ggsave_wrapper(
  plot = firth_adjusted_forest_plot_l3,
  filename = "./code_outputs/results/figs/body_comp/merkel_cell_carcinoma/firths_logistic_regression_l3_fully_adjusted_forest_plot.pdf",
  width = 2.0,
  height = 3.3
)
