# ============================================================================
# Produces: Figure 2
# (a figure listed without a panel letter contributes to a composite
#  panel whose letter could not be attributed unambiguously)
# ============================================================================
####################
# Figure 2, panels F-H: one 2x2 contingency table per cancer type with a
# significant myopenia x myosteatosis association (colorectal, lung and
# pancreatic adenocarcinoma).
#
# Each panel is rendered as a 2x2 tile grid (condition-present/absent on each
# axis), cell fill = % of that cancer type's cohort, cell label = n (%).
#
# The annotation (log OR, 95% CI, padj) comes from the adjusted logistic
# regression in myopenia_myosteatosis_by_cancer_type_mixed_model.R
# (myosteatosis ~ myopenia + age + sex + stage + ancestry + smoking, fit
# within each cancer type), rather than a raw Fisher's exact test on the 2x2
# table: that model adjusts for confounders and tests this cancer type
# alongside every other major cancer type (hence the BH-adjusted p, "padj",
# shown here). All three panels share one teal fill -- same pair, same
# encoding -- so a per-panel accent would wrongly imply they measure
# different things.
#
# Cohort and abnormal-physiology definitions (myopenia, myosteatosis) are
# copied from abnormal_physiology.R, so panel E's prevalence bars and these
# panels are read from an identical set of per-scan calls (one row per scan,
# matching the rest of Figure 2).
#
# @author: Ethan Tse
# @date: August 17th, 2026
####################


source("./src/00_shared/load_libraries.R")

body_comp_data <- read_tsv("./code_outputs/stable_data/msk/body_comp/metadata/non_clean_master_body_comp_annotated_metadata.tsv") %>%
  filter(cancer_type_detailed != "Merkel Cell Carcinoma") %>%
  # recalibrate to pre-contrast values (consistent with abnormal_physiology.R)
  mutate(
    l3total_muscle_area_precontrast = 1.0628 * l3total_muscle_area,
    l3total_muscle_median_hu_precontrast = 0.9930 * l3total_muscle_median_hu - 7.5291,
    liver_median_hu_precontrast = case_when(
      (0.7434 * liver_median_hu - 12.387) < 80 ~ (0.7434 * liver_median_hu) - 12.387,
      TRUE ~ (0.21254 * liver_median_hu) + 31.508,
    ),
  ) %>%
  mutate(SMI = l3total_muscle_area_precontrast / (height_at_scan_cm / 100)^2) %>%
  mutate(myopenia = case_when(
    # Height is missing for ~8% of scans, so SMI is NA for them. Without
    # this guard the final `TRUE ~ 0` silently calls every one of them
    # "not myopenic", inflating the denominator and biasing prevalence down.
    is.na(SMI) ~ NA_real_,
    (gender == "Male") & (SMI < 50) ~ 1,
    (gender == "Female") & (SMI < 39) ~ 1,
    TRUE ~ 0
  )) %>%
  # myosteatosis is defined relative to the muscle radiodensity distribution,
  # as the sex-specific lowest quartile of pre-contrast muscle median HU,
  # computed cohort-wide (not per cancer type) -- consistent with every other
  # script that derives this flag
  mutate(myosteatosis = case_when(
    # Same NA guard as myopenia: an unmeasured muscle HU is unknown
    # myosteatosis status, not absence of myosteatosis.
    is.na(l3total_muscle_median_hu_precontrast) ~ NA_real_,
    (gender == "Male") & l3total_muscle_median_hu_precontrast < quantile(l3total_muscle_median_hu_precontrast[gender == "Male"], 0.25, na.rm = TRUE) ~ 1,
    (gender == "Female") & l3total_muscle_median_hu_precontrast < quantile(l3total_muscle_median_hu_precontrast[gender == "Female"], 0.25, na.rm = TRUE) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(liver_steatosis = ifelse(liver_median_hu_precontrast < 57, 1, 0))

out_dir <- "./code_outputs/results/figs/body_comp/abnormal_physiology/contingency_tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Format an OR (95% CI) + p-value subtitle. ASCII hyphen, not an en dash --
# the base PDF device's default font falls back to substituting "-" for
# U+2013 anyway (with a warning) since it's not in its single-byte encoding.
format_or_subtitle <- function(or_est, ci_low, ci_high, p_val, p_label_prefix = "p") {
  or_label <- sprintf("OR = %.2f (95%% CI %.2f-%.2f)", or_est, ci_low, ci_high)
  # Small p-values are reported as the actual number in scientific notation
  # rather than floored to "< 0.001". The CRC myosteatosis x liver-steatosis
  # co-occurrence (panel G) is padj = 3.2e-26, and "< 0.001" threw away 23
  # orders of magnitude of what is the strongest co-occurrence in the sweep.
  # The %.3f branch still covers everything down to 0.001 (panel F, 0.005).
  # %.3f floored the myopenia x myosteatosis padj (0.00480) to "0.005" while
  # the text quotes 0.0048 -- two different-looking numbers for one value.
  # Three branches: scientific below 1e-3, 4 d.p. below 0.01, 3 d.p. above.
  p_label <- dplyr::case_when(
    p_val < 0.001 ~ sprintf("%s = %.1e", p_label_prefix, p_val),
    p_val < 0.01 ~ sprintf("%s = %.4f", p_label_prefix, p_val),
    TRUE ~ sprintf("%s = %.3f", p_label_prefix, p_val)
  )
  paste0(or_label, ", ", p_label)
}

# Same annotation on the LOG odds scale, for the panels whose estimate comes
# from the adjusted per-cancer-type logistic regression -- that model is fit
# on the log odds scale and its results table carries both scales side by side
# (estimate/conf.low/conf.high in log odds, odds_ratio/or_conf_* exponentiated),
# so this reads the log-scale columns rather than transforming the OR here.
#
# Two lines, split estimate / interval+p: at 5pt the one-line form is 1.86in,
# wider than the 1.4in these panels now get in Figure 2 (and a subtitle wider
# than its cell bleeds into the neighbouring panel rather than clipping).
# The split puts 0.64in on the first line and 1.20in on the second.
# " to " separates the interval bounds, not the "-" the OR form uses: both
# bounds are negative here, and "-1.03--0.25" is unreadable.
format_logor_subtitle <- function(est, ci_low, ci_high, p_val, p_label_prefix = "padj") {
  p_label <- dplyr::case_when(
    p_val < 0.001 ~ sprintf("%s = %.1e", p_label_prefix, p_val),
    p_val < 0.01 ~ sprintf("%s = %.4f", p_label_prefix, p_val),
    TRUE ~ sprintf("%s = %.3f", p_label_prefix, p_val)
  )
  sprintf("Adj. log(OR) = %.2f\n(95%% CI %.2f to %.2f), %s", est, ci_low, ci_high, p_label)
}

# Build a 2x2 contingency-table tile plot for var1 x var2 within `data`. The
# cell tiles/counts always come from the raw 2x2 table; the subtitle
# annotation defaults to an unadjusted two-sided Fisher's exact test on that
# same table, but can be overridden with `subtitle` (e.g. an adjusted
# mixed-model OR/CI/p pulled from an already-fit model) when one is
# available for this specific pair -- see format_or_subtitle().
# fill_color: the single accent hue used for the white -> fill_color cell
# scale (chosen per-panel below so the two panels are visually distinct at a
# glance while both using the same white-to-saturated encoding of cell %).
build_contingency_plot <- function(data, var1, var2, var1_label, var2_label,
                                    title, fill_color, subtitle = NULL) {
  d <- data %>%
    filter(!is.na(.data[[var1]]), !is.na(.data[[var2]]))

  if (is.null(subtitle)) {
    fisher_result <- fisher.test(table(d[[var1]], d[[var2]]))
    subtitle <- format_or_subtitle(
      unname(fisher_result$estimate),
      fisher_result$conf.int[1], fisher_result$conf.int[2],
      fisher_result$p.value
    )
  }

  cell_counts <- as.data.frame(table(d[[var1]], d[[var2]])) %>%
    rename(v1 = Var1, v2 = Var2, n = Freq) %>%
    mutate(
      pct = n / sum(n),
      v1_status = factor(ifelse(v1 == "1", "Yes", "No"), levels = c("Yes", "No")),
      v2_status = factor(ifelse(v2 == "1", "Yes", "No"), levels = c("No", "Yes"))
    )

  ggplot(cell_counts, aes(x = v2_status, y = v1_status)) +
    geom_tile(aes(fill = pct), color = "black", linewidth = 0.2) +
    geom_text(aes(label = sprintf("n = %d\n(%.1f%%)", n, pct * 100)), size = 2.2) +
    scale_fill_gradient(low = "white", high = fill_color, limits = c(0, 1), guide = "none") +
    labs(
      title = title,
      x = var2_label,
      y = var1_label,
      subtitle = subtitle
    ) +
    theme_std() +
    theme(
      plot.title = element_text(size = 6, hjust = 0.5, lineheight = 1.1),
      plot.subtitle = element_text(size = 5, hjust = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.1),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(2, 2, 2, 2, "mm")
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_fixed()
}

build_and_save <- function(data, var1, var2, var1_label, var2_label, title,
                            fill_color, out_name, subtitle = NULL) {
  p <- build_contingency_plot(data, var1, var2, var1_label, var2_label, title, fill_color, subtitle = subtitle)
  saveRDS(p, file.path(out_dir, paste0(out_name, ".rds")))
  ggsave_wrapper(
    plot = p,
    filename = file.path(out_dir, paste0(out_name, ".pdf")),
    width = 3,
    height = 3.3,
    extension = "pdf"
  )
  invisible(p)
}


# Panels F-H's log(OR)/CI/padj: pulled from the per-cancer-type adjusted
# logistic regressions (myopenia_myosteatosis_by_cancer_type_mixed_model.R),
# which must be run first. Reported padj is BH-adjusted across all major
# cancer types tested there, not just the one being plotted.
mixed_model_dir <- "./code_outputs/results/tables/body_comp/abnormal_physiology"

# One 2x2 per significant myopenia x myosteatosis result, not just the
# pancreatic one. The text reports three (colorectal, lung, pancreatic) and
# previously only pancreatic had a panel, so the two others were numbers in
# a sentence with nothing to look at. All three share the same teal fill
# (the myopenia-side "Lower" arrow color from forest_arrows(), styles.R):
# same pair, same encoding, so a per-panel accent would wrongly imply the
# panels measure different things.
myopenia_myosteatosis_cancers <- c(
  "Colorectal Adenocarcinoma", "Lung Adenocarcinoma", "Pancreatic Adenocarcinoma"
)

build_myopenia_myosteatosis_panel <- function(cancer_type, model_file, out_suffix = "",
                                               data = body_comp_data) {
  model_result <- read_tsv(file.path(mixed_model_dir, model_file), show_col_types = FALSE) %>%
    filter(cancer_type_detailed == !!cancer_type)

  build_and_save(
    data %>% filter(cancer_type_detailed == !!cancer_type),
    "myopenia", "myosteatosis",
    var1_label = "Myopenia", var2_label = "Myosteatosis",
    title = paste0("Myopenia x Myosteatosis\n(", cancer_type, ")"),
    fill_color = "#2c7fb8",
    out_name = paste0(
      "myopenia_myosteatosis_",
      tolower(unname(cancer_codes[cancer_type])), "_contingency", out_suffix
    ),
    subtitle = format_logor_subtitle(
      model_result$estimate,
      model_result$conf.low, model_result$conf.high,
      model_result$p_adj_bh
    )
  )
}

myopenia_myosteatosis_panels <- lapply(
  myopenia_myosteatosis_cancers,
  build_myopenia_myosteatosis_panel,
  model_file = "myopenia_myosteatosis_by_cancer_type_mixed_model.tsv"
)
names(myopenia_myosteatosis_panels) <- myopenia_myosteatosis_cancers

# out_name derives from cancer_codes, so the pancreatic panel keeps writing
# to myopenia_myosteatosis_pdac_contingency.rds -- the path figs.qmd has
# read since before the colorectal and lung panels existed.
