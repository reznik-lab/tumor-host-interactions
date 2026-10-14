# ============================================================================
# Produces: Supplementary Figure 5 (G, H)
# ============================================================================
################
# Supplementary Figure 5 g/h: EGFR-mutant LUAD body-composition panels
#
# Builds the two EGFR-mutant forests that sit side by side in Supplementary
# Figure 5, as a matched pair:
#   g  PI3K/AKT/mTOR defined INCLUDING STK11
#   h  PI3K/AKT/mTOR defined EXCLUDING STK11
# Both carry the same three model arms -- BMI-at-scan adjustment, the
# treatment-proximal restriction, and the two together -- so the pair is a
# like-for-like grid: one column per arm, one panel per pathway definition.
# They are built together so the pair shares one visual language (same row
# order, same arm wording, same column count) because the paragraph in
# text.qmd is arguing that the two definitions agree, and that only reads if
# the panels are directly comparable.
#
# The "Unadjusted" arm h used to open with is dropped: with g showing only
# adjusted arms the two panels had different column counts, and the
# unadjusted fit is the one estimate here that no sentence in text.qmd
# rests on.
#
# Panel g is a narrowed companion to Plot 4b of lung_pik3ca_bodycomp_analysis.R
# (pi3k_forest_luad_l3_clean_scan_arms.rds), which shows the two treatment-proximal arms
# across all three LUAD driver subtypes. That full three-subtype version is
# deliberately left in place and is still written by that script -- nothing
# here replaces it. The manuscript never cites a number from its KRAS-mutant or
# KRAS/EGFR-wild-type facets: the only treatment-proximal claim in text.qmd is
# about EGFR-mutant LUAD, and the "nothing in the other two subtypes" control
# it does make is cited to main Figure 5h, which fits those subtypes in the far
# larger all-scans cohort. The treatment-proximal cohort is small (KRAS-mut n = 63,
# KRAS/EGFR WT n = 140), so its versions of that control add width, not evidence.
#
# Panel h is a relabelled companion to the forest in
# lung_pik3ca_bodycomp_stk11_exclusion_sensitivity.R
# (pi3k_forest_luad_egfr_stk11_excluded_arms.rds), which is likewise left in
# place. The only changes are the dropped Unadjusted arm and the house x-axis
# wording.
#
# Naming the driver subtype is part of the point of rebuilding both: once g is
# restricted to EGFR-mutant, neither panel otherwise says which subtype it is
# (the old g carried it in its column strips, h never carried it at all). It
# was a right-hand facet strip until the strip column was judged to cost more
# width than it bought at this panel size; it now rides on the x-axis title,
# which both panels need anyway to name their pathway definition.
#
# No model is fit here. Per the repo convention this reads the already-computed
# results tables and only changes presentation; the estimates are identical to
# those in the two panels above.
#
# @author: Ethan Tse
# @date: October 5th, 2026
################

source("./src/00_shared/load_libraries.R")
library(ggtext)

tbl_dir <- "./code_outputs/results/tables/body_comp/pi3k"
fig_dir <- "./code_outputs/results/figs/body_comp/pi3k"

# Rendered by axis.title.x = element_markdown() below. Kept as one constant so
# the two panels cannot drift apart.
EGFR_LABEL <- "*EGFR*<sup>mut</sup> LUAD"

L3_TISSUE <- "^l3(sat_area|vat_area|total_muscle_area)$"

# Arm wording shared by both panels, ordered as the Results sentence reads
# them: each restriction layered on in turn. Both source tables also hold an
# "unadjusted" fit, which PLOT_MODELS drops before arm_label() is ever
# reached -- it is dropped by the filter rather than left out of the
# case_when, which would silently map those rows to NA arms instead.
ARM_LEVELS <- c("+ BMI at scan", "Tx-proximal", "Tx-proximal\n+ BMI at scan")
PLOT_MODELS <- c("bmi_adjusted", "clean_scan", "clean_scan_bmi_adjusted")
arm_label <- function(model) {
  factor(
    case_when(
      model == "bmi_adjusted" ~ "+ BMI at scan",
      model == "clean_scan" ~ "Tx-proximal",
      model == "clean_scan_bmi_adjusted" ~ "Tx-proximal\n+ BMI at scan"
    ),
    levels = ARM_LEVELS
  )
}

prep_panel <- function(df) {
  df %>%
    filter(grepl(L3_TISSUE, var), model %in% PLOT_MODELS) %>%
    mutate(
      arm = arm_label(model),
      metric_label = case_when(
        grepl("total_muscle_area", var) ~ "SKM",
        grepl("vat_area", var) ~ "VAT",
        grepl("sat_area", var) ~ "SAT"
      ),
      # reversed so rows read top-to-bottom SKM -> VAT -> SAT (y pos 1 = bottom)
      metric_label = factor(metric_label, levels = rev(c("SKM", "VAT", "SAT"))),
      significance = sig_level(p.value, padj)
    ) %>%
    droplevels()
}

# Visual language mirrors make_pi3k_forest() in lung_pik3ca_bodycomp_analysis.R
# (same row order, same significance colours, same stripped rows) so these read
# as the same analysis as the main-figure LUAD forests.
make_egfr_forest <- function(df, xlab) {
  ggplot(df, aes(x = estimate, y = metric_label)) +
    ggstats::geom_stripped_rows() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1) +
    geom_linerange(
      aes(xmin = conf.low, xmax = conf.high, color = significance),
      lineend = "round", linewidth = 0.1
    ) +
    geom_point(size = 2, shape = 16, aes(colour = significance)) +
    scale_colour_manual(
      name = NULL,
      # same keys and colours as every other significance forest in the figure
      values = SIG_COLOURS,
      drop = FALSE
    ) +
    significance_key(size = 2) +
    # One column per model arm, no row facet: the driver subtype is named in
    # the x-axis title instead of in a right-hand strip, which at this panel
    # size (~2.6in for three facets) cost more width than it was worth.
    facet_grid(cols = vars(arm)) +
    labs(x = xlab, y = NULL, fill = NULL) +
    theme_std() +
    theme(
      legend.position = "bottom",
      axis.ticks.y = element_blank(),
      # Strips must be element_text: gridtext reads the leading "+" of
      # "+ BMI at scan" as a markdown bullet list and dies on the <ul> it
      # generates. Only the axis title, which holds the italic/superscript
      # gene name, gets element_markdown.
      strip.text.x = element_text(size = 6),
      axis.title.x = element_markdown(size = 6, lineheight = 1.1),
      panel.spacing.x = unit(1, "pt")
    )
}

# ── Panel g: BMI / treatment-proximal arms, STK11 included ────────────────────────────
clean_scan_arms_egfr <- read_tsv(
  file.path(tbl_dir, "pi3k_bodycomp_firth_luad_subtypes.tsv"),
  show_col_types = FALSE
) %>%
  filter(cancer_type_detailed == "EGFR-mut") %>%
  prep_panel()

stopifnot(nrow(clean_scan_arms_egfr) == 9)

pi3k_forest_l3_clean_arms_egfr <- make_egfr_forest(
  clean_scan_arms_egfr,
  # The arm strips already say "Tx-proximal", so the axis title carries what the
  # strips do not: the subtype both panels are restricted to, and the pathway
  # definition that is the one thing separating this panel from h.
  xlab = paste0("Adj. log(OR) per SD<br>", EGFR_LABEL, ", STK11 included")
)

# ── Panel h: same three arms, STK11 excluded ───────────────────────────────────
stk11_excluded_arms_egfr <- read_tsv(
  file.path(tbl_dir, "pi3k_bodycomp_firth_luad_egfr_stk11_excluded.tsv"),
  show_col_types = FALSE
) %>%
  prep_panel()

stopifnot(nrow(stk11_excluded_arms_egfr) == 9)

pi3k_forest_l3_stk11_excluded_egfr <- make_egfr_forest(
  stk11_excluded_arms_egfr,
  xlab = paste0("Adj. log(OR) per SD<br>", EGFR_LABEL, ", STK11 excluded")
)

create_dir_if_needed(fig_dir)
ggsave_wrapper(
  plot = pi3k_forest_l3_clean_arms_egfr,
  filename = file.path(fig_dir, "pi3k_forest_luad_l3_clean_scan_arms_egfr.pdf"),
  width = 2.6,
  height = 1.6
)
ggsave_wrapper(
  plot = pi3k_forest_l3_stk11_excluded_egfr,
  filename = file.path(fig_dir, "pi3k_forest_luad_egfr_stk11_excluded_arms_labelled.pdf"),
  width = 3.5,
  height = 1.6
)
