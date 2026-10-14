# ============================================================================
# Produces: Supplementary Figure 2 (F); Supplementary Table 15
# ============================================================================
####################
# Forest plot: adjusted myopenia -> myosteatosis association within each
# cancer type, with its treatment-naive sensitivity re-fit.
#
# Reads the per-cancer-type myopenia coefficients from
# myopenia_myosteatosis_by_cancer_type_mixed_model.R (a separate adjusted
# logistic regression within each cancer type, adjusted for age at scan, sex,
# stage, ancestry and smoking) for both the full cohort
# and the treatment-naive ("clean scan") subset, and renders them as one
# forest with the two cohorts dodged per cancer type.
#
# This is the panel for the primary analysis behind "we identified three
# statistically significant negative associations ... colorectal, lung, and
# pancreatic adenocarcinomas": only the pancreatic cell of that model was
# plotted anywhere (Figure 2's contingency table), and the clean-scan re-fit
# was computed but never shown. Both arms now have a home.
#
# Plotted on the log(OR) scale rather than the OR scale deliberately: the
# model's native estimate is a difference of logits, Supplementary Figure 2H
# shows the same contrast as an (unadjusted, log2) odds ratio, and the
# distinction between the two panels is easier to hold onto if they do not
# share an axis scale.
#
# Cohort colors match myopenia_one_vs_rest_forest_plots.R so the two
# forests in Supplementary Figure 2 read as the same encoding.
#
# @author: Ethan Tse
# @date: October 5th, 2026
####################

source("./src/00_shared/load_libraries.R")

results_dir <- "./code_outputs/results/tables/body_comp/abnormal_physiology"
out_dir <- "./code_outputs/results/figs/body_comp/abnormal_physiology"
create_dir_if_needed(paste0(out_dir, "/"))

cohort_colors <- c("All Scans" = "#4C72B0", "Clean Scan Subset" = "#DD8452")

read_arm <- function(file, cohort) {
  read_tsv(file.path(results_dir, file), show_col_types = FALSE) %>%
    mutate(cohort = cohort)
}

combined_results <- bind_rows(
  read_arm("myopenia_myosteatosis_by_cancer_type_mixed_model.tsv", "All Scans"),
  read_arm("myopenia_myosteatosis_by_cancer_type_mixed_model_clean_scan_subset.tsv", "Clean Scan Subset")
) %>%
  mutate(
    cancer_code = unname(cancer_codes[cancer_type_detailed]),
    cohort = factor(cohort, levels = names(cohort_colors)),
    # One star = nominally significant, two = survives BH correction, matching
    # Supplementary Figure 2e.
    sig_label = case_when(
      p_adj_bh < 0.05 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Order by the full-cohort estimate so both arms of a cancer type share a
# row. The clean-scan model re-derives its own eligibility threshold (>= 10
# patients in each myopenia level), so it covers 8 of the 12 cancer types --
# the other 4 rows carry an "All Scans" point only, which is itself the
# informative thing to show.
order_key <- combined_results %>%
  filter(cohort == "All Scans") %>%
  select(cancer_code, order_estimate = estimate)

combined_results <- combined_results %>%
  left_join(order_key, by = "cancer_code") %>%
  mutate(
    order_estimate = coalesce(order_estimate, estimate),
    cancer_code = forcats::fct_reorder(cancer_code, order_estimate)
  )

glmm_forest <- combined_results %>%
  ggplot(aes(x = estimate, y = cancer_code)) +
  ggstats::geom_stripped_rows() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1) +
  geom_linerange(
    aes(xmin = conf.low, xmax = conf.high, color = cohort),
    position = position_dodge(width = 0.6),
    lineend = "round", linewidth = 0.4
  ) +
  geom_point(
    # Colour alone separates the two arms (same encoding as the one-vs-rest
    # forest beside it in Supplementary Figure 2).
    aes(color = cohort),
    position = position_dodge(width = 0.6),
    size = 1.6
  ) +
  geom_text(
    aes(x = conf.high, label = sig_label, group = cohort),
    position = position_dodge(width = 0.6),
    hjust = -0.3, size = 2, show.legend = FALSE
  ) +
  scale_colour_manual(name = NULL, values = cohort_colors, drop = FALSE) +
  # Stars sit at the upper CI bound, so the default expansion clips the ones
  # on rows whose interval reaches furthest right.
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.12))) +
  # Just the scale: the plotted estimate is the myopenia coefficient, i.e. an
  # odds ratio between myopenic and non-myopenic patients, which "log(OR) of
  # myosteatosis given myopenia" named as a conditional odds instead. The
  # contrast belongs in the caption.
  labs(x = "Adjusted log(OR)", y = NULL) +
  theme_std() +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.box.just = "center",
    axis.ticks.y = element_blank()
  )

ggsave_wrapper(
  plot = glmm_forest,
  filename = file.path(out_dir, "myopenia_myosteatosis_glmm_forest.pdf"),
  width = 3.2,
  height = 3.0
)
