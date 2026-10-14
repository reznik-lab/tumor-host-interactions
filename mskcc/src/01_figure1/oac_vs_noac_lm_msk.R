# ============================================================================
# Produces: Figure 1 (H); Supplementary Figure 1 (E); Supplementary Table 9
# ============================================================================
#####################
# Formally test if OACs have higher BMI than non-OACs, adjusting for confounders, using linear regression.
# This is a one-sided test, where the alternative hypothesis is that OACs have higher BMI than non-OACs.
#
# MSK-IMPACT cohort.
#
# @author: Ethan Tse
# @date: 2026-04-22
#####################

source("./src/00_shared/load_libraries.R")
library(ggtext)

msk_metadata <- read_tsv("./code_outputs/stable_data/msk/metadata/impact_patient_and_sample_metadata_first_surgery.tsv") %>%
  janitor::clean_names() %>%
  select(
    tumor_code = oncotree_code,
    # Deliberately the EARLIEST BMI, not bmi_nearest_result: this panel is read
    # against the All of Us trajectories in Figure 1c-f, whose exposure is each
    # individual's first recorded BMI. Every other MSK BMI analysis in the
    # manuscript uses the near-acquisition BMI -- this one stays comparable to
    # AoU instead.
    bmi = bmi_first_result,
    prior_med_to_msk,
    current_age_deid,
    gender,
    ancestry_label,
    stage_highest_recorded,
    # stage_cdm_derived is the model covariate (study-wide convention);
    # stage_cdm_derived_granular is retained only for the stage IV cohort filter
    stage_cdm_derived,
    stage_cdm_derived_granular,
    smoking_predictions_3_classes,
    primary_site,
    gene_panel,
    cancer_type_detailed, cancer_type
  ) %>%
  mutate(cancer_type_detailed = gsub("/.*", "", cancer_type_detailed)) %>%
  mutate(tumor_code = if_else(tumor_code %in% c("COAD", "READ", "COADREAD"), "COADREAD", tumor_code)) %>%
  filter(gene_panel %in% c("IMPACT341", "IMPACT410", "IMPACT468", "IMPACT505"))


adult_tumor_codes <- msk_metadata %>%
  group_by(tumor_code) %>%
  summarise(median_age = median(current_age_deid, na.rm = TRUE)) %>%
  filter(median_age >= 8) %>%
  pull(tumor_code)

msk_metadata <- msk_metadata %>%
  filter(tumor_code %in% adult_tumor_codes)

# only test tumor types with at least 100 patients, but keep all samples as background
tumor_codes_to_test <- msk_metadata %>%
  group_by(tumor_code) %>%
  filter(n() >= 100) %>%
  pull(tumor_code) %>%
  unique()

# only plot tumor types with at least 400 patients, to avoid noisy underpowered estimates in the forest plots
tumor_codes_to_plot <- msk_metadata %>%
  group_by(tumor_code) %>%
  filter(n() >= 400) %>%
  pull(tumor_code) %>%
  unique()

# Map each tumor_code to its most common primary_site for downstream classification
tumor_site_map <- msk_metadata %>%
  filter(!is.na(primary_site)) %>%
  count(tumor_code, primary_site) %>%
  group_by(tumor_code) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(tumor_code, primary_site)

# Map each tumor_code to its most common cancer_type_detailed for plot labels
tumor_label_map <- msk_metadata %>%
  filter(!is.na(cancer_type_detailed)) %>%
  count(tumor_code, cancer_type_detailed) %>%
  group_by(tumor_code) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(tumor_code, cancer_type_detailed)

# Sufficient-evidence OAC sites per NCI/IARC evidence table
is_obesity_associated_site <- function(primary_site) {
  grepl(
    paste0(
      "esophag|gastric|stomach|",
      "colon|rect|cecum|appendix|",
      "liver|gallbladder|pancrea|",
      "breast|uter|endometri|ovar|fallopian|",
      "kidney|renal|",
      "mening|thyroid|myeloma|bone marrow"
    ),
    tolower(trimws(as.character(primary_site)))
  )
}

# Tag each sample's tumor_code as obesity-associated or not, so background pools
tumor_obesity_map <- tumor_site_map %>%
  mutate(obesity_associated = is_obesity_associated_site(primary_site)) %>%
  select(tumor_code, obesity_associated)

msk_metadata <- msk_metadata %>%
  left_join(tumor_obesity_map, by = "tumor_code")

# One adjusted linear model per cancer type: the target type against every other
# cancer type pooled, fit on the whole cohort with only the binary indicator
# redefined each time.
#
# The target-vs-other effect is read straight off the `tumor_vs_rest`
# coefficient. This used to go through emmeans(pairwise ~ tumor_vs_rest), which
# for a two-level factor in a model with no interactions returns exactly the
# coefficient -- the covariate adjustment is identical for both levels and
# cancels in the difference, so estimate, SE and p agreed to machine precision.
# "other" is the reference level so the coefficient is already oriented as
# (target - other), matching the contrast emmeans used to report.
fit_one_vs_all <- function(metadata, tumor, formula_vars) {
  metadata <- metadata %>%
    mutate(
      tumor_vs_rest = factor(
        if_else(tumor_code == tumor, tumor, "other"),
        levels = c("other", tumor)
      )
    )

  formula <- as.formula(paste("bmi ~", "tumor_vs_rest +", formula_vars))

  model <- tryCatch(
    lm(formula = formula, data = metadata),
    error = function(e) NULL
  )

  if (is.null(model)) {
    return(tibble())
  }

  term <- paste0("tumor_vs_rest", tumor)
  coefs <- summary(model)$coefficients
  if (!term %in% rownames(coefs)) {
    return(tibble())
  }

  estimate <- coefs[term, "Estimate"]
  se <- coefs[term, "Std. Error"]

  tibble(
    tumor_code = tumor,
    estimate = estimate,
    # Wald interval on the normal approximation, kept rather than confint()'s
    # t interval so these match the published estimates exactly.
    confint_lower = estimate - 1.96 * se,
    confint_upper = estimate + 1.96 * se,
    p.value = coefs[term, "Pr(>|t|)"]
  )
}

# Test every tumor type with >= 100 patients (tumor_codes_to_test), using all
# samples as background; results are later restricted to tumor_codes_to_plot
# (>= 400 patients) for display, so BH correction is computed over this full
# >= 100 test set before that restriction is applied.
msk_one_vs_all <- purrr::map_dfr(
  tumor_codes_to_test,
  ~ fit_one_vs_all(
    msk_metadata, .x,
    "current_age_deid + gender + ancestry_label + stage_cdm_derived + smoking_predictions_3_classes + prior_med_to_msk"
  )
)

# Cancer types to emphasize via y-axis label colour (not tied to significance/obesity status)
highlight_cancers <- c("Stomach Adenocarcinoma", "Pancreatic Adenocarcinoma")
highlight_color <- "#D37295"

make_cancer_label <- function(cancer_type_detailed, cancer) {
  label <- if_else(
    cancer_type_detailed %in% highlight_cancers,
    paste0("<span style='color:", highlight_color, "'>", cancer_type_detailed, "</span>"),
    as.character(cancer_type_detailed)
  )
  factor(label, levels = label[order(as.integer(cancer))])
}

# Row background category: which side of the null (and how significant) each
# cancer's BMI estimate falls on, used to shade each forest-plot row.
direction_significance_levels <- c("Significantly larger", "Not significantly different", "Significantly smaller")
classify_direction_significance <- function(padj, estimate) {
  factor(
    case_when(
      padj < 0.05 & estimate > 0 ~ "Significantly larger",
      padj < 0.05 & estimate < 0 ~ "Significantly smaller",
      TRUE ~ "Not significantly different"
    ),
    levels = direction_significance_levels
  )
}

## Create a forest plot of BMI estimate differences (tumor vs other)
# BH correction is computed over every tested tumor type (tumor_codes_to_test,
# n = 109) BEFORE restricting to the tumor_codes_to_plot subset (n = 400+
# patients) -- correcting after that restriction would understate the true
# number of tests performed and inflate significance.
forest_df <- msk_one_vs_all %>%
  ungroup() %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  filter(tumor_code %in% tumor_codes_to_plot) %>%
  left_join(tumor_label_map, by = "tumor_code") %>%
  left_join(tumor_site_map, by = "tumor_code") %>%
  filter(!is.na(cancer_type_detailed)) %>%
  mutate(
    obesity_associated = is_obesity_associated_site(primary_site),
    obesity_group = factor(
      if_else(obesity_associated, "Obesity associated", "Not obesity associated"),
      levels = c("Obesity associated", "Not obesity associated")
    ),
    cancer = factor(cancer_type_detailed, levels = cancer_type_detailed[order(estimate)]),
    cancer_label = make_cancer_label(cancer_type_detailed, cancer),
    Significance = case_when(
      padj < 0.05 ~ "padj < 0.05",
      TRUE ~ "NS"
    ),
    direction_significance = classify_direction_significance(padj, estimate)
  )

# Facet vertically by obesity-association status instead of colouring by it,
# so the two groups are visually separated while the x-axis scale stays shared.

# Shared forest-plot styling, matching pi3k_bodycomp_forest_plots.R's forest_theme
forest_theme <- list(
  # Row background shaded by which side of the null the estimate falls on
  # (and whether that's significant), in place of plain zebra striping, so
  # "significantly larger/smaller vs. not significant" reads at a glance.
  geom_tile(
    aes(x = 0, y = cancer_label, fill = direction_significance),
    width = Inf, height = 1
  ),
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1),
  geom_linerange(
    aes(xmin = confint_lower, xmax = confint_upper, color = Significance),
    lineend = "round",
    linewidth = 0.1
  ),
  geom_point(aes(color = Significance), size = 1.5, shape = 16),
  scale_colour_manual(
    name = NULL,
    # red matches Figure 5D's significance colour ("#ad3f3f")
    values = c("padj < 0.05" = "#ad3f3f", "NS" = "black"),
    # legend key symbols bumped well past the in-plot point size (1.5) --
    # at that size they read as barely-visible dots once shrunk to print
    # scale -- and stacked vertically instead of ggplot's default
    # side-by-side row.
    guide = guide_legend(override.aes = list(size = 3), direction = "vertical")
  ),
  scale_fill_manual(
    name = NULL,
    values = c(
      "Significantly larger" = "#FADBD8", # light red
      "Not significantly different" = "#E5E5E5", # light grey
      "Significantly smaller" = "#D6EAF8" # light blue
    ),
    guide = "none"
  ),
  scale_x_continuous(n.breaks = 8),
  facet_wrap(~obesity_group, ncol = 1, scales = "free_y", space = "free_y"),
  theme_std(),
  theme(
    legend.position = "none",
    axis.text.y = element_markdown(),
    axis.ticks.y = element_blank()
  )
)

forest_plot <- forest_df %>%
  ggplot(aes(x = estimate, y = cancer_label)) +
  forest_theme +
  labs(x = stringr::str_wrap("Adj. BMI difference vs. all other cancers", width = 30), y = NULL, color = NULL)

ggsave_wrapper("./code_outputs/results/figs/bmi/descriptive_statistics/bmi_one_vs_all_cancer_comparision.pdf", forest_plot,
  width = 3.9, height = 3.5
)

#####################
# Repeat the analysis restricted to stage 4 cancers only.
# Stage is dropped as a covariate here since the cohort is restricted to stage 4.
#####################

stage4_metadata <- msk_metadata %>%
  filter(grepl("4.0", stage_cdm_derived_granular, ignore.case = TRUE))

# only test tumor types with at least 100 patients in the stage 4 subset
stage4_tumor_codes_to_test <- stage4_metadata %>%
  group_by(tumor_code) %>%
  filter(n() >= 100) %>%
  pull(tumor_code) %>%
  unique()

# Unlike the main panel, the stage-4 sensitivity plots everything it tests. The
# >= 400 display threshold is calibrated to the full cohort; inside the stage-4
# subset (n = 23,649) it leaves only 8 of the 40 tested types standing, which is
# too thin a slice to show that the main result survives the restriction.
stage4_tumor_codes_to_plot <- stage4_tumor_codes_to_test

stage4_one_vs_all <- purrr::map_dfr(
  stage4_tumor_codes_to_test,
  ~ fit_one_vs_all(
    stage4_metadata, .x,
    "current_age_deid + gender + ancestry_label + smoking_predictions_3_classes + prior_med_to_msk"
  )
)

stage4_forest_df <- stage4_one_vs_all %>%
  ungroup() %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  filter(tumor_code %in% stage4_tumor_codes_to_plot) %>%
  left_join(tumor_label_map, by = "tumor_code") %>%
  left_join(tumor_site_map, by = "tumor_code") %>%
  filter(!is.na(cancer_type_detailed)) %>%
  mutate(
    obesity_associated = is_obesity_associated_site(primary_site),
    obesity_group = factor(
      if_else(obesity_associated, "Obesity associated", "Not obesity associated"),
      levels = c("Obesity associated", "Not obesity associated")
    ),
    cancer = factor(cancer_type_detailed, levels = cancer_type_detailed[order(estimate)]),
    cancer_label = make_cancer_label(cancer_type_detailed, cancer),
    Significance = case_when(
      padj < 0.05 ~ "padj < 0.05",
      TRUE ~ "NS"
    ),
    direction_significance = classify_direction_significance(padj, estimate)
  )

stage4_forest_plot <- stage4_forest_df %>%
  ggplot(aes(x = estimate, y = cancer_label)) +
  forest_theme +
  labs(
    x = "Adj. BMI difference vs. all other cancers",
    y = NULL, color = NULL,
    title = "Stage 4 cancers only"
  )

# Taller than the main panel (3.5in for 34 rows) because this one now shows all
# 40 tested types; same ~0.1in per row so the two read at the same density.
ggsave_wrapper("./code_outputs/results/figs/bmi/descriptive_statistics/bmi_one_vs_all_cancer_comparision_stage4.pdf", stage4_forest_plot,
  width = 3.9, height = 4.2
)
