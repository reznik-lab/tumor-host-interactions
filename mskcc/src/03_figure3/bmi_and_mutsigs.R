# ============================================================================
# Produces: Figure 3 (A); Supplementary Table 17
# ============================================================================
##################
# Analyze deepsig mutation signatures
# Tool: https://github.com/mskcc/DeepSig
#
# DeepSig uses a deep learning model to call mutational signatures for a targeted
# sequencing sample like IMPACT.
#
#
# @author: Ethan Tse
# @date: 2026-05-01
##################

source("./src/00_shared/load_libraries.R")
library(ggtext)

long_mutsig_df <- read_tsv("./code_outputs/stable_data/msk/impact/mutation_signatures/deepsig/deepsig_exposures_long_format.tsv") %>%
  # group colorectal
  mutate(CANCER_TYPE_DETAILED = case_when(
    CANCER_TYPE_DETAILED %in% c("Colorectal Adenocarcinoma", "Rectal Adenocarcinoma", "Colon Adenocarcinoma") ~ "Colorectal Adenocarcinoma",
    TRUE ~ CANCER_TYPE_DETAILED
  ))

# do a regression for each signature, each disease
df_res <- long_mutsig_df %>%
  group_by(CANCER_TYPE_DETAILED) %>%
  filter(n_distinct(sid) >= 100) %>% # filter to cancer types with at least 100 samples
  ungroup() %>%
  mutate(weighted_exposure = exposure * M) %>%
  # drop_na() %>%
  group_by(CANCER_TYPE_DETAILED, signature) %>%
  group_modify(~ {
    # express BMI in SD units so estimates are comparable across cancer types
    .x <- .x %>% mutate(BMI_NEAREST_RESULT = as.numeric(scale(BMI_NEAREST_RESULT)))

    # build formula depending on whether GENDER has variation
    base_vars <- c("BMI_NEAREST_RESULT", "CURRENT_AGE_DEID", "STAGE_CDM_DERIVED", "ANCESTRY_LABEL", "SMOKING_PREDICTIONS_3_CLASSES", "PRIOR_MED_TO_MSK")
    if (length(unique(na.omit(.x$GENDER))) == 1) {
      formula <- as.formula(paste0("weighted_exposure ~ ", paste(base_vars, collapse = " + ")))
    } else {
      formula <- as.formula(paste0("weighted_exposure ~ ", paste(c("GENDER", base_vars), collapse = " + ")))
    }

    model <- tryCatch(lm(formula = formula, data = .x, na.action = "na.omit"), error = function(e) NULL)
    if (is.null(model)) {
      return(tibble(
        term = "BMI_NEAREST_RESULT",
        estimate = NA,
        std.error = NA,
        statistic = NA,
        p.value = NA
      ))
    }

    td <- broom::tidy(model, conf.int = TRUE)
    bmi_row <- td %>% dplyr::filter(term == "BMI_NEAREST_RESULT")
  }) %>%
  ungroup() %>%
  mutate(
    padj = p.adjust(p.value, method = "BH")
  )

# save results table
create_dir_if_needed("./code_outputs/results/tables/bmi/mutagenesis/mutation_signatures/deepsig/")
write_tsv(df_res, "./code_outputs/results/tables/bmi/mutagenesis/mutation_signatures/deepsig/mutsig_linear_regression_bmi_results.tsv")

# simply visualize with a volcano plot first
mutsig_volcano <- df_res %>%
  mutate(
    signif = ifelse(padj < 0.05, "padj < 0.05", "NS"),
    cancer_code = cancer_codes[CANCER_TYPE_DETAILED]
  ) %>%
  ggplot(aes(x = estimate, y = -log10(p.value))) +
  geom_point(aes(colour = signif), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1) +
  labs(
    x = "Adj. estimate per SD BMI increase",
    y = expression(-log[10](italic(p) * "-value")),
    color = "Significance"
  ) +
  theme_std() +
  geom_text_repel(
    data = function(df) df %>% filter(padj < 0.05),
    aes(label = paste0(CANCER_TYPE_DETAILED, ": ", signature)),
    box.padding = 0.25,
    max.overlaps = 15,
    show.legend = FALSE,
    segment.color = "black",
    min.segment.length = 0,
    force = 50,
    size = 2.5,
    # Matches the leader-line weight of the gene-level volcano this panel is
    # read against in the manuscript (Figure 4a,
    # src/04_figure4/snvs_bmi_association_volcano.R); the
    # GeomTextRepel default of 0.5 is five times heavier.
    segment.size = 0.1
  ) +
  scale_color_manual(
    name = NULL,
    # orange matches Figure 1C's CASE_COLORS (Cancer = "#D85A30")
    values = c("NS" = "darkgrey", "padj < 0.05" = "#D85A30"),
    breaks = c("padj < 0.05", "NS"),
    labels = c(
      "padj < 0.05" = "padj < 0.05",
      "NS" = "NS"
    )
  ) +
  theme(
    legend.position = "bottom",
    legend.justification = "center"
  ) +
  scale_x_continuous(n.breaks = 10)

ggsave_wrapper("./code_outputs/results/figs/bmi/mutagenesis/mutation_signatures/deepsig/mutsig_linear_regression_bmi_volcano.pdf",
  mutsig_volcano,
  width = 2.5,
  height = 2.5
)

# specifically for UCEC:
#   SBS10 (POLE-associated) -> adjust for the presence of POLE mutations
#   SBS1 (clock-like, accumulates under MMR deficiency) -> adjust for MMR
#     deficiency (MSH2 / MLH1 / MSH6 / PMS2 mutations or MSI-High status)
mutation_matrix <- read_tsv("./code_outputs/stable_data/msk/impact/snvs/matrix/panel_corrected/bmi_oncokb_annotated_first_sample_oncogenic_mutation_mtx.tsv") %>%
  dplyr::select(
    SAMPLE_ID, PATIENT_ID, CANCER_TYPE_DETAILED,
    POLE, MSH2, MLH1, MSH6, PMS2, MSI_TYPE
  ) %>%
  filter(CANCER_TYPE_DETAILED == "Uterine Endometrioid Carcinoma") %>%
  # composite MMR-deficiency indicator: any MMR gene mutation or MSI-High
  mutate(
    MMR = as.integer(
      MSI_TYPE == "Instable"
    ),
    MMR = tidyr::replace_na(MMR, 0L)
  )

ucec_mutsigs <- long_mutsig_df %>%
  inner_join(mutation_matrix, by = c(sid = "SAMPLE_ID", "CANCER_TYPE_DETAILED"))

base_covariates <- "BMI_NEAREST_RESULT + CURRENT_AGE_DEID + STAGE_CDM_DERIVED + ANCESTRY_LABEL + SMOKING_PREDICTIONS_3_CLASSES + PRIOR_MED_TO_MSK"

# ---------------------------------------------------------------
# Combined forest plot faceted by signature: the BMI coefficient
# (effect of BMI on weighted exposure) with vs. without adjusting
# for POLE mutation status.
# ---------------------------------------------------------------
sig_adjustment_models <- list(
  SBS10 = c("w/o *POLE*" = "", "w/ *POLE*" = "POLE"),
  SBS1  = c("w/o *POLE*" = "", "w/ *POLE*" = "POLE")
)

# level order is bottom-to-top (ggplot plots the first level at the
# bottom), so this puts "w/o POLE" above "w/ POLE" on the y axis;
# shared across facets since free_y only shows the levels actually
# present in each panel; POLE is a gene symbol and is italicized via
# markdown for the y axis
model_label_levels <- c("w/ *POLE*", "w/o *POLE*")

# facets should read SBS10 first, then SBS1 (not alphabetical)
signature_levels <- c("SBS10", "SBS1")

bmi_by_model_df <- purrr::imap_dfr(sig_adjustment_models, function(adj_terms, sig) {
  sig_data <- ucec_mutsigs %>%
    filter(signature == sig) %>%
    mutate(weighted_exposure = exposure * M) %>%
    # express BMI in SD units so estimates are comparable across signatures
    mutate(BMI_NEAREST_RESULT = as.numeric(scale(BMI_NEAREST_RESULT)))

  purrr::imap_dfr(adj_terms, function(adj_term, model_label) {
    covariates <- if (nzchar(adj_term)) paste(base_covariates, "+", adj_term) else base_covariates
    model <- lm(as.formula(paste("weighted_exposure ~", covariates)), data = sig_data, na.action = "na.omit")

    broom::tidy(model, conf.int = TRUE) %>%
      dplyr::filter(term == "BMI_NEAREST_RESULT") %>%
      mutate(signature = sig, model_label = model_label)
  })
})

bmi_by_model_forest <- bmi_by_model_df %>%
  mutate(
    model_label = factor(model_label, levels = model_label_levels),
    signature = factor(signature, levels = signature_levels),
    significance = ifelse(p.value < 0.05, "p < 0.05", "NS")
  ) %>%
  ggplot(aes(x = estimate, y = model_label)) +
  ggstats::geom_stripped_rows() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1) +
  geom_linerange(
    aes(xmin = estimate - 1.96 * std.error, xmax = estimate + 1.96 * std.error, color = significance),
    lineend = "round",
    linewidth = 0.1
  ) +
  geom_point(size = 2, shape = 16, aes(colour = significance)) +
  scale_colour_manual(
    name = NULL,
    # red matches Figure 5D's significance colour ("#ad3f3f")
    values = c("p < 0.05" = "#ad3f3f", "NS" = "black")
  ) +
  facet_wrap(~signature, ncol = 1, scales = "free_y") +
  labs(
    x = "Adj. change in weighted exposure per SD BMI",
    y = "Linear model specification"
  ) +
  theme_std() +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    axis.text.y = element_markdown(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 6),
    panel.spacing.y = unit(1, "pt")
  )

ggsave_wrapper(
  "./code_outputs/results/figs/bmi/mutagenesis/mutation_signatures/deepsig/ucec_bmi_adjustment_models_forest_plot.pdf",
  bmi_by_model_forest,
  width = 2.5, height = 2.5
)


# ---------------------------------------------------------------
# Same estimates, as a ggplot "table" panel plotted side-by-side with
# the forest plot. It reuses the exact same y aesthetic (model_label)
# and facet_wrap(~signature) as bmi_by_model_forest, so the rows and
# the geom_stripped_rows() shading line up automatically once the two
# panels are combined with cowplot::plot_grid(align = "hv").
# ---------------------------------------------------------------
bmi_by_model_table_plot <- bmi_by_model_df %>%
  mutate(
    model_label = factor(model_label, levels = model_label_levels),
    signature = factor(signature, levels = signature_levels),
    estimate_label = sprintf("%.2f", estimate),
    ci_label = sprintf("(%.2f, %.2f)", conf.low, conf.high),
    p_label = ifelse(p.value < 0.001, "< 0.001", formatC(p.value, format = "f", digits = 3)),
    p_signif = p.value < 0.05
  ) %>%
  ggplot(aes(y = model_label)) +
  ggstats::geom_stripped_rows() +
  geom_text(aes(x = 1, label = estimate_label), size = 2) +
  geom_text(aes(x = 2.8, label = ci_label), size = 2) +
  # colour matches bmi_by_model_forest's scale_colour_manual (red/black by
  # significance), so the p-value column reads consistently with the forest
  # plot's point/CI colouring instead of introducing a second bold-based cue
  geom_text(aes(x = 4.4, label = p_label, colour = ifelse(p_signif, "p < 0.05", "NS")), size = 2) +
  # red matches Figure 5D's significance colour ("#ad3f3f")
  scale_colour_manual(values = c("p < 0.05" = "#ad3f3f", "NS" = "black"), guide = "none") +
  scale_x_continuous(
    # lower limit trimmed from the original 0.2 -- that much of it was dead
    # space to the left of the BMI Estimate column, between the table panel's
    # left edge and where the text actually starts (widening this panel's
    # share of its column relative to how much of it the text uses), which
    # read as extra distance between the forest plot and the table
    limits = c(0.5, 5),
    breaks = c(1, 2.8, 4.4),
    labels = c("BMI Estimate", "95% CI", "*p*-value"),
    position = "top",
    expand = c(0, 0)
  ) +
  facet_wrap(~signature, ncol = 1, scales = "free_y") +
  theme_std() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    # axis.line (parent) doesn't override the axis.line.x/.y set explicitly
    # in theme_std(), so blank each one directly
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x.bottom = element_blank(),
    axis.text.x.top = ggtext::element_markdown(size = 6),
    # keep the strip the same height as the forest plot's, just invisible
    strip.text = element_text(size = 6, colour = NA),
    panel.spacing.y = unit(1, "pt"),
    plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0)
  )

# saved on its own (.pdf + .rds) so it can be assembled with
# bmi_by_model_forest (fig_3b) via cowplot::plot_grid in figs.qmd
ggsave_wrapper(
  "./code_outputs/results/figs/bmi/mutagenesis/mutation_signatures/deepsig/ucec_bmi_adjustment_models_table_plot.pdf",
  bmi_by_model_table_plot,
  width = 2.5, height = 2.5
)
