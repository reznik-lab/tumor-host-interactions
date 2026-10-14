# ============================================================================
# Produces: Figure 2 (A); Supplementary Figure 2 (C)
# ============================================================================
##################
# Visualize one-vs-rest results: is the average patient of a given cancer type
# physiologically different from the average other cancer patient, despite
# covariate adjustment?
#
# Same layout as emmeans_heatmap.R, but reads the one-vs-rest contrasts
# (bodycomp_vs_background_one_vs_rest.R) instead of the emmeans "eff" grand-mean
# contrasts.
#
# Plots two heatmaps: the full cohort, and the clean_scan sensitivity cohort
# (from sensitivity_cohorts in bodycomp_vs_background_one_vs_rest.R).
#
# @author: Ethan Tse
# @date: June 25th, 2026
##################

source("./src/00_shared/load_libraries.R")
library(ggnewscale)

tables_dir <- "./code_outputs/results/tables/body_comp/bodycomp_vs_background/lm_w_emmeans/"
figs_dir <- "./code_outputs/results/figs/body_comp/bodycomp_vs_background/"

colorbar_guide <- function(order = 2) {
  guide_colorbar(
    title.position = "top", title.hjust = 0.5, nbin = 50,
    barwidth = unit(1.2, "cm"), barheight = unit(0.15, "cm"),
    order = order
  )
}

# Build and save a one-vs-rest heatmap from a contrasts table.
build_one_vs_rest_heatmap <- function(lm_emmeans_fit, output_filename) {
  lm_emmeans_fit <- lm_emmeans_fit %>%
    mutate(Significance = case_when(
      padj < 0.05 ~ "padj < 0.05",
      p.value < 0.05 ~ "Nominal p < 0.05",
      TRUE ~ "NS"
    )) %>%
    mutate(cancer_type = gsub(" effect", "", cancer_type)) %>%
    mutate(cancer_type = gsub("[()]", "", cancer_type)) %>%
    mutate(cancer = cancer_codes[as.character(cancer_type)]) %>%
    # Strip the "/..." suffix only after the cancer_codes lookup, which is keyed
    # on the full names.
    mutate(cancer_type = trimws(gsub("/.*", "", cancer_type))) %>%
    mutate(body_compartment = gsub(".*\\.", "", body_compartment)) %>%
    mutate(
      metric_group = factor(
        case_when(
          body_compartment %in% c("l3vat_area", "l3sat_area", "l3total_muscle_area") ~ "Area",
          body_compartment %in% c("liver_volume", "kidney_volume", "pancreas_volume", "spleen_volume") ~ "Volume",
          body_compartment %in% c(
            "l3vat_median", "l3sat_median", "l3total_muscle_median_hu",
            "liver_median_hu", "kidney_median_hu", "pancreas_median_hu", "spleen_median_hu"
          ) ~ "Density"
        ),
        levels = c("Area", "Volume", "Density")
      ),
      # Facet grouping: L3 anatomical tissues (area + density) vs whole organs
      # (volume + density).
      tissue_group = factor(
        case_when(
          body_compartment %in% c(
            "l3vat_area", "l3sat_area", "l3total_muscle_area",
            "l3vat_median", "l3sat_median", "l3total_muscle_median_hu"
          ) ~ "Tissues",
          body_compartment %in% c(
            "liver_volume", "kidney_volume", "pancreas_volume", "spleen_volume",
            "liver_median_hu", "kidney_median_hu", "pancreas_median_hu", "spleen_median_hu"
          ) ~ "Organs"
        ),
        levels = c("Tissues", "Organs")
      )
    )

  # Prepare matrix for clustering: rows = body compartments, columns = cancer types.
  # Cluster (and therefore order) cancers using all metrics.
  estimate_matrix <- lm_emmeans_fit %>%
    select(body_compartment, cancer_type, estimate) %>%
    tidyr::pivot_wider(names_from = cancer_type, values_from = estimate, values_fn = mean) %>%
    column_to_rownames("body_compartment") %>%
    as.matrix()

  # Perform hierarchical clustering only on columns (cancer types)
  cancer_dist <- stats::dist(t(estimate_matrix), method = "euclidean")
  cancer_clust <- stats::hclust(cancer_dist, method = "complete")

  # Orient the dendrogram so that the Area effect runs highest -> lowest where the
  # clustering structure permits (rotate branches using mean Area estimate).
  area_score <- lm_emmeans_fit %>%
    filter(metric_group == "Area") %>%
    group_by(cancer_type) %>%
    summarise(area_mean = mean(estimate, na.rm = TRUE), .groups = "drop")

  cancer_wts <- -area_score$area_mean[match(cancer_clust$labels, area_score$cancer_type)]
  cancer_dend <- stats::reorder(stats::as.dendrogram(cancer_clust), wts = cancer_wts, agglo.FUN = mean)
  cancer_order <- labels(cancer_dend)

  # Use clustered order in heatmap
  heatmap_data <- lm_emmeans_fit %>%
    mutate(
      cancer_type = factor(cancer_type, levels = cancer_order),
      body_compartment = factor(body_compartment,
        levels = c(
          # Tissues facet: VAT/SAT/SKM area then density
          "l3vat_area", "l3sat_area", "l3total_muscle_area",
          "l3vat_median", "l3sat_median", "l3total_muscle_median_hu",
          # Organs facet: organ volume then density
          "liver_volume", "kidney_volume", "pancreas_volume", "spleen_volume",
          "liver_median_hu", "kidney_median_hu", "pancreas_median_hu", "spleen_median_hu"
        ),
        labels = c(
          "VAT", "SAT", "SKM",
          "VAT ", "SAT ", "SKM ",
          "Liver", "Kidney", "Pancreas", "Spleen",
          "Liver ", "Kidney ", "Pancreas ", "Spleen "
        )
      )
    )

  # Area estimates are asymmetric around 0; scale each side of the midpoint
  # independently so the negative end saturates instead of washing out.
  area_range <- range(heatmap_data$estimate[heatmap_data$metric_group == "Area"], na.rm = TRUE)

  heatmap_plot <- heatmap_data %>%
    ggplot(aes(x = body_compartment, y = cancer_type)) +
    geom_tile(
      data = ~ filter(.x, metric_group == "Area"),
      aes(fill = estimate), color = "black", width = 0.98, height = 0.98
    ) +
    scale_fill_gradientn(
      colours = c("#5c91a7", "white", "#a75c5c"),
      values = scales::rescale(c(area_range[1], 0, area_range[2])),
      limits = area_range, name = "Diff (cm²)",
      breaks = scales::breaks_pretty(n = 3),
      guide = colorbar_guide(order = 2)
    ) +
    new_scale_fill() +
    geom_tile(
      data = ~ filter(.x, metric_group == "Volume"),
      aes(fill = estimate), color = "black", width = 0.98, height = 0.98
    ) +
    scale_fill_gradient2(
      low = "#5c91a7", mid = "white", high = "#a75c5c",
      midpoint = 0, name = "Diff (cm³)",
      breaks = scales::breaks_pretty(n = 3),
      guide = colorbar_guide(order = 3)
    ) +
    new_scale_fill() +
    geom_tile(
      data = ~ filter(.x, metric_group == "Density"),
      aes(fill = estimate), color = "black", width = 0.98, height = 0.98
    ) +
    scale_fill_gradient2(
      low = "#5c91a7", mid = "white", high = "#a75c5c",
      midpoint = 0, name = "Diff (HU)",
      breaks = scales::breaks_pretty(n = 3),
      guide = colorbar_guide(order = 4)
    ) +
    geom_point(
      data = ~ dplyr::filter(.x, Significance != "NS"),
      aes(shape = Significance),
      color = "black", size = 0.8
    ) +
    scale_shape_manual(
      values = c("padj < 0.05" = 8, "Nominal p < 0.05" = 1),
      name = NULL,
      na.translate = FALSE
    ) +
    facet_grid(. ~ tissue_group + metric_group,
      scales = "free_x", space = "free_x",
      labeller = labeller(metric_group = c(
        Area = "Area (cm²)", Volume = "Volume (cm³)", Density = "Density (HU)"
      ))
    ) +
    labs(x = NULL, y = NULL) +
    theme_std() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "top",
      legend.justification = "center",
      legend.box = "horizontal",
      legend.title = element_text(size = 5),
      legend.text = element_text(size = 5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.1),
      panel.grid = element_blank(),
      plot.margin = margin(2, 2, 2, 2, "mm")
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    guides(shape = guide_legend(title = "", title.position = "top", title.hjust = 0.5, order = 1))

  ggsave_wrapper(
    plot = heatmap_plot,
    filename = paste0(figs_dir, output_filename),
    width = 5.1,
    height = 4,
    extension = "pdf"
  )

  heatmap_plot
}

read_contrasts <- function(filename) {
  read_tsv(paste0(tables_dir, filename)) %>%
    filter(grepl("area|volume|median|_hu", body_compartment, ignore.case = TRUE))
}

# Full cohort
full_cohort_lm_emmeans <- read_contrasts("bodycomp_vs_background_population_one_vs_rest_contrasts.tsv")

full_cohort_lm_emmeans %>%
  filter(padj < 0.05) %>%
  summarize(n = n_distinct(cancer_type))

build_one_vs_rest_heatmap(
  full_cohort_lm_emmeans,
  "bodycomp_vs_background_population_one_vs_rest_heatmap_full_cohort.pdf"
)

# Clean-scan sensitivity cohort
clean_scan_lm_emmeans <- read_contrasts("bodycomp_vs_background_population_one_vs_rest_contrasts_clean_scan.tsv")

clean_scan_lm_emmeans %>%
  filter(padj < 0.05) %>%
  summarize(n = n_distinct(cancer_type))

build_one_vs_rest_heatmap(
  clean_scan_lm_emmeans,
  "bodycomp_vs_background_population_one_vs_rest_heatmap_clean_scan_subset.pdf"
)
