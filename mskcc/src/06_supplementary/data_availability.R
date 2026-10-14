# ============================================================================
# Produces: Supplementary Figure 2 (A)
# ============================================================================
###########
# Data availability for CT scans
#
#
# @author: Ethan Tse
# @date: 2026-04-15
###########

source("./src/00_shared/load_libraries.R")

master_bodycomp_meta <- read_tsv("./code_outputs/stable_data/msk/body_comp/metadata/non_clean_master_body_comp_annotated_metadata.tsv")

# The cancer-type grouping is fixed on the full cohort: any type with >= 50
# patients keeps its own label, everything else rolls up to "Other". Recomputing
# the <50 rollup inside a smaller subset would move types that clear 50 overall
# but not within the subset (e.g. Breast Invasive Ductal Carcinoma) silently into
# "Other", so the grouping is computed once, here.
cancer_group_levels <- master_bodycomp_meta %>%
  mutate(cancer_type_detailed = gsub("/.*", "", cancer_type_detailed)) %>%
  count(cancer_type_detailed) %>%
  filter(n >= 50) %>%
  pull(cancer_type_detailed)

# barplot showing the number of CT scans for each cancer type
#
# Merkel cell carcinoma is deliberately NOT dropped here. It used to be, and
# because the filter ran after the "Other" rollup its 65 patients vanished from
# the panel entirely rather than folding into "Other", leaving 14 named types
# against the text's 15. MCC is held out of the Figure 2a one-vs-rest model
# (bodycomp_vs_background_one_vs_rest.R), not out of the cohort.
sample_size_barplot <- function(df) {
  df %>%
    mutate(
      cancer_type_detailed = gsub("/.*", "", cancer_type_detailed),
      cancer_type_detailed = ifelse(
        cancer_type_detailed %in% cancer_group_levels, cancer_type_detailed, "Other"
      )
    ) %>%
    count(cancer_type_detailed, name = "num_scans") %>%
    mutate(total_scans = num_scans) %>%
    arrange(desc(total_scans)) %>%
    ggplot(aes(
      x = reorder(cancer_type_detailed, total_scans),
      y = num_scans,
      fill = cancer_type_detailed
    )) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.7) +
    geom_text(aes(y = total_scans, label = total_scans), hjust = -0.2, size = 2) +
    scale_fill_manual(values = paletteer_d("ggthemes::Tableau_20")) +
    labs(
      x = NULL,
      y = "Number of Patients"
    ) +
    theme_std() +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.17)), breaks = function(x) pretty(x, n = 7)) +
    theme(legend.position = "none", axis.ticks.y = element_blank())
}

master_bodycomp_meta_sample_size <- sample_size_barplot(master_bodycomp_meta)

ggsave_wrapper(
  plot = master_bodycomp_meta_sample_size,
  filename = "./code_outputs/results/figs/body_comp/descriptive_overview/sample_size/master_bodycomp_meta_sample_size.pdf",
  width = 3.3,
  height = 2
)

