# ============================================================================
# Produces: Figure 5 (B); Supplementary Figure 5 (B)
# ============================================================================
################
# PI3K/AKT/mTOR mutation composition barplot — PDAC, KIRC, COADREAD
#
# For each cancer type, shows what fraction of PI3K/AKT/mTOR-altered samples
# carry each gene mutation. Bars are proportional to patient count (not raw
# mutation count) to avoid double-counting patients with multiple hits.
#
# @author: Ethan Tse
# @date: June 2026
################

source("./src/00_shared/load_libraries.R")
source("./src/00_shared/obesity_constants.R")
library(tidytext)

pathway_data <- read_tsv("./code_outputs/stable_data/msk/impact/pathways/sample_level_mutation_w_pathway_annotations.tsv")
mutation_data <- read_tsv("./code_outputs/stable_data/msk/impact/snvs/dataframes/bmi_oncokb_annotated_first_sample_oncogenic_mutations.tsv") %>%
  mutate(
    CANCER_TYPE_DETAILED = case_when(
      CANCER_TYPE_DETAILED %in% c("Colorectal Adenocarcinoma", "Rectal Adenocarcinoma", "Colon Adenocarcinoma") ~ "Colorectal Adenocarcinoma",
      TRUE ~ CANCER_TYPE_DETAILED
    ),
    MUTATION_EFFECT = case_when(
      MUTATION_EFFECT %in% c("Gain-of-function", "Likely Gain-of-function") ~ "Gain-of-function",
      MUTATION_EFFECT %in% c("Loss-of-function", "Likely Loss-of-function") ~ "Loss-of-function",
      MUTATION_EFFECT %in% c("Switch-of-function", "Likely Switch-of-function") ~ "Switch-of-function",
      TRUE ~ "Other/Unknown"
    )
  )

metadata <- read_tsv("./code_outputs/stable_data/msk/metadata/impact_patient_and_sample_metadata_first_surgery.tsv")

target_cancers <- c("Pancreatic Adenocarcinoma", "Renal Clear Cell Carcinoma", "Colorectal Adenocarcinoma")

# ── BMI cohort ─────────────────────────────────────────────────────────────────

bmi_cohort <- metadata %>%
  left_join(
    pathway_data %>% select(sample_id, pi3k_akt_mtor_pathway),
    by = "sample_id"
  ) %>%
  filter(cancer_type_detailed %in% target_cancers) %>%
  filter(!is.na(pi3k_akt_mtor_pathway))

# ── Compute per-gene mutation proportions within each cancer type ──────────────
# Denominator: unique PI3K-altered patients per cancer type.
# Numerator:   unique patients among those carrying each gene x mutation effect.

pi3k_altered_patients <- bmi_cohort %>%
  filter(pi3k_akt_mtor_pathway == 1) %>%
  distinct(patient_id, cancer_type_detailed)

n_altered_by_cancer <- pi3k_altered_patients %>%
  count(cancer_type_detailed, name = "n_altered")

mutation_composition <- mutation_data %>%
  filter(
    CANCER_TYPE_DETAILED %in% target_cancers,
    tolower(Hugo_Symbol) %in% pi3k_akt_mtor_pathway
  ) %>%
  semi_join(pi3k_altered_patients, by = c("PATIENT_ID" = "patient_id")) %>%
  group_by(CANCER_TYPE_DETAILED, Hugo_Symbol, MUTATION_EFFECT) %>%
  summarize(n_patients = n_distinct(PATIENT_ID), .groups = "drop") %>%
  left_join(n_altered_by_cancer, by = c("CANCER_TYPE_DETAILED" = "cancer_type_detailed")) %>%
  mutate(
    pct_of_altered = n_patients / n_altered,
    cancer_code    = cancer_codes[CANCER_TYPE_DETAILED]
  )

# Order genes by total proportion within each cancer type; build facet labels
mutation_composition <- mutation_composition %>%
  group_by(CANCER_TYPE_DETAILED, Hugo_Symbol) %>%
  mutate(gene_total_pct = sum(pct_of_altered)) %>%
  ungroup()

facet_levels <- n_altered_by_cancer %>%
  filter(cancer_type_detailed %in% target_cancers) %>%
  arrange(match(cancer_type_detailed, target_cancers)) %>%
  mutate(facet_label = paste0(cancer_codes[cancer_type_detailed], "\n(n=", n_altered, ")")) %>%
  pull(facet_label)

mutation_composition <- mutation_composition %>%
  mutate(
    facet_label = factor(
      paste0(cancer_code, "\n(n=", n_altered, ")"),
      levels = facet_levels
    ),
    Hugo_Symbol = reorder_within(Hugo_Symbol, gene_total_pct, facet_label)
  )

# ── Plot ───────────────────────────────────────────────────────────────────────

pi3k_mutation_composition_plot <- mutation_composition %>%
  ggplot(aes(x = Hugo_Symbol, y = pct_of_altered, fill = MUTATION_EFFECT)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.85, width = 0.7) +
  scale_x_reordered() +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.08))
  ) +
  coord_flip() +
  facet_wrap(~facet_label, scales = "free", ncol = 3) +
  labs(
    x    = NULL,
    y    = "% of PI3K/AKT/mTOR-altered samples",
    fill = "Mutation effect"
  ) +
  scale_fill_manual(values = c(
    "Gain-of-function"   = "#D55E00",
    "Loss-of-function"   = "#0072B2",
    "Switch-of-function" = "#009E73",
    "Other/Unknown"      = "grey70"
  )) +
  theme_std() +
  theme(
    legend.position = "top",
    strip.text      = element_text(face = "bold"),
    axis.text.y     = element_text(face = "italic"),
    panel.spacing   = unit(0.8, "lines")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

create_dir_if_needed("./code_outputs/results/figs/bmi_pathways/mutation_composition")

ggsave_wrapper(
  plot     = pi3k_mutation_composition_plot,
  filename = "./code_outputs/results/figs/bmi_pathways/mutation_composition/pi3k_mutation_composition_pdac_kirc_coadread.pdf",
  width    = 8,
  height   = 4
)

# ── Stacked barplot: top-5 genes by alteration share per cancer ────────────────
# Same design as pi3k_bodycomp_analysis.R: one 100%-stacked bar per cancer,
# top 4 PI3K/AKT/mTOR genes coloured and the rest collapsed to "Other". Counts
# raw mutations among PI3K-altered patients.
#
# Top 4, not top 5: the 5th-ranked gene is AKT1 in both CRC (2.9%) and PDAC
# (2.2%), a band too thin to hold its own label at this panel's size -- the
# label overflowed into the segment above it. Folding it into "Other" keeps
# every drawn segment labellable and leaves the bars summing to 100%.
pi3k_gene_stacked_data <- mutation_data %>%
  semi_join(pi3k_altered_patients, by = c("PATIENT_ID" = "patient_id")) %>%
  filter(
    CANCER_TYPE_DETAILED %in% target_cancers,
    tolower(Hugo_Symbol) %in% pi3k_akt_mtor_pathway
  ) %>%
  group_by(CANCER_TYPE_DETAILED, Hugo_Symbol) %>%
  summarize(n_mutations = n(), .groups = "drop") %>%
  group_by(CANCER_TYPE_DETAILED) %>%
  mutate(
    rank = rank(-n_mutations, ties.method = "first"),
    gene_group = if_else(rank <= 4, Hugo_Symbol, "Other")
  ) %>%
  group_by(CANCER_TYPE_DETAILED, gene_group) %>%
  summarize(n_mutations = sum(n_mutations), .groups = "drop") %>%
  group_by(CANCER_TYPE_DETAILED) %>%
  mutate(
    proportion_of_total = n_mutations / sum(n_mutations),
    cancer_code         = cancer_codes[CANCER_TYPE_DETAILED]
  ) %>%
  ungroup()

# global gene ordering by total mutations across cancers (consistent stack order /
# colours); "Other" forced to the last level = top of each stack
gene_order <- pi3k_gene_stacked_data %>%
  filter(gene_group != "Other") %>%
  group_by(gene_group) %>%
  summarize(total = sum(n_mutations), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(gene_group)

pi3k_gene_stacked_data <- pi3k_gene_stacked_data %>%
  mutate(gene_group = factor(gene_group, levels = c(gene_order, "Other")))

# fill map: each gene a fixed Tableau 20 colour, "Other" grey. Tableau 20 is a
# paired ramp (each hue followed by its light twin), so neighbouring genes in
# the rank order land on two shades of one hue -- the stacked segments are
# separated by the white hairline below rather than by hue alone.
tableau20 <- as.character(paletteer::paletteer_d("ggthemes::Tableau_20"))
stopifnot(length(gene_order) <= length(tableau20))
pi3k_gene_fill_values <- c(
  setNames(tableau20[seq_along(gene_order)], gene_order),
  Other = "grey70"
)

# Tableau 20 mixes dark hues with their light twins, so a single hardcoded text
# colour is unreadable on one half of the palette. Ink is taken from each
# fill's luminance instead.
pi3k_label_ink <- function(hex) {
  rgb_lin <- (grDevices::col2rgb(hex) / 255)^2.2
  lum <- 0.2126 * rgb_lin[1, ] + 0.7152 * rgb_lin[2, ] + 0.0722 * rgb_lin[3, ]
  ifelse(lum < 0.25, "white", "grey10")
}

# Label text and ink are materialised as columns rather than computed inside
# aes(): figs.qmd readRDS()s this object and re-renders it in its own session,
# where pi3k_label_ink()/pi3k_gene_fill_values do not exist, and an aes() that
# called them errored at print time.
pi3k_gene_stacked_data <- pi3k_gene_stacked_data %>%
  mutate(
    segment_label = as.character(gene_group),
    label_ink = pi3k_label_ink(pi3k_gene_fill_values[as.character(gene_group)])
  )

pi3k_x_labels <- setNames(
  stringr::str_wrap(target_cancers, width = 12),
  target_cancers
)

pi3k_gene_stacked_barplot <- pi3k_gene_stacked_data %>%
  ggplot(aes(x = CANCER_TYPE_DETAILED, y = proportion_of_total, fill = gene_group)) +
  # colour/linewidth draw a hairline surface gap between stacked segments;
  # alpha is gone, which was washing every fill toward white and costing the
  # separation the palette is chosen to provide.
  geom_col(
    position = position_fill(reverse = TRUE), width = 0.9,
    colour = "white", linewidth = 0.3
  ) +
  # group = gene_group is REQUIRED. Without it this layer's grouping is
  # inferred from its own discrete aesthetics (colour/fontface), which stacks
  # in a different order than the geom_col above and silently offsets every
  # label by one segment.
  #
  geom_text(
    aes(
      label = segment_label,
      colour = label_ink,
      group = gene_group,
      fontface = ifelse(gene_group == "Other", "plain", "italic")
    ),
    position = position_fill(vjust = 0.5, reverse = TRUE),
    size = 1.75, lineheight = 0.95
  ) +
  scale_colour_identity() +
  scale_fill_manual(values = pi3k_gene_fill_values) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0), labels = pi3k_x_labels) +
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = "% of PI3K/AKT/MTOR\nAlterations",
    fill = "PI3K/AKT/MTOR Gene"
  ) +
  theme_std() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 2, r = 2, b = 14, l = 14)
  )

ggsave_wrapper(
  plot     = pi3k_gene_stacked_barplot,
  filename = "./code_outputs/results/figs/bmi_pathways/mutation_composition/pi3k_gene_stacked_top5_pdac_kirc_coadread.pdf",
  width    = 1.6,
  height   = 2.5
)
