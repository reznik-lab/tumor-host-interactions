# ============================================================================
# Produces: Supplementary Figure 3 (B, C)
# ============================================================================
##################
# Validate the DeepSig mutational signature calls against known etiologies
#
# Key plots showing the method captures biologically-correct signatures:
#   1) SBS4  (tobacco)  in Lung (NSCLC)        -> enriched in smokers
#   2) SBS10 (POLE)     in Endometrial &
#      Colorectal                             -> enriched in POLE-mutant tumors
#   3) MMR signatures (SBS6, SBS14, SBS15, SBS26, SBS44) in Colorectal &
#      Endometrial                             -> enriched in MSI-H tumors
#
# Each plot shows the proportion of samples in which the signature is detected
# (exposure > 0), stratified by the relevant clinical/genomic feature, with a
# Fisher's exact test p-value.
#
# NOTE: DeepSig fits a different signature panel per cancer type, so not every
# signature is fit in every cancer type (e.g. SBS6/SBS15 are not fit in
# Endometrial). We only score samples for which a signature was actually fit.
#
# @author: Ethan Tse
# @date: 2026-06-29
##################

source("./src/00_shared/load_libraries.R")

long_mutsig_df <- read_tsv("./code_outputs/stable_data/msk/impact/mutation_signatures/deepsig/deepsig_exposures_long_format.tsv") %>%
  mutate(
    detected = exposure > 0,
    CANCER_TYPE_DETAILED = ifelse(
      CANCER_TYPE_DETAILED %in% c("Colorectal Adenocarcinoma", "Colon Adenocarcinoma", "Rectal Adenocarcinoma"),
      "Colorectal Adenocarcinoma", CANCER_TYPE_DETAILED
    )
  )

# POLE mutation status and MSI status from the oncogenic mutation matrix
mutation_status <- read_tsv("./code_outputs/stable_data/msk/impact/snvs/matrix/panel_corrected/bmi_oncokb_annotated_first_sample_oncogenic_mutation_mtx.tsv") %>%
  dplyr::select(SAMPLE_ID, POLE, MSI_TYPE)

# cancer (sub)types used for the plots below, at CANCER_TYPE_DETAILED resolution
crc_type <- "Colorectal Adenocarcinoma"
ucec_type <- "Uterine Endometrioid Carcinoma"
nsclc_type <- "Lung Adenocarcinoma"

# nicer signature labels (signature + mutational process)
signature_labels <- c(
  SBS1  = "SBS1 (Clock-like)",
  SBS4  = "SBS4 (Tobacco)",
  SBS5  = "SBS5 (Clock-like)",
  SBS10 = "SBS10 (POLE)",
  SBS6  = "SBS6 (MMR)",
  SBS14 = "SBS14 (MMR/POLE)",
  SBS15 = "SBS15 (MMR)",
  SBS26 = "SBS26 (MMR)",
  SBS44 = "SBS44 (MMR)"
)

# ---- helpers ---------------------------------------------------------------

# Fisher's exact test p-value for detected (TRUE/FALSE) vs a 2-level group
fisher_p <- function(detected, group) {
  tab <- table(group, detected)
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(NA_real_)
  }
  fisher.test(tab)$p.value
}

# Report the exact p-value rather than truncating everything below 0.001 to
# "p < 0.001": these validation panels exist to show how strong the known
# signature-phenotype associations are, and a shared "p < 0.001" floor hides
# the ~140 orders of magnitude separating SBS44/MSI-H from SBS26/MSI-H.
# Small values are written in scientific notation ("p = 6.6e-58"); figs.qmd's
# relabel_p() renders that as italic P = 6.6 x 10^-58 at presentation time.
# Fisher's exact p underflows to exactly 0 for the strongest comparisons
# (CRC SBS6), which is reported against the double-precision floor instead.
pval_label <- function(p) {
  fmt <- function(x) {
    if (is.na(x)) {
      return("p = NA")
    }
    if (x == 0) {
      return("p < 1e-300")
    }
    if (x >= 1e-3) {
      return(paste0("p = ", format(signif(x, 2), trim = TRUE)))
    }
    # Two significant figures in explicit mantissa/exponent form, so the
    # mantissa is always written to the same precision ("6.0e-72", not
    # format()'s "6e-72") and the exponent is never zero-padded ("5.0e-4",
    # not "5e-04").
    expo <- floor(log10(x))
    mant <- signif(x / 10^expo, 2)
    if (mant >= 10) {
      mant <- mant / 10
      expo <- expo + 1
    }
    paste0("p = ", formatC(mant, format = "f", digits = 1), "e", expo)
  }
  vapply(p, fmt, character(1), USE.NAMES = FALSE)
}

# two-tone palette (low vs high group)
group_pal <- c("#9ecae1", "#c44c4c")

fig_dir <- "./code_outputs/results/figs/bmi/mutagenesis/mutation_signatures/deepsig/common_signature_validation/"
create_dir_if_needed(fig_dir)

# ===========================================================================
# Plot 1: SBS4 (tobacco) in NSCLC, smokers vs never-smokers
# ===========================================================================
sbs4_lung <- long_mutsig_df %>%
  filter(
    CANCER_TYPE_DETAILED == nsclc_type,
    signature == "SBS4",
    SMOKING_PREDICTIONS_3_CLASSES %in% c("Never", "Former/Current Smoker")
  ) %>%
  mutate(group = factor(SMOKING_PREDICTIONS_3_CLASSES, levels = c("Never", "Former/Current Smoker")))

sbs4_summary <- sbs4_lung %>%
  group_by(group) %>%
  summarise(n_detected = sum(detected), n_total = dplyr::n(), frac = n_detected / n_total, .groups = "drop")

sbs4_p <- fisher_p(sbs4_lung$detected, sbs4_lung$group)

plot1 <- sbs4_summary %>%
  ggplot(aes(x = group, y = frac, fill = group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(scales::percent(frac, accuracy = 1), "\n(n=", n_total, ")")),
    vjust = -0.3, size = 2.5
  ) +
  scale_fill_manual(values = group_pal) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.18))) +
  theme_std() +
  theme(legend.position = "none") +
  labs(
    title = "SBS4 (Tobacco) in Lung Adenocarcinoma",
    subtitle = pval_label(sbs4_p),
    x = NULL,
    y = "Samples with signature detected"
  )

ggsave_wrapper(paste0(fig_dir, "sbs4_lung_by_smoking.pdf"), plot1, width = 2.5, height = 2.8)

# ===========================================================================
# Plot 2: SBS10 (POLE) in Endometrial & Colorectal, POLE-mutant vs POLE-WT
# ===========================================================================
sbs10_df <- long_mutsig_df %>%
  filter(
    CANCER_TYPE_DETAILED %in% c(ucec_type, crc_type),
    signature == "SBS10"
  ) %>%
  inner_join(mutation_status, by = c(sid = "SAMPLE_ID")) %>%
  filter(!is.na(POLE)) %>%
  mutate(group = factor(ifelse(POLE == 1, "POLE mutant", "POLE WT"), levels = c("POLE WT", "POLE mutant")))

sbs10_summary <- sbs10_df %>%
  group_by(CANCER_TYPE_DETAILED, group) %>%
  summarise(n_detected = sum(detected), n_total = dplyr::n(), frac = n_detected / n_total, .groups = "drop")

sbs10_p <- sbs10_df %>%
  group_by(CANCER_TYPE_DETAILED) %>%
  summarise(p = fisher_p(detected, group), .groups = "drop") %>%
  mutate(p_label = pval_label(p))

plot2 <- sbs10_summary %>%
  ggplot(aes(x = group, y = frac, fill = group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(scales::percent(frac, accuracy = 1), "\n(n=", n_total, ")")),
    vjust = -0.3, size = 2.5
  ) +
  geom_text(
    data = sbs10_p, aes(x = 1.5, y = Inf, label = p_label),
    inherit.aes = FALSE, vjust = 1.3, size = 2.5
  ) +
  facet_wrap(vars(CANCER_TYPE_DETAILED), scales = "free_x") +
  scale_fill_manual(values = group_pal) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.25))) +
  theme_std() +
  theme(legend.position = "none") +
  labs(
    title = "SBS10 (POLE) by cancer type",
    x = NULL,
    y = "Samples with signature detected"
  )

ggsave_wrapper(paste0(fig_dir, "sbs10_by_pole.pdf"), plot2, width = 3, height = 2.8)

# ===========================================================================
# Plot 3: MMR signatures in Colorectal & Endometrial, MSI-H vs MSS
# ===========================================================================
mmr_sigs <- c("SBS6", "SBS14", "SBS15", "SBS26", "SBS44")

mmr_df <- long_mutsig_df %>%
  filter(
    CANCER_TYPE_DETAILED %in% c(ucec_type, crc_type),
    signature %in% mmr_sigs
  ) %>%
  inner_join(mutation_status, by = c(sid = "SAMPLE_ID")) %>%
  mutate(group = dplyr::case_when(
    MSI_TYPE == "Instable" ~ "MSI-H",
    MSI_TYPE == "Stable" ~ "MSS",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group)) %>%
  mutate(
    group = factor(group, levels = c("MSS", "MSI-H")),
    signature = factor(signature, levels = mmr_sigs)
  )

mmr_summary <- mmr_df %>%
  group_by(CANCER_TYPE_DETAILED, signature, group) %>%
  summarise(n_detected = sum(detected), n_total = dplyr::n(), frac = n_detected / n_total, .groups = "drop")

mmr_p <- mmr_df %>%
  group_by(CANCER_TYPE_DETAILED, signature) %>%
  summarise(p = fisher_p(detected, group), .groups = "drop") %>%
  mutate(p_label = pval_label(p))

plot3 <- mmr_summary %>%
  ggplot(aes(x = group, y = frac, fill = group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = scales::percent(frac, accuracy = 1)), vjust = -0.3, size = 2) +
  geom_text(
    data = mmr_p, aes(x = 1.5, y = Inf, label = p_label),
    inherit.aes = FALSE, vjust = 1.3, size = 2
  ) +
  facet_wrap(vars(CANCER_TYPE_DETAILED, signature),
    scales = "free",
    ncol = 4,
    labeller = labeller(signature = signature_labels)
  ) +
  scale_fill_manual(values = group_pal) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.25))) +
  theme_std() +
  theme(legend.position = "none") +
  labs(
    title = NULL,
    x = NULL,
    y = "Samples with signature detected"
  )


ggsave_wrapper(paste0(fig_dir, "mmr_sigs_by_msi.pdf"), plot3, width = 3.5, height = 3.8)
