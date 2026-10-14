#!/usr/bin/env bash
# Runs every analysis in dependency order. Run from the repository root.
#
# Inputs, neither of which ships with the code:
#   data/                       restricted patient-level data (symlink it here)
#   code_outputs/stable_data/   the pre-assembled cohort metadata tables
# Both are read-only. All output is written under code_outputs/.
#
# Environment: R 4.4.3 and the packages pinned in env/environment.yml.
#   mamba env create -f env/environment.yml && mamba activate tumor-host-interactions
#
# Ordering: the Tier 0 analyses read only the two inputs above, so they can run
# in any order among themselves. Each Tier 1 analysis reads a table written by a
# Tier 0 analysis, so every Tier 0 analysis must finish first. Running a Tier 1
# script on its own, before its producer, fails on the missing input.
set -euo pipefail

run() { echo "== $1"; Rscript "$1"; }

# ── Tier 0 (25 analyses, no intra-bundle dependencies) ──
run src/06_supplementary/msk_cohort_characteristics_table.R
run src/01_figure1/oac_vs_noac_lm_msk.R
run src/06_supplementary/pathways_logistic_regression.R
run src/06_supplementary/pathways_logistic_regression_sex_stratified.R
run src/05_figure5/mutation_composition.R
run src/02_figure2/abnormal_physiology.R
run src/02_figure2/myopenia_myosteatosis_by_cancer_type_mixed_model.R
run src/06_supplementary/myopenia_one_vs_rest.R
run src/06_supplementary/data_availability.R
run src/06_supplementary/lung_pik3ca_bmi_analysis.R
run src/05_figure5/lung_pik3ca_bodycomp_analysis.R
run src/06_supplementary/lung_pik3ca_bodycomp_stk11_exclusion_sensitivity.R
run src/02_figure2/bodycomp_vs_background_one_vs_rest.R
run src/06_supplementary/organ_involvement_sensitivity_forest_plots.R
run src/05_figure5/pathway_enrichment_by_bodycomp.R
run src/05_figure5/pi3k_bodycomp_analysis.R
run src/04_figure4/bodycomp_associations.R
run src/04_figure4/melanoama_mcc_comparison.R
run src/04_figure4/pathology_vs_genomic_viral_status_contingency.R
run src/04_figure4/scan_bodycomp_barplot.R
run src/06_supplementary/summary_statistics.R
run src/03_figure3/bmi_and_mutsigs.R
run src/06_supplementary/bodycomp_and_mutsigs.R
run src/06_supplementary/smoking_and_mutsigs.R
run src/04_figure4/snvs_bmi_logistic_regression.R

# ── Tier 1 (10 analyses, each reads a Tier 0 output) ──
run src/05_figure5/pathways_bmi_association_volcano.R
run src/02_figure2/fig2ef_contingency_tables.R
run src/06_supplementary/myopenia_myosteatosis_glmm_forest.R
run src/06_supplementary/myopenia_one_vs_rest_forest_plots.R
run src/06_supplementary/lung_pik3ca_clean_scan_arms_egfr_panel.R
run src/02_figure2/emmeans_heatmap_one_vs_rest.R
run src/06_supplementary/pathway_enrichment_by_bodycomp_tx_naive_bmi.R
run src/05_figure5/pi3k_bodycomp_forest_plots.R
run src/06_supplementary/pi3k_bodycomp_treatment_naive_sensitivity.R
run src/04_figure4/snvs_bmi_association_volcano.R
