# Tumor-host interactions — analysis code

Code accompanying *Records of pre-diagnostic tumor-host interactions in a
pan-cancer atlas of body composition*.

The repository holds the statistical analyses behind the paper's figures and
tables. Each script reads the study's cohort tables, fits its models, and writes
its results under `code_outputs/`: one `.pdf` plus a sibling `.rds` for every
figure panel, and a `.tsv` for every result table. Those files are what the
manuscript was assembled from.

Patient-level data is not distributed here — see **Data**.

## Repository layout

| Path | Contents |
|---|---|
| `src/` | 42 R scripts, grouped by the figure they serve: 35 analyses plus 7 shared helpers they `source()`. |
| `run_all.sh` | Runs all 35 analyses in dependency order. |
| `env/` | Conda environment pinning R 4.4.3 and the 28 R packages the scripts use. |

## Requirements

R 4.4.3 and 28 CRAN packages, pinned in `env/environment.yml`:

```bash
mamba env create -f env/environment.yml
mamba activate tumor-host-interactions
```

## Data

**Patient-level data (`data/`)** — MSK-IMPACT clinical, genomic and body
composition tables. These are restricted; availability is described in the
manuscript's data availability statement. No analysis runs without them.


## Running

Scripts use paths relative to the repository root, so run them from there. Put
the input trees in place first:

```bash
ln -s /path/to/data data            # patient-level data
mkdir -p code_outputs/stable_data   # then add the cohort metadata tables
```

`data/` and `code_outputs/stable_data/` are read-only inputs; nothing writes to
them. All output lands in `code_outputs/`.

The chosen-example-scan tables are needed only to build the figures, not to run
any analysis here, so `run_all.sh` does not require them; drop them into
`code_outputs/results/tables/body_comp/scan_selection/` if you are rebuilding
the figures.

```bash
# all 35 analyses, in dependency order
bash run_all.sh

# or one analysis on its own
Rscript src/02_figure2/abnormal_physiology.R
```

`run_all.sh` runs the analyses in two tiers: 25 that depend only on the input
tables, then 10 that read a table an earlier script wrote. Each script's header
names the figure and table its outputs appear in.

## What each script produces

**`01_figure1/`** — Figure 1

| Script | Output appears in |
|---|---|
| `oac_vs_noac_lm_msk.R` | Figure 1, Supp Fig 1, Supp Table 9 |

**`02_figure2/`** — Figure 2

| Script | Output appears in |
|---|---|
| `abnormal_physiology.R` | Figure 2, Supp Table 13 |
| `bodycomp_vs_background_one_vs_rest.R` | Figure 2, Supp Fig 2, Supp Table 11 |
| `emmeans_heatmap_one_vs_rest.R` | Figure 2, Supp Fig 2 |
| `fig2ef_contingency_tables.R` | Figure 2 |
| `myopenia_myosteatosis_by_cancer_type_mixed_model.R` | Supp Fig 2, Supp Table 15 |

**`03_figure3/`** — Figure 3

| Script | Output appears in |
|---|---|
| `bmi_and_mutsigs.R` | Figure 3, Supp Table 17 |

**`04_figure4/`** — Figure 4

| Script | Output appears in |
|---|---|
| `bodycomp_associations.R` | Figure 4, Supp Fig 4, Supp Table 20 |
| `melanoama_mcc_comparison.R` | Figure 4, Supp Fig 4 |
| `pathology_vs_genomic_viral_status_contingency.R` | Figure 4 |
| `scan_bodycomp_barplot.R` | Figure 4 |
| `snvs_bmi_association_volcano.R` | Figure 4, Supp Table 18 |
| `snvs_bmi_logistic_regression.R` | Figure 4, Supp Table 18 |

**`05_figure5/`** — Figure 5

| Script | Output appears in |
|---|---|
| `lung_pik3ca_bodycomp_analysis.R` | Figure 5, Supp Fig 5, Supp Table 24 |
| `mutation_composition.R` | Figure 5, Supp Fig 5 |
| `pathway_enrichment_by_bodycomp.R` | Figure 5, Supp Fig 5, Supp Table 22, Supp Table 23 |
| `pathways_bmi_association_volcano.R` | Figure 5 |
| `pi3k_bodycomp_analysis.R` | Figure 5, Supp Fig 5, Supp Table 23 |
| `pi3k_bodycomp_forest_plots.R` | Figure 5, Supp Table 23 |

**`06_supplementary/`** — Supplementary figures and tables

| Script | Output appears in |
|---|---|
| `bodycomp_and_mutsigs.R` | Supp Fig 3 |
| `data_availability.R` | Supp Fig 2 |
| `lung_pik3ca_bmi_analysis.R` | Supp Fig 5, Supp Table 24 |
| `lung_pik3ca_bodycomp_stk11_exclusion_sensitivity.R` | Supp Fig 5, Supp Table 24 |
| `lung_pik3ca_clean_scan_arms_egfr_panel.R` | Supp Fig 5 |
| `msk_cohort_characteristics_table.R` | Supp Table 8 |
| `myopenia_myosteatosis_glmm_forest.R` | Supp Fig 2, Supp Table 15 |
| `myopenia_one_vs_rest.R` | Supp Fig 2, Supp Table 14 |
| `myopenia_one_vs_rest_forest_plots.R` | Supp Fig 2, Supp Table 14 |
| `organ_involvement_sensitivity_forest_plots.R` | Supp Fig 2, Supp Table 12 |
| `pathway_enrichment_by_bodycomp_tx_naive_bmi.R` | Supp Fig 5, Supp Table 22 |
| `pathways_logistic_regression.R` | Supp Table 22 |
| `pathways_logistic_regression_sex_stratified.R` | Supp Table 22 |
| `pi3k_bodycomp_treatment_naive_sensitivity.R` | Supp Fig 5, Supp Table 23 |
| `smoking_and_mutsigs.R` | Supp Fig 3 |
| `summary_statistics.R` | Supp Table 19 |

Supplementary Tables 1–7 come from the All of Us summaries, and Tables 10, 16,
21 and 25 are drawn straight from the cohort metadata, the pathway gene lists
and the gene-symbol mapping, so no analysis script here stands behind them.

## Shared helpers

Sourced by the analyses, never run on their own — running one directly fails by
design, since each expects the calling script's environment.

| Script (under `src/00_shared/`) | Role |
|---|---|
| `gene_symbol_harmonization.R` | synonymous MSK-IMPACT gene symbols collapsed to one name; this mapping is also Supplementary Table 25 |
| `load_libraries.R` | attaches the shared library set |
| `logistic_regression_helpers.R` | gene × BMI logistic-regression fitting and result writing |
| `obesity_constants.R` | the PI3K/AKT/mTOR gene set every pathway-alteration call in the bundle is made against |
| `pathway_constants.R` | oncogenic pathway gene-set definitions; these lists are also Supplementary Table 21 |
| `styles.R` | shared ggplot theme, significance colouring and forest-plot arrow grobs |
| `utils.R` | `ggsave_wrapper()` — saves each panel as a `.pdf` and a sibling `.rds` of the plot object |

## Scope

The scripts regenerate the panels and result tables. The code that composed
those panels into the manuscript's figures and typeset its supplementary tables
is not included.

| | |
|---|---|
| Figures 2, 3, 4, 5 | every computed panel; the example-CT-scan panels are drawn from supplied tables (see **Data**) |
| Supplementary Figures 2, 3, 4, 5 | every computed panel; same exception |
| Figure 1 | panel 1h; 1b–1g come from All of Us |
| Supplementary Figure 1 | panel 1e; the rest come from All of Us |
| Supplementary tables | the underlying result tables for the MSK tables |