# Records of pre-diagnostic tumor-host interactions in a pan-cancer atlas of body composition

Analysis code for the **All of Us** arm of:

> Ethan Tse, Sonia Boscenco, Neil Neumann, Hyung Jun Woo, Mark T. A. Donoghue,
> Christopher A. Barker, John Garrett, Perry Pickhardt, Nathaniel Swinburne,
> Ed Reznik. *Records of pre-diagnostic tumor-host interactions in a pan-cancer
> atlas of body composition.*

This repository covers the longitudinal weight-trajectory analysis of the All of
Us cohort: cohort construction and propensity score matching (34,378 cases and
153,311 matched controls), pre-diagnostic cachexia calling, BMI trajectories at
landmark timepoints, and the weight-class shift analysis — Figure 1A–F and their
supplementary figures and tables.

## Getting access to the All of Us workbench

Everything runs on the **All of Us Researcher Workbench** (Controlled Tier,
C2025Q4R6, cdrv9). You will need to create an account and register for the controlled data tier.

## Reporting results

Results in the manuscript were released under the All of Us data dissemination
policy, which requires >20 individuals per reported category.

## Layout

```
scripts/run_all.sh       the driver for the analysis — start here
scripts/
  00_query/              BigQuery: the cohort definitions and their exports
    run_queries.sh         the driver for the exports; NOT part of run_all.sh
    bq_export.R            shared BigQuery -> Cloud Storage -> workspace transport
    cancer_patients/       one script per exported domain, cases
    control_individuals/   one script per exported domain, control pool
    fig1a_exclusion_cascade.R   the Figure 1a exclusion counts; run on its own
  00_ingest/             assembling the exports into analysis tables
  01_cohort/             cohort preprocessing, one script per arm
  02_matching/           propensity score matching and its quality control
  03_cachexia/           cachexia episode calling and prevalence
  04_trajectories/       BMI trajectories and the landmark mixed-effects models
  05_weight_class/       weight-class shift before diagnosis
  06_tables/             cohort summary tables
src/setup.R              shared libraries, plotting theme, and the pipeline-wide
                         definitions: BMI usability, the minimum-cases-per-type
                         floor, the landmark windows, the tissue-group and
                         display-label maps, and the .tsv + .xlsx table writer
env.yml                  the conda environment
```

Every script finds the repository from its own path, so none of them depend on
the working directory, and every one of them sources `src/setup.R`.

### Where the data lives

No data is in this repository, and `.gitignore` is written so that none can be
added. The scripts read and write a workspace tree **two levels above the
repository** — `<repo>/../../workspace` — holding

```
cdrv9/                          the BigQuery exports, per cohort
pre_diagnosis_bmi_analysis/code_outputs_scripts/
  stable_data/                  the cohort and matched tables
  temp/                         intermediates
  cachexia_episodes/            the episode caller's output
  results/figs/, results/tables/  the manuscript artifacts
  logs/                         one log per step, written by the drivers
```

Override the root with `$AOU_WORKSPACE`. Every script honours it.

## Running it

```bash
./scripts/run_all.sh --list        # the steps of this phase and what each produces
./scripts/run_all.sh --dry-run     # print the commands, run nothing
./scripts/run_all.sh pre           # stages 00–02: everything before the Python caller
./scripts/run_all.sh python        # the cachexia episode caller
./scripts/run_all.sh post          # stages 03–06: everything after it
./scripts/run_all.sh               # all three, stopping at the first failure
```

`--keep-going` carries on past a failure and lists what failed at the end.
Interpreters default to the Workbench locations and fall back to `$PATH`;
override with `$RSCRIPT_BIN` and `$PYTHON_BIN`. Each step logs to
`<workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/logs/`.

The run splits at the Python step because `postprocess_cachexia_calls.R` reads a
pinned episode date rather than the newest file 

### The BigQuery exports are not part of `run_all.sh`

Stage 00.0 exports both cohorts out of BigQuery into Cloud Storage and mirrors
the shards down into the workspace. It costs query and storage, it needs
Workbench credentials, and its output changes only when the cohort definition
does, so it is run by hand:

```bash
./scripts/00_query/run_queries.sh --list       # the domains and what they feed
./scripts/00_query/run_queries.sh              # both cohorts, every domain
./scripts/00_query/run_queries.sh cancer       # the cases arm only
./scripts/00_query/run_queries.sh control      # the control pool only
```

Flags passed through to each domain script: `--no-query` (download what is
already in the bucket), `--no-download` (extract to Cloud Storage only),
`--date=YYYYMMDD`, `--dry-run`. The date stamp is part of the output path and is
computed **once per run**, so a run that straddles midnight cannot leave half
the domains under one date and half under the next.

`00_query/fig1a_exclusion_cascade.R` is in neither driver. It recomputes the
Figure 1a exclusion counts straight from the CDR and the dataset index, and it
takes `--print-sql`.

### You must set your own workspace identifiers

The bucket, the billing project and the CDR dataset named in `scripts/00_query/`
are **placeholders** — `YOUR-CANCER-COHORT-BUCKET`, `YOUR-CONTROL-COHORT-BUCKET`,
`YOUR-BILLING-PROJECT`, `YOUR-PROJECT.C2025Q4R6_index`. The originals identified
the workspace this analysis ran in and are not usable from another one. Set
`BUCKET`, `BILLING` and `QUERY_DATASET` per cohort at the top of each domain
script, or export `$EXPORT_BUCKET` and `$AOU_BILLING_PROJECT`. `QUERY_DATASET`
must name your workspace's copy of the Controlled Tier CDR; the manuscript used
release C2025Q4R6. `fig1a_exclusion_cascade.R` takes its two datasets from
`$AOU_CDR` and `$AOU_INDEX` instead. `00_ingest/pull_aou_enrollment_date.R` uses
`$GOOGLE_PROJECT` and `$WORKSPACE_CDR`, which the Workbench sets for you.

The dataset index is a snapshot and can be garbage-collected. If its `T_ENT_*`
tables are gone, re-export the dataset from the Workbench and point `INDEX` at
the new stamp; the concept sets are unaffected.

## What produces what

| Artifact | Script |
| --- | --- |
| Figure 1a (counts) | `00_query/fig1a_exclusion_cascade.R` |
| Figure 1b | `03_cachexia/cachexia_frequency.R` |
| Figure 1c–f, trajectories | `04_trajectories/analysis.R` |
| Figure 1c–f, landmark boxplots | `04_trajectories/mixed_effects_modelling.R` |
| Supp Fig 1a, Supp Table 3 | `06_tables/overview_of_cancers_in_cohort.R` |
| Supp Fig 1b, Supp Table 4 | `02_matching/psm_quality_control.R` |
| Supp Fig 1c, Supp Table 5 | `03_cachexia/cachexia_frequency_by_threshold.R` |
| Supp Fig 1d, Supp Table 7 | `05_weight_class/weight_cat_dist_at_timepoints.R` |
| Supp Table 1 | `06_tables/supplementary_table.R` |
| Supp Table 2 | `02_matching/supplementary_psm_quantification.R` |
| Supp Table 6 | `04_trajectories/supp_landmarks_all_types.R` |

Figure 1a itself is a schematic and is not drawn by any script here; only its
counts are computed. Every other script in the pipeline produces no manuscript
artifact — it says so in its own header — and exists to feed the ones above.
Tables are written twice, as `.tsv` and `.xlsx`, by one writer in `src/setup.R`.

## Reproducing the environment

`env.yml` is the exact environment the analysis ran in — R 4.5.3 and Python
3.14, with explicit builds — and its `name:` is deliberately blank, so pick one:

```bash
mamba env create -f env.yml -n aou
```

The Workbench images put this environment at `/opt/conda/envs/packages`, which
is where `src/setup.R` prepends its R library path and where the drivers look
for `Rscript` and `python3` first. Elsewhere, activate the environment, or point
`$RSCRIPT_BIN` and `$PYTHON_BIN` at it.
