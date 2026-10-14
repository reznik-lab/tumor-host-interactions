#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 06.1 -- cases per cancer type in the matched cohort.
##
## Produces: Supplementary Figure 1a -- matched cases per cancer type.
##           Supplementary Table 3 -- the cancer type consolidation map.
##
## Reads:
##   <workspace>/.../stable_data/psm/solid_heme/full_cohort/aggregate/psm_matched_cohort_all_cancer_types.tsv
## Writes:
##   .../results/tables/data_overview/cancer_type_mapping.tsv (+ .xlsx)
##       [Supplementary Table 3]
##   .../results/figs/data_overview/match_cancer_type_breakdown.pdf (+ .rds)
##       [Supplementary Figure 1a]
## ---------------------------------------------------------------------------

## --- locations -------------------------------------------------------------
## Found from this script's own path, so it runs from any directory.
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (!startsWith(.f, "/")) .f <- file.path(Sys.getenv("PWD", getwd()), .f)
REPO <- normalizePath(file.path(dirname(.f), "..", ".."), mustWork = TRUE)

## Data root: two levels above the repo. Override with $AOU_WORKSPACE.
WORKSPACE <- normalizePath(Sys.getenv(
  "AOU_WORKSPACE", unset = file.path(REPO, "..", "..", "workspace")),
  mustWork = TRUE)

source(file.path(REPO, "src", "setup.R"))


## --- read ------------------------------------------------------------------

cancers_cohort <- read_tsv(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/psm_matched_usable_data/solid_heme/full_cohort/all_cancers_geq1_before_diagnosis_after_psm.tsv"))


## --- cancer_type -> match_cancer_type mapping ------------------------------

cancer_type_mapping_with_counts <- cancers_cohort %>%
  select(cancer_type, match_cancer_type) %>%
  distinct() %>%
  arrange(match_cancer_type)

print(cancer_type_mapping_with_counts)

mapping_path <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/results/tables/data_overview/cancer_type_mapping.tsv")

dir.create(dirname(mapping_path), recursive = TRUE, showWarnings = FALSE)
write_table(cancer_type_mapping_with_counts, mapping_path)


## --- Supp Fig 1B -----------------------------------------------------------

## cancers_cohort is the case arm only, so this count is a count of cases.
type_counts <- cancers_cohort %>%
  select(person_id, match_cancer_type) %>%
  distinct() %>%
  count(match_cancer_type, name = "n")

## Report what the floor removes, so the dropped types are never silent.
message(sprintf(
  "Supp Fig 1B: %d of %d cancer types and %d of %d cases clear MIN_CASES_PER_TYPE = %d.",
  sum(type_counts$n >= MIN_CASES_PER_TYPE), nrow(type_counts),
  sum(type_counts$n[type_counts$n >= MIN_CASES_PER_TYPE]), sum(type_counts$n),
  MIN_CASES_PER_TYPE))
if (any(type_counts$n < MIN_CASES_PER_TYPE)) {
  dropped <- type_counts %>% filter(n < MIN_CASES_PER_TYPE) %>% arrange(desc(n))
  message("  excluded: ", paste(sprintf("%s (%d)", dropped$match_cancer_type, dropped$n),
                                collapse = ", "))
}

## The bars' own numbers, written out: the figure is a supplementary figure, so
## the counts behind it are a supplementary table rather than something to be
## read off the bar labels. Raw and display label side by side, because the
## raw value is the key every other table joins on and the display label is
## what the axis says.
cases_per_type <- type_counts %>%
  filter(n >= MIN_CASES_PER_TYPE) %>%
  transmute(match_cancer_type,
            cancer_type_label = display_type_label(match_cancer_type),
            n_cases = n) %>%
  arrange(desc(n_cases))

## Axis labels go through display_type_label() (src/setup.R), the one mapping
## every figure in the pipeline uses, rather than printing the raw
## cancer_type_grouped keys. The mapping is injective (asserted at source()
## time), so reordering on n and then labelling cannot collapse two bars.
match_cancer_type_breakdown <- cases_per_type %>%
  mutate(cancer_type_label = fct_reorder(cancer_type_label, n_cases)) %>%
  ggplot(aes(x = cancer_type_label, y = n_cases)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(n_cases, big.mark = ",", trim = TRUE)),
            hjust = -0.15, size = 2) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                     labels = function(x) format(x, big.mark = ",", trim = TRUE)) +
  labs(
    x = NULL,
    y = "Number of patients",
    subtitle = sprintf("%s matched cases across %d cancer types (>= %d cases each)",
                       format(sum(cases_per_type$n_cases), big.mark = ","),
                       nrow(cases_per_type), MIN_CASES_PER_TYPE)
  ) +
  theme_std() +
  theme(legend.position = "none")

ggsave_wrapper(file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/results/figs/data_overview/match_cancer_type_breakdown.pdf"),
  match_cancer_type_breakdown,
  width = 5, height = 5)
