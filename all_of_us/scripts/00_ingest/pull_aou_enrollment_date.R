#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 00.2 -- pull approximate All of Us enrollment date from BigQuery.
##
## Produces no manuscript figure or table: supplies the enrollment date the
## propensity score match uses as a covariate.
##
## Reads:
##   BigQuery ($WORKSPACE_CDR) -- Researcher Workbench only; needs
##   $GOOGLE_PROJECT and $WORKSPACE_CDR set.
## Writes:
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/aou_enrollment_date.tsv
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
library(bigrquery)


## --- pull ------------------------------------------------------------------

download_data <- function(query) {
  tb <- bq_project_query(Sys.getenv("GOOGLE_PROJECT"), query = query,
                         default_dataset = Sys.getenv("WORKSPACE_CDR"))
  bq_table_download(tb)
}

DATASET <- Sys.getenv("WORKSPACE_CDR")

df <- download_data(str_glue("SELECT DISTINCT
person_id,
MIN(observation_date) AS primary_consent_date
FROM `{DATASET}.concept`
JOIN `{DATASET}.concept_ancestor` on concept_id = ancestor_concept_id
JOIN `{DATASET}.observation` on descendant_concept_id = observation_concept_id
WHERE concept_name = 'Consent PII' AND concept_class_id = 'Module'
GROUP BY 1"))

print(head(df))


## --- write -----------------------------------------------------------------

out_path <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/aou_enrollment_date.tsv")

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write_tsv(df, out_path)

message("wrote ", nrow(df), " rows -> ", out_path)
