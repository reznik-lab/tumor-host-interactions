#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 00.0 -- export the `procedureOccurrence` domain for the cancer_patients cohort.
##
## Produces no manuscript figure or table: this is one of the BigQuery
## exports the case arm is built from.
##
## Reads:
##   BigQuery `YOUR-PROJECT.C2025Q4R6` -- Researcher Workbench only.
## Writes:
##   gs://YOUR-CANCER-COHORT-BUCKET/bq_exports/procedureOccurrence/<YYYYMMDD>/procedureOccurrence/procedureOccurrence_*.csv
##   <workspace>/cdrv9/cancer_patients/aou_cancer_patients_data/bq_exports/procedureOccurrence/<YYYYMMDD>/procedureOccurrence/procedureOccurrence_*.csv
##     -- the second only when the shards are mirrored down, which is the
##        default; --no-download stops at the bucket.
##
## This costs BigQuery query and storage, so it is NOT part of run_all.sh.
## See scripts/00_query/run_queries.sh for the flags and the usual invocation.
## ---------------------------------------------------------------------------

## --- locations -------------------------------------------------------------
## Found from this script's own path, so it runs from any directory.
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (!startsWith(.f, "/")) .f <- file.path(Sys.getenv("PWD", getwd()), .f)
REPO <- normalizePath(file.path(dirname(.f), "..", "..", ".."), mustWork = TRUE)

source(file.path(REPO, "src", "setup.R"))
source(file.path(REPO, "scripts", "00_query", "bq_export.R"))


## --- where this cohort's export goes ---------------------------------------
## WORKSPACE-SPECIFIC, AND REPLACED FOR PUBLICATION. These named the workspace
## the analysis ran in and are not usable from another one. Set them to your
## own, or override from the environment: $EXPORT_BUCKET, $AOU_BILLING_PROJECT,
## $AOU_EXPORT_DATE (see bq_export.R). QUERY_DATASET must name your workspace's
## copy of the Controlled Tier CDR; the manuscript used release C2025Q4R6.

DOMAIN        <- "procedureOccurrence"
BUCKET        <- "gs://YOUR-CANCER-COHORT-BUCKET"
BILLING       <- "YOUR-BILLING-PROJECT"
QUERY_DATASET <- "YOUR-PROJECT.C2025Q4R6_index"
LOCAL_SUBDIR  <- "cdrv9/cancer_patients/aou_cancer_patients_data"


## --- the cohort, as the Workbench generated it -----------------------------
procedureOccurrence_sql = r"(

    SELECT
        st.modifier_concept_id,
        dt0.concept_name AS T_DISP_modifier_concept_name,
        st.modifier_source_value,
        st.person_id,
        st.procedure_concept_id,
        st.procedure_datetime,
        st.procedure_source_concept_id,
        st.procedure_source_value,
        st.procedure_type_concept_id,
        dt1.concept_name AS T_DISP_procedure_type_concept_name,
        st.quantity,
        dt2.concept_code AS T_DISP_source_concept_code,
        dt2.concept_name AS T_DISP_source_concept_name,
        dt2.vocabulary_id AS T_DISP_source_vocabulary,
        dt3.concept_code AS T_DISP_standard_concept_code,
        dt3.concept_name AS T_DISP_standard_concept_name,
        dt3.vocabulary_id AS T_DISP_standard_vocabulary,
        st.visit_occurrence_id 
    FROM
        `YOUR-PROJECT.C2025Q4R6`.procedure_occurrence AS st 
    LEFT JOIN
        `YOUR-PROJECT.C2025Q4R6`.concept AS dt0 
            ON dt0.concept_id = st.modifier_concept_id 
    LEFT JOIN
        `YOUR-PROJECT.C2025Q4R6`.concept AS dt1 
            ON dt1.concept_id = st.procedure_type_concept_id 
    LEFT JOIN
        `YOUR-PROJECT.C2025Q4R6`.concept AS dt2 
            ON dt2.concept_id = st.procedure_source_concept_id 
    LEFT JOIN
        `YOUR-PROJECT.C2025Q4R6`.concept AS dt3 
            ON dt3.concept_id = st.procedure_concept_id 
    WHERE
        st.procedure_occurrence_id IN (SELECT
            id 
        FROM
            `@CDR@`.T_ESA_procedureOccurrence_person_id 
        WHERE
            (((person_id IN (SELECT
                person_id AS primary_id 
            FROM
                `@CDR@`.T_ENT_conditionOccurrence 
            WHERE
                (condition_concept_id IN (SELECT
                    descendant 
                FROM
                    `@CDR@`.T_HAD_conditionConcept_default 
                WHERE
                    ancestor IN (439392, 317510, 432571, 437233, 4175485) 
                UNION
                ALL SELECT
                    439392 
                UNION
                ALL SELECT
                    317510 
                UNION
                ALL SELECT
                    432571 
                UNION
                ALL SELECT
                    437233 
                UNION
                ALL SELECT
                    4175485))))) 
                AND (NOT ((person_id IN (SELECT
                    person_id AS primary_id 
                FROM
                    `@CDR@`.T_ENT_ingredientOccurrence 
                WHERE
                    drug_concept_id IN (SELECT
                        descendant 
                    FROM
                        `@CDR@`.T_HAD_ingredientConcept_default 
                    WHERE
                        ancestor = 1123618 
                    UNION
                    ALL SELECT
                        1123618))) 
                    OR (person_id IN (SELECT
                        person_id AS primary_id 
                    FROM
                        `@CDR@`.T_ENT_procedureOccurrence 
                    WHERE
                        procedure_concept_id IN (SELECT
                            descendant 
                        FROM
                            `@CDR@`.T_HAD_procedureConcept_default 
                        WHERE
                            ancestor = 4326683 
                        UNION
                        ALL SELECT
                            4326683))))) 
                        AND (person_id IN (SELECT
                            person_id AS primary_id 
                        FROM
                            `@CDR@`.T_ENT_measurementOccurrence 
                        WHERE
                            measurement_concept_id = 3038553))) 
                        AND (procedure_concept_id IN (SELECT
                            descendant 
                        FROM
                            `@CDR@`.T_HAD_procedureConcept_default 
                        WHERE
                            ancestor = 4326683 
                        UNION
                        ALL SELECT
                            4326683))))"


## --- run -------------------------------------------------------------------
export_bq_domain(
  domain        = DOMAIN,
  sql           = gsub("@CDR@", QUERY_DATASET, procedureOccurrence_sql, fixed = TRUE),
  bucket        = BUCKET,
  billing       = BILLING,
  query_dataset = QUERY_DATASET,
  local_subdir  = LOCAL_SUBDIR)
