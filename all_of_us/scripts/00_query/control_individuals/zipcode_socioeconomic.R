#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 00.0 -- export the `zipcodeSocioeconomic` domain for the control_individuals cohort.
##
## Produces no manuscript figure or table: this is one of the BigQuery
## exports the control pool is built from.
##
## Reads:
##   BigQuery `YOUR-PROJECT.C2025Q4R6` -- Researcher Workbench only.
## Writes:
##   gs://YOUR-CONTROL-COHORT-BUCKET/bq_exports/zipcodeSocioeconomic/<YYYYMMDD>/zipcodeSocioeconomic/zipcodeSocioeconomic_*.csv
##   <workspace>/cdrv9/control_individuals/aou_control_individuals_cohort_data/bq_exports/zipcodeSocioeconomic/<YYYYMMDD>/zipcodeSocioeconomic/zipcodeSocioeconomic_*.csv
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

DOMAIN        <- "zipcodeSocioeconomic"
BUCKET        <- "gs://YOUR-CONTROL-COHORT-BUCKET"
BILLING       <- "YOUR-BILLING-PROJECT"
QUERY_DATASET <- "YOUR-PROJECT.C2025Q4R6_index"
LOCAL_SUBDIR  <- "cdrv9/control_individuals/aou_control_individuals_cohort_data"


## --- the cohort, as the Workbench generated it -----------------------------
zipcodeSocioeconomic_sql = r"(

    SELECT
        american_community_survey_year,
        assisted_income,
        deprivation_index,
        high_school_education,
        id,
        median_income,
        no_health_insurance,
        observation_datetime_str,
        observation_datetime_utc,
        person_id,
        poverty,
        vacant_housing,
        zip_code 
    FROM
        `@CDR@`.T_ENT_zipcodeSocioeconomic 
    WHERE
        (
            person_id IN (SELECT
                person_id AS primary_id 
            FROM
                `@CDR@`.T_ENT_measurementOccurrence 
            WHERE
                measurement_concept_id = 3038553)
        ) 
        AND (
            NOT ((person_id IN (SELECT
                person_id AS primary_id 
            FROM
                `@CDR@`.T_ENT_conditionOccurrence 
            WHERE
                condition_concept_id IN (SELECT
                    descendant 
                FROM
                    `@CDR@`.T_HAD_conditionConcept_default 
                WHERE
                    ancestor = 439392 
                UNION
                ALL SELECT
                    439392))) 
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
                        4326683))) 
                    OR (person_id IN (SELECT
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
                            `@CDR@`.T_ENT_conditionOccurrence 
                        WHERE
                            condition_concept_id IN (SELECT
                                descendant 
                            FROM
                                `@CDR@`.T_HAD_conditionConcept_default 
                            WHERE
                                ancestor = 432571 
                            UNION
                            ALL SELECT
                                432571))) 
                            OR (person_id IN (SELECT
                                person_id AS primary_id 
                            FROM
                                `@CDR@`.T_ENT_conditionOccurrence 
                            WHERE
                                condition_concept_id IN (SELECT
                                    descendant 
                                FROM
                                    `@CDR@`.T_HAD_conditionConcept_default 
                                WHERE
                                    ancestor = 437233 
                                UNION
                                ALL SELECT
                                    437233))) 
                                OR (person_id IN (SELECT
                                    person_id AS primary_id 
                                FROM
                                    `@CDR@`.T_ENT_conditionOccurrence 
                                WHERE
                                    condition_concept_id IN (SELECT
                                        descendant 
                                    FROM
                                        `@CDR@`.T_HAD_conditionConcept_default 
                                    WHERE
                                        ancestor = 317510 
                                    UNION
                                    ALL SELECT
                                        317510))) 
                                    OR (person_id IN (SELECT
                                        person_id AS primary_id 
                                    FROM
                                        `@CDR@`.T_ENT_conditionOccurrence 
                                    WHERE
                                        condition_concept_id IN (SELECT
                                            descendant 
                                        FROM
                                            `@CDR@`.T_HAD_conditionConcept_default 
                                        WHERE
                                            ancestor = 4175485 
                                        UNION
                                        ALL SELECT
                                            4175485))) 
                                        OR (person_id IN (SELECT
                                            person_id AS primary_id 
                                        FROM
                                            `@CDR@`.T_ENT_conditionOccurrence 
                                        WHERE
                                            condition_concept_id IN (SELECT
                                                descendant 
                                            FROM
                                                `@CDR@`.T_HAD_conditionConcept_default 
                                            WHERE
                                                ancestor = 443392 
                                            UNION
                                            ALL SELECT
                                                443392))))
                                        ))"


## --- run -------------------------------------------------------------------
export_bq_domain(
  domain        = DOMAIN,
  sql           = gsub("@CDR@", QUERY_DATASET, zipcodeSocioeconomic_sql, fixed = TRUE),
  bucket        = BUCKET,
  billing       = BILLING,
  query_dataset = QUERY_DATASET,
  local_subdir  = LOCAL_SUBDIR)
