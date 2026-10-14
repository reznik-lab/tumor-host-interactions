#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Figure 1A -- the exclusion counts, recomputed so they compose.
##
## Produces: the exclusion counts drawn on Figure 1a, which is a schematic and
##           is not itself generated here.
##
## Reads:
##   BigQuery `YOUR-PROJECT.C2025Q4R6` (CDR) and
##            `YOUR-PROJECT.C2025Q4R6_index` (dataset index)
## Writes (under .../results/tables/cohort_flow/), each as .tsv and .xlsx:
##   fig1a_criteria_crosstab.tsv      -- the 32-row cross-tabulation
##   fig1a_exclusion_cascade.tsv      -- marginal vs conditional, per step
## ---------------------------------------------------------------------------

## --- locations -------------------------------------------------------------
## Found from this script's own path, so it runs from any directory.
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (!startsWith(.f, "/")) .f <- file.path(Sys.getenv("PWD", getwd()), .f)
REPO <- normalizePath(file.path(dirname(.f), "..", ".."), mustWork = TRUE)

WORKSPACE <- normalizePath(Sys.getenv(
  "AOU_WORKSPACE", unset = file.path(REPO, "..", "..", "workspace")),
  mustWork = FALSE)

source(file.path(REPO, "src", "setup.R"))

PRINT_SQL_ONLY <- "--print-sql" %in% commandArgs(TRUE)


## --- the datasets ----------------------------------------------------------
## Both come from the Workbench dataset export; see any script under
## scripts/00_query/<cohort>/. CDR is the curated data repository, INDEX the
## per-dataset snapshot holding the T_ENT_* entity tables and the T_HAD_*
## concept-ancestor rollups the cohort SQL expands through.
##
## INDEX IS A SNAPSHOT AND CAN BE COLLECTED. If the T_ENT_* tables are gone,
## re-export the dataset from the Workbench and update INDEX to the new stamp;
## the concept sets below are unaffected.

CDR     <- Sys.getenv("AOU_CDR",     unset = "YOUR-PROJECT.C2025Q4R6")
INDEX   <- Sys.getenv("AOU_INDEX",   unset = "YOUR-PROJECT.C2025Q4R6_index")
BILLING <- Sys.getenv("AOU_BILLING_PROJECT", unset = "YOUR-BILLING-PROJECT")


## --- the criteria, as the cohort SQL defines them --------------------------
## Lifted from the two cohort definitions (scripts/00_query/
## cancer_patients/person.R and control_individuals/person.R). Each concept id
## is expanded to its descendants through the matching T_HAD_* table and then
## unioned with itself, exactly as the generated SQL does -- an ancestor is not
## its own descendant in those tables, so dropping the union would silently
## lose everyone coded at the top-level concept.
##
## BMI HAS NO VALUE GATE HERE, deliberately. The cohort SQL requires only that
## a BMI row exists; the [5, 80] range in filter_usable_bmi() (src/setup.R) is
## applied later, in R. So "no EHR BMI" is "no BMI row at all", and a person
## whose only BMI is implausible passes this box and is dropped downstream.
## That is what the published figure counts, and changing it here would make
## this script disagree with the panel for a second reason.

BMI_CONCEPT       <- 3038553L
GLP1_ANCESTOR     <- 1123618L
BARIATRIC_ANCESTOR<- 4326683L

## The case arm's inclusion: primary malignant neoplasm and the haematological
## malignancies grouped with it.
PRIMARY_MALIGNANT <- c(439392L, 317510L, 432571L, 437233L, 4175485L)

## The control arm's exclusion: the five above PLUS 443392 (malignant
## neoplastic disease), so a control has no malignancy of any kind on record.
## PRIMARY_MALIGNANT is a subset of this by construction -- asserted below once
## the counts are back, since a violation would mean the concept sets have
## drifted apart.
ANY_MALIGNANT     <- c(PRIMARY_MALIGNANT, 443392L)

#' One criterion, as a pair of SQL CTEs: its expanded concept set, then the
#' people with an occurrence of any concept in that set.
#'
#' TWO CTES, AND A JOIN, NOT ONE CTE WITH AN `IN (SELECT ...)`. The obvious
#' spelling -- `WHERE drug_concept_id IN (SELECT descendant FROM <rollup>
#' WHERE ancestor IN (...))` -- is the spelling the generated cohort SQL uses,
#' and BigQuery rejects it here with
#'
#'   Correlated subqueries that reference other tables are not supported
#'   unless they can be de-correlated, such as by transforming them into an
#'   efficient JOIN.  [invalidQuery]
#'
#' The cohort SQL gets away with it because its subqueries sit inside a single
#' `person_id IN (SELECT id FROM T_ENT_person WHERE ...)` filter that the
#' planner can flatten. Nested two deep inside a WITH, as they are here, it
#' cannot. Hoisting the concept set to its own CTE and joining is the
#' transformation the error message asks for, and it changes nothing about
#' which rows come back.
#'
#' The ancestor is unioned in with its own descendants because an ancestor is
#' not its own descendant in the T_HAD_* tables -- dropping it would silently
#' lose everyone coded at the top-level concept. UNION DISTINCT rather than
#' UNION ALL so a concept listed as both cannot be matched twice; the join
#' below would otherwise duplicate person rows before DISTINCT collapses them.
#'
#' @param name      CTE name for the person set
#' @param table     T_ENT_* table holding the occurrences
#' @param id_col    the concept id column on that table
#' @param rollup    T_HAD_* table mapping ancestor -> descendant
#' @param ancestors integer concept ids to expand
criterion_ctes <- function(name, table, id_col, rollup, ancestors) {
  anc <- paste(ancestors, collapse = ", ")
  sprintf(
"  %s_concepts AS (
    SELECT descendant AS concept_id
    FROM `%s`.%s
    WHERE ancestor IN (%s)
    UNION DISTINCT
    SELECT concept_id FROM UNNEST([%s]) AS concept_id
  ),
  %s AS (
    SELECT DISTINCT o.person_id
    FROM `%s`.%s AS o
    JOIN %s_concepts AS c ON o.%s = c.concept_id
  )",
    name, INDEX, rollup, anc, anc,
    name, INDEX, table, name, id_col)
}

## ON, not USING, on every join below: `USING (person_id)` merges the two
## columns into one, after which `m.person_id IS NOT NULL` -- the whole point
## of the join -- is no longer a reference BigQuery will accept.
sql <- sprintf(
"WITH
  base AS (
    SELECT person_id FROM `%s`.person
  ),
  has_bmi AS (
    SELECT DISTINCT person_id
    FROM `%s`.T_ENT_measurementOccurrence
    WHERE measurement_concept_id = %d
  ),
%s,
%s,
%s,
%s,
  flags AS (
    SELECT
      b.person_id,
      m.person_id  IS NOT NULL AS has_bmi,
      g.person_id  IS NOT NULL AS glp1,
      s.person_id  IS NOT NULL AS bariatric,
      pm.person_id IS NOT NULL AS primary_malignant,
      am.person_id IS NOT NULL AS any_malignant
    FROM base b
    LEFT JOIN has_bmi           m  ON  m.person_id = b.person_id
    LEFT JOIN glp1              g  ON  g.person_id = b.person_id
    LEFT JOIN bariatric         s  ON  s.person_id = b.person_id
    LEFT JOIN primary_malignant pm ON pm.person_id = b.person_id
    LEFT JOIN any_malignant     am ON am.person_id = b.person_id
  )
SELECT has_bmi, glp1, bariatric, primary_malignant, any_malignant,
       COUNT(*) AS n
FROM flags
GROUP BY 1, 2, 3, 4, 5
ORDER BY 1, 2, 3, 4, 5",
  CDR, INDEX, BMI_CONCEPT,
  criterion_ctes("glp1", "T_ENT_ingredientOccurrence", "drug_concept_id",
                 "T_HAD_ingredientConcept_default", GLP1_ANCESTOR),
  criterion_ctes("bariatric", "T_ENT_procedureOccurrence", "procedure_concept_id",
                 "T_HAD_procedureConcept_default", BARIATRIC_ANCESTOR),
  criterion_ctes("primary_malignant", "T_ENT_conditionOccurrence", "condition_concept_id",
                 "T_HAD_conditionConcept_default", PRIMARY_MALIGNANT),
  criterion_ctes("any_malignant", "T_ENT_conditionOccurrence", "condition_concept_id",
                 "T_HAD_conditionConcept_default", ANY_MALIGNANT))

if (PRINT_SQL_ONLY) { cat(sql, "\n"); quit(status = 0) }


## --- run -------------------------------------------------------------------
library(bigrquery)

message("billing : ", BILLING)
message("cdr     : ", CDR)
message("index   : ", INDEX)
t0 <- Sys.time()
cells <- bq_table_download(bq_project_query(BILLING, sql))
message("queried in ", round(difftime(Sys.time(), t0, units = "secs")), "s; ",
        nrow(cells), " occupied cells of 32")

## A person with a primary malignant neoplasm must also satisfy the broader
## any-malignant set, since ANY_MALIGNANT contains PRIMARY_MALIGNANT. If this
## trips, the two concept sets have drifted and every count below is suspect.
stopifnot(
  "primary_malignant is not a subset of any_malignant -- concept sets drifted" =
    sum(cells$n[cells$primary_malignant & !cells$any_malignant]) == 0)


## --- the cascade -----------------------------------------------------------
## Everything below is arithmetic on `cells`. No further querying, and no
## non-standard evaluation: the five flags are pulled out as plain logical
## vectors over the cross-tabulation's rows, so each population is an ordinary
## boolean expression and N() sums the people in it.

bmi  <- cells$has_bmi
glp1 <- cells$glp1
bar  <- cells$bariatric
pm   <- cells$primary_malignant
am   <- cells$any_malignant

N <- function(sel) sum(cells$n[sel])

N_BASE <- sum(cells$n)

## THE TRUNK. Both arms branch from here, so neither is "after" the other.
trunk <- bmi & !glp1 & !bar

## Each row: the criterion, and the population still standing when it is
## applied. n_marginal ignores that population -- it is the criterion against
## the whole cohort, which is what the Cohort Builder reports and what the
## panel drew. Order is Figure 1A's; see the header.
cascade <- bind_rows(
  tibble(step = "All of Us participants",
         n_before = NA_integer_, n_marginal = NA_integer_,
         n_excluded_here = NA_integer_, n_after = N_BASE),

  tibble(step = "Excluded: no EHR BMI",
         n_before        = N_BASE,
         n_marginal      = N(!bmi),
         n_excluded_here = N(!bmi),
         n_after         = N(bmi)),

  tibble(step = "Excluded: received GLP-1 agonists",
         n_before        = N(bmi),
         n_marginal      = N(glp1),
         n_excluded_here = N(bmi & glp1),
         n_after         = N(bmi & !glp1)),

  tibble(step = "Excluded: received bariatric surgery",
         n_before        = N(bmi & !glp1),
         n_marginal      = N(bar),
         n_excluded_here = N(bmi & !glp1 & bar),
         n_after         = N(trunk)),

  tibble(step = "Excluded: any malignant neoplasm (control arm)",
         n_before        = N(trunk),
         n_marginal      = N(am),
         n_excluded_here = N(trunk & am),
         n_after         = N(trunk & !am)),

  tibble(step = "Included: primary malignant neoplasm (case arm)",
         n_before        = N(trunk),
         n_marginal      = N(pm),
         ## An inclusion, not an exclusion: what it removes is everyone
         ## WITHOUT one, which is why n_marginal and n_excluded_here are not
         ## comparable on this row the way they are on the others.
         n_excluded_here = N(trunk & !pm),
         n_after         = N(trunk & pm))
) %>%
  mutate(overlap_absorbed = n_marginal - n_excluded_here)

## The conditional column composes. Assert it rather than trusting the eye:
## this is the entire claim the script exists to make. The two arms are
## excluded from the chain because they branch rather than follow.
chain <- cascade %>% filter(!grepl("arm\\)$", step))
stopifnot(
  "conditional counts do not compose -- the cascade is wrong" =
    all(chain$n_after[-1] == head(chain$n_after, -1) - chain$n_excluded_here[-1]))


## --- against the panel -----------------------------------------------------
## Transcribed from the Figure 1A artwork. Update these if the panel is
## redrawn, and nothing else -- the computation above is independent of them.
##
## The control endpoint is drawn TWICE on the panel with DIFFERENT values:
## 488,528 bottom-left and 488,527 in both per-cancer-type control pools. Both
## are carried below; at most one can be right.

PUBLISHED <- tibble::tribble(
  ~step,                                            ~published,
  "All of Us participants",                             747028,
  "Excluded: no EHR BMI",                               144744,
  "Excluded: received GLP-1 agonists",                   34856,
  "Excluded: received bariatric surgery",                  298,
  "Excluded: any malignant neoplasm (control arm)",       88376,
  "Included: primary malignant neoplasm (case arm)",      63141
)
PUBLISHED_CONTROLS <- c(488528L, 488527L)

check <- cascade %>%
  left_join(PUBLISHED, by = "step") %>%
  mutate(
    reading = dplyr::case_when(
      step == "All of Us participants" ~
        ifelse(published == n_after, "base MATCHES", "base MISMATCH"),
      is.na(published)                              ~ "no published value",
      published == n_marginal & published == n_excluded_here ~
        "marginal = conditional (criterion does not overlap the others)",
      published == n_marginal                       ~ "panel drew the MARGINAL count",
      published == n_excluded_here                  ~ "panel drew the CONDITIONAL count",
      TRUE                                          ~ "MATCHES NEITHER -- investigate"
    )
  ) %>%
  select(step, published, n_marginal, n_excluded_here, overlap_absorbed,
         n_before, n_after, reading)

cat("\n")
print(as.data.frame(check), row.names = FALSE)

n_controls <- N(trunk & !am)
n_cases_in <- N(trunk & pm)

message(sprintf(
"\nEndpoints -- these are set intersections and do NOT depend on the order:
  controls entering PSM                    %s  (panel: %s)
  cases with a primary malignant neoplasm  %s  (panel: 63,141)
     -- the case arm's '>= 1 pre-diagnosis BMI' box is NOT computed here.",
  scales::comma(n_controls),
  paste(scales::comma(PUBLISHED_CONTROLS), collapse = " / "),
  scales::comma(n_cases_in)))

## The gap the panel's arithmetic opens, stated as a number rather than left
## for a reader to discover by adding up the boxes.
drawn_chain <- PUBLISHED$published[PUBLISHED$step == "All of Us participants"] -
  sum(PUBLISHED$published[grepl("^Excluded", PUBLISHED$step)])
message(sprintf(
"\nSubtracting the panel's own boxes in the panel's own order gives %s, against
%s computed. The difference of %s is people counted in more than one box.",
  scales::comma(drawn_chain), scales::comma(n_controls),
  scales::comma(abs(n_controls - drawn_chain))))


## --- where the double-counting sits ----------------------------------------
## So the gap is attributable rather than merely acknowledged. `both` is the
## CDR-wide intersection; `both_with_bmi` is the same within the population
## that survives the first box, which is the one that matters for the cascade.

## The output columns are criterion_a / criterion_b, NOT a / b. tibble()
## evaluates its arguments in order and lets later ones see earlier ones, so
## `tibble(a = label_a, both = N(a & b))` binds `a` to the character column
## just created rather than to the logical argument, and dies with
## "invalid 'type' (character) of argument" -- after everything above has
## printed and before anything is written.
overlap_pair <- function(label_a, a, label_b, b) {
  tibble(criterion_a = label_a, criterion_b = label_b,
         both = N(a & b), both_with_bmi = N(bmi & a & b))
}

pairs <- bind_rows(
  overlap_pair("GLP-1",     glp1, "any malignant neoplasm",     am),
  overlap_pair("bariatric", bar,  "any malignant neoplasm",     am),
  overlap_pair("GLP-1",     glp1, "bariatric",                  bar),
  overlap_pair("GLP-1",     glp1, "primary malignant neoplasm", pm),
  overlap_pair("bariatric", bar,  "primary malignant neoplasm", pm)
)

cat("\nPairwise overlaps between the exclusion criteria:\n")
print(as.data.frame(pairs), row.names = FALSE)


## --- write -----------------------------------------------------------------
out_dir <- file.path(WORKSPACE,
  "pre_diagnosis_bmi_analysis/code_outputs_scripts/results/tables/cohort_flow")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_table(cells, file.path(out_dir, "fig1a_criteria_crosstab"))
write_table(check, file.path(out_dir, "fig1a_exclusion_cascade"))

message("\nWrote fig1a_criteria_crosstab and fig1a_exclusion_cascade to ", out_dir)
message("Put the n_excluded_here column on the figure, or keep the panel's ",
        "marginal counts and say in the caption that the categories overlap ",
        "and do not sum.")
