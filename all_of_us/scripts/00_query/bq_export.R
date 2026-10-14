## ---------------------------------------------------------------------------
## Stage 00.0 -- the BigQuery -> Cloud Storage -> workspace transport.
##
## Produces no manuscript figure or table: shared transport, sourced by every
## domain script under scripts/00_query/<cohort>/.
##
## ---------------------------------------------------------------------------

library(bigrquery)

stopifnot("bq_export.R needs REPO, set in the caller's locations block" =
            exists("REPO") && nzchar(REPO))

## Data root: two levels above the repo. Override with $AOU_WORKSPACE.
##
## mustWork = FALSE, unlike the other stages: the query half of this script
## needs only BigQuery and the bucket, so an export can be run from a machine
## that has no workspace tree. The download half creates what it needs.
AOU_WORKSPACE <- normalizePath(Sys.getenv(
  "AOU_WORKSPACE", unset = file.path(REPO, "..", "..", "workspace")),
  mustWork = FALSE)


## --- command line ----------------------------------------------------------
## Parsed once, at source time, so every domain script takes the same flags.

.usage <- function() {
  cat(
"Usage: Rscript scripts/00_query/<cohort>/<domain>.R [options]

  --no-query        skip BigQuery; just download what is already in the bucket
  --no-download     export to Cloud Storage only; leave the workspace alone
  --date=YYYYMMDD   use this stamp instead of today's (see below)
  --bucket=gs://... export somewhere other than the cohort's own bucket
  --dry-run         print what would run, touch nothing
  -h, --help        this

Environment: $EXPORT_BUCKET, $AOU_BILLING_PROJECT, $AOU_EXPORT_DATE,
             $AOU_WORKSPACE. Flags win over all of them.

The date stamp is part of the output path, so exports on different days sit
side by side and a re-export on the same day overwrites. --date is how you
re-download an older export, or add a domain to one already under way without
splitting it across two stamps.
")
}

.args <- commandArgs(trailingOnly = TRUE)
if (any(.args %in% c("-h", "--help"))) { .usage(); quit(status = 0) }

.flag <- function(name) name %in% .args
.opt  <- function(name, default = "") {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit[length(hit)])
}

.unknown <- setdiff(
  grep("^-", .args, value = TRUE),
  c("--no-query", "--no-download", "--dry-run", "-h", "--help",
    grep("^--(date|bucket)=", .args, value = TRUE)))
if (length(.unknown)) {
  message("unknown argument(s): ", paste(.unknown, collapse = " "))
  .usage()
  quit(status = 2)
}

DO_QUERY    <- !.flag("--no-query")
DO_DOWNLOAD <- !.flag("--no-download")
DRY_RUN     <-  .flag("--dry-run")

if (!DO_QUERY && !DO_DOWNLOAD) {
  stop("--no-query and --no-download together leave nothing to do", call. = FALSE)
}

## Today, unless pinned.
EXPORT_DATE <- .opt("date", Sys.getenv("AOU_EXPORT_DATE",
                                       unset = format(Sys.Date(), "%Y%m%d")))
if (!grepl("^[0-9]{8}$", EXPORT_DATE)) {
  stop("--date/$AOU_EXPORT_DATE must be YYYYMMDD, got: ", EXPORT_DATE, call. = FALSE)
}


## ---------------------------------------------------------------------------
## export_bq_domain(domain, sql, bucket, billing, query_dataset, local_subdir)
##
## Runs one domain's query, extracts the result to
##
##   <bucket>/bq_exports/<domain>/<EXPORT_DATE>/<domain>/<domain>_*.csv
##
## -- the layout the export writes and the layout scripts/00_ingest/ reads --
## and then mirrors those shards to the same relative path under
## <workspace>/<local_subdir>/, which is where stage 00.1 and 00.3 look.
##
## BigQuery shards the extract, so the number of CSV parts is not fixed and
## only the first carries a header. Everything downstream globs the directory
## rather than naming a part, with the one exception noted in
## scripts/00_ingest/all_pulled_cancers.R, which names part 000000000000.
## ---------------------------------------------------------------------------
export_bq_domain <- function(domain, sql, bucket, billing, query_dataset,
                             local_subdir) {

  bucket  <- .opt("bucket", Sys.getenv("EXPORT_BUCKET", unset = bucket))
  billing <- Sys.getenv("AOU_BILLING_PROJECT", unset = billing)
  if (!grepl("^gs://", bucket)) {
    stop("bucket must start with gs://, got: ", bucket, call. = FALSE)
  }

  ## Repeating <domain> is not a slip: the dated directory holds a directory
  ## per export rather than loose shards.
  gcs_dir  <- paste(bucket, "bq_exports", domain, EXPORT_DATE, domain, sep = "/")
  gcs_glob <- paste0(gcs_dir, "/", domain, "_*.csv")

  message("domain    : ", domain)
  message("date      : ", EXPORT_DATE)
  message("bucket    : ", gcs_glob)
  message("billing   : ", billing)
  message("dataset   : ", query_dataset)

  local_dir <- file.path(AOU_WORKSPACE, local_subdir, "bq_exports",
                         domain, EXPORT_DATE, domain)
  if (DO_DOWNLOAD) message("workspace : ", local_dir)

  ## Checked before the query rather than after it: an extract of
  ## measurementOccurrence is minutes of billed BigQuery, and finding out only
  ## then that the shards cannot be fetched is a bad trade for one Sys.which().
  if (DO_DOWNLOAD && !DRY_RUN && Sys.which("gsutil") == "") {
    stop("gsutil is not on $PATH, so the shards could be extracted but not\n",
         "  mirrored into the workspace. Re-run with --no-download to export to\n",
         "  the bucket only, or put gsutil on $PATH.", call. = FALSE)
  }

  if (DRY_RUN) {
    if (DO_QUERY)    message("[dry-run] would query and extract to ", gcs_glob)
    if (DO_DOWNLOAD) message("[dry-run] would copy ", gcs_glob, " -> ", local_dir)
    return(invisible(NULL))
  }

  ## --- query and extract ---------------------------------------------------
  ## bq_dataset_query() lands the result in a temporary table in query_dataset,
  ## and bq_table_save() extracts it.
  ## Re-running on the same date overwrites that date's shards.
  if (DO_QUERY) {
    t0 <- Sys.time()
    bq_table_save(
      bq_dataset_query(query_dataset, sql, billing = billing),
      gcs_glob,
      destination_format = "CSV")
    message("extracted in ", round(difftime(Sys.time(), t0, units = "secs")), "s")
  } else {
    message("--no-query: using whatever is already at ", gcs_glob)
  }

  ## --- mirror into the workspace -------------------------------------------
  ## Nothing downstream wants the shards in this session's memory -- stage 00.1
  ## and 00.3 read them off disk -- so they are copied down instead.
  if (DO_DOWNLOAD) {
    dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
    status <- system2("gsutil", c("-m", "cp", shQuote(gcs_glob), shQuote(local_dir)))
    if (status != 0) {
      stop("gsutil cp exited ", status, " -- nothing was mirrored to ", local_dir,
           ".\n  The extract itself succeeded, so the shards are in the bucket:\n  ",
           gcs_glob, "\n  Re-run with --no-query to fetch them without querying again.",
           call. = FALSE)
    }

    shards <- list.files(local_dir, pattern = paste0("^", domain, "_.*\\.csv$"),
                         full.names = TRUE)
    message("mirrored ", length(shards), " shard(s), ",
            round(sum(file.size(shards)) / 2^20, 1), " MiB -> ", local_dir)
  }

  invisible(NULL)
}
