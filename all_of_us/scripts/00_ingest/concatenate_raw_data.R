#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Stage 00.1 -- concatenate the per-theme BigQuery exports.
##
## Produces no manuscript figure or table: assembles the BigQuery exports into
## the tables every later stage reads.
##
## Reads:
##   <workspace>/cdrv9/control_individuals/aou_control_individuals_cohort_data/
##   <workspace>/cdrv9/cancer_patients/aou_cancer_patients_data/
## Writes:
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/control_individuals/raw_concatenated/<theme>.tsv
##   <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/raw_concatenated/<theme>.tsv
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


## ---------------------------------------------------------------------------
## concatenate_by_theme(input_dir, output_dir): lists files in input_dir,
## groups them by "theme" (the data_<theme> segment of the filename, e.g.
## sleepLevel, activitySummary, bmi...), concatenates all files sharing a theme
## (assumed to share the same columns), and writes one combined TSV per theme
## to output_dir.
##
## Expected filename pattern:
##   <timestamp>_<id>_data_<theme>-<part_number>.csv.gz
##   e.g. 20260709_015747_1617134116_data_sleepLevel-000000000000.csv.gz
##
## .ipynb files (and anything not matching the pattern above) are ignored.
## ---------------------------------------------------------------------------
concatenate_by_theme <- function(input_dir, output_dir = input_dir) {

  # Regex to identify data files and capture the theme name
  # Group 1 = theme (e.g. "sleepLevel")
  file_pattern     <- "^\\d{8}_\\d{6}_\\d+_data_(.+)-\\d+\\.csv\\.gz$"
  file_match_regex <- "\\d{8}_\\d{6}_\\d+_data_.+-\\d+\\.csv\\.gz$"

  if (is.null(output_dir) || output_dir == "") {
    output_dir <- input_dir
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  message("Scanning: ", normalizePath(input_dir))
  message("Output to: ", normalizePath(output_dir))

  all_files <- list.files(input_dir, full.names = FALSE)

  # Explicitly drop notebooks, then keep only files matching the data pattern
  candidate_files <- all_files[!str_detect(all_files, "\\.ipynb$")]
  candidate_files <- candidate_files[!str_detect(candidate_files, regex("sleep", ignore_case = TRUE))]
  data_files      <- candidate_files[str_detect(candidate_files, file_match_regex)]

  if (length(data_files) == 0) {
    warning("No files matching the expected pattern were found in: ", input_dir)
    return(invisible(NULL))
  }

  # Extract the theme for each file
  themes <- str_match(data_files, file_pattern)[, 2]

  file_df <- tibble(
    file  = data_files,
    theme = themes
  ) %>%
    filter(!is.na(theme))

  message(
    "Found ", nrow(file_df), " data files across ",
    n_distinct(file_df$theme), " themes:"
  )
  print(count(file_df, theme, name = "n_files"), n = Inf)

  # --------------------- Concatenate per theme --------------------------
  walk(unique(file_df$theme), function(th) {

    files_for_theme <- file_df %>%
      filter(theme == th) %>%
      pull(file) %>%
      file.path(input_dir, .)

    message("\n[", th, "] reading ", length(files_for_theme), " file(s)...")

    combined <- map_dfr(files_for_theme, function(f) {
      read_csv(f, col_types = cols(.default = col_character()),
               show_col_types = FALSE, progress = FALSE)
    })
    out_path <- file.path(output_dir, paste0(th, ".tsv"))
    write_tsv(combined, out_path)

    message("[", th, "] wrote ", nrow(combined), " rows -> ", out_path)
  })

  message("\nDone. Combined TSVs written to: ", normalizePath(output_dir))
  invisible(unique(file_df$theme))
}


## --- run -------------------------------------------------------------------

# concatenate control data
concatenate_by_theme(
  file.path(WORKSPACE, "cdrv9/control_individuals/aou_control_individuals_cohort_data"),
  file.path(WORKSPACE, "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/control_individuals/raw_concatenated")
)

# concatenate cancer data
concatenate_by_theme(
  file.path(WORKSPACE, "cdrv9/cancer_patients/aou_cancer_patients_data"),
  file.path(WORKSPACE, "pre_diagnosis_bmi_analysis/code_outputs_scripts/stable_data/cancer_patients/raw_concatenated")
)
