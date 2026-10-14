#!/usr/bin/env bash
#
# Produces no manuscript figure or table itself: runs the BigQuery exports
# the cohorts are built from.
# Export every domain of both cohorts out of BigQuery.
#
#   ./scripts/00_query/run_queries.sh             # both cohorts, every domain
#   ./scripts/00_query/run_queries.sh cancer      # the cases arm only
#   ./scripts/00_query/run_queries.sh control     # the control pool only
#   ./scripts/00_query/run_queries.sh --list      # show the domains and exit
#   ./scripts/00_query/run_queries.sh --dry-run   # print what would run
#
# This is stage 00.0 -- it produces the CSVs that scripts/00_ingest/ reads, and
# it is the only stage that talks to BigQuery for bulk data. It is NOT in
# run_all.sh: it costs query and storage, it needs Workbench credentials, and
# the exports change only when the cohort definition does. Run it by hand when
# you have re-exported the dataset from the Workbench, then run_all.sh pre.
#
# ONE DATE STAMP FOR THE WHOLE RUN
#
#   Each domain lands in a directory stamped with the export date, so exports
#   from different days sit side by side. The stamp is computed ONCE here and
#   passed to every script, so a run that straddles midnight -- these are big
#   extracts -- cannot leave half the domains under one date and half under the
#   next. Pin it yourself with --date=YYYYMMDD to add a domain to an export
#   already on disk, or to re-download one without re-querying.
#
# FLAGS PASSED THROUGH TO EACH DOMAIN SCRIPT
#
#   --no-query      download what is already in the bucket; run no SQL
#   --no-download   extract to Cloud Storage only; leave the workspace alone
#   --dry-run       print what each script would do, touch nothing
#   --date=YYYYMMDD use this stamp
#
# Domains run one at a time. They are independent BigQuery jobs and would
# parallelise, but each extract is large and a failure part-way through a
# parallel run is far harder to resume from than re-running one named script.
#
# Interpreter: override with $RSCRIPT_BIN. Data root: $AOU_WORKSPACE.

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

# ---------------------------------------------------------------------------
# The domains, per cohort. One entry per generated script; regenerate those
# Add a new domain's script to the list here.
# ---------------------------------------------------------------------------
CANCER=(
  "scripts/00_query/cancer_patients/person.R"
  "scripts/00_query/cancer_patients/measurement_occurrence.R"
  "scripts/00_query/cancer_patients/condition_occurrence.R"
  "scripts/00_query/cancer_patients/procedure_occurrence.R"
  "scripts/00_query/cancer_patients/ingredient_occurrence.R"
  "scripts/00_query/cancer_patients/zipcode_socioeconomic.R"
)

CONTROL=(
  "scripts/00_query/control_individuals/person.R"
  "scripts/00_query/control_individuals/measurement_occurrence.R"
  "scripts/00_query/control_individuals/condition_occurrence.R"
  "scripts/00_query/control_individuals/procedure_occurrence.R"
  "scripts/00_query/control_individuals/ingredient_occurrence.R"
  "scripts/00_query/control_individuals/zipcode_socioeconomic.R"
)

# What the rest of the pipeline does with each domain, for --list.
label_for() {
  case "$1" in
    */person.R)                   echo "demographics -- PSM covariates (stage 01)" ;;
    */measurement_occurrence.R)   echo "BMI measurements -- every downstream analysis" ;;
    */condition_occurrence.R)     echo "diagnoses -- first_cancer_diagnosis_dates (stage 00.3)" ;;
    */procedure_occurrence.R)     echo "procedures -- cohort definition only" ;;
    */ingredient_occurrence.R)    echo "drug exposures -- cohort definition only" ;;
    */zipcode_socioeconomic.R)    echo "deprivation_index -- PSM covariate (unused, see README)" ;;
    *)                            echo "" ;;
  esac
}

# ---------------------------------------------------------------------------
COHORT="all"; LIST=0; KEEP_GOING=0; PASSTHRU=(); DATE_OPT=""
for arg in "$@"; do
  case "$arg" in
    cancer|cancers|cases)  COHORT="cancer" ;;
    control|controls)      COHORT="control" ;;
    all|both)              COHORT="all" ;;
    --list)                LIST=1 ;;
    --keep-going)          KEEP_GOING=1 ;;
    --date=*)              DATE_OPT="${arg#--date=}" ;;
    --no-query|--no-download|--dry-run) PASSTHRU+=("$arg") ;;
    -h|--help)             sed -n '3,40p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown argument: $arg (try --help)" >&2; exit 2 ;;
  esac
done

case "$COHORT" in
  cancer)  STEPS=("${CANCER[@]}") ;;
  control) STEPS=("${CONTROL[@]}") ;;
  all)     STEPS=("${CANCER[@]}" "${CONTROL[@]}") ;;
esac

# The lists above are hand-ordered -- person and the BMI measurements first,
# the exploratory sleep domains last -- rather than globbed, so a long run
# produces what the pipeline needs soonest. That ordering is worth keeping, but
# it means a newly added domain script would be skipped silently.
# So: say so, loudly, and keep going.
for f in "$REPO"/scripts/00_query/*/*.R; do
  rel="${f#$REPO/}"
  [[ "$rel" == "scripts/00_query/bq_export.R" ]] && continue
  listed=0
  for s in "${CANCER[@]}" "${CONTROL[@]}"; do
    [[ "$s" == "$rel" ]] && { listed=1; break; }
  done
  if [[ $listed -eq 0 ]]; then
    echo "WARNING: $rel exists but is not in this script's CANCER/CONTROL list," >&2
    echo "         so it will NOT be exported. Add it." >&2
  fi
done

if [[ $LIST -eq 1 ]]; then
  printf '%-3s %-60s %s\n' "#" "script" "feeds"
  i=0
  for s in "${STEPS[@]}"; do
    i=$((i + 1)); printf '%-3s %-60s %s\n' "$i" "$s" "$(label_for "$s")"
  done
  exit 0
fi

# ---------------------------------------------------------------------------
# Interpreter: the Workbench location, falling back to whatever is on PATH.
# ---------------------------------------------------------------------------
if [[ -z "${RSCRIPT_BIN:-}" ]]; then
  for c in /opt/conda/envs/packages/bin/Rscript "$(command -v Rscript || true)"; do
    [[ -x "${c:-}" ]] && { RSCRIPT_BIN="$c"; break; }
  done
fi
[[ -x "${RSCRIPT_BIN:-}" ]] || { echo "no Rscript found; set \$RSCRIPT_BIN" >&2; exit 2; }
export RSCRIPT_BIN

# One stamp for every domain in this run -- see the header.
EXPORT_DATE="${DATE_OPT:-${AOU_EXPORT_DATE:-$(date +%Y%m%d)}}"
if [[ ! "$EXPORT_DATE" =~ ^[0-9]{8}$ ]]; then
  echo "--date must be YYYYMMDD, got: $EXPORT_DATE" >&2; exit 2
fi
export AOU_EXPORT_DATE="$EXPORT_DATE"

# Logs go beside run_all.sh's, when there is a workspace to put them in. With
# --no-download there need not be one, so a missing tree streams to the console
# rather than stopping the run.
WORKSPACE="${AOU_WORKSPACE:-$REPO/../../workspace}"
LOG_DIR=""
if [[ -d "$WORKSPACE" ]]; then
  WORKSPACE="$(cd "$WORKSPACE" && pwd)"
  LOG_DIR="$WORKSPACE/pre_diagnosis_bmi_analysis/code_outputs_scripts/logs"
  mkdir -p "$LOG_DIR" 2>/dev/null || LOG_DIR=""
fi

echo "repo      : $REPO"
echo "workspace : ${WORKSPACE:-<none>}"
echo "Rscript   : $RSCRIPT_BIN"
echo "cohort    : $COHORT (${#STEPS[@]} domains)"
echo "date      : $EXPORT_DATE"
echo "options   : ${PASSTHRU[*]:-<none>}"
if [[ -n "$LOG_DIR" ]]; then echo "logs      : $LOG_DIR"
else echo "logs      : <none; output goes to the console>"; fi
echo

# ---------------------------------------------------------------------------
started=$(date +%s)
i=0
FAILED=()
for s in "${STEPS[@]}"; do
  i=$((i + 1))
  name="$(basename "$s" .R)"
  cohort_dir="$(basename "$(dirname "$s")")"

  printf '[%2d/%2d] %-56s ' "$i" "${#STEPS[@]}" "${cohort_dir}/${name}"
  t0=$(date +%s)

  if [[ -n "$LOG_DIR" ]]; then
    log="$LOG_DIR/00_query_${cohort_dir}_${name}_${EXPORT_DATE}.log"
    "$RSCRIPT_BIN" "$REPO/$s" ${PASSTHRU[@]+"${PASSTHRU[@]}"} > "$log" 2>&1
    rc=$?
  else
    echo
    "$RSCRIPT_BIN" "$REPO/$s" ${PASSTHRU[@]+"${PASSTHRU[@]}"}
    rc=$?
    log=""
  fi

  if [[ $rc -eq 0 ]]; then
    printf 'ok    %5ds\n' "$(( $(date +%s) - t0 ))"
  else
    printf 'FAIL  %5ds  (exit %d)\n' "$(( $(date +%s) - t0 ))" "$rc"
    FAILED+=("$s")
    if [[ $KEEP_GOING -eq 0 ]]; then
      if [[ -n "$log" ]]; then
        echo
        echo "--- last 25 lines of $log ---"
        tail -25 "$log"
        echo "-----------------------------"
      fi
      echo
      echo "stopped at domain $i. Fix, then re-run just that one:"
      echo "  AOU_EXPORT_DATE=$EXPORT_DATE $RSCRIPT_BIN $REPO/$s"
      echo "Or carry on past failures:  $0 $COHORT --keep-going"
      exit 1
    fi
    [[ -n "$log" ]] && echo "         see $log"
  fi
done

total=$(( $(date +%s) - started ))
echo
if [[ ${#FAILED[@]} -gt 0 ]]; then
  echo "$COHORT finished with ${#FAILED[@]} of ${#STEPS[@]} domains FAILED, in $((total / 60))m $((total % 60))s"
  echo
  for f in "${FAILED[@]}"; do echo "  FAILED  $f"; done
  echo
  echo "Re-run a failed domain on this run's stamp:  AOU_EXPORT_DATE=$EXPORT_DATE $RSCRIPT_BIN $REPO/<script>"
  exit 1
fi
echo "$COHORT complete: ${#STEPS[@]} domains in $((total / 60))m $((total % 60))s, stamped $EXPORT_DATE"
echo
echo "Next:  the stamp is part of the path, so point stage 00 at it --"
echo "       scripts/00_ingest/all_pulled_cancers.R names a date in its"
echo "       conditionOccurrence path. Then:  ./scripts/run_all.sh pre"
