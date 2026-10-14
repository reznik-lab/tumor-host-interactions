#!/usr/bin/env bash
#
# Produces no manuscript figure or table itself: runs the scripts that do,
# in dependency order.
# Run the analysis scripts in dependency order.
#
#   ./scripts/run_all.sh            # everything, stopping at the first failure
#   ./scripts/run_all.sh pre        # stages 00-02: everything BEFORE the Python caller
#   ./scripts/run_all.sh python     # just the Python cachexia episode caller
#   ./scripts/run_all.sh post       # stages 03-06: everything AFTER it
#   ./scripts/run_all.sh --list     # show the steps and exit
#   ./scripts/run_all.sh --dry-run  # print the commands, run nothing
#
#   ./scripts/run_all.sh --keep-going
#       Don't stop at the first failure. The 03_cachexia chain cannot run as
#       written (see README.md "Known issues"), but 04/05/06 depend only on
#       02_matching/post_psm_filter_ctrls, so this still produces 5 of the 7
#       manuscript artifacts in one pass and lists what failed at the end.
#
# WHY THE SPLIT AT THE PYTHON STEP
#
#   calling_cachexia_episodes.py writes
#       df_episodes_all_precomp_WL10_<TODAY>.csv
#   but postprocess_cachexia_calls.R reads a hard-coded
#       df_episodes_all_precomp_WL10_20260929.csv
#
#   So `post` will not pick up a fresh `python` run until you either edit the
#   date in postprocess_cachexia_calls.R or rename the file. Run `pre`, then
#   `python`, sort the filename out, then `post`. See README.md.
#
# WHERE THE INPUTS COME FROM
#
#   Stage 00.0, scripts/00_query/, is NOT run from here: it exports the cohorts
#   out of BigQuery, which costs query and storage and needs Workbench
#   credentials, and its outputs change only when the cohort definition does.
#   Run it by hand after a fresh dataset export, then come back to `pre`:
#
#       ./scripts/00_query/run_queries.sh --list
#       ./scripts/00_query/run_queries.sh
#
# Each step logs to <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/logs/.
# Interpreters: override with $RSCRIPT_BIN / $PYTHON_BIN. Data root: $AOU_WORKSPACE.

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ---------------------------------------------------------------------------
# Steps, in dependency order. post_psm_filter_ctrls is the hub: everything
# after it reads its tables.
# ---------------------------------------------------------------------------
PRE=(
  "scripts/00_ingest/concatenate_raw_data.R"
  "scripts/00_ingest/pull_aou_enrollment_date.R"
  "scripts/00_ingest/all_pulled_cancers.R"
  "scripts/01_cohort/preprocess_cancer_data.R"
  "scripts/01_cohort/preprocess_control_data.R"
  "scripts/02_matching/psm.R"
  "scripts/02_matching/psm_quality_control.R"
  "scripts/02_matching/post_psm_filter_ctrls.R"
  "scripts/02_matching/supplementary_psm_quantification.R"
)

PYTHON_STEP=(
  "scripts/03_cachexia/calling_cachexia_episodes.py"
)

POST=(
  "scripts/03_cachexia/postprocess_cachexia_calls.R"
  "scripts/03_cachexia/cachexia_frequency.R"
  "scripts/03_cachexia/cachexia_frequency_by_threshold.R"
  "scripts/04_trajectories/analysis.R"
  "scripts/04_trajectories/setup_me_modelling_df.R"
  "scripts/04_trajectories/mixed_effects_modelling.R"
  "scripts/04_trajectories/supp_landmarks_all_types.R"
  "scripts/05_weight_class/weight_cat_dist_at_timepoints.R"
  "scripts/06_tables/overview_of_cancers_in_cohort.R"
  "scripts/06_tables/supplementary_table.R"
)

# Manuscript artifact produced, for --list.
label_for() {
  case "$1" in
    *cachexia_frequency_by_threshold.R) echo "Supp Fig 1c, Supp Table 5" ;;
    *cachexia_frequency.R) echo "Figure 1b" ;;
    *04_trajectories/analysis.R) echo "Figure 1c-f (trajectories)" ;;
    *mixed_effects_modelling.R) echo "Figure 1c-f (boxplots)" ;;
    *supp_landmarks_all_types.R) echo "Supp Table 6" ;;
    *weight_cat_dist_at_timepoints.R) echo "Supp Fig 1d, Supp Table 7" ;;
    *overview_of_cancers_in_cohort.R) echo "Supp Fig 1a, Supp Table 3" ;;
    *psm_quality_control.R) echo "Supp Fig 1b, Supp Table 4" ;;
    *supplementary_psm_quantification.R) echo "Supp Table 2" ;;
    *supplementary_table.R) echo "Supp Table 1" ;;
    *)                                echo "" ;;
  esac
}

# ---------------------------------------------------------------------------
PHASE="all"; DRY=0; LIST=0; KEEP_GOING=0
for arg in "$@"; do
  case "$arg" in
    pre|python|post|all) PHASE="$arg" ;;
    --list)              LIST=1 ;;
    --dry-run)           DRY=1 ;;
    --keep-going)        KEEP_GOING=1 ;;
    -h|--help)           sed -n '3,40p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown argument: $arg (try --help)" >&2; exit 2 ;;
  esac
done

case "$PHASE" in
  pre)    STEPS=("${PRE[@]}") ;;
  python) STEPS=("${PYTHON_STEP[@]}") ;;
  post)   STEPS=("${POST[@]}") ;;
  all)    STEPS=("${PRE[@]}" "${PYTHON_STEP[@]}" "${POST[@]}") ;;
esac

if [[ $LIST -eq 1 ]]; then
  printf '%-3s %-56s %s\n' "#" "script" "produces"
  i=0
  for s in "${STEPS[@]}"; do
    i=$((i + 1)); printf '%-3s %-56s %s\n' "$i" "$s" "$(label_for "$s")"
  done
  exit 0
fi

# ---------------------------------------------------------------------------
# Interpreters: the Workbench locations, falling back to whatever is on PATH.
# ---------------------------------------------------------------------------
if [[ -z "${RSCRIPT_BIN:-}" ]]; then
  for c in /opt/conda/envs/packages/bin/Rscript "$(command -v Rscript || true)"; do
    [[ -x "${c:-}" ]] && { RSCRIPT_BIN="$c"; break; }
  done
fi
if [[ -z "${PYTHON_BIN:-}" ]]; then
  for c in /opt/conda/envs/packages/bin/python3.14 "$(command -v python3 || true)"; do
    [[ -x "${c:-}" ]] && { PYTHON_BIN="$c"; break; }
  done
fi
[[ -x "${RSCRIPT_BIN:-}" ]] || { echo "no Rscript found; set \$RSCRIPT_BIN" >&2; exit 2; }
export RSCRIPT_BIN PYTHON_BIN

# Data root, derived exactly as the scripts derive it.
WORKSPACE="${AOU_WORKSPACE:-$REPO/../../workspace}"
if [[ ! -d "$WORKSPACE" ]]; then
  echo "workspace not found: $WORKSPACE" >&2
  echo "Set \$AOU_WORKSPACE to the directory holding cdrv9/ and pre_diagnosis_bmi_analysis/." >&2
  exit 2
fi
WORKSPACE="$(cd "$WORKSPACE" && pwd)"
LOG_DIR="$WORKSPACE/pre_diagnosis_bmi_analysis/code_outputs_scripts/logs"

if [[ $DRY -eq 0 ]] && ! mkdir -p "$LOG_DIR" 2>/dev/null; then
  echo "cannot create the log directory: $LOG_DIR" >&2
  exit 2
fi

echo "repo      : $REPO"
echo "workspace : $WORKSPACE"
echo "Rscript   : $RSCRIPT_BIN"
echo "python    : ${PYTHON_BIN:-<none>}"
echo "phase     : $PHASE (${#STEPS[@]} steps)"
[[ $DRY -eq 0 ]] && echo "logs      : $LOG_DIR"
echo

# ---------------------------------------------------------------------------
started=$(date +%s)
i=0
FAILED=()
for s in "${STEPS[@]}"; do
  i=$((i + 1))
  name="$(basename "$s")"; name="${name%.*}"
  stage="$(basename "$(dirname "$s")")"
  log="$LOG_DIR/$(printf '%02d' "$i")_${stage}_${name}.log"

  case "$s" in
    *.R)  cmd=("$RSCRIPT_BIN" "$REPO/$s") ;;   # absolute: scripts must not depend on cwd
    *.py) cmd=("$PYTHON_BIN"  "$REPO/$s") ;;
  esac

  if [[ $DRY -eq 1 ]]; then
    printf '[%2d/%2d] would run: %s\n' "$i" "${#STEPS[@]}" "${cmd[*]}"
    continue
  fi

  printf '[%2d/%2d] %-52s ' "$i" "${#STEPS[@]}" "$s"
  t0=$(date +%s)
  if "${cmd[@]}" > "$log" 2>&1; then
    printf 'ok    %5ds\n' "$(( $(date +%s) - t0 ))"
  else
    rc=$?
    printf 'FAIL  %5ds  (exit %d)\n' "$(( $(date +%s) - t0 ))" "$rc"
    FAILED+=("$s")
    if [[ $KEEP_GOING -eq 0 ]]; then
      echo
      echo "--- last 25 lines of $log ---"
      tail -25 "$log"
      echo "-----------------------------"
      echo
      echo "stopped at step $i. Fix, then re-run this phase:  ./scripts/run_all.sh $PHASE"
      echo "Or run the rest anyway:                          ./scripts/run_all.sh $PHASE --keep-going"
      exit 1
    fi
    echo "         see $log"
  fi
done

[[ $DRY -eq 1 ]] && exit 0

total=$(( $(date +%s) - started ))
echo
if [[ ${#FAILED[@]} -gt 0 ]]; then
  echo "$PHASE finished with ${#FAILED[@]} of ${#STEPS[@]} steps FAILED, in $((total / 60))m $((total % 60))s"
  echo
  for f in "${FAILED[@]}"; do echo "  FAILED  $f"; done
  echo
  echo "Logs are in $LOG_DIR"
  exit 1
fi
echo "$PHASE complete: ${#STEPS[@]} steps in $((total / 60))m $((total % 60))s"
if [[ "$PHASE" == "pre" ]]; then
  echo
  echo "Next:  ./scripts/run_all.sh python"
elif [[ "$PHASE" == "python" ]]; then
  echo
  echo "Next:  point postprocess_cachexia_calls.R at the file just written --"
  ls -1 "$WORKSPACE/pre_diagnosis_bmi_analysis/code_outputs_scripts/cachexia_episodes/processed_data"/df_episodes_all_precomp_*.csv 2>/dev/null | tail -3 | sed 's/^/           /'
  echo "       then:  ./scripts/run_all.sh post"
fi
