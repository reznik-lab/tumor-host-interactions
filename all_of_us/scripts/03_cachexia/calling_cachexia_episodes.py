#!/usr/bin/env python3
"""
Stage 03.1 -- smooth BMI trajectories and identify cachectic episodes.

Produces no manuscript figure or table: calls the weight-loss episodes that
Figure 1b and Supplementary Figure 1c are counted from.

The EWMA episode caller is adapted from a prior MSK codebase; its original
header follows.

    Created on Fri Feb  6 06:35:46 2026
    @author: castilv

The MSK driver code that accompanied these functions -- the MRN-keyed loaders,
the pre-edema QC pass and the edema filter -- does not run here and is not
reproduced. The function definitions are kept as-is, including the ones this
run does not call (see "unused in this run" notes).

Reads .../psm_matched_usable_data/solid_heme/full_cohort/all_{cancers,controls}_all_bmi_measurements_geq4_eligible.tsv
-- the paths post_psm_filter_ctrls.R writes.

SPEED (no change to what is computed)
The episode caller is vectorised and evaluates the three weight-loss thresholds
in a single pass over people; --jobs N spreads that pass across processes.
identify_cachexia_episodes_2 and merge_episodes below are untouched and remain
the definition of an episode -- call_episodes_fast() reproduces them exactly,
edge cases included (see its comments). --selfcheck N re-runs those original
functions on the first N people and asserts the frames match, so the
equivalence can be confirmed against the real cohort at any time:

    python3 scripts/03_cachexia/calling_cachexia_episodes.py --jobs 4 --selfcheck 500

!! DATE-STAMPED OUTPUT !!
Output filenames embed today's date (STAMP). postprocess_cachexia_calls.R reads
a hard-coded df_episodes_all_precomp_WL10_20260929.csv, so a fresh run here will
not be picked up by that script without editing the date.

Writes:
  <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/cachexia_episodes/processed_data/smoothed_bmi_all_patients_<STAMP>_alpha0.2.csv
  <workspace>/pre_diagnosis_bmi_analysis/code_outputs_scripts/cachexia_episodes/processed_data/df_episodes_all_precomp_<WL label>_<STAMP>.csv
"""

import argparse
import os
from datetime import datetime
from multiprocessing import get_context

import numpy as np
import pandas as pd
from tqdm import tqdm


# --- options ----------------------------------------------------------------
# --jobs only splits the episode pass across processes; it does not alter what
# is computed. It matches PSM_CORES in psm.R: unset means one worker per core
# bar one, --jobs 1 is the serial run, and $AOU_N_JOBS sets it without editing
# run_all.sh's bare `python3 <script>` invocation.
_parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
_parser.add_argument("--jobs", type=int, default=None, metavar="N",
                     help="processes for the episode pass. Default: one per "
                          "core bar one ($AOU_N_JOBS overrides). 1 = serial.")
_parser.add_argument("--selfcheck", type=int, default=0, metavar="N",
                     help="before the real run, re-run the original "
                          "identify_cachexia_episodes_2 on the first N people "
                          "and assert the rewrite matches it exactly.")
ARGS = _parser.parse_args()
if ARGS.jobs is None:
    _env = os.environ.get("AOU_N_JOBS", "").strip()
    ARGS.jobs = int(_env) if _env.isdigit() and int(_env) > 0 else (os.cpu_count() or 2) - 1
ARGS.jobs = max(1, ARGS.jobs)


# --- locations --------------------------------------------------------------
# Found from this script's own path, so it runs from any directory.
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

# Data root: two levels above the repo. Override with $AOU_WORKSPACE.
WORKSPACE = os.path.abspath(
    os.environ.get("AOU_WORKSPACE") or os.path.join(REPO, "..", "..", "workspace")
)

ANALYSIS = os.path.join(WORKSPACE, "pre_diagnosis_bmi_analysis")


# ---------------- Functions ------------------------------------------------
# ------------- Smoothing ----------------------
def smooth_bmi_ewma(df, smooth_col, alpha):
    # Unused in this run: the driver below does the same EWMA via
    # groupby().transform(), which survives the pandas 2.2 change to how
    # groupby().apply() handles the grouping column. Kept for reference.
    # Smooth the BMI values with an exponentially weighted moving average set to 0.2
    df['smoothed_BMI'] = df[smooth_col].ewm(alpha=alpha, ignore_na=True).mean()
    return df


# ------------- Recovery Episode Identification ----------------------
# Unused in this run: identify_cachexia_episodes_2 is called with recovery=False.
def identify_recovery_episodes(patient_data, merged_episodes_df, time_col, bmi_col):
    # Initialize new columns for recovery day, date, and BMI
    merged_episodes_df['recovery_day'] = None
    merged_episodes_df['recovery_date'] = None
    merged_episodes_df['recovery_bmi'] = None

    merged_episodes_df['recovery_smoothed_bmi'] = None
    if not merged_episodes_df['start_day'].isna().all():
        for i in range(merged_episodes_df.shape[0]):
            t0 = merged_episodes_df['end_day'][i]
            if i < merged_episodes_df.shape[0] - 1:
                t1 = merged_episodes_df['start_day'][i + 1]
            else:
                t1 = patient_data[time_col].max()

            # Get the BMI at the end of the episode
            end_index = patient_data.loc[patient_data[time_col] == t0].index[0]
            b0 = patient_data[bmi_col][end_index]

            # Detect recovery: first time BMI exceeds the episode end BMI by 5%
            recovery_data = patient_data[(patient_data[time_col] >= t0) & (patient_data[time_col] <= t1)
                                         & (patient_data[bmi_col] > b0 + np.log(1.05))]

            if not recovery_data.empty:
                recovery_day = recovery_data[time_col].min()
                recovery_date = recovery_data.loc[recovery_data[time_col] == recovery_day, 'datetime'].values[0]
                recovery_bmi = recovery_data.loc[recovery_data[time_col] == recovery_day, bmi_col].values[0]
                recovery_smoothed_bmi = recovery_data.loc[recovery_data[time_col] == recovery_day, 'smoothed_BMI'].values[0]

                merged_episodes_df.at[i, 'recovery_day'] = recovery_day
                merged_episodes_df.at[i, 'recovery_date'] = recovery_date
                merged_episodes_df.at[i, 'recovery_bmi'] = recovery_bmi
                merged_episodes_df.at[i, 'recovery_smoothed_bmi'] = recovery_smoothed_bmi

    return merged_episodes_df


# ------------- Cachexia Episode Identification ----------------------
# Unused in this run: identify_cachexia_episodes_2, below, is the one called.
def identify_cachexia_episodes(df, time_col, bmi_col, recovery=True, wl_frac=0.05):
    '''
    Detect cachexia episodes, including start, onset, end, and recovery information.
    wl_frac: 0.05 (5%), 0.10 (10%), 0.15 (15%)
    NOTE: keeps the 180-day window exactly as in your current code.
    '''
    results = []

    df.reset_index(drop=True, inplace=True)

    wl_log = np.log(1 - wl_frac)

    for i in range(df.shape[0]):
        t0 = df[time_col][i]            # Start day of observation
        b0 = df[bmi_col][i]             # Baseline log-smoothed BMI on start day
        b0_raw = df['smoothed_BMI'][i]
        current_date = df['datetime'][i]
        current_bmi = df['bmi'][i]

        # KEEP 180 DAYS exactly
        df_window = df[(df[time_col] > t0) & (df[time_col] <= t0 + 180)]

        # ONLY CHANGE: threshold is parameterized
        df_filtered = df_window[df_window[bmi_col] < b0 + wl_log]

        if not df_filtered.empty:
            onset_day = df_filtered[time_col].min()
            onset_date = df_filtered.loc[df_filtered[time_col] == onset_day, 'datetime'].values[0] if not df_filtered.loc[df_filtered[time_col] == onset_day].empty else None
            onset_bmi = df_filtered.loc[df_filtered[time_col] == onset_day, 'smoothed_BMI'].values[0] if not df_filtered.loc[df_filtered[time_col] == onset_day].empty else None

            max_t = df_filtered[time_col].max()
            end_date = df_filtered.loc[df_filtered[time_col] == max_t, 'datetime'].values[0] if not df_filtered.loc[df_filtered[time_col] == max_t].empty else None
            end_bmi = df_filtered.loc[df_filtered[time_col] == max_t, 'smoothed_BMI'].values[0] if not df_filtered.loc[df_filtered[time_col] == max_t].empty else None

            results.append({
                'start_day': t0,
                'start_date': current_date,
                'start_bmi': current_bmi,
                'onset_day': onset_day,
                'onset_date': onset_date,
                'onset_bmi': onset_bmi,
                'end_day': max_t,
                'end_date': end_date,
                'end_bmi': end_bmi
            })

    if len(results) == 0:
        results.append({
            'start_day': None,
            'start_date': None,
            'start_bmi': None,
            'onset_day': None,
            'onset_date': None,
            'onset_bmi': None,
            'end_day': None,
            'end_date': None,
            'end_bmi': None
        })

    merged_episodes_df = merge_episodes(pd.DataFrame(results), start_col='start_day', end_col='end_day')

    if recovery:
        merged_episodes_df = identify_recovery_episodes(df, merged_episodes_df, time_col, bmi_col)

    return merged_episodes_df


# ------------- Cachexia Episode Identification 5% , 2% if BMI<20----------------------
def identify_cachexia_episodes_2(df, time_col, bmi_col, recovery=True, wl_frac=0.05, bmi20_rule=False):
    '''
    Detect cachexia episodes, including start, onset, end, and recovery information.
    wl_frac: 0.05 (5%), 0.10 (10%), 0.15 (15%)
    bmi20_rule: if True, use 2% WL threshold when baseline BMI < 20 (else wl_frac).
    NOTE: keeps the 180-day window exactly as in your current code.
    '''
    results = []
    df.reset_index(drop=True, inplace=True)

    wl_log_default = np.log(1 - wl_frac)   # e.g., log(0.95)
    wl_log_lowbmi = np.log(0.98)           # 2% WL threshold for BMI<20

    for i in range(df.shape[0]):
        t0 = df[time_col][i]
        b0 = df[bmi_col][i]                    # baseline log-smoothed BMI
        baseline_bmi = df["smoothed_BMI"][i]   # baseline BMI in linear space (smoothed)

        current_date = df["datetime"][i]
        current_bmi = df["bmi"][i]

        df_window = df[(df[time_col] > t0) & (df[time_col] <= t0 + 180)]

        wl_log_i = wl_log_lowbmi if (bmi20_rule and baseline_bmi < 20) else wl_log_default
        df_filtered = df_window[df_window[bmi_col] < b0 + wl_log_i]

        if not df_filtered.empty:
            onset_day = df_filtered[time_col].min()
            onset_date = df_filtered.loc[df_filtered[time_col] == onset_day, "datetime"].values[0] if not df_filtered.loc[df_filtered[time_col] == onset_day].empty else None
            onset_bmi = df_filtered.loc[df_filtered[time_col] == onset_day, "smoothed_BMI"].values[0] if not df_filtered.loc[df_filtered[time_col] == onset_day].empty else None

            max_t = df_filtered[time_col].max()
            end_date = df_filtered.loc[df_filtered[time_col] == max_t, "datetime"].values[0] if not df_filtered.loc[df_filtered[time_col] == max_t].empty else None
            end_bmi = df_filtered.loc[df_filtered[time_col] == max_t, "smoothed_BMI"].values[0] if not df_filtered.loc[df_filtered[time_col] == max_t].empty else None

            results.append({
                "start_day": t0,
                "start_date": current_date,
                "start_bmi": baseline_bmi,  # changed from current_bmi, as the others are using the smoothed units
                "onset_day": onset_day,
                "onset_date": onset_date,
                "onset_bmi": onset_bmi,
                "end_day": max_t,
                "end_date": end_date,
                "end_bmi": end_bmi
            })

    if len(results) == 0:
        results.append({
            "start_day": None,
            "start_date": None,
            "start_bmi": None,
            "onset_day": None,
            "onset_date": None,
            "onset_bmi": None,
            "end_day": None,
            "end_date": None,
            "end_bmi": None
        })

    merged_episodes_df = merge_episodes(pd.DataFrame(results), start_col="start_day", end_col="end_day")

    if recovery:
        merged_episodes_df = identify_recovery_episodes(df, merged_episodes_df, time_col, bmi_col)

    return merged_episodes_df


# ------------- Merging Episodes ----------------------
def merge_episodes(df, start_col, end_col):
    '''
    Merges overlapping cachexia episodes, keeping start, onset, end, and associated dates/BMI.
    '''
    merged_episodes = []
    current_episode = None
    df = df.dropna(subset=[start_col, end_col]).copy()

    for i in range(df.shape[0]):
        if current_episode is None:
            # Initialize the first episode
            current_episode = df.iloc[i].to_dict()
        elif df[start_col][i] <= current_episode['end_day']:
            # Merge episodes if overlapping
            current_episode['end_day'] = max(current_episode['end_day'], df[end_col][i])
            current_episode['end_date'] = df.loc[i, 'end_date']
            current_episode['end_bmi'] = df.loc[i, 'end_bmi']
        else:
            # No overlap, append the previous episode and start a new one
            merged_episodes.append(current_episode)
            current_episode = df.iloc[i].to_dict()

    # Add the last episode
    if current_episode is not None:
        merged_episodes.append(current_episode)

    return pd.DataFrame(merged_episodes)


# ------------- Vectorised episode caller (speed only) ----------------------
# identify_cachexia_episodes_2 + merge_episodes, per person, without the
# per-row pandas indexing. Same calls, same rows, same order. The driver loop
# this replaces did three things that dominated the runtime:
#
#   1. re-scanned the whole measurement table once per person
#      (smoothed_measurements[smoothed_measurements["person_id"] == pid]),
#   2. grew the output with pd.concat inside the loop, recopying it each time,
#   3. repeated both for each of the three weight-loss thresholds.
#
# Here the frame is sliced once into per-person views, every threshold is
# evaluated in the same pass, and the output is concatenated once at the end.

EPISODE_COLS = ["start_day", "start_date", "start_bmi",
                "onset_day", "onset_date", "onset_bmi",
                "end_day", "end_date", "end_bmi"]

WL_LOG_LOWBMI = np.log(0.98)   # the bmi20_rule threshold, as in the original


def _episodes_for_person(t, logb, sbmi, wl_log_default, bmi20_rule):
    """Merged episodes for one person, as positions into that person's arrays.

    Returns (start_pos, onset_pos, end_day_value, end_pos) or None when the
    person has no episode -- which is what the original contributes too: its
    all-None placeholder row is dropped by merge_episodes' dropna, leaving a
    frame with no rows.
    """
    n = t.shape[0]
    # The original's window is value-based -- (t > t0) & (t <= t0 + 180) -- so
    # rows sharing t0's day are excluded even where they sort after it. The day
    # column is sorted within a person, so searchsorted reproduces that exactly.
    # NaN days (from NaT measurement_datetime) sort last and fall outside every
    # window, matching the original, whose comparisons are all False for NaN.
    lo = np.searchsorted(t, t, side="right")
    hi = np.searchsorted(t, t + 180, side="right")

    s_pos, o_pos, e_pos = [], [], []
    for i in range(n):
        j0, j1 = lo[i], hi[i]
        if j1 <= j0:
            continue
        thr = logb[i] + (WL_LOG_LOWBMI if (bmi20_rule and sbmi[i] < 20) else wl_log_default)
        rel = np.flatnonzero(logb[j0:j1] < thr)
        if rel.size == 0:
            continue
        idx = j0 + rel
        tv = t[idx]
        # onset: df_filtered is sorted, so .min() is its first row.
        # end:   the original takes .values[0] of the rows matching max_t, i.e.
        #        the FIRST row on that day, not the last -- which differs from
        #        the last only when several measurements share max_t.
        s_pos.append(i)
        o_pos.append(idx[0])
        e_pos.append(idx[np.searchsorted(tv, tv[-1], side="left")])

    if not s_pos:
        return None

    # merge_episodes: overlapping [start_day, end_day] runs collapse into one.
    # end_date/end_bmi are taken from the merging row unconditionally, even when
    # its end_day is below the running max, so end_date can disagree with
    # end_day. That is the original's behaviour and is kept.
    ms, mo, med, mep = [], [], [], []
    cs = co = ce_day = ce_pos = None
    for k in range(len(s_pos)):
        if cs is None:
            cs, co, ce_day, ce_pos = s_pos[k], o_pos[k], t[e_pos[k]], e_pos[k]
        elif t[s_pos[k]] <= ce_day:
            ce_day = max(ce_day, t[e_pos[k]])
            ce_pos = e_pos[k]
        else:
            ms.append(cs); mo.append(co); med.append(ce_day); mep.append(ce_pos)
            cs, co, ce_day, ce_pos = s_pos[k], o_pos[k], t[e_pos[k]], e_pos[k]
    ms.append(cs); mo.append(co); med.append(ce_day); mep.append(ce_pos)
    return np.array(ms), np.array(mo), np.array(med), np.array(mep)


# Set in the parent before the pool is created; fork gives the workers these
# arrays without pickling the frame.
_SHARED = {}


def _episodes_for_chunk(bounds_slice):
    """Every threshold, for one contiguous block of people."""
    T, LB, SB = _SHARED["T"], _SHARED["LB"], _SHARED["SB"]
    pid, bounds, wl_runs = _SHARED["pid"], _SHARED["bounds"], _SHARED["wl_runs"]
    k0, k1 = bounds_slice
    out = {}
    for label, wl_frac in wl_runs.items():
        wl_log_default = np.log(1 - wl_frac)
        S, O, ED, EP, P = [], [], [], [], []
        any_empty = False
        first_had_episode = None
        for k in range(k0, k1):
            a, b = bounds[k], bounds[k + 1]
            r = _episodes_for_person(T[a:b], LB[a:b], SB[a:b], wl_log_default, True)
            if first_had_episode is None:
                first_had_episode = r is not None
            if r is None:
                any_empty = True
                continue
            ms, mo, med, mep = r
            S.append(a + ms); O.append(a + mo); ED.append(med); EP.append(a + mep)
            P.append(np.full(ms.shape[0], pid[a]))
        out[label] = {
            "S": np.concatenate(S) if S else None,
            "O": np.concatenate(O) if S else None,
            "ED": np.concatenate(ED) if S else None,
            "EP": np.concatenate(EP) if S else None,
            "P": np.concatenate(P) if S else None,
            "any_empty": any_empty,
            "first_had_episode": bool(first_had_episode),
        }
    return out


def call_episodes_fast(smoothed_measurements, wl_runs, time_col, bmi_log_col,
                       n_jobs=1):
    """{label: episode DataFrame}, identical to the original driver's output."""
    sm = smoothed_measurements
    pid = sm["person_id"].to_numpy()
    # Already sorted by (person_id, datetime) above, so each person is a
    # contiguous run and the run order matches the original's .unique() order.
    bounds = np.flatnonzero(np.r_[True, pid[1:] != pid[:-1], True])
    n_people = len(bounds) - 1

    _SHARED.update(T=sm[time_col].to_numpy(), LB=sm[bmi_log_col].to_numpy(),
                   SB=sm["smoothed_BMI"].to_numpy(), pid=pid, bounds=bounds,
                   wl_runs=wl_runs)
    T, SB = _SHARED["T"], _SHARED["SB"]
    DT = sm["datetime"].to_numpy()

    # A few chunks per worker, so people with long trajectories spread out.
    n_chunks = max(1, min(n_people, n_jobs * 4))
    edges = np.linspace(0, n_people, n_chunks + 1).astype(int)
    tasks = [(int(edges[c]), int(edges[c + 1])) for c in range(n_chunks)
             if edges[c + 1] > edges[c]]

    if n_jobs > 1 and len(tasks) > 1:
        # fork: workers inherit _SHARED instead of pickling the frame.
        with get_context("fork").Pool(n_jobs) as pool:
            chunk_results = list(tqdm(pool.imap(_episodes_for_chunk, tasks),
                                      total=len(tasks), desc="Calling episodes"))
    else:
        chunk_results = [_episodes_for_chunk(t)
                         for t in tqdm(tasks, desc="Calling episodes")]

    out = {}
    for label in wl_runs:
        per_chunk = [cr[label] for cr in chunk_results]
        any_person_empty = any(c["any_empty"] for c in per_chunk)
        first_person_had_episode = per_chunk[0]["first_had_episode"]
        filled = [c for c in per_chunk if c["S"] is not None]

        # Chunks are contiguous and in person order, so concatenating them in
        # order reproduces the original's row order exactly.
        if filled:
            S = np.concatenate([c["S"] for c in filled])
            O = np.concatenate([c["O"] for c in filled])
            ED = np.concatenate([c["ED"] for c in filled])
            EP = np.concatenate([c["EP"] for c in filled])
            df = pd.DataFrame({
                "start_day": T[S],  "start_date": DT[S],  "start_bmi": SB[S],
                "onset_day": T[O],  "onset_date": DT[O],  "onset_bmi": SB[O],
                "end_day": ED,      "end_date": DT[EP],   "end_bmi": SB[EP],
            })
            df.insert(0, "person_id", np.concatenate([c["P"] for c in filled]))
            df["wl_threshold"] = label

            # A person with no episodes hands the original's accumulator a frame
            # carrying no episode columns, and the concat upcasts the day columns
            # int64 -> float64 to make room for the implied NaNs. Reproduced
            # under the same condition, so a cohort in which everyone has an
            # episode still comes out int64, as it does today.
            if any_person_empty:
                for c in ("start_day", "onset_day", "end_day"):
                    df[c] = df[c].astype("float64")

            # Column order is likewise decided by the FIRST person: an empty
            # first person seeds the accumulator with person_id/wl_threshold,
            # pushing wl_threshold ahead of the episode columns.
            if not first_person_had_episode:
                df = df[["person_id", "wl_threshold"] + EPISODE_COLS]
        else:
            # Nobody had an episode: the original's accumulator only ever saw
            # column-less frames, so it ends with person_id/wl_threshold and no
            # episode columns at all. Degenerate, but reproduced.
            df = pd.DataFrame({c: pd.Series(dtype="object")
                               for c in ("person_id", "wl_threshold")})
        out[label] = df
    return out


# ------------- Edema QC ----------------------
# Unused in this run: the edema driver code does not run here. The function
# definitions are kept for completeness.
def identify_edema_windows_log(df_mrn, time_col, log_bmi_col,
                               window_days=30, up_frac=0.05, return_frac=0.02):
    """
    Detect edema-like reversible *increase* windows :
      baseline -> >= baseline + log(1+up_frac) -> <= baseline + log(1+return_frac)
    within [baseline_day, baseline_day + window_days].

    Returns DataFrame with: ed_start_day, ed_peak_day, ed_end_day
    """
    df_mrn = df_mrn[[time_col, log_bmi_col]].dropna().sort_values(time_col).reset_index(drop=True)
    if df_mrn.shape[0] < 3:
        return pd.DataFrame(columns=["ed_start_day", "ed_peak_day", "ed_end_day"])

    t = df_mrn[time_col].to_numpy()
    x = df_mrn[log_bmi_col].to_numpy()

    up_log = np.log(1.0 + up_frac)
    ret_log = np.log(1.0 + return_frac)

    out = []
    n = len(t)

    for i in range(n - 2):
        baseline_day = t[i]
        baseline_log = x[i]
        t_end = baseline_day + window_days

        j_end = np.searchsorted(t, t_end, side="right")
        if j_end <= i + 1:
            continue

        win_t = t[i + 1:j_end]
        win_x = x[i + 1:j_end]

        inc_idx_rel = np.where(win_x >= (baseline_log + up_log))[0]
        if inc_idx_rel.size == 0:
            continue

        peak_rel = inc_idx_rel[0]
        peak_day = win_t[peak_rel]

        after_x = win_x[peak_rel + 1:]
        after_t = win_t[peak_rel + 1:]
        if after_x.size == 0:
            continue

        ret_idx_rel = np.where(after_x <= (baseline_log + ret_log))[0]
        if ret_idx_rel.size == 0:
            continue

        end_day = after_t[ret_idx_rel[0]]
        out.append((baseline_day, peak_day, end_day))

    if len(out) == 0:
        return pd.DataFrame(columns=["ed_start_day", "ed_peak_day", "ed_end_day"])

    return pd.DataFrame(out, columns=["ed_start_day", "ed_peak_day", "ed_end_day"])


def merge_overlapping_windows(ed):
    """Merge overlapping edema windows per MRN (ed_start_day/ed_end_day)."""
    if ed.empty:
        return ed
    ed = ed.sort_values("ed_start_day").reset_index(drop=True)
    merged = []
    cur_s = ed.loc[0, "ed_start_day"]
    cur_e = ed.loc[0, "ed_end_day"]

    for k in range(1, ed.shape[0]):
        s = ed.loc[k, "ed_start_day"]
        e = ed.loc[k, "ed_end_day"]
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e

    merged.append((cur_s, cur_e))
    return pd.DataFrame(merged, columns=["ed_start_day", "ed_end_day"])


def flag_episode_edema_overlap(ep_mrn, ed_mrn):
    """Add boolean edema_overlap for episodes that overlap any edema window."""
    ep_mrn = ep_mrn.copy()
    ep_mrn["edema_overlap"] = False

    if ed_mrn.empty or ep_mrn.empty:
        return ep_mrn

    for idx in ep_mrn.index:
        s = ep_mrn.at[idx, "start_day"]
        e = ep_mrn.at[idx, "end_day"]
        if pd.isna(s) or pd.isna(e):
            continue
        ep_mrn.at[idx, "edema_overlap"] = bool(((ed_mrn["ed_start_day"] <= e) & (ed_mrn["ed_end_day"] >= s)).any())

    return ep_mrn


# --- Load longitudinal BMI measurements -------------------------------------

# the filtered, matched cohort cancer patient BMIs
cancer_measurements = pd.read_csv(
    os.path.join(ANALYSIS,
                 "code_outputs_scripts/stable_data/psm_matched_usable_data/solid_heme/full_cohort/all_cancers_all_bmi_measurements_geq4_eligible.tsv"),
    delimiter="\t")
cancer_measurements["cohort"] = "cancer"

# the filtered, matched cohort non-cancer matched patient BMIs
control_measurements = pd.read_csv(
    os.path.join(ANALYSIS,
                 "code_outputs_scripts/stable_data/psm_matched_usable_data/solid_heme/full_cohort/all_controls_all_bmi_measurements_geq4_eligible.tsv"),
    delimiter="\t")
control_measurements["cohort"] = "control"

# combine the two to create one dataframe to process
measurements = pd.concat([cancer_measurements, control_measurements])

# select columns
measurements_minimal = measurements[["person_id",
                                     "cancer_type",
                                     "measurement_datetime",
                                     "value_as_number"]]
measurements_minimal = measurements_minimal.rename(columns={"value_as_number": "bmi"})
measurements_minimal = measurements_minimal.rename(columns={"measurement_datetime": "datetime"})

# coerce datatypes
measurements_minimal["datetime"] = pd.to_datetime(measurements_minimal["datetime"], errors="coerce")
measurements_minimal["bmi"] = pd.to_numeric(measurements_minimal["bmi"], errors="coerce")

# print statistics
print(measurements_minimal["person_id"].nunique())
print(measurements_minimal["cancer_type"].unique())

# define the start date, which is the minimum.
# get each person's earliest measurement date (anchor)
anchor_dates = measurements_minimal.groupby("person_id")["datetime"].transform("min")

# compute days since anchor
measurements_minimal["days_since_anchor"] = (
    measurements_minimal["datetime"] - anchor_dates
).dt.days

# de-dupe per MRN-day (take lowest BMI)
measurements_minimal = (measurements_minimal.sort_values(["person_id", "datetime", "bmi"])
                        .drop_duplicates(subset=["person_id", "datetime"], keep="first")
                        .sort_values(["person_id", "datetime"])
                        .reset_index(drop=True))


# --- constants --------------------------------------------------------------

STAMP = datetime.now().strftime("%Y%m%d")
# All three thresholds run.
WL_RUNS = {"WL5": 0.05, "WL10": 0.10, "WL15": 0.15}
TIME_COL = "days_since_anchor"
BMI_LOG_COL = "log_smoothed_BMI"

OUT_DIR = os.path.join(ANALYSIS, "code_outputs_scripts/cachexia_episodes/processed_data")
os.makedirs(OUT_DIR, exist_ok=True)


# --- smooth BMIs ------------------------------------------------------------

# Set constants for smoothing
alpha = 0.2
smooth_col = "bmi"
# Smooth within person. Uses groupby().transform() rather than
# groupby().apply(smooth_bmi_ewma, ...): since pandas 2.2 the grouping column is
# excluded from the frame handed to apply(), so person_id was dropped from the
# result -- and .reset_index(drop=True) discarded it from the index too, leaving
# no person_id column at all. transform() keeps every column, and the frame is
# already sorted by (person_id, datetime) above, so the EWMA still sees each
# person's measurements in date order.
smoothed_measurements = measurements_minimal.copy()
smoothed_measurements["smoothed_BMI"] = (
    smoothed_measurements.groupby("person_id")[smooth_col]
    .transform(lambda s: s.ewm(alpha=alpha, ignore_na=True).mean())
)

smoothed_measurements[BMI_LOG_COL] = np.log(smoothed_measurements["smoothed_BMI"])

smoothed_bmi_file = os.path.join(OUT_DIR, f"smoothed_bmi_all_patients_{STAMP}_alpha{alpha}.csv")
smoothed_measurements.to_csv(smoothed_bmi_file, index=False)
print("Saved:", smoothed_bmi_file)


# --- detect episodes --------------------------------------------------------

PERSON_IDS = smoothed_measurements["person_id"].unique()

# One pass over people covers all three thresholds. --jobs splits that pass
# across processes; --selfcheck re-runs the original functions on the first N
# people and asserts the frames match, which is how the rewrite is validated
# against the real cohort.
if ARGS.selfcheck:
    head_ids = PERSON_IDS[:ARGS.selfcheck]
    head = smoothed_measurements[
        smoothed_measurements["person_id"].isin(head_ids)].reset_index(drop=True)

    reference = {}
    for label, wl_frac in WL_RUNS.items():
        acc = pd.DataFrame()
        for person_id in tqdm(head_ids, desc=f"selfcheck reference {label}"):
            df = head[head["person_id"] == person_id].copy()
            df.reset_index(drop=True, inplace=True)
            merged = identify_cachexia_episodes_2(
                df, time_col=TIME_COL, bmi_col=BMI_LOG_COL,
                recovery=False, wl_frac=wl_frac, bmi20_rule=True)
            merged["person_id"] = person_id
            merged["wl_threshold"] = label
            acc = pd.concat([acc, merged], ignore_index=True)
        cols = ["person_id"] + [c for c in acc.columns if c != "person_id"]
        reference[label] = acc[cols]

    rewritten = call_episodes_fast(head, WL_RUNS, TIME_COL, BMI_LOG_COL,
                                   n_jobs=ARGS.jobs)
    for label in WL_RUNS:
        pd.testing.assert_frame_equal(reference[label], rewritten[label],
                                      check_dtype=True, check_exact=True)
        print(f"selfcheck {label}: identical "
              f"({len(reference[label])} episode rows over {len(head_ids)} people)")
    print("selfcheck passed on", len(head_ids), "people")

episodes_by_label = call_episodes_fast(smoothed_measurements, WL_RUNS,
                                       TIME_COL, BMI_LOG_COL, n_jobs=ARGS.jobs)

for label in WL_RUNS:

    df_episodes_all = episodes_by_label[label]

    out_fp = os.path.join(OUT_DIR, f"df_episodes_all_precomp_{label}_{STAMP}.csv")
    df_episodes_all.to_csv(out_fp, index=False)
    print("Saved:", out_fp)
    print("Rows:", df_episodes_all.shape[0], "| Unique person_id:", df_episodes_all["person_id"].nunique())

    n_total = df_episodes_all["person_id"].nunique()
    n_any = df_episodes_all.dropna(subset=["start_day"])["person_id"].nunique()
    print(f"Incidence (>=1 episode, pre-QC): {n_any:,}/{n_total:,} ({100 * n_any / n_total:.2f}%)")
