# pytest tests for sweep_calprotectin_cor_larry.py
# Run with:  pytest -q
import math
import numpy as np
import pandas as pd
import importlib.util
from pathlib import Path
import pytest

# ---- load the script as a module (adjust path if needed)
SCRIPT = Path(__file__).resolve().parents[1] / "strainscape" / "sweep_calprotectin_cor_larry.py"
spec = importlib.util.spec_from_file_location("sweepmod", str(SCRIPT))
sweepmod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sweepmod)

def as_weeks(n):
    return np.arange(n, dtype=float)

# ---------- WINDOWING: BASIC INCREASE WITH SMALL DIPS ----------

def test_mutation_increase_single_window_with_tiny_drawdowns():
    times = as_weeks(12)
    # 0.0 -> 0.85 overall rise, with tiny dips (< 0.6 amp, <= 3 points)
    vals = np.array([0.0, 0.05, 0.10, 0.14, 0.20, 0.18, 0.24, 0.30, 0.45, 0.44, 0.60, 0.85])
    windows = sweepmod.extract_piecewise_windows(
        times=times,
        values=vals,
        series_type="mutation",
        amp_max=0.60,          # your tolerance
        len_max=3,
        plateau_length=5,
        plateau_threshold=0.01
    )
    assert len(windows) >= 1, "Expected at least one window for a clear increase"
    w = windows[0]
    # Direction should be increase; start at t=0, end at last point
    assert w['direction'] == 'increase', "Trend misdetected (check start-index / slope logic)"
    assert math.isclose(w['start_time'], 0.0)
    assert math.isclose(w['end_time'], 11.0)
    assert w['amplitude'] >= 0.60  # should pass amplitude gate

# ---------- STOP/SPLIT ON BIG COUNTER-MOVE (AMP VIOLATION) ----------

def test_mutation_split_on_amp_violation():
    times = as_weeks(10)
    # climb to 0.7 then a large drop to 0.0 (> 0.6 amp) -> should split
    vals = np.array([0.0, 0.2, 0.4, 0.7, 0.68, 0.66, 0.10, 0.05, 0.02, 0.0])
    windows = sweepmod.extract_piecewise_windows(
        times, vals, "mutation", amp_max=0.60, len_max=3, plateau_length=5, plateau_threshold=0.01
    )
    assert len(windows) >= 2, f"Expected split into two windows, got {len(windows)}"
    assert windows[0]['direction'] == 'increase'
    assert windows[1]['direction'] == 'decrease'
    # Ensure first window ends before the deep drop completes
    assert windows[0]['end_time'] <= 5.0

# ---------- PLATEAU STOP (5 points; < 1% per point) ----------

def test_mutation_plateau_stop_rule():
    times = as_weeks(12)
    vals = np.array([0.0, 0.15, 0.30, 0.45, 0.60, 0.61, 0.61, 0.615, 0.618, 0.62, 0.621, 0.622])
    windows = sweepmod.extract_piecewise_windows(
        times, vals, "mutation", amp_max=0.60, len_max=3, plateau_length=5, plateau_threshold=0.01
    )
    assert len(windows) >= 1
    w = windows[0]
    # Plateau should truncate end near index 7 (start of ~flat region)
    assert w['end_time'] <= 8.0, "Plateau not stopping the window as specified (L=5 & <1%/pt)."
    assert w['split_reason'] in {"S2", "S1", "S3"}  # ideally S2 for plateau

# ---------- MERGE RULE (< 3-week gap and small gap amplitude) ----------

def test_merge_same_direction_small_gap():
    # Two rise segments separated by a 2-week small wobble within amp tolerance
    times = np.array([0,1,2,3,4,5,6,7,8,9], dtype=float)
    vals  = np.array([0.0, 0.2, 0.4, 0.45, 0.5, 0.48, 0.49, 0.7, 0.8, 0.9])
    windows = sweepmod.extract_piecewise_windows(
        times, vals, "mutation", amp_max=0.60, len_max=3, plateau_length=5, plateau_threshold=0.01
    )
    assert len(windows) >= 1
    # Should be merged into one long increase window
    assert windows[0]['direction'] == 'increase'
    assert math.isclose(windows[0]['start_time'], 0.0)
    assert math.isclose(windows[0]['end_time'], 9.0)

# ---------- CALPROTECTIN EVENTS VIA PIECEWISE ----------

def test_calprotectin_surge_event_detected():
    # 40 -> 170 surge (130 µg/g), small dips inside tolerance
    weeks = as_weeks(8)
    cal   = np.array([40, 45, 50, 60, 95, 130, 170, 155])
    df = pd.DataFrame({'week_num': weeks, 'calprotectin': cal})
    events = sweepmod.detect_calprotectin_events(df, abs_change_threshold=50.0)
    assert len(events) >= 1, "Expected at least one calprotectin event"
    e = events[0]
    assert e['event_type'] == 'surge'
    assert e['abs_delta'] >= 50.0

@pytest.mark.xfail(strict=True, reason="Current implementation drops windows with amplitude < abs_change_threshold")
def test_calprotectin_smaller_episode_kept_if_spec_intended():
    # If you WANT to keep 30–40 µg/g episodes, this should pass after you decouple min amplitude from amp_max
    weeks = as_weeks(6)
    cal   = np.array([40, 55, 65, 75, 70, 68])  # ~35 rise then small fall
    df = pd.DataFrame({'week_num': weeks, 'calprotectin': cal})
    events = sweepmod.detect_calprotectin_events(df, abs_change_threshold=50.0)
    assert len(events) >= 1, "Spec says 'keep any window that passes rules'; amplitude < 50 shouldn't auto-drop."

# ---------- CORRELATION LOGIC (OVERLAP WITH LAG/EFFECT, IN WEEKS) ----------

def test_correlation_overlap_with_lag_effect():
    sweep_events = [{
        'sweep_type': 'increase',
        'start_time': 5.0,
        'peak_time':  8.0,
        'end_time':  10.0,
        'start_freq': 0.1, 'peak_freq': 0.6, 'end_freq': 0.7,
        'total_change': 0.6, 'window_size': 5
    }]
    cal_events = [{
        'event_type': 'surge',
        'start_time': 10.5,
        'peak_time':  11.0,
        'end_time':   12.0,
        'start_level': 100.0, 'peak_level': 160.0, 'end_level': 150.0,
        'fold_change': 1.6, 'baseline': np.nan, 'baseline_category': 'borderline',
        'peak_category': 'active', 'abs_delta': 60.0
    }]
    # lag/effect in WEEKS (your function treats them as weeks)
    corr = sweepmod.correlate_sweep_calprotectin(sweep_events, cal_events, lag_days=2, effect_days=2)
    # Event window expands to [10.5-2, 12+2] = [8.5, 14]; overlaps sweep [5,10] -> yes.
    assert len(corr['surge_correlated']) == 1

# ---------- INTEGRATION: ANALYZE MUTATIONS + CALPROTECTIN ----------

def test_analyze_mutations_end_to_end_basic():
    # mutation table with two rows; one has freq_range >= 0.6
    mutation_data = pd.DataFrame({
        'id': ['a','b'],
        'freq_range': [0.65, 0.2],
        '0': [0.0, 0.05],
        '1': [0.1, 0.04],
        '2': [0.2, 0.03],
        '3': [0.5, 0.02],
        '4': [0.7, 0.01],
        '5': [0.8, 0.00],
    }).set_index('id', drop=False)
    cal_df = pd.DataFrame({'week_num': [0,1,2,3,4,5], 'calprotectin': [40, 45, 60, 110, 170, 160]})
    traj_cols = ['0','1','2','3','4','5']  # simulate your detection
    res = sweepmod.analyze_mutation_trajectories(
        mutation_data.reset_index(drop=True),
        cal_df,
        traj_cols,
        min_sweep_range=0.6,
        lag_days=2,
        effect_days=2
    )
    assert not res.empty, "Should produce at least one correlated sweep record"
    # Only the first mutation row should contribute rows
    assert set(res['sweep_type'].unique()) <= {'increase','decrease'}

# ---------- COLUMN NAME NUMERIC DETECTION (FYI) ----------

def test_numeric_column_names_scientific_notation():
    # This mirrors your heuristic; update if you change it in script.
    cols = ['0', '1.0', '1e-3', '2.85714e-05', '-3', '+4', '.5', 'week', 'foo']
    df = pd.DataFrame(columns=cols)
    # emulate your detection
    traj_cols = [c for c in df.columns
                 if c.replace('.', '').replace('e-', '').replace('e+', '')
                   .replace('-', '').replace('+', '').replace('E-', '').replace('E+', '')
                   .isdigit()]
    # .5 will be excluded; scientific notation should be included by your heuristic
    assert '2.85714e-05' in traj_cols
    assert '1e-3' in traj_cols
    assert '.5' not in traj_cols  # FYI, this is excluded by your current logic

# ---------------- 1) TWO-POINT EPISODES (must be counted) ----------------

def test_two_point_calprotectin_episode_kept():
    weeks = np.array([0.0, 1.0])
    cal   = np.array([40.0, 95.0])  # +55 µg/g ≥ 50
    df = pd.DataFrame({'week_num': weeks, 'calprotectin': cal})
    events = sweepmod.detect_calprotectin_events(df, abs_change_threshold=50.0)
    assert len(events) == 1
    e = events[0]
    assert e['event_type'] == 'surge'
    assert math.isclose(e['start_time'], 0.0) and math.isclose(e['end_time'], 1.0)

def test_two_point_mutation_sweep_kept():
    # one mutation row; freq_range gate passes (0.7)
    mut = pd.DataFrame([{'id':'m1','freq_range':0.7,'0':0.0,'1':0.7}])
    cal = pd.DataFrame({'week_num':[0.0,1.0], 'calprotectin':[40.0,120.0]})
    traj_cols = ['0','1']
    res = sweepmod.analyze_mutation_trajectories(mut, cal, traj_cols,
                                                 min_sweep_range=0.6,
                                                 lag_days=0, effect_days=0)
    assert not res.empty
    assert set(res['sweep_type']) == {'increase'}

# --------------- 2) NEGATIVE RELATIONSHIP (surge vs decrease) ------------

def test_negative_relationship_sign():
    # mutation decreases while calprotectin surges; overlapping windows
    sweep_events = [{
        'sweep_type': 'decrease',
        'start_time': 2.0, 'peak_time': 2.0, 'end_time': 6.0,
        'start_freq': 0.8, 'peak_freq': 0.8, 'end_freq': 0.1,
        'total_change': -0.7, 'window_size': 4
    }]
    cal_events = [{
        'event_type': 'surge',
        'start_time': 4.0, 'peak_time': 5.0, 'end_time': 7.0,
        'start_level': 60.0, 'peak_level': 160.0, 'end_level': 150.0,
        'fold_change': 2.5, 'baseline': np.nan,
        'baseline_category': 'borderline', 'peak_category': 'active',
        'abs_delta': 90.0
    }]
    corr = sweepmod.correlate_sweep_calprotectin(sweep_events, cal_events,
                                                 lag_days=0, effect_days=0)
    assert len(corr['surge_correlated']) == 1
    rec = corr['surge_correlated'][0]
    # requires code patch to add 'relationship_sign'
    assert rec.get('relationship_sign') == 'negative'

# ----------- 3) BOUNDARY CONDITIONS ON THRESHOLDS & DURATIONS -------------

def test_amp_exactly_at_tolerance_is_allowed():
    # tolerance boundary: exactly 0.60 for mutation
    times = as_weeks(3); vals = np.array([0.0, 0.6, 0.6])
    windows = sweepmod.extract_piecewise_windows(times, vals, "mutation",
                                                 amp_max=0.60, len_max=3,
                                                 plateau_length=5, plateau_threshold=0.01)
    # Either one longer window or a 2-pt window; but should not be dropped for amplitude
    assert any(w['amplitude'] >= 0.60 for w in windows)

def test_amp_just_under_tolerance_is_not_forced_to_split():
    times = as_weeks(4); vals = np.array([0.0, 0.20, 0.39, 0.40])  # max counter-move < 0.60
    windows = sweepmod.extract_piecewise_windows(times, vals, "mutation",
                                                 amp_max=0.30, len_max=3,
                                                 plateau_length=5, plateau_threshold=0.01)
    assert len(windows) >= 1

# -------------------- 4) MERGE RULE GAP HANDLING --------------------------

def test_merge_gap_less_than_threshold_merges():
    times = np.array([0,1,2,3,4,5,6], float)
    vals  = np.array([0.0,0.2,0.4, 0.35,0.36, 0.7,0.9])  # 2-point gap
    ws = sweepmod.extract_piecewise_windows(times, vals, "mutation",
                                            amp_max=0.60, len_max=3,
                                            plateau_length=5, plateau_threshold=0.01)
    assert len(ws) >= 1 and ws[0]['direction'] == 'increase'

def test_merge_gap_equal_threshold_does_not_merge():
    # if your code uses < 3.0, a 3.0-week gap must not merge
    # craft times so the gap is exactly 3.0
    times = np.array([0,1,2,5,6,7], float)  # gap from 2->5 is 3 weeks
    vals  = np.array([0.0,0.2,0.4, 0.5,0.7,0.9])
    ws = sweepmod.extract_piecewise_windows(times, vals, "mutation",
                                            amp_max=0.30, len_max=3,
                                            plateau_length=5, plateau_threshold=0.01)
    # Expect two windows if strict '< 3.0' is enforced
    assert len(ws) >= 2

# --------------------- 5) MULTIPLE EVENTS REPORTING -----------------------

def test_report_multiple_events_nonoverlapping():
    # Two non-overlapping calprotectin surges; one mutation increase spanning both → multiple
    cal_events = [
        {'event_type':'surge','start_time':2.0,'peak_time':3.0,'end_time':4.0,
         'start_level':60,'peak_level':140,'end_level':120,'fold_change':2.3,
         'baseline':np.nan,'baseline_category':'borderline','peak_category':'active','abs_delta':80.0},
        {'event_type':'drop','start_time':8.0,'peak_time':9.0,'end_time':10.0,
         'start_level':160,'peak_level':80,'end_level':90,'fold_change':0.5,
         'baseline':np.nan,'baseline_category':'active','peak_category':'borderline','abs_delta':70.0}
    ]
    sweep_events = [{
        'sweep_type':'increase','start_time':1.0,'peak_time':9.5,'end_time':11.0,
        'start_freq':0.1,'peak_freq':0.9,'end_freq':0.85,'total_change':0.75,'window_size':10
    }]
    corr = sweepmod.correlate_sweep_calprotectin(sweep_events, cal_events, lag_days=0, effect_days=0)
    assert len(corr['multiple_events']) == 1
    me = corr['multiple_events'][0]
    assert len(me['calprotectin_events']) == 2
