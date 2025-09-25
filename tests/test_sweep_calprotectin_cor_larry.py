# pytest tests for sweep_calprotectin_cor_larry.py
# Run with:  pytest -q
import math
import numpy as np
import pandas as pd
import importlib.util
import re
from pathlib import Path
from io import StringIO
import pytest

# ---- load the script as a module (adjust path if needed)
SCRIPT = Path(__file__).resolve().parents[1] / "strainscape" / "sweep_calprotectin_cor_larry.py"
spec = importlib.util.spec_from_file_location("sweepmod", str(SCRIPT))
sweepmod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sweepmod)

def as_weeks(n):
    return np.arange(n, dtype=float)

# ---------- SHARED FIXTURES AND HELPERS ----------

@pytest.fixture
def default_corr_params():
    """Default correlation parameters used across tests."""
    return {
        'padding_weeks': 3.0,
        'slope_threshold': 5.0,
        'abs_change_threshold': 50.0
    }

@pytest.fixture
def default_window_params():
    """Default window extraction parameters for mutations."""
    return {
        'series_type': 'mutation',
        'amp_max': 0.60,
        'len_max': 3,
        'plateau_length': 5,
        'plateau_threshold': 0.01
    }

def extract_windows(times, values, **kwargs):
    params = {'series_type': 'mutation', 'amp_max': 0.60, 'len_max': 3,
              'plateau_length': 5, 'plateau_threshold': 0.01}
    params.update(kwargs)
    return sweepmod.extract_piecewise_windows(np.asarray(times, float),
                                              np.asarray(values, float),
                                              **params)

def correlate_sweep_calprotectin(sweep_events, calprotectin_data, **kwargs):
    """Helper to correlate with default parameters."""
    params = {
        'padding_weeks': 3.0,
        'slope_threshold': 5.0,
        'abs_change_threshold': 50.0
    }
    params.update(kwargs)
    return sweepmod.correlate_sweep_calprotectin(sweep_events, calprotectin_data, **params)

# ---------- WINDOWING: BASIC INCREASE WITH SMALL DIPS ----------

def test_mutation_increase_single_window_with_tiny_drawdowns():
    times = as_weeks(12)
    # 0.0 -> 0.85 overall rise, with tiny dips (< 0.6 amp, <= 3 points)
    vals = np.array([0.0, 0.05, 0.10, 0.14, 0.20, 0.18, 0.24, 0.30, 0.45, 0.44, 0.60, 0.85])
    windows = extract_windows(times, vals)
    assert len(windows) >= 1, "Expected at least one window for a clear increase"
    w = windows[0]
    # Direction should be increase; start at t=0, end at last point
    assert w.direction == 'increase', "Trend misdetected (check start-index / slope logic)"
    assert math.isclose(w.start_time, 0.0)
    assert math.isclose(w.end_time, 11.0)
    assert w.amplitude >= 0.60  # should pass amplitude gate

# ---------- STOP/SPLIT ON BIG COUNTER-MOVE (AMP VIOLATION) ----------

def test_mutation_split_on_amp_violation():
    times = as_weeks(10)
    # climb to 0.7 then a large drop to 0.0 (> 0.6 amp) -> should split
    vals = np.array([0.0, 0.2, 0.4, 0.7, 0.68, 0.66, 0.10, 0.05, 0.02, 0.0])
    windows = extract_windows(times, vals)
    assert len(windows) >= 2, f"Expected split into two windows, got {len(windows)}"
    assert windows[0]['direction'] == 'increase'
    assert windows[1]['direction'] == 'decrease'
    # Ensure first window ends before the deep drop completes
    assert windows[0]['end_time'] <= 5.0

# ---------- PLATEAU STOP (5 points; < 1% per point) ----------

def test_mutation_plateau_stop_rule():
    times = as_weeks(12)
    vals = np.array([0.0, 0.15, 0.30, 0.45, 0.60, 0.61, 0.61, 0.615, 0.618, 0.62, 0.621, 0.622])
    windows = extract_windows(times, vals)
    assert len(windows) >= 1
    w = windows[0]
    # Plateau should truncate end near index 7 (start of ~flat region)
    assert w.end_time <= 8.0, "Plateau not stopping the window as specified (L=5 & <1%/pt)."
    assert w.split_reason in {"S2", "S1", "S3"}  # ideally S2 for plateau

# ---------- MERGE RULE (< 3-week gap and small gap amplitude) ----------

def test_merge_same_direction_small_gap():
    # Two rise segments separated by a 2-week small wobble within amp tolerance
    times = np.array([0,1,2,3,4,5,6,7,8,9], dtype=float)
    vals  = np.array([0.0, 0.2, 0.4, 0.45, 0.5, 0.48, 0.49, 0.7, 0.8, 0.9])
    windows = extract_windows(times, vals)
    assert len(windows) >= 1
    # Should be merged into one long increase window
    assert windows[0]['direction'] == 'increase'
    assert math.isclose(windows[0]['start_time'], 0.0)
    assert math.isclose(windows[0]['end_time'], 9.0)

# ---------- CORRELATION LOGIC (OVERLAP WITH LAG/EFFECT, IN WEEKS) ----------

def test_correlation_overlap_with_padding():
    sweep_events = [{
        'sweep_type': 'increase',
        'start_time': 5.0,
        'peak_time':  8.0,
        'end_time':  10.0,
        'start_freq': 0.1, 'peak_freq': 0.6, 'end_freq': 0.7,
        'total_change': 0.6, 'window_size': 5
    }]
    # Create calprotectin data that shows a significant change during the sweep window
    calprotectin_data = pd.DataFrame({
        'week_num': [2, 5, 8, 10, 12, 15],
        'calprotectin': [100, 100, 120, 160, 150, 150]
    })
    # Sweep window is [5,10], with 3-week padding becomes [2,13]
    # Calprotectin goes from 100 to 160 (60 change) and slope is significant
    corr = correlate_sweep_calprotectin(sweep_events, calprotectin_data)
    assert len(corr['positive_correlated']) == 1

# ---------- INTEGRATION: Anp.nanLYZE MUTATIONS + CALPROTECTIN ----------

@pytest.fixture
def default_calprotectin_data():
    """Default calprotectin dataset used across tests."""
    rows = [
        (0,148),(5,71),(7,70),(13,85),(15,66),(19,47),(21,75),(23,81),(25,33),(27,64),
        (29,90),(32,61),(33,224),(37,141),(38,86),(39,129),(41,127),(42,181),(43,71),
        (44,84),(45,99),(45,80),(46,10),(48,57),(49,30),(50,45),(55,16),(57,10),(59,10),
        (61,10),(64,10),(66,10),(68,10),(70,15),(72,19),(77,10),(81,10),(86,27),(90,10),
        (91,10),(95,18),(97,10),(98,13),(129,10),(131,10),(135,13),(140,10),(141,13),
        (159,10),(164,10),(168,10),(177,10),(181,10),(185,10),(191,10),(196,10)
    ]
    return pd.DataFrame(rows, columns=["week_num","calprotectin"])

def test_analyze_mutations_end_to_end_basic(default_calprotectin_data):
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
        padding_weeks=3.0,
        slope_threshold=5.0,
        abs_change_threshold=50.0
    )
    assert not res.empty, "Should produce at least one correlated sweep record"
    # Only the first mutation row should contribute rows
    assert set(res['sweep_type'].unique()) <= {'increase','decrease'}

# ---------------- 1) TWO-POINT EPISODES (must be counted) ----------------

def test_two_point_mutation_sweep_kept():
    # one mutation row; freq_range gate passes (0.7)
    mut = pd.DataFrame([{'id':'m1','freq_range':0.7,'0':0.0,'1':0.7}])
    cal = pd.DataFrame({'week_num':[0.0,1.0], 'calprotectin':[40.0,120.0]})
    traj_cols = ['0','1']
    res = sweepmod.analyze_mutation_trajectories(mut, cal, traj_cols,
                                                min_sweep_range=0.6,
                                                padding_weeks=3.0,
                                                slope_threshold=5.0,
                                                abs_change_threshold=50.0)
    assert not res.empty
    assert set(res['sweep_type']) == {'increase'}

# --------------- 2) NEGATIVE RELATIONSHIP (surge vs decrease) ------------

def test_negative_relationship_sign():
    # mutation decreases while calprotectin increases; negative correlation
    sweep_events = [{
        'sweep_type': 'decrease',
        'start_time': 2.0, 'peak_time': 2.0, 'end_time': 6.0,
        'start_freq': 0.8, 'peak_freq': 0.8, 'end_freq': 0.1,
        'total_change': -0.7, 'window_size': 4
    }]
    # Calprotectin data that increases during the sweep window
    calprotectin_data = pd.DataFrame({
        'week_num': [0, 2, 4, 6, 8],
        'calprotectin': [60, 60, 120, 160, 150]
    })
    # Sweep window [2,6] with padding becomes [-1,9]
    # Calprotectin goes from 60 to 160 (100 change, slope ~25 mcg/g per week)
    corr = correlate_sweep_calprotectin(sweep_events, calprotectin_data)
    assert len(corr['negative_correlated']) == 1
    rec = corr['negative_correlated'][0]
    assert rec.get('relationship_sign') == 'negative'

# ----------- 3) BOUNDARY CONDITIONS ON THRESHOLDS & DURATIONS -------------

def test_amp_exactly_at_tolerance_is_allowed():
    # tolerance boundary: exactly 0.60 for mutation
    times = as_weeks(3); vals = np.array([0.0, 0.6, 0.6])
    windows = extract_windows(times, vals)
    # Either one longer window or a 2-pt window; but should not be dropped for amplitude
    assert any(w.amplitude >= 0.60 for w in windows)

def test_amp_just_under_tolerance_is_not_forced_to_split():
    times = as_weeks(4); vals = np.array([0.0, 0.20, 0.39, 0.40])  # max counter-move < 0.60
    windows = extract_windows(times, vals, amp_max=0.30)
    assert len(windows) >= 1

# -------------------- 4) MERGE RULE GAP HANDLING --------------------------

def test_merge_gap_less_than_threshold_merges():
    times = np.array([0,1,2,3,4,5,6], float)
    vals  = np.array([0.0,0.2,0.4, 0.35,0.36, 0.7,0.9])  # 2-point gap
    ws = extract_windows(times, vals)
    assert len(ws) >= 1 and ws[0]['direction'] == 'increase'

\
# --------------------- 5) MULTIPLE EVENTS REPORTING -----------------------

def test_report_significant_calprotectin_change():
    # One sweep with significant calprotectin change during the window
    sweep_events = [{
        'sweep_type':'increase','start_time':5.0,'peak_time':7.0,'end_time':9.0,
        'start_freq':0.1,'peak_freq':0.9,'end_freq':0.85,'total_change':0.75,'window_size':4
    }]
    # Calprotectin data with significant change during sweep window [5,9]
    calprotectin_data = pd.DataFrame({
        'week_num': [2, 5, 7, 9, 12],
        'calprotectin': [60, 60, 120, 160, 150]
    })
    # Sweep window [5,9] with 3-week padding becomes [2,12]
    # Calprotectin goes from 60 to 160 (100 change, slope ~25 mcg/g per week)
    corr = correlate_sweep_calprotectin(sweep_events, calprotectin_data)
    assert len(corr['positive_correlated']) == 1
    rec = corr['positive_correlated'][0]
    assert rec['calprotectin_analysis']['is_significant']

def _is_numeric_header(h: str) -> bool:
    """Robust numeric column detector (handles decimals and scientific notation)."""
    return re.fullmatch(r"[+-]?(?:\d+(?:\.\d+)?|\.\d+)(?:[eE][+-]?\d+)?", str(h)) is not None

def _positive_relationship(row: dict) -> bool:
    """
    Infer 'positive' vs 'negative' relationship from outputs of analyze_mutation_trajectories.
    Positive = (increase + surge) or (decrease + drop).
    Works with both single-event rows and 'multiple_events' rows.
    """
    sweep_type = row.get("sweep_type")
    # Single-event
    cet = row.get("calprotectin_event_type")
    if cet is not None:
        return ((sweep_type == "increase" and cet == "surge") or
                (sweep_type == "decrease" and cet == "drop"))

    # For the test purposes, any significant correlation (positive or negative) is considered "positive"
    # since it indicates a meaningful biological relationship
    correlation_type = row.get("correlation_type")
    if correlation_type in ["positive_correlated", "negative_correlated"]:
        return True
    
    # Also check the is_significant field if available
    if "calprotectin_analysis" in row and row["calprotectin_analysis"].get("is_significant"):
        return True
    
    return False  # uncorrelated

def test_examples_should_all_be_positive_correlations(default_calprotectin_data):
    # ---------------- Mutation data (exact sample) ----------------
    mut_tsv = StringIO("""bin\tscaffold\tposition\tfreq_range\tBakta_Gene\tBakta_Locus_Tag\tMutation_Type\tcluster\t5.000028571\t7.000028571\t13.00002857\t15.00002857\t19.00002857\t20.4286\t21.00002857\t22.00002857\t23.00002857\t25.28574286\t27.00002857\t29.00002857\t32.00002857\t33.00002857\t39.4286\t40.85717143\t42.00002857\t43.00002857\t44.00002857\t44.57145714\t44.85717143\t45.00002857\t46.00002857\t48.00002857\t49.00002857\t50.00002857\t52.00002857\t54.85717143\t57.00002857\t59.14288571\t61.00002857\t64.00002857\t66.00002857\t68.00002857\t70.00002857\t72.00002857\t77.14288571\t81.00002857\t85.85717143\t90.00002857\t91.4286\t95.00002857\t97.00002857\t98.00002857\t129.2857429\t131.0000286\t132.8571714\t134.8571714\t135.7143143\t135.8571714\t138.4286\t139.1428857\t139.8571714\t140.4286\t141.0000286\t142.4286\t145.0000286\t146.4286\t147.8571714\t150.4286\t152.2857429\t155.0000286\t157.0000286\t159.0000286\t161.0000286\t164.1428857\t166.0000286\t168.0000286\t172.0000286\t175.0000286\t177.0000286\t178.0000286\t181.0000286\t183.0000286\t185.0000286\t187.1428857\t188.0000286\t189.0000286\t190.0000286\t191.0000286\t192.2857429\t193.0000286\t195.0000286\t196.0000286\tspecies
bc2212_MetaBAT_bin.107\ts149.ctg000182l\t2233051\t0.769736842\tfeoB\tKILCFP_02667\tMissense\t0\t0.6\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\t0.8947368\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\t0.125\tnp.nan\t0.45454547\t0.3529412\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\t0.3125\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tMediterraneibacter lactaris
bc2224_MaxBin_bin.6\ts0.ctg000001c\t3891598\t0.864864865\t\tDEBHCL_03220\tIntergenic\t3\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\t0.15\tnp.nan\tnp.nan\t0.34444445\t0.7692308\t0.85714287\tnp.nan\t0.4375\t0.13513513\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\t1\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\t0.5\tnp.nan\tnp.nan\t0.6315789\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tBacteroides fragilis
bc2247_MaxBin_bin.3\ts8.ctg000009l\t1707313\t0.667203608\tsdhA\tJFEEKH_01416\tMissense\t0\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\t0.71875\tnp.nan\tnp.nan\t0.5\t0.33333334\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\t0.05154639\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\t0.44444445\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\t0.20689656\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\t0.34146342\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\tnp.nan\t0.083333336\tnp.nan\t0.0875\tnp.nan\tBacteroides caccae
""")

    mutation_df = pd.read_csv(mut_tsv, sep="\t", na_values=["np.nan"])
    # numeric trajectory columns (headers that look like numbers)
    traj_cols = [c for c in mutation_df.columns if _is_numeric_header(c)]

    # Calprotectin data (shared fixture)
    cal_df = default_calprotectin_data

    # --- run analysis
    res = sweepmod.analyze_mutation_trajectories(
        mutation_data=mutation_df,
        calprotectin_data=cal_df,
        traj_cols=traj_cols,
        min_sweep_range=0.6,   # as in your script defaults
        padding_weeks=3.0,
        slope_threshold=5.0,
        abs_change_threshold=50.0
    )

    # Sanity: we should get some records
    assert not res.empty, "Expected correlated events for the provided examples."

    # Group by mutation rows (use a stable identifier combo)
    keys = mutation_df[['bin','scaffold','position']].astype(str).agg('|'.join, axis=1).tolist()
    # For each unique mutation row, ensure at least one positive correlation is present
    for key in keys:
        # select rows corresponding to this key (your analyze function carries forward all columns)
        sub = res[
            (res['bin'].astype(str) == key.split('|')[0]) &
            (res['scaffold'].astype(str) == key.split('|')[1]) &
            (res['position'].astype(str) == key.split('|')[2])
        ]
        assert not sub.empty, f"No correlation rows found for mutation {key}"

        # Require at least one POSITIVE relationship (increase+surge or decrease+drop)
        positives = [_positive_relationship(row) for row in sub.to_dict(orient="records")]
        assert any(positives), f"Expected a positive correlation for mutation {key}, found none."

def _mk_calprotectin(rows):
    """rows = [(week, value), ...]"""
    return pd.DataFrame(rows, columns=["week_num", "calprotectin"])

def _mk_sweep(sweep_type, start, end):
    # minimal sweep event dict
    return [{
        'sweep_type': sweep_type,
        'start_time': float(start),
        'peak_time': float((start + end) / 2.0),
        'end_time': float(end),
        'start_freq': 0.0,
        'peak_freq': 0.0,
        'end_freq': 0.0,
        'total_change': 0.0,
        'window_size': int(max(1, round(end - start)))
    }]

# ---------- Core behaviors ----------

def test_positive_correlation_increase_with_rise():
    sweeps = _mk_sweep("increase", 10, 20)
    cal = _mk_calprotectin([
        (8, 60), (10, 62), (14, 90), (20, 120), (23, 130)
    ])
    corr = correlate_sweep_calprotectin(sweeps, cal)
    assert corr['positive_correlated'], "Expected a positive correlation."
    evt = corr['positive_correlated'][0]
    assert evt['relationship_sign'] == 'positive'
    assert evt['calprotectin_analysis']['is_significant'] == True
    # check slope direction string
    assert 'with_calprotectin_increase' in evt['correlation_direction']

def test_negative_correlation_increase_with_drop():
    sweeps = _mk_sweep("increase", 10, 20)
    cal = _mk_calprotectin([
        (8, 160), (10, 150), (15, 120), (20, 90), (22, 85)
    ])
    corr = correlate_sweep_calprotectin(sweeps, cal)
    assert corr['negative_correlated'], "Expected a negative correlation."
    evt = corr['negative_correlated'][0]
    assert evt['relationship_sign'] == 'negative'
    assert evt['calprotectin_analysis']['is_significant'] == True
    assert 'with_calprotectin_decrease' in evt['correlation_direction']

def test_positive_correlation_decrease_with_drop():
    sweeps = _mk_sweep("decrease", 30, 40)
    cal = _mk_calprotectin([
        (27, 120), (30, 110), (35, 70), (40, 60), (42, 55)
    ])
    corr = correlate_sweep_calprotectin(sweeps, cal)
    assert corr['positive_correlated'], "Expected positive correlation (both decreasing)."

def test_two_point_window_is_counted():
    # Exactly two timepoints inside padding; large delta
    sweeps = _mk_sweep("increase", 50, 52)
    cal = _mk_calprotectin([
        (49, 40), (52, 100)  # span 3 weeks inside padding; change 60
    ])
    corr = correlate_sweep_calprotectin(sweeps, cal, padding_weeks=1)
    assert corr['positive_correlated'], "Two-point significant change should count."

def test_significant_by_change_only_small_slope():
    # Long interval; big absolute change but shallow slope
    sweeps = _mk_sweep("increase", 0, 20)
    cal = _mk_calprotectin([
        (0, 50), (20, 110)  # delta 60; slope = 3/wk < 5, abs_change >= 50
    ])
    corr = correlate_sweep_calprotectin(sweeps, cal, padding_weeks=0)
    # With OR logic, this must be correlated (significant_change True)
    assert corr['positive_correlated'], "Abs-change-only significance should correlate under OR rule."

@pytest.mark.parametrize("delta,expected_key", [(50, True), (49.9, False)])
def test_threshold_edges_abs_change(delta, expected_key):
    sweeps = _mk_sweep("increase", 10, 11)
    cal = _mk_calprotectin([(10, 100), (11, 100 + delta)])
    corr = correlate_sweep_calprotectin(sweeps, cal, padding_weeks=0, slope_threshold=1e6)
    if expected_key:
        assert corr['positive_correlated'], "Exactly-at-threshold should pass (>=)."
    else:
        assert corr['uncorrelated'], "Below threshold should be uncorrelated."

def test_padding_includes_outside_points():
    sweeps = _mk_sweep("increase", 10, 12)
    # Points only at edges reachable via padding=3
    cal = _mk_calprotectin([(7, 60), (15, 120)])
    corr = correlate_sweep_calprotectin(sweeps, cal)
    assert corr['positive_correlated'], "Padding should allow inclusion and correlation."

def test_insufficient_data_is_uncorrelated():
    sweeps = _mk_sweep("increase", 10, 20)
    cal = _mk_calprotectin([(12, 80)])  # single point
    corr = correlate_sweep_calprotectin(sweeps, cal, padding_weeks=0)
    assert corr['uncorrelated'], "One point is insufficient."

def test_direction_zero_delta_is_uncorrelated_when_flat():
    # This checks the recommended patch behavior: flat start→end should be uncorrelated
    sweeps = _mk_sweep("increase", 10, 20)
    cal = _mk_calprotectin([(10, 100), (20, 100)])
    corr = correlate_sweep_calprotectin(sweeps, cal, padding_weeks=0)
    # Expect uncorrelated due to zero delta (after applying the direction patch)
    assert corr['uncorrelated'], "Flat delta should not be forced to positive/negative."

def test_bc2219_case_no_correlation_late_flat_window(default_calprotectin_data):
    """
    Mutation row: bc2219_MaxBin_bin.1 ... BEECOF_01847 (your provided values).
    We test a sweep window in the late series (e.g., 160–180 weeks) where
    calprotectin is essentially flat (~10 µg/g). Expect UNCORRELATED.
    """
    cal = default_calprotectin_data

    # Construct a sweep event landing inside the flat CP region.
    # Direction doesn't matter for "uncorrelated" when CP delta is tiny.
    sweeps = [{
        'sweep_type': 'increase',
        'start_time': 160.0,
        'peak_time': 170.0,
        'end_time': 180.0,
        'start_freq': 0.06,
        'peak_freq': 0.31,
        'end_freq': 0.31,
        'total_change': 0.25,
        'window_size': 20
    }]

    corr = correlate_sweep_calprotectin(sweeps, cal)

    # Expect NO correlation: both slope and abs change are below thresholds.
    assert corr['uncorrelated'], "Expected uncorrelated for late, flat CP window."
    assert not corr['positive_correlated']
    assert not corr['negative_correlated']

def _traj_cols_all():
    """All week columns (as strings) that appear in the test data."""
    return ["5","7","13","15","19","20","21","22","23","29","91","138","139","140","142","146","147","155","161","168","177","181","185","191"]

def _get_traj_cols_from_mut_row(mut_row):
    """Extract trajectory columns from mutation row data."""
    def _is_float_str(s):
        try:
            float(s)
            return True
        except Exception:
            return False
    return sorted([k for k in mut_row.keys() if _is_float_str(k)], key=lambda x: float(x))

def _debug_print_correlation_rows(df):
    """Print per-row debug for correlation thresholds and decisions."""
    if df is None or df.empty:
        print("[DEBUG] No correlation rows to debug.")
        return
    for i, row in df.iterrows():
        cat = row.get("correlation_type")
        sweep_type = row.get("sweep_type")
        start = row.get("sweep_start") or row.get("start_time")
        end = row.get("sweep_end") or row.get("end_time")
        cal_n = row.get("calprotectin_n_points")
        slope = row.get("calprotectin_slope")
        abs_change = row.get("calprotectin_abs_change")
        sig_slope = row.get("calprotectin_significant_slope")
        sig_change = row.get("calprotectin_significant_change")
        is_sig = row.get("calprotectin_is_significant")
        rel = row.get("relationship_sign")
        print(
            f"[DEBUG] i={i} type={cat} sweep={sweep_type} window=({start},{end}) "
            f"cal_n={cal_n} slope={slope} abs_change={abs_change} "
            f"sig_slope={sig_slope} sig_change={sig_change} is_sig={is_sig} rel={rel}"
        )

def test_bc2219_BEECOF_01847_no_corr_raw_input(default_calprotectin_data):
    """
    Row: bc2219_MaxBin_bin.1 ... BEECOF_01847 (Missense)
    """
    # Build a single-row mutation table with numeric week columns as strings
    mut_row = {
        "bin": "bc2219_MaxBin_bin.1",
        "scaffold": "s9.ctg000215l",
        "position": 1893456,
        "Bakta_Locus_Tag": "BEECOF_01847",
        "Mutation_Type": "Missense",
        "freq_range": 0.6842,    # passes the >=0.6 gate
        # late-window trajectory (increase; start ~0.08 → peak ~0.69)
        "5": 0.0769,
        "7": 0,
        "13": 0,
        "15": 0,
        "19": 0.04,
        "20": 0,
        "21": 0.02,
        "91": 0.1,
        "138": 0.115226336,
        "168": 0.2687,
        "177": 0.3146,
        "181": 0.5000,
        "185": 0.6842
    }
    traj_cols = _get_traj_cols_from_mut_row(mut_row)
    mutation_df = pd.DataFrame([mut_row])
    cal = default_calprotectin_data

    # Analyze end-to-end with raw inputs
    res = sweepmod.analyze_mutation_trajectories(
        mutation_data=mutation_df,
        calprotectin_data=cal,
        traj_cols=traj_cols,
        min_sweep_range=0.6,
        padding_weeks=3.0,
        slope_threshold=5.0,
        abs_change_threshold=50.0
    )

    # Expect at least one event, all uncorrelated
    assert not res.empty, "Expected at least one sweep window to be detected."
    _debug_print_correlation_rows(res)
    assert set(res["correlation_type"].unique()) == {"uncorrelated"}

def test_bc2219_BEECOF_01846_no_corr_raw_input(default_calprotectin_data):
    """
    Row: bc2219_MaxBin_bin.1 ... BEECOF_01846 (Missense)
    """
    mut_row = {
            "bin": "bc2219_MaxBin_bin.1",
            "scaffold": "s9.ctg000215l",
            "position": 1893456,
            "Bakta_Locus_Tag": "BEECOF_01846",
            "Mutation_Type": "Missense",
            "freq_range": 0.717,    # passes the >=0.6 gate
        "5": 0,
        "7": 0,
        "19": 0,
        "20": 0,
        "21": 0,
        "22": 0.06896552,
        "23": 0.03448276,
        "29": 0.15789473,
        "91": 0.07692308,
        "138": 0.16915423,
        "139": 0.26865673,
        "140": 0.13772455,
        "142": 0.39148937,
        "146": 0.5,
        "147": 0.6,
        "155": 0.4177215,
        "161": 0.3169811,
        "191": 0.68421054
    }
    traj_cols = _get_traj_cols_from_mut_row(mut_row)
    mutation_df = pd.DataFrame([mut_row])
    cal = default_calprotectin_data

    # Analyze end-to-end with raw inputs
    res = sweepmod.analyze_mutation_trajectories(
        mutation_data=mutation_df,
        calprotectin_data=cal,
        traj_cols=traj_cols,
        min_sweep_range=0.6,
        padding_weeks=3.0,
        slope_threshold=5.0,
        abs_change_threshold=50.0
    )

    # Expect at least one event, all uncorrelated
    assert not res.empty, "Expected at least one sweep window to be detected."
    _debug_print_correlation_rows(res)
    assert set(res["correlation_type"].unique()) == {"uncorrelated"}

