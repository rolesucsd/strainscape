#!/usr/bin/env python3
"""
sweep_calprotectin_cor_larry.py

===============================================================================
Windowing & Correlation Rules (Final Spec)
===============================================================================

Goal
----
Segment each time series (mutation frequency or calprotectin) into
"piecewise, mostly one-direction" windows that tolerate brief counter-moves,
then correlate mutation windows to calprotectin windows.

Vocabulary
----------
- Point: one measurement in time (usually ~weekly).
- Trend: direction of change (increase or decrease)
- Temporary counter-move (a.k.a. drawdown/drawup): a brief move opposite the
  current trend.
- Window: a contiguous span of points that is monotone "with tolerance."

Series-specific thresholds
--------------------------
We apply the same rules to both series, but with different amplitudes.

• MUTATION FREQUENCY (0–1 scale; volatile)
  - Allowable temporary counter-move amplitude (δ_amp_max_mut): 0.60
  - Allowable temporary counter-move length    (δ_len_max):     3 points (≈ 3 wks)
  - Plateau stop test:
      — No net gain in trend direction for L = 5 points
      — AND average per-point gain < 0.01 (1 percentage point)
  - Sweep call: no minimum duration; "net change" gate is handled elsewhere
    (e.g., freq_range ≥ 0.60).

• CALPROTECTIN (µg/g)
  - Allowable temporary counter-move amplitude (δ_amp_max_cal): 50 µg/g
  - Allowable temporary counter-move length    (δ_len_max):     3 points (≈ 3 wks)
  - Plateau stop test:
      — No net gain in trend direction for L = 5 points
      — AND average per-point gain < 10 µg/g
  - Biological context retained (for reporting, not for splitting):
      Normal < 50; Borderline 50–120; Active ≥ 120 µg/g.

Piecewise window extraction (applies to BOTH series)
----------------------------------------------------
Given a time-ordered series after light smoothing:

1) Seed / start:
   Scan left→right. Start a new window at a local extremum where the robust
   slope (e.g., over last 2–3 points) establishes a trend sign.

2) Grow (monotone-with-tolerance):
   While growing an INCREASE window, maintain the running peak value P.
   Allow counter-moves that satisfy BOTH:
     (a) instantaneous amplitude against trend ≤ δ_amp_max (series-specific), and
     (b) consecutive length of that counter-move ≤ δ_len_max (3 points).
   Update P whenever a new high is observed. (Mirror with trough for DECREASE.)

3) Stop & split (start next window at the reversal extremum) when ANY holds:
   S1. Counter-move violates tolerance: amplitude > δ_amp_max OR lasts > δ_len_max.
   S2. Plateau: over the last L = 5 points, net gain in the trend direction ≤ 0
       AND the average per-point gain < series threshold (0.01 for mutations,
       10 µg/g for calprotectin).
   S3. Persistent reversal: robust slope opposite to the window’s direction for
       ≥ 3 consecutive points (use if you prefer a stricter reversal criterion).

4) Minimum to KEEP a window:
   - MUTATIONS: rely on your existing net-change gate (e.g., freq_range ≥ 0.60).
     No explicit minimum duration imposed here.
   - CALPROTECTIN: keep any window that passes the rules above; amplitude naturally
     clears ≥ 50 µg/g when counter-moves are limited to ≤ 50 µg/g.

5) Merge (prevent over-segmentation):
   If two adjacent windows have the SAME trend direction and the gap between
   them is < 3 weeks AND the gap's counter-move amplitude stays within the
   series-specific tolerance (≤ δ_amp_max), MERGE them into one longer window.

Correlation between mutation and calprotectin windows
-----------------------------------------------------
- Time units: weeks for everything (lags already converted by caller).
- For each calprotectin window, define an expanded “influence” interval by
  adding a pre-event lag and post-event effect (both in weeks).
- A mutation window is considered correlated if its [start,end] interval
  overlaps the expanded calprotectin interval.
- If multiple calprotectin windows overlap a mutation window and are mostly
  non-overlapping with each other, REPORT THEM ALL (multi-event evidence).
  When you need a single link, choose the calprotectin window whose midpoint
  is closest to the mutation window midpoint; break ties by larger amplitude.

Reporting
---------
For EVERY kept window, record:
  series_type (mutation/calprotectin),
  direction (increase/decrease),
  start_time, end_time, peak_time (extreme within the window),
  start_value, end_value, peak_value,
  amplitude = |end − start|,
  length_weeks = end − start,
  split_reason (S1/S2/S3), merged_from (IDs) if applicable.

Notes & Rationale
-----------------
- High δ_amp_max for mutations (0.60) is intentional to absorb volatility and
  preserve long sweeps; the plateau rule (L = 5 with 1%/pt) prevents windows
  from drifting forever with tiny gains.
- Calprotectin uses absolute units (50 µg/g tolerance) so windows reflect
  clinically meaningful swings; 10 µg/g per-point plateau is ~20% of the
  50 µg/g absolute threshold.
- The merge rule (< 4 weeks gap) captures short interruptions without losing
  a biologically continuous episode.
- Multi-event reporting is ON: if both series go up together, then weeks later
  both go down together, both episodes are reported—useful corroborative signal.
===============================================================================
"""

import argparse
import os
import pandas as pd
import numpy as np
import logging
from typing import Tuple, List, Dict

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# -----------------------------
# Piecewise window extraction (both mutation and calprotectin)
# -----------------------------

def calculate_3point_slope(values: np.ndarray, idx: int) -> float:
    """Calculate 3-point slope at index idx."""
    if idx < 1 or idx >= len(values) - 1:
        return 0.0
    return (values[idx + 1] - values[idx - 1]) / 2.0

def extract_piecewise_windows(
    times: np.ndarray,
    values: np.ndarray,
    series_type: str,
    amp_max: float,
    len_max: int = 3,
    plateau_length: int = 5,
    plateau_threshold: float = 0.01,
) -> List[Dict]:
    """
    Extract piecewise windows following the specification.
    
    Args:
        times: Time array (weeks)
        values: Value array
        series_type: 'mutation' or 'calprotectin'
        amp_max: Maximum counter-move amplitude tolerance
        len_max: Maximum counter-move length (points)
        plateau_length: Points to check for plateau (L=5)
        plateau_threshold: Per-point gain threshold for plateau detection
    """
    if len(values) < 2:
        return []
    
    windows = []
    n = len(values)
    i = 0
    
    while i < n - 1:
        # Start from current position - be more permissive about starting windows
        start_idx = i
        
        # If we're at the very beginning, start there regardless of slope
        if start_idx == 0:
            pass  # Start from beginning
        else:
            # For other positions, look for a clear trend
            while start_idx < n - 1:
                if start_idx < 1 or start_idx >= n - 1:
                    break
                slope = calculate_3point_slope(values, start_idx)
                if abs(slope) > 1e-6:  # Found a clear trend
                    break
                start_idx += 1
            
            # If we couldn't find a clear trend, but we're near the end, start anyway
            if start_idx >= n - 1 and i < n - 2:
                start_idx = i  # Start from current position
        
        if start_idx >= n - 1:
            break
            
        # Determine initial direction
        if start_idx == 0:
            # At the beginning, determine direction by looking at the first few points
            if len(values) >= 2:
                is_increase = values[1] > values[0]
            else:
                is_increase = True  # Default to increase
        else:
            slope = calculate_3point_slope(values, start_idx)
            is_increase = slope > 0
        
        # Grow window
        end_idx = start_idx
        peak_val = values[start_idx]
        peak_idx = start_idx
        counter_move_length = 0
        
        for j in range(start_idx + 1, n):
            current_val = values[j]
            
            # Check for large time gaps (should split windows)
            if j > 0:
                time_gap = times[j] - times[j-1]
                if time_gap >= 3.0:  # Large time gap detected
                    break
            
            # Check if this point continues the trend
            if is_increase and current_val >= peak_val:
                # New high - continue trend
                peak_val = current_val
                peak_idx = j
                end_idx = j
                counter_move_length = 0
            elif not is_increase and current_val <= peak_val:
                # New low - continue trend
                peak_val = current_val
                peak_idx = j
                end_idx = j
                counter_move_length = 0
            else:
                # Counter-move detected
                counter_move_length += 1
                
                # Check amplitude tolerance (one consecutive point)
                if is_increase:
                    amp_violation = (peak_val - current_val) > amp_max
                else:
                    amp_violation = (current_val - peak_val) > amp_max
                
                # Check length tolerance
                len_violation = counter_move_length > len_max
                
                # Check if this counter-move would lead to a large drop
                # If the current counter-move is significant, stop the window
                if is_increase and (peak_val - current_val) > amp_max * 0.8:  # 80% of threshold
                    # Significant counter-move detected, stop window
                    break
                elif not is_increase and (current_val - peak_val) > amp_max * 0.8:  # 80% of threshold
                    # Significant counter-move detected, stop window
                    break
                
                if amp_violation or len_violation:
                    # Stop window here
                    break
                else:
                    # Continue with counter-move
                    end_idx = j
        
        # Check plateau condition (only when about to stop naturally, not due to violations)
        net_gain = 0  # Initialize for split_reason logic
        # Only check plateau if we didn't stop due to amplitude/length violation
        if end_idx - start_idx >= plateau_length and counter_move_length <= len_max:
            recent_points = values[end_idx - plateau_length + 1:end_idx + 1]
            if is_increase:
                net_gain = recent_points[-1] - recent_points[0]
                avg_gain = net_gain / plateau_length
            else:
                net_gain = recent_points[0] - recent_points[-1]
                avg_gain = net_gain / plateau_length
            
            if avg_gain < plateau_threshold:
                # Plateau detected - stop window
                end_idx = end_idx - plateau_length + 1
        
        # Create window if it has meaningful length AND meets amplitude threshold
        # Use peak value for amplitude calculation, not end value
        window_amplitude = abs(peak_val - values[start_idx])
        if end_idx > start_idx and window_amplitude >= amp_max:
            # Determine split reason
            if counter_move_length > len_max:
                split_reason = 'S1'
            elif end_idx - start_idx >= plateau_length and net_gain <= 0:
                split_reason = 'S2'
            else:
                split_reason = 'S3'
                
            window = {
                'series_type': series_type,
                'direction': 'increase' if is_increase else 'decrease',
                'start_time': float(times[start_idx]),
                'end_time': float(times[end_idx]),
                'peak_time': float(times[peak_idx]),
                'start_value': float(values[start_idx]),
                'end_value': float(values[end_idx]),
                'peak_value': float(peak_val),
                'amplitude': window_amplitude,
                'length_weeks': float(times[end_idx] - times[start_idx]),
                'split_reason': split_reason,
                'merged_from': []
            }
            windows.append(window)
        
        # Start next window from the peak of the current window, not the end
        i = peak_idx + 1
    
    # Merge adjacent windows with same direction and < 3 week gap
    merged_windows = []
    i = 0
    while i < len(windows):
        current = windows[i]
        j = i + 1
        
        # Look for windows to merge
        while j < len(windows):
            next_window = windows[j]
            
            # Check if mergeable
            gap = next_window['start_time'] - current['end_time']
            same_direction = current['direction'] == next_window['direction']
            small_gap = gap < 3.0
            
            if same_direction and small_gap:
                # Check if gap amplitude is within tolerance
                gap_start_idx = np.where(times == current['end_time'])[0][0]
                gap_end_idx = np.where(times == next_window['start_time'])[0][0]
                gap_values = values[gap_start_idx:gap_end_idx + 1]
                
                if len(gap_values) > 0:
                    gap_amp = max(gap_values) - min(gap_values)
                    if gap_amp <= amp_max:
                        # Calculate what the merged amplitude would be
                        merged_amplitude = abs(next_window['end_value'] - current['start_value'])
                        
                        # Only merge if the merged amplitude meets the threshold
                        if merged_amplitude >= amp_max:
                            # Merge windows
                            current['end_time'] = next_window['end_time']
                            current['end_value'] = next_window['end_value']
                            current['amplitude'] = merged_amplitude
                            current['length_weeks'] = current['end_time'] - current['start_time']
                            current['merged_from'].append(j)
                            
                            # Update peak if next window has higher peak
                            if (current['direction'] == 'increase' and next_window['peak_value'] > current['peak_value']) or \
                               (current['direction'] == 'decrease' and next_window['peak_value'] < current['peak_value']):
                                current['peak_time'] = next_window['peak_time']
                                current['peak_value'] = next_window['peak_value']
                            
                            j += 1
                        else:
                            # Merged amplitude would be too small, stop merging
                            break
                    else:
                        break
                else:
                    break
            else:
                break
        
        merged_windows.append(current)
        i = j
    
    return merged_windows

def get_calprotectin_category(value: float, normal_upper: float, active_lower: float) -> str:
    if value < normal_upper:
        return 'normal'
    if value < active_lower:
        return 'borderline'
    return 'active'

def detect_calprotectin_events(
    calprotectin_data: pd.DataFrame,
    surge_threshold: float = 2.0,  # retained but unused in pairwise mode
    drop_threshold: float = 0.5,   # retained but unused in pairwise mode
    window_size: int = 2,          # retained for logging context only
    normal_upper: float = 50.0,
    active_lower: float = 120.0,
    abs_change_threshold: float = 50.0,
) -> List[Dict]:
    """
    Detect calprotectin surges and drops with biological interpretation.

    Biological categories (mcg/g):
      - Normal:     < normal_upper (default < 50)
      - Borderline: [normal_upper, active_lower) (default 50–120)
      - Active:     >= active_lower (default >= 120)

    Pairwise window criteria (no baseline): for any pair of weeks (i<j):
      - Absolute change: |value_j - value_i| >= abs_change_threshold, OR
      - Category change: category(value_i) != category(value_j).
    Direction: surge if value_j > value_i, drop if value_j < value_i.

    Args:
        calprotectin_data: DataFrame with 'week_num' and 'calprotectin' columns.
        surge_threshold: Optional fold increase to flag surge.
        drop_threshold: Optional fold decrease to flag drop.
        window_size: Half-window for rolling median baseline (2*window_size span).
        normal_upper: Upper bound of normal range (mcg/g).
        active_lower: Lower bound of active inflammation (mcg/g).
        abs_change_threshold: Absolute change threshold (mcg/g) to flag event.

    Returns:
        List of dicts, one per event window, including start/end times/levels,
        categories at both ends, direction, and absolute delta.
    """
    data = calprotectin_data.sort_values('week_num')
    
    if len(data) < 2:
        return []
    
    calprotectin = data['calprotectin'].values
    weeks = data['week_num'].values
    
    logger.info(
        "Calprotectin piecewise detection: normal<%.1f, borderline[%.1f, %.1f), active>=%.1f; amp_max=%.1f",
        normal_upper, normal_upper, active_lower, active_lower, abs_change_threshold
    )

    # Extract windows using piecewise algorithm
    windows = extract_piecewise_windows(
        times=weeks,
        values=calprotectin,
        series_type='calprotectin',
        amp_max=abs_change_threshold,
        len_max=3,
        plateau_length=5,
        plateau_threshold=10.0  # 10 mcg/g per point for calprotectin
    )
    
    # Convert to event format for compatibility
    events = []
    for window in windows:
        events.append({
            'event_type': 'surge' if window['direction'] == 'increase' else 'drop',
            'start_time': window['start_time'],
            'peak_time': window['peak_time'],
            'end_time': window['end_time'],
            'start_level': window['start_value'],
            'peak_level': window['peak_value'],
            'end_level': window['end_value'],
            'fold_change': (window['end_value'] / window['start_value']) if window['start_value'] != 0 else np.nan,
            'baseline': np.nan,
            'baseline_category': get_calprotectin_category(window['start_value'], normal_upper, active_lower),
            'peak_category': get_calprotectin_category(window['peak_value'], normal_upper, active_lower),
            'abs_delta': window['amplitude']
        })

    logger.info("Detected %d calprotectin events", len(events))
    
    return events

def correlate_sweep_calprotectin(
    sweep_events: List[Dict],
    calprotectin_events: List[Dict],
    lag_days: int = 2,
    effect_days: int = 2,
) -> Dict:
    """
    Associate sweep windows with calprotectin event windows.

    A calprotectin event window is expanded by a permissive margin to account for
    anticipated biological lag before sweeps (lag_days, in weeks) and extended downstream
    effects after the event (effect_days, in weeks). Any overlap is considered a positive
    association; directionality is annotated (surge vs drop, increase vs decrease).

    Args:
        sweep_events: List of inferred sweep windows and metrics for one row.
        calprotectin_events: Global list of detected calprotectin events.
        lag_days: Pre-event margin for permissive matching (weeks).
        effect_days: Post-event margin for permissive matching (weeks).

    Returns:
        Dict keyed by correlation type containing matched sweep-to-event associations.
    """
    correlations = {
        'surge_correlated': [],
        'drop_correlated': [],
        'uncorrelated': [],
        'multiple_events': []
    }
    
    for sweep in sweep_events:
        sweep_start = sweep['start_time']
        sweep_end = sweep['end_time']
        
        # Find overlapping calprotectin events
        overlapping_events = []
        
        for cal_event in calprotectin_events:
            # Create permissive window for calprotectin event
            event_start = cal_event['start_time'] - lag_days
            event_end = cal_event['end_time'] + effect_days
            
            # Check for overlap
            if (sweep_start <= event_end and sweep_end >= event_start):
                overlapping_events.append(cal_event)
        
        # Classify the sweep (any correlation is positive signal)
        if len(overlapping_events) == 0:
            correlations['uncorrelated'].append({
                'sweep': sweep,
                'calprotectin_events': []
            })
        elif len(overlapping_events) == 1:
            cal_event = overlapping_events[0]
            if cal_event['event_type'] == 'surge':
                # Determine relationship sign: positive if both increase or both decrease, negative otherwise
                if (sweep['sweep_type'] == 'increase' and cal_event['event_type'] == 'surge') or \
                   (sweep['sweep_type'] == 'decrease' and cal_event['event_type'] == 'drop'):
                    relationship_sign = 'positive'
                else:
                    relationship_sign = 'negative'
                
                correlations['surge_correlated'].append({
                    'sweep': sweep,
                    'calprotectin_event': cal_event,
                    'correlation_direction': f"sweep_{sweep['sweep_type']}_with_calprotectin_{cal_event['event_type']}",
                    'relationship_sign': relationship_sign
                })
            else:
                # Determine relationship sign: positive if both increase or both decrease, negative otherwise
                if (sweep['sweep_type'] == 'increase' and cal_event['event_type'] == 'drop') or \
                   (sweep['sweep_type'] == 'decrease' and cal_event['event_type'] == 'surge'):
                    relationship_sign = 'positive'
                else:
                    relationship_sign = 'negative'
                
                correlations['drop_correlated'].append({
                    'sweep': sweep,
                    'calprotectin_event': cal_event,
                    'correlation_direction': f"sweep_{sweep['sweep_type']}_with_calprotectin_{cal_event['event_type']}",
                    'relationship_sign': relationship_sign
                })
        else:
            correlations['multiple_events'].append({
                'sweep': sweep,
                'calprotectin_events': overlapping_events,
                'correlation_direction': f"sweep_{sweep['sweep_type']}_with_multiple_calprotectin_events"
            })
    
    return correlations

def analyze_mutation_trajectories(
    mutation_data: pd.DataFrame,
    calprotectin_data: pd.DataFrame,
    traj_cols: List[str],
    min_sweep_range: float = 0.6,
    surge_threshold: float = 2.0,
    drop_threshold: float = 0.5,
    lag_days: int = 2,
    effect_days: int = 2,
) -> pd.DataFrame:
    """
    For each mutation row, if `freq_range >= min_sweep_range`, infer a single
    sweep window and direction from trajectory min/max, then correlate it with
    detected calprotectin events using permissive windows.

    Notes on direction inference:
      - If the minimum occurs earlier than the maximum, we call it an increase.
      - Otherwise, we call it a decrease.

    Returns a DataFrame of annotated rows preserving all original columns plus
    added sweep/correlation annotations.
    """
    logger.info(f"Analyzing {len(mutation_data)} mutations for sweep events")
    
    # Filter mutation data first to only include rows with freq_range >= min_sweep_range
    filtered_mutations = mutation_data[
        mutation_data['freq_range'].notna() & 
        (mutation_data['freq_range'] >= min_sweep_range)
    ].copy()
    
    logger.info(f"Filtered to {len(filtered_mutations)} mutations with freq_range >= {min_sweep_range}")
    
    # Detect calprotectin events once
    calprotectin_events = detect_calprotectin_events(
        calprotectin_data,
        surge_threshold=surge_threshold,
        drop_threshold=drop_threshold,
    )
    logger.info(f"Detected {len(calprotectin_events)} calprotectin events")
    
    results = []
    num_processed = 0
    
    for idx, row in filtered_mutations.iterrows():
        
        # Extract trajectory and normalize time index to numeric weeks
        trajectory = row[traj_cols].dropna()
        if len(trajectory) < 2:
            continue
        week_index = pd.to_numeric(trajectory.index, errors='coerce')
        valid_mask = ~np.isnan(week_index)
        if not np.any(valid_mask):
            continue
        traj = pd.Series(trajectory.values[valid_mask], index=week_index[valid_mask]).sort_index()
        if len(traj) < 2:
            continue
        
        # Extract piecewise windows from mutation trajectory
        times = np.array(traj.index, dtype=float)
        values = np.array(traj.values, dtype=float)
        
        mutation_windows = extract_piecewise_windows(
            times=times,
            values=values,
            series_type='mutation',
            amp_max=min_sweep_range,  # Use freq_range threshold as amplitude tolerance
            len_max=3,
            plateau_length=5,
            plateau_threshold=0.01  # 1% per point for mutations
        )
        
        if not mutation_windows:
            continue
        
        # Convert to sweep events format for compatibility
        sweep_events = []
        for window in mutation_windows:
            sweep_events.append({
                'sweep_type': window['direction'],
                'start_time': window['start_time'],
                'peak_time': window['peak_time'],
                'end_time': window['end_time'],
                'start_freq': window['start_value'],
                'peak_freq': window['peak_value'],
                'end_freq': window['end_value'],
                'total_change': window['end_value'] - window['start_value'],
                'window_size': int(max(1, round(window['length_weeks'])))
            })
        
        # Correlate each sweep event with calprotectin events
        correlations = correlate_sweep_calprotectin(sweep_events, calprotectin_events, lag_days, effect_days)
        
        # Add results for each sweep event
        for correlation_type, events in correlations.items():
            for event_data in events:
                sweep = event_data['sweep']
                
                # Start with all original columns from the input data
                result = row.to_dict()
                
                # Add new correlation analysis columns
                result.update({
                    'correlation_type': correlation_type,
                    'sweep_type': sweep['sweep_type'],
                    'sweep_start_time': sweep['start_time'],
                    'sweep_peak_time': sweep['peak_time'],
                    'sweep_end_time': sweep['end_time'],
                    'sweep_window_start': sweep['start_time'],
                    'sweep_window_end': sweep['end_time'],
                    'sweep_window_span_weeks': float(max(0.0, sweep['end_time'] - sweep['start_time'])),
                    'sweep_window_str': f"{sweep['start_time']}-{sweep['end_time']}",
                    'sweep_start_freq': sweep['start_freq'],
                    'sweep_peak_freq': sweep['peak_freq'],
                    'sweep_end_freq': sweep['end_freq'],
                    'sweep_change': sweep['total_change'],
                    'sweep_window_size': sweep['window_size']
                })
                
                # Add correlation direction info
                if 'correlation_direction' in event_data:
                    result['correlation_direction'] = event_data['correlation_direction']
                
                # Add calprotectin event info if available
                if 'calprotectin_event' in event_data:
                    cal_event = event_data['calprotectin_event']
                    result.update({
                        'calprotectin_event_type': cal_event['event_type'],
                        'calprotectin_start_time': cal_event['start_time'],
                        'calprotectin_peak_time': cal_event['peak_time'],
                        'calprotectin_end_time': cal_event['end_time'],
                        'calprotectin_fold_change': cal_event['fold_change'],
                        'calprotectin_baseline': cal_event['baseline']
                    })
                elif 'calprotectin_events' in event_data:
                    # Multiple events case
                    cal_events = event_data['calprotectin_events']
                    result.update({
                        'calprotectin_event_types': [e['event_type'] for e in cal_events],
                        'calprotectin_fold_changes': [e['fold_change'] for e in cal_events],
                        'num_calprotectin_events': len(cal_events)
                    })
                
                results.append(result)
                num_processed += 1

    logger.info(f"Processed {num_processed} sweep events from {len(filtered_mutations)} filtered mutations")
    
    return pd.DataFrame(results)

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Correlate mutation sweep events (freq_range-gated) with biologically "
            "interpreted calprotectin events using permissive temporal windows."
        )
    )
    parser.add_argument("mutation_file", help="analyzed_mutation_types.tsv file")
    parser.add_argument("calprotectin_file", help="CSV with week_num and calprotectin columns")
    parser.add_argument("output_file", help="Output TSV file")
    parser.add_argument("--lag-days", type=int, default=2, help="Weeks before calprotectin event to consider")
    parser.add_argument("--effect-days", type=int, default=2, help="Weeks after calprotectin event to consider")
    parser.add_argument("--surge-threshold", type=float, default=2.0, help="Fold increase for calprotectin surge")
    parser.add_argument("--drop-threshold", type=float, default=0.5, help="Fold decrease for calprotectin drop")
    parser.add_argument("--min-sweep-change", type=float, default=0.6, help="Minimum frequency change (increase or decrease) for sweep")
    parser.add_argument("--min-peak-freq", type=float, default=0.1, help="Minimum peak frequency for sweep")
    
    args = parser.parse_args()
    
    logger.info("=== Starting Sweep-Calprotectin Correlation Analysis ===")
    logger.info(f"Parameters (weeks): lag={args.lag_days}, effect={args.effect_days}")
    logger.info(f"Thresholds: surge={args.surge_threshold}, drop={args.drop_threshold}")
    logger.info("Biology: normal<50, borderline[50,120), active>=120; abs_change>=50 mcg/g")
    
    # Load data
    logger.info(f"Loading mutation data from {args.mutation_file}")
    mutation_data = pd.read_csv(args.mutation_file, sep='\t')
    
    logger.info(f"Loading calprotectin data from {args.calprotectin_file}")
    # Auto-detect delimiter for calprotectin file (handles .csv/.tsv/.txt)
    try:
        calprotectin_data = pd.read_csv(args.calprotectin_file, sep=None, engine='python')
    except Exception:
        # Fallbacks: try tab then comma
        try:
            calprotectin_data = pd.read_csv(args.calprotectin_file, sep='\t')
        except Exception:
            calprotectin_data = pd.read_csv(args.calprotectin_file, sep=',')
    
    # Normalize expected column names
    lower_cols = {c.lower(): c for c in calprotectin_data.columns}
    if 'week_num' not in lower_cols:
        # try alternatives
        for cand in ['week', 'weeknumber', 'week_num', 'weekid']:
            if cand in lower_cols:
                lower_cols['week_num'] = lower_cols[cand]
                break
    if 'calprotectin' not in lower_cols:
        for cand in ['fecal_calprotectin', 'fcal', 'calprot', 'calprotectin']:
            if cand in lower_cols:
                lower_cols['calprotectin'] = lower_cols[cand]
                break
    
    if 'week_num' not in lower_cols or 'calprotectin' not in lower_cols:
        raise ValueError("Calprotectin file must contain columns for week_num and calprotectin (case-insensitive)")
    
    calprotectin_data = calprotectin_data.rename(columns={
        lower_cols['week_num']: 'week_num',
        lower_cols['calprotectin']: 'calprotectin'
    })[['week_num','calprotectin']]
    
    # Detect trajectory columns
    traj_cols = [c for c in mutation_data.columns if c.replace('.', '').replace('e-', '').replace('e+', '').replace('-', '').replace('+', '').replace('E-', '').replace('E+', '').isdigit()]
    logger.info(f"Detected {len(traj_cols)} trajectory columns")
    
    # Analyze correlations
    results = analyze_mutation_trajectories(
        mutation_data, 
        calprotectin_data, 
        traj_cols,
        min_sweep_range=args.min_sweep_change,
        surge_threshold=args.surge_threshold,
        drop_threshold=args.drop_threshold,
        lag_days=args.lag_days,
        effect_days=args.effect_days
    )
    
    # Ensure output directory exists
    out_dir = os.path.dirname(args.output_file)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    
    # Save results
    results.to_csv(args.output_file, sep='\t', index=False)
    logger.info(f"Saved {len(results)} correlation results to {args.output_file}")
    
    # Print summary
    if len(results) > 0:
        summary = results['correlation_type'].value_counts()
        logger.info("Correlation Summary:")
        for corr_type, count in summary.items():
            logger.info(f"  {corr_type}: {count} sweep events")
        
        # Print sweep type summary
        sweep_summary = results['sweep_type'].value_counts()
        logger.info("Sweep Type Summary:")
        for sweep_type, count in sweep_summary.items():
            logger.info(f"  {sweep_type}: {count} sweep events")
        
        # Print correlation direction summary
        if 'correlation_direction' in results.columns:
            corr_dir_summary = results['correlation_direction'].value_counts()
            logger.info("Correlation Direction Summary:")
            for corr_dir, count in corr_dir_summary.items():
                logger.info(f"  {corr_dir}: {count} events")
    else:
        logger.info("No sweep events detected in the data")

if __name__ == "__main__":
    main()
