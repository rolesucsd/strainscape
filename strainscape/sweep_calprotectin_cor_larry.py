#!/usr/bin/env python3
"""
sweep_calprotectin_cor_larry.py

Correlate mutation sweep events with calprotectin dynamics during sweep timepoints.
Uses direct slope and absolute change analysis of calprotectin data within sweep windows.

===============================================================================
New Correlation Approach (v3.0)
===============================================================================

Goal
----
For each mutation sweep event, extract calprotectin data within the sweep timepoints
(± padding), calculate slope and absolute change, and determine if there's a significant
association based on biologically meaningful thresholds.

Key Algorithm Changes (v3.0)
----------------------------
- Direct calprotectin analysis: Extract data during sweep timepoints ± padding_weeks
- Slope-based correlation: Calculate calprotectin slope (mcg/g per week) during sweep
- Absolute change correlation: Calculate total calprotectin change during sweep
- Threshold-based significance: Either absolute slope (|slope|) or absolute change must meet thresholds
- Direction classification: Positive (same direction) vs Negative (opposite direction)

Mutation Sweep Detection (Unchanged)
------------------------------------
- Hard amplitude threshold: ALL mutation windows MUST have sweep_change ≥ 0.6
- Initial direction: Compare values[1] vs values[0] at start_idx=0
- Amplitude calculation: Use |peak_val − values[start_idx]|
- Plateau detection: avg_gain < plateau_threshold (mutations: 0.01 per point)
- Window continuation: Start next window from peak_idx + 1
- Merge guard: Only merge windows that maintain amplitude threshold

Calprotectin Analysis (New)
---------------------------
• Extraction window: calprotectin data during [sweep_start − padding_weeks, sweep_end + padding_weeks]
• Slope (mcg/g per week): (end_value − start_value) / (end_week − start_week)
• Absolute change (mcg/g): |end_value − start_value|
• Significance thresholds (defaults):
  - slope_threshold: |slope| ≥ 5.0 mcg/g/week
  - abs_change_threshold: |Δ| ≥ 50.0 mcg/g
  - is_significant = (|slope| ≥ slope_threshold) OR (|Δ| ≥ abs_change_threshold)
• Correlation classification:
  - positive_correlated: sweep direction matches calprotectin direction (both rise or both fall)
  - negative_correlated: sweep direction opposes calprotectin direction (one rises, other falls)
  - uncorrelated: insufficient points (<2) or not significant

Output Columns
--------------
• All original mutation data columns preserved
• New correlation columns:
  - correlation_type: positive_correlated/negative_correlated/uncorrelated
  - relationship_sign: positive/negative (for correlated events)
  - calprotectin_slope: mcg/g per week during sweep window
  - calprotectin_abs_change: mcg/g absolute change during sweep window
  - calprotectin_significant_slope: Boolean (slope >= threshold)
  - calprotectin_significant_change: Boolean (abs_change >= threshold)
  - calprotectin_is_significant: Boolean (either slope or change significant)
----------------------------------------------------
Given a time-ordered series (raw data, no smoothing):

1) Seed / start:
   Scan left→right. Start a new window at a local extremum where the robust
   slope (e.g., over last 2–3 points) establishes a trend sign.
   SPECIAL CASE: At start_idx=0, determine initial direction by comparing
   values[1] vs values[0] directly (3-point slope returns 0.0 at index 0).

2) Grow (monotone-with-tolerance):
   While growing an INCREASE window, maintain the running peak value P.
   Allow counter-moves that satisfy BOTH:
     (a) instantaneous amplitude against trend ≤ δ_amp_max (series-specific), and
     (b) consecutive length of that counter-move ≤ δ_len_max (3 points).
   Update P whenever a new high is observed. (Mirror with trough for DECREASE.)

3) Stop & split (start next window at the reversal extremum) when ANY holds:
   S1. Counter-move violates tolerance: amplitude > δ_amp_max OR lasts > δ_len_max.
   S2. Plateau: average per-point gain < series threshold (0.01 for mutations,
       10 µg/g for calprotectin).
   S4. Early counter-move detection: significant counter-move > 0.8 * amp_max
       detected during window growth (prevents continuing through large drops).

4) Minimum to KEEP a window:
   - MUTATIONS: HARD THRESHOLD - sweep_change MUST be >= 0.6 (amplitude from start to peak)
   - CALPROTECTIN: keep any window that passes the rules above; amplitude naturally
     clears ≥ 50 µg/g when counter-moves are limited to ≤ 50 µg/g.

5) Merge (prevent over-segmentation):
   If two adjacent windows have the SAME trend direction and the gap between
   them is < 3 weeks AND the gap's counter-move amplitude stays within the
   series-specific tolerance (≤ δ_amp_max), MERGE them into one longer window.
   CRITICAL: Only merge if the merged window maintains the amplitude threshold
   (>= 0.6 for mutations, >= 50 µg/g for calprotectin).

6) Window continuation:
   After creating a window, start the next window from peak_idx + 1 (not end_idx + 1)
   to allow detection of subsequent events (e.g., decrease after increase).

7) Adjacent timepoint gap cap:
   A window may span more than 100 weeks overall, but if any consecutive
   timepoints within the candidate window are separated by ≥ 100 weeks,
   cancel that window (do not keep it). This prevents single-step jumps from
   stitching together unrelated epochs.

Correlation analysis (direct, windowed)
---------------------------------------
- Time units: weeks for everything.
- For each mutation sweep window, analyze calprotectin within the padded
  interval [sweep_start − padding_weeks, sweep_end + padding_weeks].
- Compute slope and absolute change; apply thresholds as above to set
  is_significant, then classify positive/negative by relative directions.
  If not significant or <2 points, classify as uncorrelated.

Reporting
---------
For EVERY kept window, record:
  series_type (mutation/calprotectin),
  direction (increase/decrease),
  start_time, end_time, peak_time (extreme within the window),
  start_value, end_value, peak_value,
  amplitude = |peak_value − start_value| (corrected calculation),
  length_weeks = end − start,
  correlation_type (uncorrelated/surge_correlated/drop_correlated/multiple_events),
  relationship_sign ('positive' or 'negative' correlation direction),
  split_reason (S1/S2/S4), merged_from (IDs) if applicable.

Notes & Rationale
-----------------
- Mutation amp tolerance (δ_amp_max = 0.60) absorbs small drawdowns while the
  plateau rule (avg_gain < 0.01 per point) prevents drifting on noise.
- Calprotectin uses absolute units (50 mcg/g) so windows reflect clinically
  meaningful swings; 10 mcg/g per-point plateau is ~20% of that absolute threshold.
- Merge rule (< 3 weeks gap) captures short interruptions without losing a
  biologically continuous episode; merged amplitude must still meet the threshold.
- Hard amplitude threshold (0.6 for mutations) ensures detected sweep events
  meet biological significance.
- Amplitude uses peak_value − start_value (not end − start) to reflect the
  sweep magnitude.
- Starting next window from peak_idx + 1 enables detecting multiple events
  (e.g., increase followed by decrease).
- relationship_sign records correlation directionality explicitly.
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
    
    # Filter out NA values by creating masks
    valid_mask = ~np.isnan(values)
    if not np.any(valid_mask):
        return []
    
    # Apply masks to get clean data
    times = times[valid_mask]
    values = values[valid_mask]
    
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
        
        cancelled_due_to_gap = False
        for j in range(start_idx + 1, n):
            current_val = values[j]
            # Cancel window if any adjacent timepoint gap within it is ≥ 100 weeks
            if times[j] - times[j - 1] >= 100.0:
                cancelled_due_to_gap = True
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
        if not cancelled_due_to_gap and end_idx > start_idx and window_amplitude >= amp_max:
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
        # But ensure we advance to avoid infinite loops
        i = max(peak_idx, i + 1)
    
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
    calprotectin_data: pd.DataFrame,
    padding_weeks: float = 3.0,
    slope_threshold: float = 5.0,  # mcg/g per week
    abs_change_threshold: float = 50.0,  # mcg/g absolute change
    min_sweep_range: float = 0.6,
) -> Dict:
    """
    Correlate sweep events with calprotectin dynamics during the sweep timepoints.
    
    For each sweep event:
    1. Extract calprotectin data within sweep timepoints ± padding_weeks
    2. Calculate slope and absolute change of calprotectin during this period
    3. Apply thresholds to determine if association is significant
    4. Classify as positive/negative/uncorrelated based on sweep direction vs calprotectin direction
    
    Args:
        sweep_events: List of sweep event dictionaries
        calprotectin_data: DataFrame with week_num and calprotectin columns
        padding_weeks: Weeks to pad around sweep timepoints
        slope_threshold: Minimum slope (mcg/g per week) for significant association
        abs_change_threshold: Minimum absolute change (mcg/g) for significant association
    
    Returns:
        Dictionary with correlation classifications
    """
    correlations = {
        'positive_correlated': [],
        'negative_correlated': [],
        'uncorrelated': []
    }
    
    for sweep in sweep_events:
        sweep_start = sweep['start_time']
        sweep_end = sweep['end_time']

        # Use the mutation's active change segment (start -> peak) for calprotectin analysis
        # Fall back to end_time if peak_time is unavailable
        mutation_change_end = sweep.get('peak_time', sweep_end)
        # Ensure ordering (some edge cases may have peak before start due to sparse data)
        if mutation_change_end < sweep_start:
            mutation_change_end = sweep_end

        # Active slope segment: recent monotonic ramp into the peak
        # This captures where the mutation actually changes (not the full sweep span)
        mutation_change_start = sweep_start
        if 'local_times' in sweep and 'local_values' in sweep:
            lt = np.array(sweep['local_times'], dtype=float)
            lv = np.array(sweep['local_values'], dtype=float)
            try:
                # indices for start and peak within local window
                start_idx = int(np.argmin(np.abs(lt - sweep_start)))
                peak_idx = int(np.argmin(np.abs(lt - mutation_change_end)))
                if start_idx > peak_idx:
                    start_idx, peak_idx = peak_idx, start_idx
                # walk backwards from peak while monotonic in sweep direction
                i = max(start_idx, peak_idx - 1)
                if sweep.get('sweep_type') == 'increase':
                    while i >= start_idx and lv[i] <= lv[i + 1]:
                        i -= 1
                else:
                    while i >= start_idx and lv[i] >= lv[i + 1]:
                        i -= 1
                mutation_change_start = lt[max(start_idx, i + 1)]
            except Exception:
                pass
        
        # Define analysis window with padding around the active mutation change segment
        window_start = mutation_change_start - padding_weeks
        window_end = mutation_change_end + padding_weeks
        
        # Extract calprotectin data within the window
        window_data = calprotectin_data[
            (calprotectin_data['week_num'] >= window_start) & 
            (calprotectin_data['week_num'] <= window_end)
        ].copy()
        
        if len(window_data) < 2:
            # Not enough data points for analysis
            correlations['uncorrelated'].append({
                'sweep': sweep,
                'calprotectin_analysis': {
                    'n_points': len(window_data),
                    'slope': None,
                    'abs_change': None,
                    'start_value': None,
                    'end_value': None,
                    'reason': 'insufficient_data'
                }
            })
            continue
        
        # Sort by time
        window_data = window_data.sort_values('week_num')
        
        # Calculate calprotectin dynamics
        # Find min and max values within the window (regardless of temporal position)
        min_value = window_data['calprotectin'].min()
        max_value = window_data['calprotectin'].max()
        abs_change = abs(max_value - min_value)
        
        # Calculate slope using first and last values for temporal trend
        start_value = window_data['calprotectin'].iloc[0]
        end_value = window_data['calprotectin'].iloc[-1]
        time_diff = window_data['week_num'].iloc[-1] - window_data['week_num'].iloc[0]
        if time_diff > 0:
            slope = (end_value - start_value) / time_diff
        else:
            slope = 0.0
        
        # Determine if association is significant
        significant_slope = abs(slope) >= slope_threshold
        significant_change = abs_change >= abs_change_threshold
        is_significant = significant_slope or significant_change
        
        calprotectin_analysis = {
            'n_points': len(window_data),
            'slope': slope,
            'abs_change': abs_change,
            'start_value': start_value,
            'end_value': end_value,
            'min_value': min_value,
            'max_value': max_value,
            'significant_slope': significant_slope,
            'significant_change': significant_change,
            'is_significant': is_significant
        }
        
        if not is_significant:
            correlations['uncorrelated'].append({
                'sweep': sweep,
                'calprotectin_analysis': calprotectin_analysis
            })
        else:
            # Determine correlation direction
            # Positive: sweep and calprotectin change in same direction
            # Negative: sweep and calprotectin change in opposite directions
            
            sweep_direction = 1 if sweep['sweep_type'] == 'increase' else -1
            calprotectin_direction = 1 if slope > 0 else -1
            
            if sweep_direction * calprotectin_direction > 0:
                correlation_type = 'positive_correlated'
                relationship_sign = 'positive'
            else:
                correlation_type = 'negative_correlated'
                relationship_sign = 'negative'
            
            correlations[correlation_type].append({
                'sweep': sweep,
                'calprotectin_analysis': calprotectin_analysis,
                'correlation_direction': f"sweep_{sweep['sweep_type']}_with_calprotectin_{'increase' if slope > 0 else 'decrease'}",
                'relationship_sign': relationship_sign
            })
    
    return correlations

def analyze_mutation_trajectories(
    mutation_data: pd.DataFrame,
    calprotectin_data: pd.DataFrame,
    traj_cols: List[str],
    min_sweep_range: float = 0.6,
    padding_weeks: float = 3.0,
    slope_threshold: float = 5.0,
    abs_change_threshold: float = 50.0,
) -> pd.DataFrame:
    """
    For each mutation row, if `freq_range >= min_sweep_range`, extract sweep events
    using piecewise window detection, then correlate each sweep with calprotectin
    dynamics during the sweep timepoints ± padding.

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
            # Collect local trajectory within this window for later narrowing
            window_mask = (times >= window['start_time']) & (times <= window['end_time'])
            local_times = times[window_mask]
            local_values = values[window_mask]

            sweep_events.append({
                'sweep_type': window['direction'],
                'start_time': window['start_time'],
                'peak_time': window['peak_time'],
                'end_time': window['end_time'],
                'start_freq': window['start_value'],
                'peak_freq': window['peak_value'],
                'end_freq': window['end_value'],
                'total_change': window['end_value'] - window['start_value'],
                'window_size': int(max(1, round(window['length_weeks']))),
                'local_times': local_times.tolist(),
                'local_values': local_values.tolist()
            })
        
        # Correlate each sweep event with calprotectin dynamics
        correlations = correlate_sweep_calprotectin(
            sweep_events, 
            calprotectin_data, 
            padding_weeks=padding_weeks,
            slope_threshold=slope_threshold,
            abs_change_threshold=abs_change_threshold,
            min_sweep_range=min_sweep_range
        )
        
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
                if 'relationship_sign' in event_data:
                    result['relationship_sign'] = event_data['relationship_sign']
                
                # Add calprotectin analysis info
                if 'calprotectin_analysis' in event_data:
                    cal_analysis = event_data['calprotectin_analysis']
                    result.update({
                        'calprotectin_n_points': cal_analysis.get('n_points', 0),
                        'calprotectin_slope': cal_analysis.get('slope'),
                        'calprotectin_abs_change': cal_analysis.get('abs_change'),
                        'calprotectin_start_value': cal_analysis.get('start_value'),
                        'calprotectin_end_value': cal_analysis.get('end_value'),
                        'calprotectin_min_value': cal_analysis.get('min_value'),
                        'calprotectin_max_value': cal_analysis.get('max_value'),
                        'calprotectin_significant_slope': cal_analysis.get('significant_slope', False),
                        'calprotectin_significant_change': cal_analysis.get('significant_change', False),
                        'calprotectin_is_significant': cal_analysis.get('is_significant', False)
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
    parser.add_argument("--padding-weeks", type=float, default=3.0, help="Weeks to pad around sweep timepoints for calprotectin analysis")
    parser.add_argument("--slope-threshold", type=float, default=5.0, help="Minimum slope (mcg/g per week) for significant calprotectin association")
    parser.add_argument("--abs-change-threshold", type=float, default=50.0, help="Minimum absolute change (mcg/g) for significant calprotectin association")
    parser.add_argument("--min-sweep-change", type=float, default=0.6, help="Minimum frequency change (increase or decrease) for sweep")
    parser.add_argument("--min-peak-freq", type=float, default=0.1, help="Minimum peak frequency for sweep")
    
    args = parser.parse_args()
    
    logger.info("=== Starting Sweep-Calprotectin Correlation Analysis ===")
    logger.info(f"Parameters: padding={args.padding_weeks} weeks")
    logger.info(f"Thresholds: slope>={args.slope_threshold} mcg/g/week, abs_change>={args.abs_change_threshold} mcg/g")
    logger.info("Analysis: Extract calprotectin dynamics during sweep timepoints ± padding")
    
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
        padding_weeks=args.padding_weeks,
        slope_threshold=args.slope_threshold,
        abs_change_threshold=args.abs_change_threshold
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
