#!/usr/bin/env python3
"""
sweep_calprotectin_cor_larry.py

Purpose
-------
Rapidly correlate mutation sweep events with host inflammation dynamics measured
by fecal calprotectin. The script keeps all input columns from the mutation file
and appends correlation annotations per detected sweep.

High-level flow
---------------
1) Input mutation table is scanned for timepoint trajectory columns and for an
   existing per-row sweep flag via `freq_range` (range across the trajectory).
   We treat rows with `freq_range >= min_sweep_range` as sweeps.
2) For each swept row, we infer sweep direction (increase/decrease) and window
   using the min/max of the trajectory and their times.
3) Calprotectin time series is parsed with robust delimiter detection and then
   calprotectin events (surges/drops) are identified using biologically grounded
   thresholds and a rolling baseline.
4) A permissive temporal overlap check associates each sweep with zero, one,
   or multiple calprotectin events, recording the directionality context.

Biological thresholds (calprotectin)
------------------------------------
Values in mcg/g are categorized as follows:
- Normal: < 50
- Borderline: 50–120
- Active: >= 120

An event is called when either a category transition occurs relative to the
rolling baseline category, or when the absolute deviation from baseline is
>= 50 mcg/g. Legacy fold thresholds may also trigger events but are not required.

Outputs
-------
Appends to each matched sweep row:
- correlation_type: uncorrelated / surge_correlated / drop_correlated / multiple_events
- sweep_type: increase / decrease
- sweep window metrics (start/end/peak times and frequencies)
- correlation_direction (e.g., sweep_increase_with_calprotectin_surge)
- calprotectin event details where applicable
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
# Calprotectin event detection (pairwise windows)
# -----------------------------

# -----------------------------
# Sweep event extraction (trajectory-level)
# -----------------------------
def extract_pairwise_sweeps(
    traj: pd.Series,
    min_abs_change_freq: float = 0.6,
) -> List[Dict]:
    """
    Generate all sweep windows from any pair of timepoints (i<j) whose
    absolute frequency change >= min_abs_change_freq. No preference for
    shorter windows; return all windows.

    Direction is increase if freq[j] > freq[i], decrease if freq[j] < freq[i].
    """
    if len(traj) < 2:
        return []

    times = np.asarray(traj.index, dtype=float)
    vals = np.asarray(traj.values, dtype=float)
    n = len(vals)
    sweeps: List[Dict] = []
    for i in range(n - 1):
        for j in range(i + 1, n):
            vi, vj = float(vals[i]), float(vals[j])
            if abs(vj - vi) >= float(min_abs_change_freq):
                ti, tj = float(times[i]), float(times[j])
                if vj > vi:
                    sweep_type = 'increase'
                elif vj < vi:
                    sweep_type = 'decrease'
                else:
                    sweep_type = 'flat'
                sweeps.append({
                    'sweep_type': sweep_type,
                    'start_time': ti,
                    'peak_time': tj,
                    'end_time': tj,
                    'start_freq': vi,
                    'peak_freq': vj,
                    'end_freq': vj,
                    'total_change': vj - vi,
                    'window_size': int(max(1, round(tj - ti + 1)))
                })
    return sweeps

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
    events = []
    data = calprotectin_data.sort_values('week_num')
    
    if len(data) < window_size * 2:
        return events
    
    calprotectin = data['calprotectin'].values
    weeks = data['week_num'].values
    
    # Calculate rolling baseline
    logger.info(
        "Calprotectin pairwise detection: normal<%.1f, borderline[%.1f, %.1f), active>=%.1f; abs_change>=%.1f",
        normal_upper, normal_upper, active_lower, active_lower, abs_change_threshold
    )

    events: List[Dict] = []
    for i in range(len(calprotectin) - 1):
        vi = float(calprotectin[i]); wi = float(weeks[i])
        ci = get_calprotectin_category(vi, normal_upper, active_lower)
        for j in range(i+1, len(calprotectin)):
            vj = float(calprotectin[j]); wj = float(weeks[j])
            cj = get_calprotectin_category(vj, normal_upper, active_lower)
            abs_delta = abs(vj - vi)
            category_change = (ci != cj)
            if abs_delta >= abs_change_threshold or category_change:
                event_type = 'surge' if vj > vi else ('drop' if vj < vi else 'flat')
                events.append({
                    'event_type': event_type,
                    'start_time': wi,
                    'peak_time': wj,
                    'end_time': wj,
                    'start_level': vi,
                    'peak_level': vj,
                    'end_level': vj,
                    'fold_change': (vj / vi) if vi not in (0.0,) else np.nan,
                    'baseline': np.nan,
                    'baseline_category': ci,
                    'peak_category': cj,
                    'abs_delta': abs_delta
                })

    logger.info("Detected %d calprotectin events", len(events))
    
    return events

def correlate_sweep_calprotectin(
    sweep_events: List[Dict],
    calprotectin_events: List[Dict],
    lag_days: int = 7,
    effect_days: int = 14,
) -> Dict:
    """
    Associate sweep windows with calprotectin event windows.

    A calprotectin event window is expanded by a permissive margin to account for
    anticipated biological lag before sweeps (lag_days) and extended downstream
    effects after the event (effect_days). Any overlap is considered a positive
    association; directionality is annotated (surge vs drop, increase vs decrease).

    Args:
        sweep_events: List of inferred sweep windows and metrics for one row.
        calprotectin_events: Global list of detected calprotectin events.
        lag_days: Pre-event margin for permissive matching.
        effect_days: Post-event margin for permissive matching.

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
                correlations['surge_correlated'].append({
                    'sweep': sweep,
                    'calprotectin_event': cal_event,
                    'correlation_direction': f"sweep_{sweep['sweep_type']}_with_calprotectin_{cal_event['event_type']}"
                })
            else:
                correlations['drop_correlated'].append({
                    'sweep': sweep,
                    'calprotectin_event': cal_event,
                    'correlation_direction': f"sweep_{sweep['sweep_type']}_with_calprotectin_{cal_event['event_type']}"
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
    lag_days: int = 7,
    effect_days: int = 14,
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
    
    # Detect calprotectin events once
    calprotectin_events = detect_calprotectin_events(
        calprotectin_data,
        surge_threshold=surge_threshold,
        drop_threshold=drop_threshold,
    )
    
    results = []
    
    num_gated_out = 0
    num_processed = 0
    for idx, row in mutation_data.iterrows():
        # Gate by existing sweep annotation using freq_range
        freq_range_val = row.get('freq_range', np.nan)
        if pd.isna(freq_range_val) or float(freq_range_val) < float(min_sweep_range):
            num_gated_out += 1
            continue
        
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
        
        # Extract all pairwise sweep windows that pass threshold
        sweep_events = extract_pairwise_sweeps(traj, min_abs_change_freq=float(min_sweep_range))
        if not sweep_events:
            continue
        
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

    logger.info(
        "Sweep gating summary: kept %d (freq_range>=%.2f), filtered %d",
        num_processed, min_sweep_range, num_gated_out
    )
    
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
    parser.add_argument("--lag-days", type=int, default=7, help="Days before calprotectin event to consider")
    parser.add_argument("--effect-days", type=int, default=14, help="Days after calprotectin event to consider")
    parser.add_argument("--surge-threshold", type=float, default=2.0, help="Fold increase for calprotectin surge")
    parser.add_argument("--drop-threshold", type=float, default=0.5, help="Fold decrease for calprotectin drop")
    parser.add_argument("--min-sweep-change", type=float, default=0.6, help="Minimum frequency change (increase or decrease) for sweep")
    parser.add_argument("--min-peak-freq", type=float, default=0.1, help="Minimum peak frequency for sweep")
    
    args = parser.parse_args()
    
    logger.info("=== Starting Sweep-Calprotectin Correlation Analysis ===")
    logger.info(f"Parameters: lag_days={args.lag_days}, effect_days={args.effect_days}")
    logger.info(f"Thresholds: surge={args.surge_threshold}, drop={args.drop_threshold}")
    logger.info(
        "Biology: normal<50, borderline[50,120), active>=120; abs_change>=50 mcg/g; lag=%dd, effect=%dd",
        args.lag_days,
        args.effect_days,
    )
    
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
