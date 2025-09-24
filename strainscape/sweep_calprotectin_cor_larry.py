#!/usr/bin/env python3
"""
sweep_calprotectin_cor_larry.py

Fast correlation analysis between mutation sweep events and calprotectin dynamics.
- Identifies sweep events from mutation trajectories
- Detects calprotectin surges/drops from metadata
- Correlates sweep timing with calprotectin events using permissive windows
- Accounts for biological lag times and extended effects
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

# (detect_sweep_events and merge_overlapping_sweeps were removed; we now gate by freq_range)

def detect_calprotectin_events(calprotectin_data: pd.DataFrame,
                              surge_threshold: float = 2.0,
                              drop_threshold: float = 0.5,
                              window_size: int = 2) -> List[Dict]:
    """
    Detect calprotectin surges and drops.
    
    Args:
        calprotectin_data: DataFrame with 'week_num' and 'calprotectin' columns
        surge_threshold: Fold increase to define a surge
        drop_threshold: Fold decrease to define a drop
        window_size: Window for calculating baseline
    """
    events = []
    data = calprotectin_data.sort_values('week_num')
    
    if len(data) < window_size * 2:
        return events
    
    calprotectin = data['calprotectin'].values
    weeks = data['week_num'].values
    
    # Calculate rolling baseline
    baseline = pd.Series(calprotectin, index=weeks).rolling(window=window_size*2, center=True).median().values
    
    for i in range(window_size, len(calprotectin) - window_size):
        current_val = calprotectin[i]
        current_baseline = baseline[i]
        
        if pd.isna(current_baseline) or current_baseline == 0:
            continue
            
        fold_change = current_val / current_baseline
        
        if fold_change >= surge_threshold:
            # Find the extent of the surge
            start_idx = i
            end_idx = i
            
            # Find start
            for j in range(i-1, -1, -1):
                if calprotectin[j] <= current_baseline * 1.2:  # 20% above baseline
                    start_idx = j
                    break
            
            # Find end
            for j in range(i+1, len(calprotectin)):
                if calprotectin[j] <= current_baseline * 1.2:
                    end_idx = j
                    break
            
            events.append({
                'event_type': 'surge',
                'start_time': weeks[start_idx],
                'peak_time': weeks[i],
                'end_time': weeks[end_idx],
                'start_level': calprotectin[start_idx],
                'peak_level': current_val,
                'end_level': calprotectin[end_idx],
                'fold_change': fold_change,
                'baseline': current_baseline
            })
            
        elif fold_change <= drop_threshold:
            # Find the extent of the drop
            start_idx = i
            end_idx = i
            
            # Find start
            for j in range(i-1, -1, -1):
                if calprotectin[j] >= current_baseline * 0.8:  # 20% below baseline
                    start_idx = j
                    break
            
            # Find end
            for j in range(i+1, len(calprotectin)):
                if calprotectin[j] >= current_baseline * 0.8:
                    end_idx = j
                    break
            
            events.append({
                'event_type': 'drop',
                'start_time': weeks[start_idx],
                'peak_time': weeks[i],
                'end_time': weeks[end_idx],
                'start_level': calprotectin[start_idx],
                'peak_level': current_val,
                'end_level': calprotectin[end_idx],
                'fold_change': fold_change,
                'baseline': current_baseline
            })
    
    return events

def correlate_sweep_calprotectin(sweep_events: List[Dict],
                                calprotectin_events: List[Dict],
                                lag_days: int = 7,
                                effect_days: int = 14) -> Dict:
    """
    Correlate sweep events with calprotectin events using permissive windows.
    
    Args:
        lag_days: Days before calprotectin event to consider (biological lag)
        effect_days: Days after calprotectin event to consider (extended effects)
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

def analyze_mutation_trajectories(mutation_data: pd.DataFrame,
                                calprotectin_data: pd.DataFrame,
                                traj_cols: List[str],
                                min_sweep_range: float = 0.6,
                                surge_threshold: float = 2.0,
                                drop_threshold: float = 0.5,
                                lag_days: int = 7,
                                effect_days: int = 14) -> pd.DataFrame:
    """
    Analyze all mutation trajectories for sweep events and correlate with calprotectin.
    """
    logger.info(f"Analyzing {len(mutation_data)} mutations for sweep events")
    
    # Detect calprotectin events once
    calprotectin_events = detect_calprotectin_events(calprotectin_data, surge_threshold, drop_threshold)
    logger.info(f"Detected {len(calprotectin_events)} calprotectin events")
    
    results = []
    
    for idx, row in mutation_data.iterrows():
        # Gate by existing sweep annotation using freq_range
        freq_range_val = row.get('freq_range', np.nan)
        if pd.isna(freq_range_val) or float(freq_range_val) < float(min_sweep_range):
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
        
        # Infer sweep direction and window from min/max
        min_idx = int(np.argmin(traj.values))
        max_idx = int(np.argmax(traj.values))
        t_min = float(traj.index[min_idx])
        t_max = float(traj.index[max_idx])
        v_min = float(traj.values[min_idx])
        v_max = float(traj.values[max_idx])
        
        if t_min < t_max and v_max - v_min > 0:
            sweep_type = 'increase'
            start_time, end_time = t_min, t_max
            start_freq, end_freq = v_min, v_max
            peak_time = t_max
            peak_freq = v_max
            total_change = v_max - v_min
        else:
            sweep_type = 'decrease'
            start_time, end_time = t_max, t_min
            start_freq, end_freq = v_max, v_min
            peak_time = t_min
            peak_freq = v_min
            total_change = v_min - v_max  # negative
        
        sweep_events = [{
            'sweep_type': sweep_type,
            'start_time': start_time,
            'peak_time': peak_time,
            'end_time': end_time,
            'start_freq': start_freq,
            'peak_freq': peak_freq,
            'end_freq': end_freq,
            'total_change': total_change,
            'window_size': int(round(end_time - start_time + 1)) if end_time >= start_time else 1
        }]
        
        # Correlate with calprotectin events
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
    
    return pd.DataFrame(results)

def main():
    parser = argparse.ArgumentParser(description="Correlate mutation sweep events with calprotectin dynamics")
    parser.add_argument("mutation_file", help="analyzed_mutation_types.tsv file")
    parser.add_argument("calprotectin_file", help="CSV with week_num and calprotectin columns")
    parser.add_argument("output_file", help="Output TSV file")
    parser.add_argument("--lag-days", type=int, default=7, help="Days before calprotectin event to consider")
    parser.add_argument("--effect-days", type=int, default=14, help="Days after calprotectin event to consider")
    parser.add_argument("--surge-threshold", type=float, default=2.0, help="Fold increase for calprotectin surge")
    parser.add_argument("--drop-threshold", type=float, default=0.5, help="Fold decrease for calprotectin drop")
    parser.add_argument("--min-sweep-change", type=float, default=0.2, help="Minimum frequency change (increase or decrease) for sweep")
    parser.add_argument("--min-peak-freq", type=float, default=0.1, help="Minimum peak frequency for sweep")
    
    args = parser.parse_args()
    
    logger.info("=== Starting Sweep-Calprotectin Correlation Analysis ===")
    logger.info(f"Parameters: lag_days={args.lag_days}, effect_days={args.effect_days}")
    logger.info(f"Thresholds: surge={args.surge_threshold}, drop={args.drop_threshold}")
    
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
