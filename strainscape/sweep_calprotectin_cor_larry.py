#!/usr/bin/env python3
"""
sweep_calprotectin_cor_larry.py  (refactor v3.0)

Purpose
-------
Correlate mutation sweep events with calprotectin dynamics during the sweep's
active change segment, using (a) slope (µg/g/week) and (b) absolute change (µg/g)
within a padded window. Classification: positive / negative / uncorrelated.

Logic preservation note
-----------------------
This refactor preserves your original rules:

- Mutation windows:
  * piecewise monotone-with-tolerance growth from a start extremum
  * hard amplitude threshold (>= min_sweep_change, default 0.6)
  * plateau rule: avg per-point gain < 0.01 (mutations) over L=5
  * stop on counter-move amplitude > δ_amp_max or length > 3
  * early counter-move stop at > 0.8 * δ_amp_max
  * cancel window if any adjacent gap >= 100 weeks
  * start next window from peak_idx + 1
  * merge adjacent same-direction windows if gap < 3 weeks AND gap amplitude
    ≤ δ_amp_max AND merged amplitude still meets threshold
  * amplitude uses peak_value − start_value (not end − start)

- Calprotectin:
  * same piecewise extraction (amp_max = 50 µg/g; plateau 10 µg/g/point)
  * categories helper retained
  * correlation computed on the calprotectin series *inside*
    [mutation_active_start − padding, mutation_peak + padding]
    (falls back to end if needed)

- Correlation significance:
  * slope threshold |slope| ≥ 5 µg/g/week OR
  * absolute change threshold |Δ| ≥ 50 µg/g
  * direction: positive if sweep direction matches sign(slope), else negative
  * if < 2 points or not significant ⇒ uncorrelated
"""

from __future__ import annotations

import argparse
import logging
import os
from dataclasses import dataclass, asdict
from typing import List, Dict, Sequence, Optional, Tuple

import numpy as np
import pandas as pd

# -------------------------
# Logging
# -------------------------

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger("sweep_calprotectin_cor_larry")


# -------------------------
# Config & Data classes
# -------------------------

@dataclass(frozen=True)
class SeriesRules:
    """Tolerance/plateau rules used by piecewise extraction."""
    amp_max: float                 # maximum counter-move amplitude tolerated within a window
    len_max: int = 3               # maximum consecutive counter-move length (points)
    plateau_length: int = 5        # number of points to compute plateau avg gain/loss
    plateau_threshold: float = 0.01  # per-point gain threshold (mutations); use 10.0 for calprotectin

@dataclass
class Window:
    """A contiguous monotone-with-tolerance segment."""
    series_type: str               # 'mutation' or 'calprotectin'
    direction: str                 # 'increase' or 'decrease'
    start_time: float
    end_time: float
    peak_time: float
    start_value: float
    end_value: float
    peak_value: float
    amplitude: float               # |peak_value - start_value|
    length_weeks: float
    split_reason: str              # 'S1' (violation), 'S2' (plateau), 'S3' (other)
    merged_from: List[int]         # indices merged into this window (bookkeeping only)

@dataclass
class SweepEvent:
    """Mutation sweep window in correlation-friendly form."""
    sweep_type: str                # 'increase'/'decrease'
    start_time: float
    peak_time: float
    end_time: float
    start_freq: float
    peak_freq: float
    end_freq: float
    total_change: float
    window_size: int
    local_times: List[float]
    local_values: List[float]

@dataclass
class CalprotectinAnalysis:
    n_points: int
    slope: Optional[float]
    abs_change: Optional[float]
    start_value: Optional[float]
    end_value: Optional[float]
    min_value: Optional[float]
    max_value: Optional[float]
    significant_slope: bool
    significant_change: bool
    is_significant: bool

@dataclass(frozen=True)
class CorrelationThresholds:
    slope_threshold: float = 5.0     # µg/g/week
    abs_change_threshold: float = 50.0  # µg/g

# -------------------------
# Utilities
# -------------------------

def _three_point_slope(values: np.ndarray, idx: int) -> float:
    """Robust 3-point slope centered at idx (0 at edges)."""
    if idx < 1 or idx >= len(values) - 1:
        return 0.0
    return (values[idx + 1] - values[idx - 1]) / 2.0

def _has_large_internal_gap(times: np.ndarray, threshold_weeks: float = 100.0) -> bool:
    """True if any adjacent gap within candidate span is ≥ threshold."""
    if len(times) < 2:
        return False
    return np.any(np.diff(times) >= threshold_weeks)

def _nearest_index(times: np.ndarray, t: float) -> int:
    """Index of nearest time to t (no exact-float equality required)."""
    pos = np.searchsorted(times, t)
    if pos == 0:
        return 0
    if pos >= len(times):
        return len(times) - 1
    # choose closer of pos-1/pos
    return pos if abs(times[pos] - t) < abs(times[pos - 1] - t) else pos - 1

def _is_numeric_header(h: str) -> bool:
    """Return True if a column name looks like a numeric week label."""
    try:
        float(h)
        return True
    except Exception:
        return False

# -------------------------
# Piecewise window extraction (shared)
# -------------------------

def extract_piecewise_windows(
    times: np.ndarray,
    values: np.ndarray,
    series_type: str,
    rules: SeriesRules,
) -> List[Window]:
    """
    Extract windows per monotone-with-tolerance rules.

    Implementation notes (logic identical to original):
    - Start at a local extremum where 3-pt slope gives a non-zero sign;
      special-case index 0 by comparing values[1] vs values[0].
    - Grow while respecting counter-move amplitude ≤ amp_max and consecutive
      length ≤ len_max; early-stop if counter-move exceeds 0.8*amp_max.
    - Plateau detection over the last `plateau_length` points using average
      per-point gain (< plateau_threshold).
    - Cancel window if any adjacent gap ≥ 100 weeks within candidate span.
    - Amplitude computed as |peak - start|.
    - Start next window from peak_idx + 1.
    - Merge adjacent same-direction windows when gap < 3 weeks AND gap amplitude
      ≤ amp_max and merged amplitude still meets amp_max.
    """
    # drop NaNs
    m = ~np.isnan(values)
    times, values = times[m], values[m]
    n = len(values)
    if n < 2:
        return []

    out: List[Window] = []
    i = 0

    while i < n - 1:
        # find a start with non-zero slope (or index 0)
        start_idx = i
        if start_idx > 0:
            while start_idx < n - 1 and abs(_three_point_slope(values, start_idx)) <= 1e-6:
                start_idx += 1
            if start_idx >= n - 1:
                break

        # direction
        is_inc = (values[1] > values[0]) if start_idx == 0 else (_three_point_slope(values, start_idx) > 0)

        end_idx = start_idx
        peak_idx = start_idx
        peak_val = values[start_idx]
        counter_len = 0
        split_reason = "S3"
        cancelled = False

        for j in range(start_idx + 1, n):
            # internal large gap check
            if times[j] - times[j - 1] >= 100.0:
                cancelled = True
                break

            v = values[j]
            if (is_inc and v >= peak_val) or ((not is_inc) and v <= peak_val):
                # trend continues; update peak
                peak_val, peak_idx, end_idx, counter_len = v, j, j, 0
                continue

            # counter-move
            counter_len += 1
            drop = (peak_val - v) if is_inc else (v - peak_val)
            amp_violation = drop > rules.amp_max
            len_violation = counter_len > rules.len_max
            early_violation = drop > 0.8 * rules.amp_max

            if early_violation or amp_violation or len_violation:
                split_reason = "S1"
                break
            else:
                end_idx = j  # tolerate small counter-move; keep extending

        # plateau rule over last plateau_length points (only if we didn't break via S1)
        if end_idx - start_idx >= rules.plateau_length and split_reason != "S1":
            seg = values[end_idx - rules.plateau_length + 1 : end_idx + 1]
            if is_inc:
                net_gain = seg[-1] - seg[0]
            else:
                net_gain = seg[0] - seg[-1]
            avg_gain = net_gain / rules.plateau_length
            if avg_gain < rules.plateau_threshold:
                end_idx = end_idx - rules.plateau_length + 1
                split_reason = "S2"

        amp = abs(peak_val - values[start_idx])
        if (not cancelled) and end_idx > start_idx and amp >= rules.amp_max:
            out.append(
                Window(
                    series_type=series_type,
                    direction="increase" if is_inc else "decrease",
                    start_time=float(times[start_idx]),
                    end_time=float(times[end_idx]),
                    peak_time=float(times[peak_idx]),
                    start_value=float(values[start_idx]),
                    end_value=float(values[end_idx]),
                    peak_value=float(peak_val),
                    amplitude=float(amp),
                    length_weeks=float(times[end_idx] - times[start_idx]),
                    split_reason=split_reason,
                    merged_from=[],
                )
            )

        # advance from peak to allow immediate reversal detection
        i = max(peak_idx, i + 1)

    # Merge adjacent windows (same direction, <3 weeks gap, gap amp ≤ amp_max, merged amp ≥ amp_max)
    if not out:
        return out

    merged: List[Window] = []
    w = list(out)  # copy
    base = w[0]
    k = 1
    while k < len(w):
        nxt = w[k]
        if base.direction == nxt.direction and (nxt.start_time - base.end_time) < 3.0:
            # estimate gap amplitude from raw series by nearest indices
            gap_start_idx = _nearest_index(times, base.end_time)
            gap_end_idx = _nearest_index(times, nxt.start_time)
            if gap_end_idx < gap_start_idx:
                gap_start_idx, gap_end_idx = gap_end_idx, gap_start_idx
            gap_vals = values[gap_start_idx : gap_end_idx + 1]
            gap_amp = float(np.max(gap_vals) - np.min(gap_vals)) if len(gap_vals) else 0.0

            if gap_amp <= rules.amp_max:
                merged_amp = abs(nxt.end_value - base.start_value)
                if merged_amp >= rules.amp_max:
                    # apply merge
                    base.end_time = nxt.end_time
                    base.end_value = nxt.end_value
                    base.length_weeks = base.end_time - base.start_time
                    base.amplitude = merged_amp
                    if (base.direction == "increase" and nxt.peak_value > base.peak_value) or \
                       (base.direction == "decrease" and nxt.peak_value < base.peak_value):
                        base.peak_time = nxt.peak_time
                        base.peak_value = nxt.peak_value
                    base.merged_from.append(k)
                    k += 1
                    continue
        # cannot merge; finalize base
        merged.append(base)
        base = nxt
        k += 1
    merged.append(base)
    return merged


# -------------------------
# Calprotectin helpers
# -------------------------

def get_calprotectin_category(value: float, normal_upper: float, active_lower: float) -> str:
    if value < normal_upper:
        return "normal"
    if value < active_lower:
        return "borderline"
    return "active"

def detect_calprotectin_events(
    calprotectin_data: pd.DataFrame,
    abs_change_threshold: float = 50.0,
    normal_upper: float = 50.0,
    active_lower: float = 120.0,
) -> List[Dict]:
    """
    Detect calprotectin windows via piecewise rules (amp_max=abs_change_threshold,
    plateau_threshold=10 µg/g/point). Returned in a compact dict form for logging/QA.
    """
    df = calprotectin_data.sort_values("week_num")
    if len(df) < 2:
        return []

    weeks = df["week_num"].to_numpy(dtype=float)
    calp  = df["calprotectin"].to_numpy(dtype=float)

    rules = SeriesRules(
        amp_max=abs_change_threshold,
        len_max=3,
        plateau_length=5,
        plateau_threshold=10.0,
    )
    windows = extract_piecewise_windows(weeks, calp, "calprotectin", rules)

    events: List[Dict] = []
    for w in windows:
        events.append(
            {
                "event_type": "surge" if w.direction == "increase" else "drop",
                "start_time": w.start_time,
                "peak_time": w.peak_time,
                "end_time": w.end_time,
                "start_level": w.start_value,
                "peak_level": w.peak_value,
                "end_level": w.end_value,
                "fold_change": (w.end_value / w.start_value) if w.start_value != 0 else np.nan,
                "baseline": np.nan,
                "baseline_category": get_calprotectin_category(w.start_value, normal_upper, active_lower),
                "peak_category": get_calprotectin_category(w.peak_value, normal_upper, active_lower),
                "abs_delta": w.amplitude,
            }
        )
    logger.info("Detected %d calprotectin events", len(events))
    return events


# -------------------------
# Correlation
# -------------------------

def _active_mutation_segment(sweep: SweepEvent) -> Tuple[float, float]:
    """Return (start, end) for the active change segment (start→peak), falling back to end if needed."""
    start = sweep.start_time
    end = sweep.peak_time if sweep.peak_time is not None else sweep.end_time
    if end < start:  # sparse edge case
        end = sweep.end_time
    return start, end

def _analyze_calprotectin_window(
    calprotectin_data: pd.DataFrame,
    t0: float,
    t1: float,
    padding_weeks: float,
    thr: CorrelationThresholds,
) -> CalprotectinAnalysis:
    """Compute slope and absolute change within [t0 - pad, t1 + pad]."""
    lo, hi = t0 - padding_weeks, t1 + padding_weeks
    window = calprotectin_data[(calprotectin_data["week_num"] >= lo) & (calprotectin_data["week_num"] <= hi)].copy()
    if len(window) < 2:
        return CalprotectinAnalysis(0, None, None, None, None, None, None, False, False, False)

    window.sort_values("week_num", inplace=True)
    vals = window["calprotectin"].to_numpy(dtype=float)
    wks  = window["week_num"].to_numpy(dtype=float)

    start_val, end_val = float(vals[0]), float(vals[-1])
    dt = float(wks[-1] - wks[0]) if wks[-1] >= wks[0] else 0.0
    slope = (end_val - start_val) / dt if dt > 0 else 0.0
    abs_change = float(np.max(vals) - np.min(vals))

    sig_slope = abs(slope) >= thr.slope_threshold
    sig_change = abs_change >= thr.abs_change_threshold
    return CalprotectinAnalysis(
        n_points=len(window),
        slope=slope,
        abs_change=abs_change,
        start_value=start_val,
        end_value=end_val,
        min_value=float(np.min(vals)),
        max_value=float(np.max(vals)),
        significant_slope=sig_slope,
        significant_change=sig_change,
        is_significant=(sig_slope or sig_change),
    )

def correlate_sweep_calprotectin(
    sweeps: List[SweepEvent],
    calprotectin_data: pd.DataFrame,
    padding_weeks: float,
    thr: CorrelationThresholds,
) -> Dict[str, List[Dict]]:
    """Return dict with keys: positive_correlated, negative_correlated, uncorrelated."""
    buckets: Dict[str, List[Dict]] = {"positive_correlated": [], "negative_correlated": [], "uncorrelated": []}

    for sw in sweeps:
        t0, t1 = _active_mutation_segment(sw)
        cal = _analyze_calprotectin_window(calprotectin_data, t0, t1, padding_weeks, thr)

        if not cal.is_significant or cal.n_points < 2:
            buckets["uncorrelated"].append({"sweep": sw.__dict__, "calprotectin_analysis": asdict(cal)})
            continue

        sweep_dir = 1 if sw.sweep_type == "increase" else -1
        cal_dir = 1 if (cal.slope or 0.0) > 0 else -1
        if sweep_dir * cal_dir > 0:
            key, sign = "positive_correlated", "positive"
        else:
            key, sign = "negative_correlated", "negative"

        buckets[key].append(
            {
                "sweep": sw.__dict__,
                "calprotectin_analysis": asdict(cal),
                "correlation_direction": f"sweep_{sw.sweep_type}_with_calprotectin_{'increase' if cal_dir > 0 else 'decrease'}",
                "relationship_sign": sign,
            }
        )
    return buckets


# -------------------------
# Top-level analysis
# -------------------------

def analyze_mutation_trajectories(
    mutation_data: pd.DataFrame,
    calprotectin_data: pd.DataFrame,
    traj_cols: Sequence[str],
    min_sweep_range: float = 0.6,
    padding_weeks: float = 3.0,
    slope_threshold: float = 5.0,
    abs_change_threshold: float = 50.0,
) -> pd.DataFrame:
    """Filter by freq_range and emit one row per sweep correlation event."""
    logger.info("Analyzing %d mutations for sweep events", len(mutation_data))

    m = mutation_data["freq_range"].notna() & (mutation_data["freq_range"] >= min_sweep_range)
    filt = mutation_data.loc[m].copy()
    logger.info("Filtered to %d mutations with freq_range >= %.3f", len(filt), min_sweep_range)

    thr = CorrelationThresholds(slope_threshold, abs_change_threshold)
    mut_rules = SeriesRules(amp_max=min_sweep_range, len_max=3, plateau_length=5, plateau_threshold=0.01)

    all_rows: List[Dict] = []
    processed = 0

    for _, row in filt.iterrows():
        # build trajectory (index are week numbers in column headers)
        traj_series = row[list(traj_cols)].dropna()
        if len(traj_series) < 2:
            continue

        # turn column names into numeric 'week' and sort
        week_idx = pd.to_numeric(traj_series.index, errors="coerce")
        valid = ~np.isnan(week_idx)
        if not valid.any():
            continue
        traj = pd.Series(traj_series.values[valid], index=week_idx[valid]).sort_index()

        times = traj.index.to_numpy(dtype=float)
        vals = traj.to_numpy(dtype=float)

        windows = extract_piecewise_windows(times, vals, "mutation", mut_rules)
        if not windows:
            continue

        # convert to SweepEvent
        sweeps: List[SweepEvent] = []
        for w in windows:
            mask = (times >= w.start_time) & (times <= w.end_time)
            local_times = times[mask].tolist()
            local_vals = vals[mask].tolist()
            sweeps.append(
                SweepEvent(
                    sweep_type=w.direction,
                    start_time=w.start_time,
                    peak_time=w.peak_time,
                    end_time=w.end_time,
                    start_freq=w.start_value,
                    peak_freq=w.peak_value,
                    end_freq=w.end_value,
                    total_change=w.end_value - w.start_value,
                    window_size=int(max(1, round(w.length_weeks))),
                    local_times=local_times,
                    local_values=local_vals,
                )
            )

        buckets = correlate_sweep_calprotectin(sweeps, calprotectin_data, padding_weeks, thr)

        # materialize rows
        for corr_type, events in buckets.items():
            for ev in events:
                sw = ev["sweep"]
                out = row.to_dict()
                out.update(
                    {
                        "correlation_type": corr_type,
                        "sweep_type": sw["sweep_type"],
                        "sweep_start_time": sw["start_time"],
                        "sweep_peak_time": sw["peak_time"],
                        "sweep_end_time": sw["end_time"],
                        "sweep_window_start": sw["start_time"],
                        "sweep_window_end": sw["end_time"],
                        "sweep_window_span_weeks": float(max(0.0, sw["end_time"] - sw["start_time"])),
                        "sweep_window_str": f"{sw['start_time']}-{sw['end_time']}",
                        "sweep_start_freq": sw["start_freq"],
                        "sweep_peak_freq": sw["peak_freq"],
                        "sweep_end_freq": sw["end_freq"],
                        "sweep_change": sw["total_change"],
                        "sweep_window_size": sw["window_size"],
                    }
                )
                cal = ev.get("calprotectin_analysis", {})
                out.update(
                    {
                        "calprotectin_n_points": cal.get("n_points", 0),
                        "calprotectin_slope": cal.get("slope"),
                        "calprotectin_abs_change": cal.get("abs_change"),
                        "calprotectin_start_value": cal.get("start_value"),
                        "calprotectin_end_value": cal.get("end_value"),
                        "calprotectin_min_value": cal.get("min_value"),
                        "calprotectin_max_value": cal.get("max_value"),
                        "calprotectin_significant_slope": cal.get("significant_slope", False),
                        "calprotectin_significant_change": cal.get("significant_change", False),
                        "calprotectin_is_significant": cal.get("is_significant", False),
                    }
                )
                if "correlation_direction" in ev:
                    out["correlation_direction"] = ev["correlation_direction"]
                if "relationship_sign" in ev:
                    out["relationship_sign"] = ev["relationship_sign"]

                all_rows.append(out)
                processed += 1

    logger.info("Processed %d sweep events from %d filtered mutations", processed, len(filt))
    return pd.DataFrame(all_rows)


# -------------------------
# CLI
# -------------------------

def _load_calprotectin(path: str) -> pd.DataFrame:
    """Load and normalize calprotectin CSV/TSV to columns: week_num, calprotectin."""
    # try python engine auto-sep; fall back to tab, then comma
    for kwargs in (dict(sep=None, engine="python"), dict(sep="\t"), dict(sep=",")):
        try:
            df = pd.read_csv(path, **kwargs)
            break
        except Exception:
            df = None
    if df is None:
        raise ValueError(f"Failed to read calprotectin file: {path}")

    # normalize columns (case-insensitive)
    cols = {c.lower(): c for c in df.columns}
    wk = next((cols[k] for k in ("week_num", "week", "weeknumber", "weekid") if k in cols), None)
    cp = next((cols[k] for k in ("calprotectin", "fecal_calprotectin", "fcal", "calprot") if k in cols), None)
    if wk is None or cp is None:
        raise ValueError("Calprotectin file must contain columns like week_num and calprotectin (case-insensitive).")

    out = df.rename(columns={wk: "week_num", cp: "calprotectin"})[["week_num", "calprotectin"]].copy()
    out["week_num"] = pd.to_numeric(out["week_num"], errors="coerce")
    out["calprotectin"] = pd.to_numeric(out["calprotectin"], errors="coerce")
    return out.dropna(subset=["week_num", "calprotectin"])

def _detect_traj_columns(mutation_df: pd.DataFrame) -> List[str]:
    """Columns whose headers look numeric serve as week-indexed trajectories."""
    traj = [c for c in mutation_df.columns if _is_numeric_header(str(c))]
    logger.info("Detected %d trajectory columns", len(traj))
    return traj

def main():
    p = argparse.ArgumentParser(
        description=(
            "Correlate mutation sweep events (freq_range-gated) with calprotectin dynamics "
            "during sweep timepoints ± padding. Logic matches v3.0 rules."
        )
    )
    p.add_argument("mutation_file", help="analyzed_mutation_types.tsv")
    p.add_argument("calprotectin_file", help="CSV/TSV with columns week_num, calprotectin")
    p.add_argument("output_file", help="Output TSV")
    p.add_argument("--padding-weeks", type=float, default=3.0, help="Padding around mutation active segment")
    p.add_argument("--slope-threshold", type=float, default=5.0, help="Calprotectin |slope| threshold (µg/g/week)")
    p.add_argument("--abs-change-threshold", type=float, default=50.0, help="Calprotectin |Δ| threshold (µg/g)")
    p.add_argument("--min-sweep-change", type=float, default=0.6, help="Mutation sweep hard amplitude threshold")
    _ = p.add_argument("--min-peak-freq", type=float, default=0.1, help="(Unused in v3.0; retained for CLI compat)")

    args = p.parse_args()

    logger.info("=== Sweep–Calprotectin Correlation (v3.0 logic, refactor) ===")
    logger.info("Params: padding=%.2f, slope>=%.2f, absΔ>=%.2f, sweep_amp>=%.2f",
                args.padding_weeks, args.slope_threshold, args.abs_change_threshold, args.min_sweep_change)

    mut = pd.read_csv(args.mutation_file, sep="\t")
    cal = _load_calprotectin(args.calprotectin_file)

    traj_cols = _detect_traj_columns(mut)

    res = analyze_mutation_trajectories(
        mut,
        cal,
        traj_cols,
        min_sweep_range=args.min_sweep_change,
        padding_weeks=args.padding_weeks,
        slope_threshold=args.slope_threshold,
        abs_change_threshold=args.abs_change_threshold,
    )

    out_dir = os.path.dirname(args.output_file)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    res.to_csv(args.output_file, sep="\t", index=False)
    logger.info("Saved %d correlation rows to %s", len(res), args.output_file)

    if len(res):
        logger.info("Correlation Summary:\n%s", res["correlation_type"].value_counts().to_string())
        logger.info("Sweep Type Summary:\n%s", res["sweep_type"].value_counts().to_string())
        if "correlation_direction" in res.columns:
            logger.info("Direction Summary:\n%s", res["correlation_direction"].value_counts().to_string())
    else:
        logger.info("No sweep events detected.")
        
if __name__ == "__main__":
    main()
