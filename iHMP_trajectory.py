import os
import numpy as np
import pandas as pd

def get_value(row):
    """
    Retrieve frequency corresponding to 'new_base'.
    """
    nuc = row['new_base']
    valid_bases = ['A', 'T', 'C', 'G']
    return row[nuc] if nuc in valid_bases and nuc in row else np.nan

def prepare_trajectory(df):
    """
    Transform multi-model output DataFrame into a wide-format trajectory table.
    """
    df = df.sort_values(['scaffold', 'position', 'week_num']).reset_index(drop=True)
    df['Type'] = df['ref_base'].astype(str) + '>' + df['new_base'].astype(str)
    df['Value'] = df.apply(get_value, axis=1)
    model_cols = ['OLS_slope','OLS_pvalue','OLS_fit']
    model_cols = [c for c in model_cols if c in df.columns]
    pivot_index = ['scaffold', 'position', 'ref_base', 'new_base', 'Type']
    pivot_index = list(dict.fromkeys(pivot_index))
    trajectory_wide = df.pivot_table(index=pivot_index, columns='week_num', values='Value', aggfunc='first').reset_index()
    df_subset = df[['scaffold', 'position', 'ref_base', 'new_base', 'Type'] + model_cols].drop_duplicates().reset_index(drop=True)
    trajectory_wide = pd.merge(trajectory_wide, df_subset, how="left", on=['scaffold', 'position', 'ref_base', 'new_base', 'Type'])
    trajectory_wide.rename(columns={'scaffold': 'Chromosome', 'position': 'Position'}, inplace=True)
    week_columns = sorted(df['week_num'].unique())
    final_cols = ['Trajectory', 'Chromosome', 'Position', 'ref_base', 'new_base', 'Type'] + model_cols + week_columns
    final_cols = [c for c in final_cols if c in trajectory_wide.columns]
    trajectory_wide = trajectory_wide[final_cols]
    return trajectory_wide
