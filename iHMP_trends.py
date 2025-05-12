import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize
from scipy.stats import f, beta, t as tdist
from joblib import Parallel, delayed
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import pdist, squareform
from scipy.stats import zscore

def fast_ols_single(x, y):
    """
    Compute OLS slope, two-sided p-value, and R² for a single SNP's trajectory,
    ignoring NaNs.
    """
    n = len(x)
    if n < 2:
        return np.nan, np.nan, np.nan
    xbar = np.mean(x)
    ybar = np.mean(y)
    ss_xx = np.sum((x - xbar)**2)
    ss_yy = np.sum((y - ybar)**2)
    ss_xy = np.sum((x - xbar)*(y - ybar))
    if ss_xx < 1e-12:
        return np.nan, np.nan, np.nan
    slope = ss_xy / ss_xx
    intercept = ybar - slope * xbar
    y_hat = slope * x + intercept
    rss = np.sum((y - y_hat)**2)
    dof = n - 2
    if dof <= 0:
        return slope, np.nan, np.nan
    s2 = rss / dof
    if s2 <= 0:
        return slope, np.nan, np.nan
    se_slope = np.sqrt(s2 / ss_xx)
    if se_slope < 1e-12:
        return slope, np.nan, np.nan
    t_slope = slope / se_slope
    p_slope = 2 * (1 - tdist.cdf(abs(t_slope), dof))
    r2 = 1 - rss / ss_yy if ss_yy > 1e-12 else np.nan
    return slope, p_slope, r2

def process_one_group(key, subdf, slope_thresh=0.0001, p_thresh=0.1, do_logistic_stepwise=True):
    """
    Process a group of SNPs (grouped by scaffold, position, ref_base) to compute
    OLS summary statistics and, if conditions are met, fit logistic and stepwise models.
    """
    scaffold, position, ref_base = key
    nucleotides = [nuc for nuc in ['A', 'T', 'C', 'G'] if nuc != ref_base]
    row_out = {
        'scaffold': scaffold,
        'position': position,
        'new_base': ref_base,
        'OLS_slope': 0,
        'OLS_pvalue': 1,
        'OLS_fit': 0
    }
    for nuc in nucleotides:
        tmp = subdf[['week_num', nuc]].dropna()
        if tmp.empty or tmp[nuc].std() == 0:
            continue
        x, y = tmp['week_num'].values, tmp[nuc].values
        slope, p_val, r2 = fast_ols_single(x, y)
        if r2 > row_out['OLS_fit']:
            row_out.update({
                'OLS_slope': slope,
                'OLS_pvalue': p_val,
                'OLS_fit': r2,
                'new_base': nuc
            })
    return row_out

def calculate_trend(df, n_jobs=1):
    """
    Group the DataFrame by scaffold, position, and ref_base; process each group in parallel.
    """
    grouped = df.groupby(['scaffold', 'position', 'ref_base'], sort=False)
    results_list = Parallel(n_jobs=n_jobs)(
        delayed(process_one_group)(key, subdf) for key, subdf in grouped
    )
    return pd.DataFrame(results_list)
