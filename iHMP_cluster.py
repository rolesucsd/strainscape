import os
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import zscore, beta
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import pdist, squareform
from sklearn.impute import KNNImputer

def compute_blockwise_corrcoef(X, block_size=1000):
    """
    Compute the Pearson correlation coefficient matrix for a 2D array X in a blockwise fashion.
    """
    n, T = X.shape
    mu = X.mean(axis=1, keepdims=True)
    X_centered = X - mu
    sigma = np.sqrt(np.sum(X_centered**2, axis=1) / (T - 1))
    corr = np.empty((n, n), dtype=float)
    for i in range(0, n, block_size):
        i_end = min(i + block_size, n)
        X_i = X_centered[i:i_end, :]
        sigma_i = sigma[i:i_end]
        for j in range(i, n, block_size):
            j_end = min(j + block_size, n)
            X_j = X_centered[j:j_end, :]
            sigma_j = sigma[j:j_end]
            cov_block = np.dot(X_i, X_j.T) / (T - 1)
            corr_block = cov_block / np.outer(sigma_i, sigma_j)
            corr[i:i_end, j:j_end] = corr_block
            if i != j:
                corr[j:j_end, i:i_end] = corr_block.T
    return corr

def cluster_snps_combined(data, timepoint_cols, correlation_threshold=0.65,
                          raw_weight=0.7, summary_weight=0.3,
                          variability_threshold=0.001):
    """
    Cluster SNPs by combining:
      (a) raw allele frequency trajectories,
      (b) summary features (slope, -log10(p_value), r_squared),
      (c) LD similarity.
    """
    timepoint_cols = [str(col) for col in timepoint_cols]
    data.columns = data.columns.astype(str)

    variability = data[timepoint_cols].astype(float).var(axis=1, ddof=1)
    data = data[variability >= variability_threshold].copy()
    print("Length after variability filtering: ", len(data))

    trajectories = data[timepoint_cols].astype(float).values
    knn_imputer = KNNImputer(n_neighbors=5, weights='uniform')
    trajectories = knn_imputer.fit_transform(trajectories)

    temporal_corr = compute_blockwise_corrcoef(trajectories)
    raw_distance = 1 - temporal_corr

    features = data[['OLS_slope', 'OLS_pvalue', 'OLS_fit']].copy()
    features_norm = features.apply(zscore)
    summary_distance = squareform(pdist(features_norm.values, metric='euclidean'))

    combined_distance = (raw_weight * raw_distance +
                         summary_weight * summary_distance)


    # Agglomerative clustering using the combined distance.
    clustering = AgglomerativeClustering(
        metric='precomputed',  # Updated from 'affinity'
        linkage='average',
        distance_threshold=1 - correlation_threshold,
        n_clusters=None
    )
    labels = clustering.fit_predict(combined_distance)
    data['strain_cluster'] = labels

    cluster_sizes = data['strain_cluster'].value_counts().to_dict()

    data['cluster_size'] = data['strain_cluster'].map(cluster_sizes)

    return data

def transition_density(x_next, x, s, dt):
    """
    Approximate transition density from x to x_next under a Wright-Fisher-like model 
    using a Beta distribution parameterization.
    """
    # Mean frequency after time dt, given selection coefficient s
    m = x + s * x * (1 - x) * dt
    # Variance-like term
    v = x * (1 - x) * dt
    if v <= 0:
        # Avoid invalid Beta parameters
        return 1e-10
    factor = (m * (1 - m)) / v - 1
    alpha_param = m * factor
    beta_param  = (1 - m) * factor
    if alpha_param <= 0 or beta_param <= 0:
        return 1e-10

    return beta.pdf(x_next, alpha_param, beta_param)

def negative_log_likelihood(s, traj, timepoints):
    """
    Compute negative log-likelihood for a cluster's mean frequency trajectory 'traj' 
    given a selection coefficient s and variable time intervals in 'timepoints'.
    
    traj       : 1D array of shape (T,) with the average frequency at each timepoint
    timepoints : 1D array of shape (T,) giving the actual times (non-uniform).
    """
    nll = 0.0
    for i in range(1, len(traj)):
        dt_i = timepoints[i] - timepoints[i - 1]  # variable interval
        p = transition_density(traj[i], traj[i - 1], s, dt_i)
        nll -= np.log(p + 1e-10)  # stable log
    return nll

def infer_selection_for_cluster(traj, timepoints):
    """
    Infer best-fit selection coefficient for one cluster's trajectory 'traj'
    given non-uniform 'timepoints'.
    
    Returns the scalar s_hat.
    """
    def objective(s):
        return negative_log_likelihood(s, traj, timepoints)

    # Minimize negative log-likelihood starting from s=0.0
    for s0 in [0, 0.05]:
        result = minimize(lambda s: negative_log_likelihood(s, traj, timepoints), x0=[s0], method='L-BFGS-B')

    return result.x[0]  # The best-fit s

def add_selection_coeff(
    df, 
    timepoint_cols,     # columns holding allele frequencies for each time
    cluster_col='strain_cluster', 
    timepoints=None
):
    """
    For each cluster in 'df[cluster_col]', average the allele frequency columns 
    specified by 'timepoint_cols', then infer a single selection coefficient 
    s for that cluster. 
    Store the result in a new column 'selection_coeff' in the original DataFrame.
    
    timepoints : array-like or None
        If None, we try to parse the time from 'timepoint_cols' by converting them to float.
        If not None, it should be an array of the same length as timepoint_cols.
    """
    # Convert the timepoint columns to strings if needed
    timepoint_cols = [str(col) for col in timepoint_cols]
    df.columns = df.columns.astype(str)

    # Build the actual timepoints array if user didn't supply one
    if timepoints is None:
        # If your columns are named like "0.0", "1.0", "4.5", parse them:
        timepoints = np.array([float(c) for c in timepoint_cols], dtype=float)
    else:
        timepoints = np.array(timepoints, dtype=float)
    
    cluster2s = {}

    # Group by the cluster column
    for cluster_label, group in df.groupby(cluster_col):
        trajectories = group[timepoint_cols].astype(float).values
        knn_imputer = KNNImputer(n_neighbors=5, weights='uniform')
        trajectories = knn_imputer.fit_transform(trajectories)

        # compute mean trajectory across timepoint_cols
        mean_traj = trajectories.mean(axis=0)
        # infer s
        s_hat = infer_selection_for_cluster(mean_traj, timepoints)
        cluster2s[cluster_label] = s_hat
    
    # Map the result back to each row
    df['selection_coeff'] = df[cluster_col].map(cluster2s)
    return df
