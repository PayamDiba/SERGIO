import pandas as pd
import numpy as np
import scanpy as sc

from tqdm import tqdm
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import pairwise_distances
from joblib import Parallel, delayed
from .edistance import edist

# TODO make etest allow for multiple controls (accept list of controls)

def etest(adata, obs_key='perturbation', obsm_key='X_pca', dist='sqeuclidean',
          control='control', alpha=0.05, runs=1000, flavor=1, n_jobs=1,
          correction_method='holm-sidak', verbose=True):
    """Performs Monte Carlo permutation test with E-distance as test statistic.
    Tests for each group of cells defined in adata.obs[obs_key] if it is significantly
    different from control based on the E-distance in adata.obsm[obsm_key] space.
    Does multiple-testing correction using per default with Holm-Sidak.

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    obs_key: `str` in adata.obs.keys() (default: `perturbation`)
        Key in adata.obs specifying the groups to consider.
    obsm_key: `str` in adata.obsm (default: `adata.obsm['X_pca']`)
        Key for embedding coordinates to use.
    dist: `str` for any distance in scipy.spatial.distance (default: `sqeuclidean`)
        Distance metric to use in embedding space.
    control: `str` (default: `'control'`)
        Defines the control group in adata.obs[obs_key] to test against.
    alpha: `float` between `0` and `1` (default: `0.05`)
        significance cut-off for the test to annotate significance.
    runs: `int` (default: `100`)
        Number of iterations for the permutation test. This is basically the resolution of the E-test p-value.
        E.g. if you choose two iterations, then the p-value can only have 3 values `(0, .5, 1)`. Lower numbers will be much faster.
        We do not recommend going lower than `100` and suggest between `100` and `10000` iterations.
    correction_method: `None` or any valid method for statsmodels.stats.multitest.multipletests (default: `'holm-sidak'`)
        Method used for multiple-testing correction, since we are testing each group in `adata.obs[obs_key]`.
    verbose: `bool` (default: `True`)
        Whether to show a progress bar iterating over all groups.

    Returns
    -------
    tab: pandas.DataFrame
        E-test results for each group in adata.obs[obs_key] with columns
        - edist: E-distance to control
        - pvalue: E-test p-value if group is different from control
        - significant: If p-value < alpha
        - pvalue_adj: Multiple-testing corrected E-test p-value
        - significant_adj: If p-value_adj < alpha
    """

    groups = pd.unique(adata.obs[obs_key])
    
    # Compute pairwise distances selectively once
    # (we need pairwise distances within each group and between each group and control)
    # Note: this could be improved further, since we compute distances within control multiple times here. Speedup likely minimal though.
    pwds = {}
    for group in groups:
        x = adata[adata.obs[obs_key].isin([group, control])].obsm[obsm_key].copy()
        pwd = pairwise_distances(x,x, metric=dist)
        pwds[group] = pwd

    # Approximate sampling from null distribution (equal distributions)
    res = []
    fct = tqdm if verbose else lambda x: x
    M = np.sum(adata.obs[obs_key]==control)
    def one_step():
        # per perturbation, shuffle with control and compute e-distance
        df = pd.DataFrame(index=groups, columns=['edist'], dtype=float)
        for group in groups:
            if group==control:
                df.loc[group] = [0]
                continue
            N = np.sum(adata.obs[obs_key]==group)
            # shuffle the labels
            labels = adata.obs[obs_key].values[adata.obs[obs_key].isin([group, control])]
            shuffled_labels = np.random.permutation(labels)

            # use precomputed pairwise distances
            sc_pwd = pwds[group]  # precomputed pairwise distances between single cells
            idx = shuffled_labels==group

            # Note that this is wrong: sc_pwd[idx, ~idx] but this is correct: sc_pwd[idx, :][:, ~idx]
            # The first produces a vector, the second a matrix (we need the matrix)
            factor = N / (N-1) if flavor==1 else 1
            factor_c = M / (M-1) if flavor==1 else 1
            delta = np.sum(sc_pwd[idx, :][:, ~idx]) / (N * M)
            sigma = np.sum(sc_pwd[idx, :][:, idx]) / (N * N) * factor
            sigma_c = np.sum(sc_pwd[~idx, :][:, ~idx]) / (M * M) * factor_c

            edistance = 2 * delta - sigma - sigma_c

            df.loc[group] = edistance
        return df.sort_index()
    res = Parallel(n_jobs=n_jobs)(delayed(one_step)() for i in fct(range(runs)))
    
    # "Sampling" from original distribution without shuffling (hypothesis)
    df_old = edist(adata, obs_key, obsm_key=obsm_key, dist=dist, verbose=False).loc[control]
    df_old.columns = ['edist']
    
    # the following is faster than the above and produces the same result
    original = []
    for group in groups:
        if group==control:
            original.append(0)
            continue
        N = np.sum(adata.obs[obs_key]==group)
        # shuffle the labels
        labels = adata.obs[obs_key].values[adata.obs[obs_key].isin([group, control])]
        
        # use precomputed pairwise distances
        sc_pwd = pwds[group]  # precomputed pairwise distances between single cells
        idx = labels==group
        
        # Note that this is wrong: sc_pwd[idx, ~idx] but this is correct: sc_pwd[idx, :][:, ~idx]
        # The first produces a vector, the second a matrix (we need the matrix)
        factor = N / (N-1) if flavor==1 else 1
        factor_c = M / (M-1) if flavor==1 else 1
        delta = np.mean(sc_pwd[idx, :][:, ~idx]) / (N * M)
        sigma = np.mean(sc_pwd[idx, :][:, idx]) / (N * N) * factor
        sigma_c = np.mean(sc_pwd[~idx, :][:, ~idx]) / (M * M) * factor
        
        edistance = 2 * delta - sigma - sigma_c
        original.append(edistance)
    df = pd.DataFrame(original, index=groups, columns=['edist'])
    df = df.sort_index()

    # Evaluate test (hypothesis vs null hypothesis)
    # count times shuffling resulted in larger e-distance
    results = np.array(pd.concat([r['edist'] - df['edist'] for r in res], axis=1) > 0, dtype=int)
    n_failures = pd.Series(np.clip(np.sum(results, axis=1), 1, np.inf), index=df.index)
    pvalues = n_failures / runs

    # Apply multiple testing correction
    significant_adj, pvalue_adj, _, _ = multipletests(pvalues.values, alpha=alpha, method=correction_method)

    # Aggregate results
    tab = pd.DataFrame({'edist': df['edist'], 'pvalue': pvalues, 
                        'significant': pvalues < alpha, 'pvalue_adj': pvalue_adj, 
                        'significant_adj': significant_adj}, index=df.index)
    return tab