
import pandas as pd
import numpy as np
import scanpy as sc

from tqdm import tqdm
from sklearn.metrics import pairwise_distances
from warnings import warn

def pairwise_pca_distances(adata, obs_key, obsm_key='X_pca', dist='sqeuclidean',
                           correction_factor=False, verbose=True):
    """Average of pairwise PCA distances between cells of each group in obs_key.
    For each pair of groups defined in adata.obs[obs_key] (e.g. perturbations)
    computes all pairwise distances between cells in adata.obsm[obsm_key] (e.g. PCA space)
    and averages them per group-pair. This results in a distance matrix between all groups.
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    obs_key: `str` in adata.obs.keys()
        Key in adata.obs specifying the groups to consider.
    obsm_key: `str` in adata.obsm (default: `adata.obsm['X_pca']`)
        Key for embedding coordinates to use.
    dist: `str` for any distance in scipy.spatial.distance (default: `sqeuclidean`)
        Distance metric to use in embedding space.
    correction_factor: `bool` (default: `False`)
        Whether make the estimator for sigma more unbiased (dividing by N-1 instead of N, similar to sample and population variance).
    verbose: `bool` (default: `True`)
        Whether to show a progress bar iterating over all groups.
    Returns
    -------
    pwd: pandas.DataFrame
        DataFrame with average of pairwise PCA distances between all groups.
    """

    if obsm_key=='X_pca' and 'X_pca' not in adata.obsm.keys():
        warn('PCA embedding not found, computing...')
        sc.pp.pca(adata)
    
    X = adata.obsm[obsm_key].copy()
    y = adata.obs[obs_key].astype(str).copy()
    groups = pd.unique(y)
    df = pd.DataFrame(index=groups, columns=groups, dtype=float)
    fct = tqdm if verbose else lambda x: x
    for i, p1 in enumerate(fct(groups)):
        x1 = X[y==p1].copy()
        N = len(x1)
        for p2 in groups[i:]:
            x2 = X[y==p2].copy()
            pwd = pairwise_distances(x1, x2, metric=dist)
            M = len(x2) if (p1!=p2) & ~correction_factor else len(x2)-1
            factor = N * M
            mean_pwd = np.sum(pwd) / factor
            df.loc[p1, p2] = mean_pwd
            df.loc[p2, p1] = mean_pwd
    df.index.name = obs_key
    df.columns.name = obs_key
    df.name = 'pairwise PCA distances'
    return df

def edist(adata, obs_key='perturbation', obsm_key='X_pca', pwd=None, 
          dist='sqeuclidean', flavor=0, verbose=True):
    """Computes the edistance to control. Accepts precomputed pwd.
    Computes the pairwise E-distances between all groups of cells defined in
    adata.obs[obs_key] (e.g. perturbations). Distances are computed in embedding
    space given by adata.obsm[obsm_key] (e.g. PCA space).
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
    correction_factor: `bool` (default: `False`)
        Whether make the estimator for sigma more unbiased (dividing by N-1 instead of N, similar to sample and population variance).
    verbose: `bool` (default: `True`)
        Whether to show a progress bar iterating over all groups.
    Returns
    -------
    estats: pandas.DataFrame
        DataFrame with pairwise E-distances between all groups.
    """
    correction_factor = flavor==1
    pwd = pairwise_pca_distances(adata, obs_key=obs_key, obsm_key=obsm_key, 
                                 dist=dist, correction_factor=correction_factor, 
                                 verbose=verbose) if pwd is None else pwd
    # derive basic statistics
    sigmas = np.diag(pwd)
    deltas = pwd
    estats = 2 * deltas - sigmas - sigmas[:, np.newaxis]
    if flavor == 2:
        sizes = np.array([np.sum(adata.obs[obs_key]==g) for g in pwd.index])
        corrections = np.outer(sizes, sizes) / np.add.outer(sizes, sizes)
        estats = estats * corrections
    if flavor == 3:
        estats = estats / deltas
    return estats

def onesided_pca_distances(adata, obs_key, control, obsm_key='X_pca', 
                           dist='sqeuclidean', correction_factor=False, 
                           verbose=True):
    """Average of pairwise PCA distances between cells of each group in obs_key with control group.
    For each group defined in adata.obs[obs_key] (e.g. perturbations)
    computes all pairwise distances between cells in adata.obsm[obsm_key] (e.g. PCA space)
    and averages them per group-control-pair. This results in a distance vector with a value for each group.
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    obs_key: `str` in adata.obs.keys()
        Key in adata.obs specifying the groups to consider.
    control: `str` of a category in adata.obs[obs_key]
        Group in obs_key for control cells.
    obsm_key: `str` in adata.obsm (default: `adata.obsm['X_pca']`)
        Key for embedding coordinates to use.
    dist: `str` for any distance in scipy.spatial.distance (default: `sqeuclidean`)
        Distance metric to use in embedding space.
    correction_factor: `bool` (default: `False`)
        Whether make the estimator for sigma more unbiased (dividing by N-1 instead of N, similar to sample and population variance).
    verbose: `bool` (default: `True`)
        Whether to show a progress bar iterating over all groups.
    Returns
    -------
    pwd: pandas.DataFrame
        DataFrame with average PCA distances to control for all groups.
    """

    if obsm_key=='X_pca' and 'X_pca' not in adata.obsm.keys():
        warn('PCA embedding not found, computing...')
        sc.pp.pca(adata)

    groups = pd.unique(adata.obs[obs_key])
    assert control in groups, f'No cells of control group "{control}" were not found in groups defined by "{obs_key}".'
    df = pd.DataFrame(index=groups, columns=['distance'], dtype=float)
    fct = tqdm if verbose else lambda x: x
    
    x1 = adata[adata.obs[obs_key]==control].obsm[obsm_key].copy()
    N = len(x1)
    for p in fct(groups):
        x2 = adata[adata.obs[obs_key]==p].obsm[obsm_key].copy()
        pwd = pairwise_distances(x1, x2, metric=dist)
        M = len(x2) if (p==control) & ~correction_factor else len(x2)-1
        factor = N * M  # Thanks to Garrett Wong for finding this bug
        mean_pwd = np.sum(pwd) / factor
        df.loc[p] = mean_pwd
    df.index.name = obs_key
    df.name = f'PCA distances to {control}'
    return df

def self_pca_distances(adata, obs_key, obsm_key='X_pca', dist='sqeuclidean', 
                       correction_factor=False, verbose=True):
    """Average of pairwise PCA distances between cells within each group in obs_key.
    For each group defined in adata.obs[obs_key] (e.g. perturbations)
    computes all pairwise distances between cells within in the space given by adata.obsm[obsm_key] (e.g. PCA space)
    and averages them per group. This results in a distance vector with a value for each group.
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    obs_key: `str` in adata.obs.keys()
        Key in adata.obs specifying the groups to consider.
    obsm_key: `str` in adata.obsm (default: `adata.obsm['X_pca']`)
        Key for embedding coordinates to use.
    dist: `str` for any distance in scipy.spatial.distance (default: `sqeuclidean`)
        Distance metric to use in embedding space.
    verbose: `bool` (default: `True`)
        Whether to show a progress bar iterating over all groups.
    Returns
    -------
    pwd: pandas.DataFrame
        DataFrame with average PCA distances to control for all groups.
    """

    if obsm_key=='X_pca' and 'X_pca' not in adata.obsm.keys():
        warn('PCA embedding not found, computing...')
        sc.pp.pca(adata)

    groups = pd.unique(adata.obs[obs_key])
    df = pd.DataFrame(index=groups, columns=['distance'], dtype=float)
    fct = tqdm if verbose else lambda x: x
    
    for p in fct(groups):
        x = adata[adata.obs[obs_key]==p].obsm[obsm_key].copy()
        pwd = pairwise_distances(x, x, metric=dist)
        N = len(x)
        factor = N**2-N if correction_factor else N**2
        mean_pwd = np.sum(pwd) / factor
        df.loc[p] = mean_pwd
    df.index.name = obs_key
    df.name = 'PCA distances within groups'
    return df

def edist_to_control(adata, obs_key='perturbation', control='control', 
                     obsm_key='X_pca', dist='sqeuclidean',  
                     flavor=0, verbose=True):
    """Computes the edistance to control.
    Computes the all E-distances between all groups of cells defined in
    adata.obs[obs_key] (e.g. perturbations) and control cells. Distances are computed in embedding
    space given by adata.obsm[obsm_key] (e.g. PCA space).
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
    flavor: `int` (default: `0`)
        Which flavor of E-distance to use.
        - `0`: E-distance no correction factor
        - `1`: E-distance with simple correction factor (N-1)
        - `2`: E-distance with correction factor NM/(N+M)
        - `3`: E-distance with correction division by delta
    verbose: `bool` (default: `True`)
        Whether to show a progress bar iterating over all groups.
    Returns
    -------
    ed_to_c: pandas.DataFrame
        DataFrame with E-distances between all groups and control group.
    """
    correction_factor = flavor==1
    deltas_to_c = onesided_pca_distances(adata, obs_key=obs_key, control=control, 
                                         obsm_key=obsm_key, dist=dist, 
                                         correction_factor=correction_factor, 
                                         verbose=verbose)
    sigmas = self_pca_distances(adata, obs_key, obsm_key=obsm_key, dist=dist, 
                                correction_factor=correction_factor, verbose=False)
    # derive basic statistics
    ed_to_c = 2 * deltas_to_c - sigmas - sigmas.loc[control]
    if flavor == 2:
        sizes = np.array([np.sum(adata.obs[obs_key]==g) for g in sigmas.index])
        corrections = (sizes * sizes[sigmas.index==control]) / (sizes + sizes[sigmas.index==control])
        ed_to_c = ed_to_c * corrections.reshape(-1,1)
    if flavor == 3:
        ed_to_c = ed_to_c / deltas_to_c
    return ed_to_c
Footer
© 2023 GitHub, Inc.
Footer navigation
Terms
Privacy
Security
Status
Docs
Contact GitHub
Pricing
API
Training
Blog
About
scPerturb/edistance.py at master · sanderlab/scPerturb