"""Within-biopsy spatial graphs and exploratory ligand/receptor association."""
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import false_discovery_control
import squidpy as sq


def validate_graph(adata, sample_key='sample'):
    graph = sparse.csr_matrix(adata.obsp['spatial_connectivities'])
    graph.eliminate_zeros()
    if graph.shape != (adata.n_obs, adata.n_obs) or not np.isfinite(graph.data).all():
        raise ValueError('Invalid spatial graph')
    if (graph.data <= 0).any() or graph.diagonal().any() or (graph != graph.T).nnz:
        raise ValueError('Spatial graph must be positive, symmetric, and without self edges')
    if adata.obs[sample_key].isna().any():
        raise ValueError('Missing sample labels')
    source, target = graph.nonzero()
    samples = adata.obs[sample_key].to_numpy()
    if np.any(samples[source] != samples[target]):
        raise ValueError('Spatial graph contains cross-sample edges')
    return graph.astype(bool).astype(np.float64)


def build_spatial_graph(adata, *, n_neighs=6, coord_type='grid'):
    if adata.obs['sample'].isna().any() or not np.isfinite(adata.obsm['spatial']).all():
        raise ValueError('Finite coordinates and complete sample labels are required')
    if n_neighs < 1 or adata.obs.groupby('sample', observed=True).size().min() <= n_neighs:
        raise ValueError('Each sample must have more spots than n_neighs')
    adata.obs['sample'] = adata.obs['sample'].astype('category')
    sq.gr.spatial_neighbors(adata, library_key='sample', coord_type=coord_type, n_neighs=n_neighs)
    graph = adata.obsp['spatial_connectivities'].maximum(adata.obsp['spatial_connectivities'].T)
    graph.setdiag(0)
    graph.eliminate_zeros()
    adata.obsp['spatial_connectivities'] = graph.astype(bool).astype(float)
    distances = adata.obsp['spatial_distances']
    adata.obsp['spatial_distances'] = distances.maximum(distances.T)
    validate_graph(adata)
    return adata


def spatial_lr(adata, pairs, *, cluster_key='leiden', n_perms=999, seed=42):
    """Mean log-expression product on directed edges, separately for each biopsy.

    The null shuffles cluster labels within each biopsy while keeping expression
    and geometry fixed. This tests exploratory label association, not signalling
    activity or a disease effect. Zero-edge combinations have score=0 and p=1.
    BH covers every sample × present-cluster pair × available LR pair in the call.
    """
    if not isinstance(n_perms, int) or n_perms < 1:
        raise ValueError('n_perms must be a positive integer')
    if adata.n_obs == 0:
        raise ValueError('LR requires at least one spot')
    if 'lognorm' not in adata.layers:
        raise ValueError('Missing lognorm layer; raw counts are not accepted')
    if adata.obs[cluster_key].isna().any():
        raise ValueError('Missing cluster labels')
    graph = validate_graph(adata)
    pairs = list(dict.fromkeys(tuple(pair) for pair in pairs))
    pairs = [(lig, rec) for lig, rec in pairs if lig in adata.var_names and rec in adata.var_names]
    if not pairs:
        raise ValueError('No complete ligand/receptor pair is present')
    genes = sorted({gene for pair in pairs for gene in pair})
    expression = adata.layers['lognorm'][:, adata.var_names.get_indexer(genes)]
    expression = expression.toarray() if sparse.issparse(expression) else np.asarray(expression)
    expression = expression.astype(np.float64)
    if not np.isfinite(expression).all() or (expression < 0).any():
        raise ValueError('LR requires finite nonnegative log-normalized expression')
    gene_index = {gene: index for index, gene in enumerate(genes)}
    ligands = [gene_index[lig] for lig, _ in pairs]
    receptors = [gene_index[rec] for _, rec in pairs]
    rng = np.random.default_rng(seed)
    tables = []
    for sample in sorted(adata.obs['sample'].unique()):
        indices = np.flatnonzero(adata.obs['sample'].to_numpy() == sample)
        labels, levels = pd.factorize(adata.obs.iloc[indices][cluster_key].astype(str), sort=True)
        n_groups = len(levels)
        source, target = graph[indices][:, indices].nonzero()
        # Compute the LR edge products once, never copy AnnData in the null loop.
        values = expression[indices]
        products = values[source][:, ligands] * values[target][:, receptors]

        def scores_for(cluster_labels):
            codes = cluster_labels[source] * n_groups + cluster_labels[target]
            counts = np.bincount(codes, minlength=n_groups ** 2)
            sums = np.column_stack([np.bincount(codes, weights=products[:, p], minlength=n_groups ** 2)
                                    for p in range(len(pairs))])
            means = np.divide(sums, counts[:, None], out=np.zeros_like(sums, dtype=float), where=counts[:, None] > 0)
            return means, counts

        observed, counts = scores_for(labels)
        exceedances = np.ones_like(observed, dtype=np.int64)
        for _ in range(n_perms):
            null, _ = scores_for(rng.permutation(labels))
            exceedances += null >= observed
        table = pd.DataFrame({
            'sample': sample,
            'source': np.repeat(np.repeat(levels, n_groups), len(pairs)),
            'target': np.repeat(np.tile(levels, n_groups), len(pairs)),
            'ligand': np.tile([lig for lig, _ in pairs], n_groups ** 2),
            'receptor': np.tile([rec for _, rec in pairs], n_groups ** 2),
            'n_edges': np.repeat(counts, len(pairs)),
            'score': observed.ravel(),
            'pvalue': (exceedances / (n_perms + 1)).ravel(),
        })
        for key in ['condition', 'rejection_type']:
            metadata = adata.obs.iloc[indices][key].unique()
            if len(metadata) != 1 or pd.isna(metadata[0]):
                raise ValueError(f'{key} must be constant and nonmissing per sample')
            table[key] = metadata[0]
        tables.append(table)
    result = pd.concat(tables, ignore_index=True)
    result['padj'] = false_discovery_control(result['pvalue'].to_numpy(), method='bh')
    result['n_perms'] = n_perms
    result['fdr_family_size'] = len(result)
    return result
