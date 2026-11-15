import copy
import gzip
import io
import json
from pathlib import Path
import tarfile

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from scipy import sparse
from scipy.stats import false_discovery_control

from src import preprocessing, utils, visualization
from src.prepare_data import prepare_archive
from src.spatial import build_spatial_graph, spatial_lr, validate_graph


def toy():
    obs = pd.DataFrame({'sample': ['s1'] * 4, 'condition': ['control'] * 4,
                        'rejection_type': ['none'] * 4, 'leiden': ['A', 'A', 'B', 'B']},
                       index=list('abcd'))
    a = ad.AnnData(sparse.csr_matrix([[10., 1., 8.], [9., 2., 8.], [1., 10., 8.], [2., 9., 8.]]),
                   obs=obs, var=pd.DataFrame(index=['L', 'R', 'other']))
    a.obsm['spatial'] = np.column_stack([np.arange(4), np.zeros(4)])
    a.obsp['spatial_connectivities'] = sparse.csr_matrix(np.ones((4, 4)) - np.eye(4))
    return preprocessing.normalize_data(a)


def test_counts_preserved_and_normalization_repeatable():
    a = toy()
    counts, lognorm = a.layers['counts'].copy(), a.layers['lognorm'].copy()
    sc.pp.scale(a)
    preprocessing.normalize_data(a)
    np.testing.assert_allclose(a.layers['counts'].toarray(), counts.toarray())
    np.testing.assert_allclose(a.layers['lognorm'].toarray(), lognorm.toarray())


def test_marker_fold_changes_use_lognorm_even_after_scaling():
    a = toy()
    sc.pp.scale(a)
    sc.tl.rank_genes_groups(a, 'leiden', method='wilcoxon', layer='lognorm', use_raw=False)
    table = sc.get.rank_genes_groups_df(a, group=None)
    assert np.isfinite(table.logfoldchanges).all()
    assert table.query("group == 'A' and names == 'L'").logfoldchanges.iloc[0] > 0


def test_spatial_lr_known_edge_product_and_geometry_dependence():
    a = toy()
    table = spatial_lr(a, [('L', 'R')], n_perms=19)
    row = table.query("source == 'A' and target == 'B'").iloc[0]
    x = a.layers['lognorm'].toarray()
    expected = np.mean([x[i, 0] * x[j, 1] for i in [0, 1] for j in [2, 3]])
    assert row.n_edges == 4
    assert row.score == pytest.approx(expected)
    a.obsp['spatial_connectivities'] = sparse.csr_matrix((4, 4))
    empty = spatial_lr(a, [('L', 'R')], n_perms=19)
    assert (empty.score == 0).all() and (empty.pvalue == 1).all()
    assert (empty.padj == 1).all() and len(empty) == len(table)


def test_depth_only_changes_do_not_change_lr():
    a, b = toy(), toy()
    b.layers['counts'] *= 3
    preprocessing.normalize_data(b)
    one = spatial_lr(a, [('L', 'R')], n_perms=19, seed=8)
    two = spatial_lr(b, [('L', 'R')], n_perms=19, seed=8)
    np.testing.assert_allclose(one.score, two.score, rtol=1e-7)


def test_fdr_full_family_determinism_and_bounds():
    a = toy()
    first = spatial_lr(a, [('L', 'R'), ('R', 'L')], n_perms=39, seed=10)
    second = spatial_lr(a, [('L', 'R'), ('R', 'L')], n_perms=39, seed=10)
    pd.testing.assert_frame_equal(first, second)
    assert len(first) == 8 and (first.fdr_family_size == 8).all()
    np.testing.assert_allclose(first.padj, false_discovery_control(first.pvalue))
    assert first.padj.between(0, 1).all()
    assert first.pvalue.min() >= 1 / 40


def test_graph_and_null_respect_sample_boundaries():
    a = toy()
    a.obs['sample'] = ['s1', 's1', 's2', 's2']
    with pytest.raises(ValueError, match='cross-sample'):
        spatial_lr(a, [('L', 'R')], n_perms=3)
    a.obsp['spatial_connectivities'] = sparse.block_diag([[[0, 1], [1, 0]]] * 2).tocsr()
    # Each biopsy has only one label. Within-biopsy permutation cannot change it.
    result = spatial_lr(a, [('L', 'R')], n_perms=19)
    assert len(result) == 2 and (result.pvalue == 1).all()
    assert (result.source == result.target).all()


def test_graph_builder_never_joins_overlapping_slides():
    coords = np.column_stack([np.arange(8) * 10., np.zeros(8)])
    a = ad.AnnData(np.ones((16, 2)), obs=pd.DataFrame({'sample': ['s1'] * 8 + ['s2'] * 8}, index=[str(i) for i in range(16)]))
    a.obsm['spatial'] = np.vstack([coords, coords + [0.1, 0]])
    build_spatial_graph(a, n_neighs=6, coord_type='generic')
    graph = validate_graph(a)
    assert graph[:8, 8:].nnz == 0 and graph[8:, :8].nnz == 0
    assert graph.nnz > 0 and (graph != graph.T).nnz == 0


def test_reference_route_is_explicitly_disabled():
    with pytest.raises(NotImplementedError, match='not validated'):
        utils.score_cell_types(toy(), method='cell2location')


def test_marker_scores_read_lognorm_without_forcing_fractions():
    rng = np.random.default_rng(33)
    counts = rng.poisson(np.linspace(2, 30, 300), size=(24, 300)).astype(float)
    a = ad.AnnData(counts, var=pd.DataFrame(index=[f'g{i}' for i in range(300)]))
    preprocessing.normalize_data(a)
    markers = {'T': ['g60', 'g65'], 'B': ['g150', 'g155']}
    utils.score_cell_types(a, markers)
    before = a.obs[['T_score', 'B_score']].copy()
    a.X = np.zeros_like(a.X)
    utils.score_cell_types(a, markers)
    pd.testing.assert_frame_equal(before, a.obs[['T_score', 'B_score']])
    assert not any(column.endswith('_prop') for column in a.obs)
    assert not np.allclose(before.sum(axis=1), 1)


def test_sample_summaries_do_not_create_disease_tests():
    a = toy()
    a.obs['T_score'] = [-2., -1., 1., 2.]
    a.obs['marker_profile'] = pd.Categorical(['0', '0', '1', '1'])
    summary = utils.sample_score_summary(a, ['T_score'])
    assert len(summary) == 1 and summary.n_spots.iloc[0] == 4
    assert not {'pvalue', 'padj', 'log2FC'} & set(summary.columns)
    frequencies = utils.profile_summary(a)
    assert frequencies.profile_0.iloc[0] == 0.5
    a.obs.loc['a', 'condition'] = 'rejection'
    with pytest.raises(ValueError, match='constant'):
        utils.sample_score_summary(a, ['T_score'])


def test_config_expansion_checkpoint_and_stale_run_protection(tmp_path, monkeypatch):
    monkeypatch.setenv('SKR_DATA_DIR', str(tmp_path / 'inputs'))
    monkeypatch.setenv('SKR_RUN_DIR', str(tmp_path / 'run'))
    settings = utils.load_settings()
    assert settings['raw_dir'] == str(tmp_path / 'inputs/data/raw')
    run = utils.begin_stage(settings, 'one.h5ad')
    utils.save_checkpoint(toy(), run / 'processed/one.h5ad')
    with pytest.raises(FileExistsError):
        utils.save_checkpoint(toy(), run / 'processed/one.h5ad')
    utils.begin_stage(settings, 'two.h5ad', 'one.h5ad')
    changed = copy.deepcopy(settings)
    changed['params']['seed'] += 1
    with pytest.raises(ValueError, match='changed'):
        utils.begin_stage(changed, 'two.h5ad', 'one.h5ad')
    with pytest.raises(FileExistsError, match='not empty'):
        utils.begin_stage(settings, 'one.h5ad')


def test_lr_plot_honors_top_n():
    import matplotlib.pyplot as plt
    result = spatial_lr(toy(), [('L', 'R')], n_perms=3)
    fig = visualization.plot_ligand_receptor(result, top_n=2)
    assert len(fig.axes[0].patches) == 2
    plt.close(fig)


def test_spatial_profile_colors_do_not_shift_when_groups_are_absent():
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgba
    a = toy()
    a.obs['sample'] = ['s1', 's1', 's2', 's2']
    a.obs['marker_profile'] = pd.Categorical(['0', '0', '1', '1'])
    a.uns['spatial'] = {sample: {'images': {'hires': np.ones((10, 10, 3))},
        'scalefactors': {'tissue_hires_scalef': 1, 'spot_diameter_fullres': 1}}
        for sample in ['s1', 's2']}
    fig = visualization.plot_spatial_clusters(a, cluster_key='marker_profile', ncols=2)
    colors = a.uns['marker_profile_colors']
    assert len(colors) == 2 and colors[0] != colors[1]
    for index in range(2):
        actual = fig.axes[index].collections[0].get_facecolors()[0]
        np.testing.assert_allclose(actual, to_rgba(colors[index]))
    plt.close(fig)


def test_geo_preparation_rejects_nonfiles_and_never_extracts_paths(tmp_path):
    archive = tmp_path / 'malicious.tar'
    with tarfile.open(archive, 'w') as tar:
        entry = tarfile.TarInfo('s1_matrix.mtx.gz')
        entry.type = tarfile.SYMTYPE
        entry.linkname = '/etc/passwd'
        tar.addfile(entry)
    with pytest.raises(ValueError, match='regular file'):
        prepare_archive(archive, tmp_path / 'raw', ['s1'])
    assert not (tmp_path / 'raw').exists()


def test_positions_and_barcode_alignment(tmp_path):
    from scipy.io import mmwrite
    for sample in ['s1', 's2']:
        root = tmp_path / sample
        matrix = root / 'filtered_feature_bc_matrix'
        matrix.mkdir(parents=True)
        (root / 'spatial').mkdir()
        mmwrite(matrix / 'matrix.mtx', sparse.coo_matrix([[1, 2], [3, 4]]))
        (matrix / 'features.tsv').write_text('g1\tG1\ng2\tG2\n')
        (matrix / 'barcodes.tsv').write_text('a\nb\n')
        # Reverse order in positions to check barcode, not positional, matching.
        (root / 'spatial/tissue_positions.csv').write_text('barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\nb,1,0,1,20,200\na,1,0,0,10,100\n')
        (root / 'spatial/scalefactors_json.json').write_text('{}')
    data = preprocessing.load_all_samples(tmp_path, [{'name': s, 'condition': 'control', 'rejection_type': 'none'} for s in ['s1', 's2']])
    assert data.obs_names.is_unique
    np.testing.assert_array_equal(data.obsm['spatial'], [[100, 10], [200, 20]] * 2)
    assert data.obs['rejection_type'].eq('none').all()
