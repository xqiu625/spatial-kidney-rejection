"""Configuration, run provenance and descriptive marker summaries."""
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import tempfile
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import yaml


def load_config(config_path):
    with Path(config_path).open() as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError(f'Expected a configuration mapping: {config_path}')
    return config


def load_settings(config_dir=None):
    config_dir = Path(config_dir or os.environ.get('SKR_CONFIG_DIR') or Path(__file__).resolve().parents[1] / 'config').resolve()
    paths = load_config(config_dir / 'paths.yaml')
    params = load_config(config_dir / 'analysis_params.yaml')
    root = config_dir.parent
    base = Path(os.environ.get('SKR_DATA_DIR', paths['base_dir'])).expanduser()
    if not base.is_absolute():
        base = (root / base).resolve()
    raw = Path(paths['data']['raw'].replace('${base_dir}', str(base))).expanduser()
    if not raw.is_absolute():
        raw = (root / raw).resolve()
    run = Path(os.environ.get('SKR_RUN_DIR', paths['run_dir'])).expanduser()
    if not run.is_absolute():
        run = (root / run).resolve()
    if '${' in str(raw):
        raise ValueError('Unresolved variable in raw data path')
    return {'raw_dir': str(raw), 'run_dir': str(run), 'samples': paths['samples'], 'params': params}


def begin_stage(settings, output_name, input_name=None):
    """Require a fresh run for stage 01, and matching configuration thereafter."""
    run = Path(settings['run_dir'])
    manifest = run / 'run.json'
    output = run / 'processed' / output_name
    if input_name is None:
        run.mkdir(parents=True, exist_ok=True)
        if any(run.iterdir()):
            raise FileExistsError(f'Run directory is not empty. Set a new SKR_RUN_DIR: {run}')
        versions = {}
        for package in ['scanpy', 'squidpy', 'anndata', 'numpy', 'pandas', 'scipy', 'scikit-learn', 'leidenalg']:
            versions[package] = importlib.metadata.version(package)
        revision = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=Path(__file__).resolve().parents[1],
                                  text=True, capture_output=True, check=False).stdout.strip()
        dirty = bool(subprocess.run(['git', 'status', '--porcelain'], cwd=Path(__file__).resolve().parents[1],
                                    text=True, capture_output=True, check=False).stdout.strip())
        source = Path(settings['raw_dir']) / 'source.json'
        manifest.write_text(json.dumps({'settings': settings, 'versions': versions,
                                       'input_archive': json.loads(source.read_text()) if source.exists() else None,
                                       'git_commit': revision, 'working_tree_dirty': dirty}, indent=2) + '\n')
        for directory in ['processed', 'tables', 'figures']:
            (run / directory).mkdir()
    else:
        previous = json.loads(manifest.read_text())
        if previous['settings'] != settings:
            raise ValueError('Configuration changed during this run. Start a fresh SKR_RUN_DIR.')
        if not (run / 'processed' / input_name).is_file():
            raise FileNotFoundError(f'Run the previous notebook first: {input_name}')
    if output.exists():
        raise FileExistsError(f'Stage already completed. Use a new SKR_RUN_DIR: {output}')
    return run


def save_checkpoint(adata, path):
    """Do not replace completed checkpoints; publish only a complete H5AD."""
    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    with tempfile.NamedTemporaryFile(dir=path.parent, suffix='.partial.h5ad', delete=False) as handle:
        temporary = Path(handle.name)
    try:
        adata.write_h5ad(temporary)
        # Hard linking fails if another writer already published this checkpoint.
        os.link(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def get_kidney_markers():
    """
    Get marker genes for kidney cell types.

    Returns
    -------
    dict
        Dictionary of cell types and their marker genes
    """
    markers = {
        'Proximal_Tubule': ['SLC34A1', 'LRP2', 'CUBN', 'SLC13A3'],
        'Loop_of_Henle': ['SLC12A1', 'UMOD', 'CLDN16'],
        'Distal_Tubule': ['SLC12A3', 'CALB1', 'TRPM6'],
        'Collecting_Duct_PC': ['AQP2', 'AQP3', 'FXYD4'],
        'Collecting_Duct_IC': ['SLC4A1', 'ATP6V1B1', 'SLC26A4'],
        'Podocyte': ['NPHS1', 'NPHS2', 'WT1', 'PODXL'],
        'Endothelial': ['PECAM1', 'VWF', 'CDH5', 'FLT1'],
        'Mesangial': ['PDGFRB', 'ITGA8', 'DES'],
        'Fibroblast': ['DCN', 'COL1A1', 'COL3A1', 'LUM'],
        'Macrophage': ['CD68', 'CD163', 'CSF1R', 'MARCO'],
        'T_Cell': ['CD3D', 'CD3E', 'CD4', 'CD8A'],
        'B_Cell': ['CD79A', 'MS4A1', 'CD19'],
        'NK_Cell': ['NKG7', 'GNLY', 'KLRD1'],
        'Plasma_Cell': ['MZB1', 'SDC1', 'JCHAIN'],
    }

    return markers


def score_cell_types(adata, markers=None, *, method='marker', seed=42):
    """Background-corrected log-expression scores; never cell fractions."""
    if method != 'marker':
        raise NotImplementedError('Only marker scoring is supported. Reference deconvolution is not validated.')
    if 'lognorm' not in adata.layers:
        raise ValueError('Run normalization first: missing lognorm layer')
    markers = get_kidney_markers() if markers is None else markers
    scored = []
    for cell_type, genes in markers.items():
        present = [gene for gene in genes if gene in adata.var_names]
        if len(present) < 2:
            warnings.warn(f'Skipping {cell_type}: fewer than two markers are present', stacklevel=2)
            continue
        name = f'{cell_type}_score'
        sc.tl.score_genes(adata, present, score_name=name, use_raw=False,
                         layer='lognorm', random_state=seed)
        scored.append(name)
    if not scored:
        raise ValueError('No cell types have enough markers to score')
    if not np.isfinite(adata.obs[scored].to_numpy()).all():
        raise ValueError('Marker scoring produced non-finite values')
    adata.uns['marker_score_columns'] = scored
    adata.uns['cell_type_method'] = 'background_corrected_marker_score_not_abundance'
    return adata


def sample_score_summary(adata, columns):
    """One row per biopsy, without spot-level disease hypothesis tests."""
    keys = ['sample', 'condition', 'rejection_type']
    for key in keys:
        if key not in adata.obs or adata.obs[key].isna().any():
            raise ValueError(f'Missing sample metadata: {key}')
    if (adata.obs.groupby('sample', observed=True)[keys[1:]].nunique() != 1).any().any():
        raise ValueError('Condition/subtype must be constant within each sample')
    groups = adata.obs.groupby(keys, observed=True)
    summary = groups[list(columns)].mean()
    summary.insert(0, 'n_spots', groups.size())
    return summary.reset_index()


def profile_summary(adata, key='marker_profile'):
    """Per-sample profile frequencies, including absent profiles as zero."""
    metadata = sample_score_summary(adata, []).set_index('sample')
    counts = pd.crosstab(adata.obs['sample'], adata.obs[key]).reindex(metadata.index, fill_value=0)
    fractions = counts.div(counts.sum(axis=1), axis=0).add_prefix('profile_')
    return metadata.join(fractions).reset_index()
