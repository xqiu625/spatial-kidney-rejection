"""
Preprocessing functions for spatial transcriptomics data
"""

import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Optional, Dict, Union


def load_visium_sample(
    sample_path: Union[str, Path],
    sample_name: Optional[str] = None
) -> sc.AnnData:
    """
    Load a single Visium sample.

    Parameters
    ----------
    sample_path : str or Path
        Path to the sample directory containing filtered_feature_bc_matrix and spatial folders
    sample_name : str, optional
        Name to assign to the sample

    Returns
    -------
    AnnData
        Loaded spatial transcriptomics data
    """
    import json
    from PIL import Image
    from scipy.io import mmread
    from scipy.sparse import csr_matrix

    sample_path = Path(sample_path)
    library_id = sample_name or sample_path.name

    # Load count matrix manually (handles both compressed and uncompressed)
    matrix_dir = sample_path / 'filtered_feature_bc_matrix'

    # Find matrix file (try both .mtx and .mtx.gz)
    mtx_path = matrix_dir / 'matrix.mtx'
    if not mtx_path.exists():
        mtx_path = matrix_dir / 'matrix.mtx.gz'

    # Read matrix
    matrix = mmread(str(mtx_path)).T.tocsr()  # transpose: genes x cells -> cells x genes

    # Find and read barcodes
    barcodes_path = matrix_dir / 'barcodes.tsv'
    if not barcodes_path.exists():
        barcodes_path = matrix_dir / 'barcodes.tsv.gz'
    barcodes = pd.read_csv(barcodes_path, header=None, sep='\t')[0].values

    # Find and read features/genes
    features_path = matrix_dir / 'features.tsv'
    if not features_path.exists():
        features_path = matrix_dir / 'features.tsv.gz'
    if not features_path.exists():
        features_path = matrix_dir / 'genes.tsv'
    if not features_path.exists():
        features_path = matrix_dir / 'genes.tsv.gz'

    features = pd.read_csv(features_path, header=None, sep='\t')
    gene_ids = features[0].values
    gene_names = features[1].values if features.shape[1] > 1 else gene_ids

    # Create AnnData
    adata = sc.AnnData(
        X=csr_matrix(matrix),
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame({'gene_ids': gene_ids, 'gene_symbols': gene_names}, index=gene_names)
    )

    adata.var_names_make_unique()
    adata.obs['sample'] = library_id

    # Load spatial data
    spatial_dir = sample_path / 'spatial'

    # Load scale factors
    scalefactors_path = spatial_dir / 'scalefactors_json.json'
    with open(scalefactors_path, 'r') as f:
        scalefactors = json.load(f)

    # Load tissue positions
    positions_path = spatial_dir / 'tissue_positions.csv'
    if positions_path.exists():
        positions = pd.read_csv(positions_path)
    else:
        positions = pd.read_csv(spatial_dir / 'tissue_positions_list.csv', header=None,
            names=['barcode', 'in_tissue', 'array_row', 'array_col',
                   'pxl_row_in_fullres', 'pxl_col_in_fullres'])
    positions = positions.set_index('barcode')
    if not positions.index.is_unique or not adata.obs_names.is_unique:
        raise ValueError(f'Duplicate barcodes in {library_id}')
    missing = adata.obs_names.difference(positions.index)
    if len(missing):
        raise ValueError(f'{library_id}: {len(missing)} matrix barcodes lack spatial coordinates')
    positions = positions.loc[adata.obs_names]
    adata.obs['barcode'] = adata.obs_names
    adata.obs['in_tissue'] = positions['in_tissue'].to_numpy()
    adata.obsm['spatial'] = positions[['pxl_col_in_fullres', 'pxl_row_in_fullres']].to_numpy(float)
    if not np.isfinite(adata.obsm['spatial']).all():
        raise ValueError(f'Non-finite coordinates in {library_id}')
    adata = adata[adata.obs['in_tissue'] == 1].copy()

    # Load images
    hires_path = spatial_dir / 'tissue_hires_image.png'
    lowres_path = spatial_dir / 'tissue_lowres_image.png'

    adata.uns['spatial'] = {library_id: {'images': {}, 'scalefactors': scalefactors}}

    if hires_path.exists():
        adata.uns['spatial'][library_id]['images']['hires'] = np.array(Image.open(hires_path))
    if lowres_path.exists():
        adata.uns['spatial'][library_id]['images']['lowres'] = np.array(Image.open(lowres_path))

    return adata


def load_all_samples(
    raw_dir: Union[str, Path],
    sample_info: Optional[List[Dict]] = None
) -> sc.AnnData:
    """
    Load and concatenate all Visium samples.

    Parameters
    ----------
    raw_dir : str or Path
        Directory containing all sample folders
    sample_info : list of dict, optional
        Sample metadata with 'name' and 'condition' keys

    Returns
    -------
    AnnData
        Concatenated spatial data
    """
    raw_dir = Path(raw_dir)

    if sample_info is None:
        # Auto-detect sample directories
        sample_dirs = sorted([d for d in raw_dir.iterdir() if d.is_dir() and d.name.startswith('GSM')])
        sample_info = [{'name': d.name, 'condition': 'unknown'} for d in sample_dirs]

    if not sample_info or len({x['name'] for x in sample_info}) != len(sample_info):
        raise ValueError('Provide at least one sample, with unique sample names')
    adatas = []
    spatial_dict = {}  # Store spatial info for all samples

    for info in sample_info:
        sample_path = raw_dir / info['name']
        if sample_path.exists():
            print(f"Loading {info['name']}...")
            adata = load_visium_sample(sample_path, info['name'])
            adata.obs['condition'] = info.get('condition', 'unknown')
            adata.obs['rejection_type'] = info.get('rejection_type', 'unknown')

            # Store spatial info before concatenation
            if 'spatial' in adata.uns:
                spatial_dict.update(adata.uns['spatial'])

            adatas.append(adata)
        else:
            raise FileNotFoundError(f'Required sample missing: {sample_path}')

    # Concatenate
    adata = sc.concat(adatas, label='sample', keys=[x['name'] for x in sample_info],
                      index_unique='::', merge='same')

    # Restore spatial info after concatenation
    adata.uns['spatial'] = spatial_dict

    print(f"Loaded {len(adatas)} samples with {adata.n_obs} total spots")
    return adata


def calculate_qc_metrics(adata: sc.AnnData) -> sc.AnnData:
    """
    Calculate QC metrics for spatial data.

    Parameters
    ----------
    adata : AnnData
        Input data

    Returns
    -------
    AnnData
        Data with QC metrics in obs
    """
    # Mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')

    # Ribosomal genes
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))

    # Hemoglobin genes
    adata.var['hb'] = adata.var_names.str.match(r'^HB(?!P)')

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt', 'ribo', 'hb'],
        percent_top=None,
        log1p=False,
        inplace=True
    )

    return adata


def filter_spots(
    adata: sc.AnnData,
    min_genes: int = 200,
    min_cells: int = 3,
    max_pct_mt: float = 20.0
) -> sc.AnnData:
    """
    Filter spots and genes based on QC metrics.

    Parameters
    ----------
    adata : AnnData
        Input data with QC metrics
    min_genes : int
        Minimum genes per spot
    min_cells : int
        Minimum spots per gene
    max_pct_mt : float
        Maximum mitochondrial percentage

    Returns
    -------
    AnnData
        Filtered data
    """
    n_spots_before = adata.n_obs
    n_genes_before = adata.n_vars

    # Filter spots
    sc.pp.filter_cells(adata, min_genes=min_genes)

    # Filter genes
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # Filter by mitochondrial content
    adata = adata[adata.obs['pct_counts_mt'] < max_pct_mt, :].copy()

    print(f"Filtered: {n_spots_before} -> {adata.n_obs} spots")
    print(f"Filtered: {n_genes_before} -> {adata.n_vars} genes")

    return adata


def normalize_data(adata: sc.AnnData, target_sum: float = 1e4) -> sc.AnnData:
    """
    Normalize and log-transform data.

    Parameters
    ----------
    adata : AnnData
        Input data
    target_sum : float
        Target sum for normalization

    Returns
    -------
    AnnData
        Normalized data
    """
    # Reuse counts when called again; never replace counts with normalized data.
    if 'counts' not in adata.layers:
        values = adata.X.data if hasattr(adata.X, 'toarray') else np.asarray(adata.X)
        if not np.isfinite(values).all() or (values < 0).any() or not np.allclose(values, np.rint(values)):
            raise ValueError('Normalization requires finite, nonnegative UMI counts')
        adata.layers['counts'] = adata.X.copy()
    adata.X = adata.layers['counts'].copy()
    if np.any(np.asarray(adata.X.sum(axis=1)).ravel() <= 0):
        raise ValueError('Remove zero-count spots before normalization')

    # Normalize
    sc.pp.normalize_total(adata, target_sum=target_sum)

    # Log transform
    sc.pp.log1p(adata)
    adata.layers['lognorm'] = adata.X.copy()

    return adata


def select_hvg(
    adata: sc.AnnData,
    n_top_genes: int = 2000,
    batch_key: Optional[str] = 'sample'
) -> sc.AnnData:
    """
    Select highly variable genes.

    Parameters
    ----------
    adata : AnnData
        Normalized data
    n_top_genes : int
        Number of HVGs to select
    batch_key : str, optional
        Key for batch correction in HVG selection

    Returns
    -------
    AnnData
        Data with HVG annotation
    """
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        batch_key=batch_key,
        subset=False
    )

    n_hvg = adata.var['highly_variable'].sum()
    print(f"Selected {n_hvg} highly variable genes")

    return adata


def run_preprocessing_pipeline(
    raw_dir: Union[str, Path],
    output_path: Union[str, Path],
    sample_info: Optional[List[Dict]] = None,
    min_genes: int = 200,
    min_cells: int = 3,
    max_pct_mt: float = 20.0,
    n_top_genes: int = 2000
) -> sc.AnnData:
    """
    Run the full preprocessing pipeline.

    Parameters
    ----------
    raw_dir : str or Path
        Directory containing raw data
    output_path : str or Path
        Path to save processed data
    sample_info : list of dict, optional
        Sample metadata
    min_genes, min_cells, max_pct_mt, n_top_genes :
        QC and HVG parameters

    Returns
    -------
    AnnData
        Processed data
    """
    print("=" * 50)
    print("SPATIAL TRANSCRIPTOMICS PREPROCESSING PIPELINE")
    print("=" * 50)

    if Path(output_path).exists():
        raise FileExistsError(f'Use a new run directory: {output_path}')

    # Load data
    print("\n1. Loading samples...")
    adata = load_all_samples(raw_dir, sample_info)

    # QC metrics
    print("\n2. Calculating QC metrics...")
    adata = calculate_qc_metrics(adata)

    # Filter
    print("\n3. Filtering spots and genes...")
    adata = filter_spots(adata, min_genes, min_cells, max_pct_mt)

    # Normalize
    print("\n4. Normalizing data...")
    adata = normalize_data(adata)

    # HVG selection
    print("\n5. Selecting highly variable genes...")
    adata = select_hvg(adata, n_top_genes)

    # Save
    print(f"\n6. Saving to {output_path}...")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(output_path)

    print("\n" + "=" * 50)
    print("PREPROCESSING COMPLETE")
    print(f"Final data: {adata.n_obs} spots x {adata.n_vars} genes")
    print("=" * 50)

    return adata
