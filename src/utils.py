"""
Utility functions for spatial transcriptomics analysis
"""

import os
import yaml
import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Any, Optional, Union
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_config(config_path: Union[str, Path]) -> Dict[str, Any]:
    """
    Load YAML configuration file.

    Parameters
    ----------
    config_path : str or Path
        Path to config file

    Returns
    -------
    dict
        Configuration dictionary
    """
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def setup_directories(base_dir: Union[str, Path]) -> Dict[str, Path]:
    """
    Create standard directory structure.

    Parameters
    ----------
    base_dir : str or Path
        Base project directory

    Returns
    -------
    dict
        Dictionary of directory paths
    """
    base_dir = Path(base_dir)

    dirs = {
        'data_raw': base_dir / 'data' / 'raw',
        'data_processed': base_dir / 'data' / 'processed',
        'data_reference': base_dir / 'data' / 'reference',
        'results': base_dir / 'results',
        'tables': base_dir / 'results' / 'tables',
        'figures': base_dir / 'results' / 'figures',
        'notebooks': base_dir / 'notebooks',
    }

    for name, path in dirs.items():
        path.mkdir(parents=True, exist_ok=True)
        logger.info(f"Created/verified directory: {path}")

    return dirs


def load_or_compute(
    path: Union[str, Path],
    compute_fn,
    *args,
    force_recompute: bool = False,
    **kwargs
) -> sc.AnnData:
    """
    Load cached result or compute and save.

    Parameters
    ----------
    path : str or Path
        Path to cached file
    compute_fn : callable
        Function to compute result
    *args :
        Positional arguments for compute_fn
    force_recompute : bool
        If True, recompute even if cache exists
    **kwargs :
        Keyword arguments for compute_fn

    Returns
    -------
    AnnData
        Loaded or computed data
    """
    path = Path(path)

    if path.exists() and not force_recompute:
        logger.info(f"Loading cached data from {path}")
        return sc.read_h5ad(path)

    logger.info(f"Computing result...")
    result = compute_fn(*args, **kwargs)

    path.parent.mkdir(parents=True, exist_ok=True)
    result.write(path)
    logger.info(f"Saved result to {path}")

    return result


def get_sample_info() -> pd.DataFrame:
    """
    Get sample information for the kidney transplant dataset.

    Returns
    -------
    DataFrame
        Sample metadata
    """
    samples = [
        {'gsm_id': 'GSM9155022', 'name': 'control_58055', 'condition': 'control', 'rejection_type': 'none'},
        {'gsm_id': 'GSM9155023', 'name': 'active_AMR_58056', 'condition': 'rejection', 'rejection_type': 'active_AMR'},
        {'gsm_id': 'GSM9155024', 'name': 'acute_TCMR_58057', 'condition': 'rejection', 'rejection_type': 'acute_TCMR'},
        {'gsm_id': 'GSM9155025', 'name': 'chronic_active_AMR_58058', 'condition': 'rejection', 'rejection_type': 'chronic_active_AMR'},
    ]

    df = pd.DataFrame(samples)
    df['full_name'] = df['gsm_id'] + '_' + df['name']

    return df


def get_kidney_markers() -> Dict[str, list]:
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


def score_cell_types(
    adata: sc.AnnData,
    markers: Optional[Dict[str, list]] = None
) -> sc.AnnData:
    """
    Score cell types based on marker genes.

    Parameters
    ----------
    adata : AnnData
        Input data
    markers : dict, optional
        Cell type markers. If None, uses kidney markers.

    Returns
    -------
    AnnData
        Data with cell type scores in obs
    """
    if markers is None:
        markers = get_kidney_markers()

    for cell_type, genes in markers.items():
        # Filter to genes present in data
        genes_present = [g for g in genes if g in adata.var_names]

        if len(genes_present) > 0:
            sc.tl.score_genes(adata, genes_present, score_name=f'{cell_type}_score')
            logger.info(f"Scored {cell_type} with {len(genes_present)}/{len(genes)} markers")
        else:
            logger.warning(f"No markers found for {cell_type}")

    return adata


def differential_expression(
    adata: sc.AnnData,
    groupby: str = 'condition',
    groups: Optional[list] = None,
    reference: str = 'control',
    method: str = 'wilcoxon'
) -> pd.DataFrame:
    """
    Perform differential expression analysis.

    Parameters
    ----------
    adata : AnnData
        Input data
    groupby : str
        Column to group by
    groups : list, optional
        Groups to test
    reference : str
        Reference group
    method : str
        Statistical test method

    Returns
    -------
    DataFrame
        Differential expression results
    """
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        groups=groups,
        reference=reference,
        method=method
    )

    # Extract results to DataFrame
    result = sc.get.rank_genes_groups_df(adata, group=None)

    return result


def export_results(
    adata: sc.AnnData,
    output_dir: Union[str, Path],
    prefix: str = 'results'
):
    """
    Export analysis results to files.

    Parameters
    ----------
    adata : AnnData
        Data with analysis results
    output_dir : str or Path
        Output directory
    prefix : str
        File prefix
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Export obs (spot metadata)
    adata.obs.to_csv(output_dir / f'{prefix}_spot_metadata.csv')
    logger.info(f"Exported spot metadata")

    # Export var (gene metadata)
    adata.var.to_csv(output_dir / f'{prefix}_gene_metadata.csv')
    logger.info(f"Exported gene metadata")

    # Export h5ad
    adata.write(output_dir / f'{prefix}_adata.h5ad')
    logger.info(f"Exported AnnData object")
