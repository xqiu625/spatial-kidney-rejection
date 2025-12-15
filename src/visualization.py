"""
Visualization functions for spatial transcriptomics data
"""

import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Optional, Union, Tuple
import warnings

warnings.filterwarnings('ignore')

# Set style
plt.style.use('default')
sns.set_palette("husl")


def plot_qc_metrics(
    adata: sc.AnnData,
    save_path: Optional[Union[str, Path]] = None,
    figsize: Tuple[int, int] = (15, 10)
) -> plt.Figure:
    """
    Plot QC metrics for spatial data.

    Parameters
    ----------
    adata : AnnData
        Data with QC metrics
    save_path : str or Path, optional
        Path to save figure
    figsize : tuple
        Figure size

    Returns
    -------
    Figure
        Matplotlib figure
    """
    fig, axes = plt.subplots(2, 3, figsize=figsize)

    # Total counts distribution
    sns.histplot(data=adata.obs, x='total_counts', hue='sample', ax=axes[0, 0])
    axes[0, 0].set_title('Total Counts per Spot')
    axes[0, 0].set_xlabel('Total Counts')

    # Genes detected
    sns.histplot(data=adata.obs, x='n_genes_by_counts', hue='sample', ax=axes[0, 1])
    axes[0, 1].set_title('Genes Detected per Spot')
    axes[0, 1].set_xlabel('Number of Genes')

    # MT percentage
    sns.histplot(data=adata.obs, x='pct_counts_mt', hue='sample', ax=axes[0, 2])
    axes[0, 2].set_title('Mitochondrial %')
    axes[0, 2].set_xlabel('% MT')

    # Violin plots by condition
    if 'condition' in adata.obs.columns:
        sc.pl.violin(adata, ['total_counts'], groupby='condition', ax=axes[1, 0], show=False)
        sc.pl.violin(adata, ['n_genes_by_counts'], groupby='condition', ax=axes[1, 1], show=False)
        sc.pl.violin(adata, ['pct_counts_mt'], groupby='condition', ax=axes[1, 2], show=False)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved to {save_path}")

    return fig


def plot_spatial_feature(
    adata: sc.AnnData,
    feature: str,
    sample: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
    title: Optional[str] = None,
    spot_size: float = 1.3,
    cmap: str = 'viridis',
    **kwargs
) -> plt.Axes:
    """
    Plot a feature on spatial coordinates.

    Parameters
    ----------
    adata : AnnData
        Spatial data
    feature : str
        Feature to plot (gene name or obs column)
    sample : str, optional
        Sample to plot (if multi-sample)
    ax : Axes, optional
        Matplotlib axes
    title : str, optional
        Plot title
    spot_size : float
        Size of spots
    cmap : str
        Colormap

    Returns
    -------
    Axes
        Plot axes
    """
    if sample is not None:
        adata_plot = adata[adata.obs['sample'] == sample].copy()
        lib_id = sample
    else:
        adata_plot = adata
        # Get library_id from uns if single sample
        lib_id = list(adata_plot.uns.get('spatial', {}).keys())[0] if 'spatial' in adata_plot.uns else None

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))

    sq.pl.spatial_scatter(
        adata_plot,
        color=feature,
        library_id=lib_id,
        size=spot_size,
        cmap=cmap,
        title=title or feature,
        ax=ax,
        **kwargs
    )

    ax.set_aspect('equal')

    return ax


def plot_spatial_multi_feature(
    adata: sc.AnnData,
    features: List[str],
    sample: Optional[str] = None,
    ncols: int = 3,
    spot_size: float = 1.3,
    figsize_per_panel: Tuple[float, float] = (5, 5),
    save_path: Optional[Union[str, Path]] = None,
    **kwargs
) -> plt.Figure:
    """
    Plot multiple features on spatial coordinates.

    Parameters
    ----------
    adata : AnnData
        Spatial data
    features : list of str
        Features to plot
    sample : str, optional
        Sample to plot
    ncols : int
        Number of columns
    spot_size : float
        Size of spots
    figsize_per_panel : tuple
        Size per panel
    save_path : str or Path, optional
        Path to save figure

    Returns
    -------
    Figure
        Matplotlib figure
    """
    nrows = int(np.ceil(len(features) / ncols))
    figsize = (figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes = axes.flatten() if nrows * ncols > 1 else [axes]

    for i, feature in enumerate(features):
        plot_spatial_feature(
            adata, feature, sample=sample,
            ax=axes[i], spot_size=spot_size, **kwargs
        )

    # Hide empty axes
    for i in range(len(features), len(axes)):
        axes[i].axis('off')

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved to {save_path}")

    return fig


def plot_spatial_clusters(
    adata: sc.AnnData,
    cluster_key: str = 'leiden',
    samples: Optional[List[str]] = None,
    ncols: int = 4,
    spot_size: float = 1.3,
    save_path: Optional[Union[str, Path]] = None,
    **kwargs
) -> plt.Figure:
    """
    Plot clusters for all samples.

    Parameters
    ----------
    adata : AnnData
        Spatial data with clusters
    cluster_key : str
        Key for cluster labels
    samples : list of str, optional
        Samples to plot
    ncols : int
        Number of columns
    spot_size : float
        Size of spots
    save_path : str or Path, optional
        Path to save figure

    Returns
    -------
    Figure
        Matplotlib figure
    """
    if samples is None:
        samples = adata.obs['sample'].unique().tolist()

    nrows = int(np.ceil(len(samples) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows))
    axes = axes.flatten() if nrows * ncols > 1 else [axes]

    for i, sample in enumerate(samples):
        adata_sample = adata[adata.obs['sample'] == sample].copy()
        sq.pl.spatial_scatter(
            adata_sample,
            color=cluster_key,
            library_id=sample,
            size=spot_size,
            title=sample,
            ax=axes[i],
            **kwargs
        )
        axes[i].set_aspect('equal')

    for i in range(len(samples), len(axes)):
        axes[i].axis('off')

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved to {save_path}")

    return fig


def plot_deconvolution_results(
    adata: sc.AnnData,
    cell_types: List[str],
    sample: Optional[str] = None,
    ncols: int = 4,
    save_path: Optional[Union[str, Path]] = None,
    **kwargs
) -> plt.Figure:
    """
    Plot cell type deconvolution results.

    Parameters
    ----------
    adata : AnnData
        Data with deconvolution results in obsm
    cell_types : list of str
        Cell types to plot
    sample : str, optional
        Sample to plot
    ncols : int
        Number of columns
    save_path : str or Path, optional
        Path to save figure

    Returns
    -------
    Figure
        Matplotlib figure
    """
    return plot_spatial_multi_feature(
        adata, cell_types, sample=sample,
        ncols=ncols, save_path=save_path,
        cmap='Reds', **kwargs
    )


def plot_nhood_enrichment(
    adata: sc.AnnData,
    cluster_key: str = 'leiden',
    save_path: Optional[Union[str, Path]] = None,
    figsize: Tuple[int, int] = (8, 8)
) -> plt.Figure:
    """
    Plot neighborhood enrichment heatmap.

    Parameters
    ----------
    adata : AnnData
        Data with neighborhood enrichment computed
    cluster_key : str
        Key for cluster labels
    save_path : str or Path, optional
        Path to save figure
    figsize : tuple
        Figure size

    Returns
    -------
    Figure
        Matplotlib figure
    """
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    sq.pl.nhood_enrichment(
        adata,
        cluster_key=cluster_key,
        ax=ax
    )

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved to {save_path}")

    return fig


def plot_ligand_receptor(
    adata: sc.AnnData,
    cluster_key: str = 'leiden',
    top_n: int = 20,
    save_path: Optional[Union[str, Path]] = None,
    figsize: Tuple[int, int] = (12, 8)
) -> plt.Figure:
    """
    Plot ligand-receptor analysis results.

    Parameters
    ----------
    adata : AnnData
        Data with L-R analysis results
    cluster_key : str
        Key for cluster labels
    top_n : int
        Number of top interactions to show
    save_path : str or Path, optional
        Path to save figure
    figsize : tuple
        Figure size

    Returns
    -------
    Figure
        Matplotlib figure
    """
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    sq.pl.ligrec(
        adata,
        cluster_key=cluster_key,
        source_groups=None,
        target_groups=None,
        ax=ax
    )

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved to {save_path}")

    return fig


def create_publication_figure(
    panels: List[plt.Figure],
    layout: Tuple[int, int],
    labels: Optional[List[str]] = None,
    figsize: Optional[Tuple[int, int]] = None,
    save_path: Optional[Union[str, Path]] = None
) -> plt.Figure:
    """
    Combine multiple panels into a publication-ready figure.

    Parameters
    ----------
    panels : list of Figure
        Panel figures
    layout : tuple
        (nrows, ncols) for layout
    labels : list of str, optional
        Panel labels (A, B, C, ...)
    figsize : tuple, optional
        Figure size
    save_path : str or Path, optional
        Path to save figure

    Returns
    -------
    Figure
        Combined figure
    """
    nrows, ncols = layout
    if figsize is None:
        figsize = (6 * ncols, 6 * nrows)

    if labels is None:
        labels = [chr(65 + i) for i in range(len(panels))]  # A, B, C, ...

    fig = plt.figure(figsize=figsize)

    for i, (panel, label) in enumerate(zip(panels, labels)):
        ax = fig.add_subplot(nrows, ncols, i + 1)
        ax.text(-0.1, 1.1, label, transform=ax.transAxes,
                fontsize=16, fontweight='bold', va='top', ha='right')

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved to {save_path}")

    return fig
