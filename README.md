# Spatial Transcriptomics Analysis of Kidney Transplant Rejection

[![Python 3.10+](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Scanpy](https://img.shields.io/badge/Scanpy-1.9+-orange.svg)](https://scanpy.readthedocs.io/)
[![Squidpy](https://img.shields.io/badge/Squidpy-1.3+-purple.svg)](https://squidpy.readthedocs.io/)
[![Paper](https://img.shields.io/badge/Paper-Frontiers_in_Immunology-red.svg)](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2025.1654741/full)

> **Advanced spatial transcriptomics pipeline for analyzing immune microenvironments in kidney allograft rejection using 10X Visium data.**

This project extends our [published research in Frontiers in Immunology](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2025.1654741/full) (Data: [GSE304669](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304669)) with advanced computational methods including spatial neighborhood analysis, cell-cell communication inference, and condition-specific ligand-receptor profiling.

---

## Key Results

| Metric | Value |
|--------|-------|
| **Spots Analyzed** | 3,431 |
| **Genes** | 18,027 |
| **Samples** | 4 (1 Control + 3 Rejection Types) |
| **Spatial Clusters** | 11 |
| **Tissue Niches** | 6 |
| **L-R Pairs Tested** | 50 kidney-relevant pairs |

---

## Analysis Pipeline

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                     SPATIAL TRANSCRIPTOMICS PIPELINE                         │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│   ┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐             │
│   │  01_QC   │───▶│02_Cluster│───▶│03_CellType│───▶│04_Niches │             │
│   │ Loading  │    │  Leiden  │    │  Scoring  │    │ Analysis │             │
│   └──────────┘    └──────────┘    └──────────┘    └──────────┘             │
│        │               │               │               │                    │
│        ▼               ▼               ▼               ▼                    │
│   Filter spots    Spatial PCA    Marker-based    Neighborhood              │
│   Normalize       UMAP/Leiden    cell scoring    enrichment                │
│   HVG selection   11 clusters    8 cell types    Co-occurrence             │
│                                                                              │
│                           ┌──────────┐                                      │
│                           │05_LigRec │                                      │
│                           │ Analysis │                                      │
│                           └──────────┘                                      │
│                                │                                            │
│                                ▼                                            │
│                    Custom L-R interaction scoring                           │
│                    Permutation-based p-values                               │
│                    Control vs Rejection comparison                          │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Technical Highlights

### Methods Implemented

| Analysis | Method | Description |
|----------|--------|-------------|
| **Spatial Clustering** | Leiden + Spatial PCA | Graph-based clustering with spatial context |
| **Cell Type Scoring** | Marker-based | Custom scoring for 8 kidney cell types |
| **Niche Identification** | Squidpy | Neighborhood enrichment & co-occurrence |
| **L-R Communication** | Custom permutation | 500 permutations for statistical significance |
| **Differential Analysis** | Log2FC comparison | Control vs rejection L-R changes |

### Cell Types Analyzed

- **Proximal Tubule** (CUBN, LRP2, SLC34A1)
- **Distal Tubule** (SLC12A3, CALB1)
- **Collecting Duct** (AQP2, AQP3)
- **Endothelial** (PECAM1, VWF, CD34)
- **Immune Cells** (PTPRC, CD68, CD3D)
- **Fibroblasts** (COL1A1, DCN, VIM)
- **Podocytes** (NPHS1, NPHS2, WT1)
- **Mesangial** (PDGFRB, ACTA2)

### Ligand-Receptor Pairs Profiled

Curated kidney-relevant interactions including:
- **Immune signaling**: CCL2-CCR2, CXCL10-CXCR3, TNF-TNFRSF1A
- **Fibrosis**: TGFB1-TGFBR1/2, PDGFB-PDGFRB
- **Complement**: C3-C3AR1, C5-C5AR1
- **Immune checkpoints**: CD274-PDCD1 (PD-L1/PD-1)
- **Antigen presentation**: HLA-A-CD8A, HLA-DRA-CD4

---

## Repository Structure

```
spatial-kidney-rejection/
│
├── notebooks/                    # Analysis notebooks (run sequentially)
│   ├── 01_data_loading_qc.ipynb       # Data loading, QC, normalization
│   ├── 02_spatial_clustering.ipynb    # Leiden clustering, UMAP
│   ├── 03_cell2location_deconv.ipynb  # Cell type scoring
│   ├── 04_spatial_niche_analysis.ipynb # Spatial neighborhoods
│   └── 05_ligand_receptor.ipynb       # L-R communication analysis
│
├── src/                          # Reusable Python modules
│   ├── preprocessing.py               # Data loading & QC functions
│   ├── visualization.py               # Spatial plotting functions
│   └── utils.py                       # Helper functions & sample info
│
├── config/                       # Configuration files
│   ├── analysis_params.yaml           # Analysis parameters
│   └── paths.yaml                     # Data paths
│
├── results/                      # Analysis outputs
│   ├── figures/                       # Generated visualizations
│   └── tables/                        # CSV results
│
├── environment.yml               # Conda environment
├── requirements.txt              # Pip dependencies
└── README.md
```

---

## Sample Visualizations

### Spatial Clustering
![Spatial Clusters](results/figures/spatial_clusters_all_samples.png)
*Leiden clusters mapped to tissue coordinates across all samples*

### Spatial Niches
![Spatial Niches](results/figures/spatial_niches_all_samples.png)
*Tissue microenvironment niches identified across samples*

### Cell Type Distribution
![Cell Types](results/figures/celltype_spatial_GSM9155022_control_58055.png)
*Marker-based cell type scoring showing spatial distribution*

### Neighborhood Enrichment
![Niche Enrichment](results/figures/nhood_enrichment_clusters.png)
*Cell type co-localization patterns in the tissue microenvironment*

### L-R Communication Network
![Communication Network](results/figures/communication_network.png)
*Cell-cell communication network based on ligand-receptor interactions*

### Top Ligand-Receptor Interactions
![Top LR Pairs](results/figures/top_LR_pairs.png)
*Ranked ligand-receptor interactions by interaction score*

---

## Quick Start

### 1. Clone Repository
```bash
git clone https://github.com/yourusername/spatial-kidney-rejection.git
cd spatial-kidney-rejection
```

### 2. Create Environment
```bash
# Using conda/mamba
mamba env create -f environment.yml
conda activate spatial-kidney

# Or using pip
pip install -r requirements.txt
```

### 3. Download Data
```bash
# Data from GEO: GSE304669
# See notebooks/01_data_loading_qc.ipynb for download instructions
```

### 4. Run Analysis
```bash
# Launch Jupyter
jupyter notebook

# Run notebooks 01-05 sequentially
```

---

## Key Findings

### Spatial Niches in Rejection

The analysis identified **6 distinct tissue niches** with differential cell type compositions between control and rejection samples:

1. **Tubular-dominant niches** - Proximal/distal tubule enriched
2. **Immune-infiltrated zones** - High immune cell density in rejection
3. **Fibrotic regions** - Fibroblast accumulation in chronic rejection
4. **Vascular niches** - Endothelial-rich areas
5. **Mixed cellular regions** - Transitional zones
6. **Glomerular regions** - Podocyte and mesangial enriched

### Differential L-R Communication

Comparing rejection vs control samples revealed:
- **Upregulated in rejection**: Inflammatory chemokine signaling (CXCL9/10-CXCR3)
- **Immune checkpoint changes**: Altered PD-L1/PD-1 axis
- **Fibrotic signaling**: Enhanced TGF-β pathway activity

---

## Technologies & Skills

| Category | Technologies |
|----------|-------------|
| **Languages** | Python 3.10+, R |
| **Spatial Analysis** | Scanpy, Squidpy, AnnData |
| **Visualization** | Matplotlib, Seaborn, NetworkX |
| **Statistics** | Permutation testing, Multiple testing correction |
| **Data Formats** | 10X Visium, H5AD, MTX |
| **Environment** | Conda, Jupyter, HPC (SLURM) |

---

## Data Source & Citation

- **Publication**: [Frontiers in Immunology (2025)](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2025.1654741/full)
- **GEO Accession**: [GSE304669](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304669)
- **Platform**: 10X Genomics Visium FFPE
- **Tissue**: Human kidney transplant biopsies
- **Conditions**: Control (n=1), Active AMR (n=1), Acute TCMR (n=1), Chronic Active AMR (n=1)

---

## Future Directions

- [ ] Integration with reference scRNA-seq for probabilistic deconvolution (Cell2location)
- [ ] Spatial-aware clustering with BayesSpace
- [ ] Trajectory analysis of rejection progression
- [ ] Deep learning-based spatial pattern recognition

---

## Author

**Xinru Qiu**
Computational Biology | Spatial Transcriptomics | Machine Learning

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- Original data from kidney transplant rejection study (GSE304669)
- Scanpy and Squidpy development teams
- 10X Genomics for Visium technology

