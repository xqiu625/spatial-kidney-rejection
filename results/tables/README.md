# Analysis Results Tables

This directory contains CSV outputs from the analysis pipeline.

## Expected Outputs

| File | Notebook | Description |
|------|----------|-------------|
| `significant_LR_pairs.csv` | 05 | Significant ligand-receptor interactions |
| `LR_comparison_control_vs_rejection.csv` | 05 | Differential L-R analysis |
| `cluster_markers.csv` | 02 | Marker genes per cluster |
| `celltype_scores.csv` | 03 | Cell type scoring results |

## Data Dictionary

### significant_LR_pairs.csv
| Column | Description |
|--------|-------------|
| ligand | Ligand gene name |
| receptor | Receptor gene name |
| source | Source cluster |
| target | Target cluster |
| score | Interaction score (ligand_expr * receptor_expr) |
| pvalue | Permutation-based p-value |

### LR_comparison_control_vs_rejection.csv
| Column | Description |
|--------|-------------|
| ligand | Ligand gene name |
| receptor | Receptor gene name |
| source | Source cluster |
| target | Target cluster |
| score_ctrl | Interaction score in control |
| score_rej | Interaction score in rejection |
| log2FC | Log2 fold change (rejection/control) |
