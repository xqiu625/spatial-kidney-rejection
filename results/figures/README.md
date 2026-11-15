# Regenerated figures

All eight PNGs come from the corrected real-data workflow. The two categorical tissue figures were rendered again from the final checkpoint after fixing cross-biopsy color consistency. See [the run manifest](../run_manifest.json) for both code revisions.

- `expression_umap.png`: expression embedding colored by Leiden cluster and biopsy.
- `expression_clusters_on_tissue.png`: expression clusters mapped to each biopsy's own coordinates.
- `marker_scores_by_sample.png`: signed background-corrected marker scores, not cell fractions.
- `marker_scores_on_tissue.png`: marker scores in the control biopsy; each panel has its own numeric color scale.
- `marker_profiles_on_tissue.png`: configured KMeans profile groups; no anatomical niche labels are asserted.
- `neighborhood_enrichment_by_sample.png`: separate within-biopsy permutation Z-scores; panels have individual color scales.
- `spatial_lr_fdr.png`: top 20 scores among the BH-filtered exploratory associations.
- `qc.png`: descriptive QC summaries; violin plots do not perform a disease test.
