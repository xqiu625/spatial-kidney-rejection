# Corrected tables

See the [table index and interpretation](../README.md). The input archive, settings, library versions and code commit are recorded in [run_manifest.json](../run_manifest.json).

LR fields: `sample`, `condition`, `rejection_type`, `source`, `target`, `ligand`, `receptor`, `n_edges`, `score`, `pvalue`, `padj`, `n_perms`, `fdr_family_size`. Source/target are expression-cluster IDs, not inferred cell types. The score is a mean log-expression product over the indicated directed spatial edges. BH adjustment uses all 8,505 tests in the published run.

Cluster marker p-values are exploratory within-dataset cluster comparisons. No table supplies a patient-level test of rejection. Marker scores are not abundance estimates.
