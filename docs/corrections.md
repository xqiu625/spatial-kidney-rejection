# Analysis validity corrections

The review baseline was `72856f7`. The original results are preserved in `results/legacy_20260724`; those outputs are superseded.

| Review finding | Change |
|---|---|
| Scaled expression used for marker fold changes | Keep counts/lognorm separate; scale only a temporary PCA matrix; marker testing explicitly uses lognorm |
| Spatial coordinates pooled across biopsies | Build a sample-separated graph and reject cross-sample edges; run enrichment/co-occurrence separately per biopsy |
| LR calculation ignored the graph | Score products on actual same-biopsy edges; no-edge combinations score zero |
| Raw-count LR depth confounding | LR reads only lognorm; regression check changes count depth without changing relative expression |
| FDR absent and wrong number of hypotheses stated | Save the full family with global BH padj and its size; use adjusted screening results in plots |
| Spots treated as biological replicates | Remove disease hypothesis tests; summarize each biopsy and subtype descriptively |
| Scores presented as cell fractions | Keep background-corrected scores; rename composition groups as marker profiles |
| Broken reference branch | Replace notebook 03 with marker scoring; explicitly reject reference methods until validated |
| Methods documentation did not describe the implementation | Rewrite the README, diagram, notebook explanations and output names |
| Hard-coded paths and unused configuration | Load existing YAML, resolve paths, support fresh run directories and record configuration/versions |
| Whole AnnData copied and all genes densified per permutation | Extract only LR genes once; reuse edge products and permute only labels |
| Empty publication helper / ignored top_n | Remove the unused empty-panel helper and old deconvolution wrapper; use top_n in the actual LR table plot |

No new Cell2location training workflow or patient-level inferential model was added. The data cannot support those claims merely through a code change.

## Validation

Validation results and the real-data run summary will be recorded after execution. A complete run must finish all five notebooks and pass graph, expression and output checks; syntax/import success alone is insufficient.
