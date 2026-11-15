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
| Category colors shifted when a biopsy lacked a group | Assign the full category palette before slicing; test that absent groups do not change colors |

No new Cell2location training workflow or patient-level inferential model was added. The data cannot support those claims merely through a code change.

## Validation

The initial correction was verified on Python 3.12 with the environment recorded in [the results manifest](../results/run_manifest.json):

- 15 regression tests passed, including cross-biopsy category-color consistency.
- All five notebooks completed on a synthetic Visium archive and passed saved-output assertions.
- All five notebooks completed on the actual GSE304669 RAW archive at analysis code commit `75354df38d1679765dea7cb3dcfcbe07611009e1` with a clean working tree.
- The real run retained 3,431 spots and 18,027 genes across four biopsies.
- All 8,414 undirected spatial edges stayed within samples; none crossed biopsies.
- All 1,100 exported marker log fold changes were finite (top 100 per each of 11 clusters).
- The full LR table contained 8,505 hypotheses. Independent output checks reproduced its BH adjustment; all p-values and adjusted values were in [0,1]. Zero-edge combinations had score=0 and p=1. The filtered CSV exactly matched the full-table filter (228 rows).
- Full H5AD reload confirmed unique spot IDs, preserved integer counts, nonnegative finite lognorm, X matching lognorm, and no *_prop columns.
- Figures were inspected for labels and rendering. Some upstream dependencies emit deprecation warnings; these are visible rather than globally suppressed.

After the complete run, the two categorical tissue figures were rendered again from its unchanged final checkpoint with the palette correction in commit `d57fde1`. This changes the displayed colors only; the numerical tables are unchanged. The complete notebook execution used `75354df`, and the 15-test regression run includes the subsequent palette fix. Both revisions are recorded in the results manifest.

The first hosted integration run exposed a serialization regression in that palette correction: Seaborn's `ColorPalette` rendered correctly but could not be written into an H5AD checkpoint. Colors are now stored as a plain list. The color regression test also saves and reloads the checkpoint, covering the previously missed path. This fix leaves the color values and numerical results unchanged.

The compact metrics and portable provenance are published under [results](../results/README.md). The original archive SHA-256 is `c798a1065a04cf967cb388f4c03fd40aa0b15421ba115b8124674f59317f211c`.

GitHub Actions runs the regression and synthetic notebook checks on every push; consult the [workflow runs](https://github.com/xqiu625/spatial-kidney-rejection/actions/workflows/tests.yml) for the status of each revision. Reference training and additional biological replication were not performed. These checks verify implementation and reproducibility, not clinical inference or experimentally confirmed communication.

## Security dependency update — November 15, 2026

Updated PyArrow to 23.0.1, setuptools to 84.1.0 and pytest to 9.0.3. Removed the unused nbconvert dependency and its exclusive dependencies; notebook execution continues through nbclient. SpatialData 0.7.2 removes the obsolete xarray-schema/pkg_resources dependency, allowing the setuptools update. Its compatible stack requires Squidpy 1.8.0, AnnData 0.12.6, NumPy 2.2.6, Zarr 3.1.3 and related lock updates. Dask and Distributed are pinned to 2026.1.1, including the dashboard security fix while respecting SpatialData's version constraints.

In a fresh Python 3.12 environment, `pip check` passed, all 15 regression tests passed, and all five notebooks completed on the synthetic Visium fixture with saved-output assertions passing. An OSV query of all 161 locked package versions returned no known vulnerability matches on November 15, 2026. GitHub Actions also checks dependency consistency with `pip check`.

Analysis source code and the published real-data outputs are unchanged. The four-biopsy analysis was not rerun for this dependency update; its original versions and validation remain in the results manifest.
