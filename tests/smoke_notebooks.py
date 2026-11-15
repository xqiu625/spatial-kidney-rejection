"""Run every public notebook on a deterministic, small synthetic Visium archive."""
import gzip
import io
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import anndata as ad
import numpy as np
import pandas as pd
from PIL import Image
from scipy import sparse
from scipy.io import mmwrite
import yaml

from src.prepare_data import prepare_archive
from src.utils import get_kidney_markers


def main():
    with tempfile.TemporaryDirectory(prefix='kidney-smoke-') as directory:
        root = Path(directory)
        rng = np.random.default_rng(51)
        pairs = pd.read_csv(ROOT / 'config/ligand_receptor_pairs.csv')
        genes = sorted({g for markers in get_kidney_markers().values() for g in markers} |
                       set(pairs.ligand) | set(pairs.receptor)) + [f'background_{i}' for i in range(400)]
        coords = np.array([(10 + 10 * x, 10 + 10 * y) for x in range(8) for y in range(8)])
        barcodes = [f'barcode_{i}' for i in range(64)]
        archive = root / 'synthetic.tar'
        with tarfile.open(archive, 'w') as tar:
            for sample in ['s1', 's2']:
                counts = rng.poisson(2, size=(64, len(genes))).astype(float)
                for group in range(4):
                    counts[group * 16:(group + 1) * 16, group * 20:(group + 1) * 20] += 15
                buffer = io.BytesIO()
                mmwrite(buffer, sparse.coo_matrix(counts.T))
                positions = pd.DataFrame({'barcode': barcodes, 'in_tissue': 1,
                    'array_row': np.repeat(np.arange(8), 8), 'array_col': np.tile(np.arange(8), 8),
                    'pxl_row_in_fullres': coords[:, 1], 'pxl_col_in_fullres': coords[:, 0]})
                image = io.BytesIO()
                Image.fromarray(np.full((100, 100, 3), 245, dtype=np.uint8)).save(image, format='PNG')
                contents = {'matrix.mtx': buffer.getvalue(),
                    'barcodes.tsv': ('\n'.join(barcodes) + '\n').encode(),
                    'features.tsv': ''.join(f'{g}\t{g}\tGene Expression\n' for g in genes).encode(),
                    'tissue_positions.csv': positions.to_csv(index=False).encode(),
                    'scalefactors_json.json': b'{"tissue_hires_scalef": 1, "tissue_lowres_scalef": 1, "spot_diameter_fullres": 6}',
                    'tissue_hires_image.png': image.getvalue(), 'tissue_lowres_image.png': image.getvalue()}
                for name, value in contents.items():
                    compressed = gzip.compress(value)
                    member = tarfile.TarInfo(f'{sample}_{name}.gz')
                    member.size = len(compressed)
                    tar.addfile(member, io.BytesIO(compressed))
        raw = root / 'raw'
        prepare_archive(archive, raw, ['s1', 's2'])
        config = root / 'config'
        config.mkdir()
        (config / 'paths.yaml').write_text(yaml.safe_dump({'base_dir': str(root), 'data': {'raw': str(raw)},
            'run_dir': str(root / 'run'), 'samples': [
                {'name': 's1', 'condition': 'control', 'rejection_type': 'none'},
                {'name': 's2', 'condition': 'rejection', 'rejection_type': 'synthetic'}]}))
        params = yaml.safe_load((ROOT / 'config/analysis_params.yaml').read_text())
        params['preprocessing'].update(min_genes=50, n_top_genes=200)
        params['clustering'].update(n_pcs=12, n_neighbors=8, top_markers=20)
        params['spatial'].update(coord_type='generic', n_neighs=4, n_perms=9, distance_bins=5)
        params['ligand_receptor']['n_perms'] = 9
        params['marker_profiles']['n_clusters'] = 4
        (config / 'analysis_params.yaml').write_text(yaml.safe_dump(params))
        run = root / 'run'
        subprocess.run([sys.executable, '-m', 'src.run_notebooks', '--run-dir', str(run),
                        '--config-dir', str(config)], cwd=ROOT, check=True,
                       env={**os.environ, 'MPLBACKEND': 'Agg', 'NUMBA_NUM_THREADS': '1'})
        final = ad.read_h5ad(run / 'processed/adata_final.h5ad')
        assert final.n_obs == 128 and final.obs_names.is_unique
        assert not any(column.endswith('_prop') for column in final.obs)
        assert set(final.layers) >= {'counts', 'lognorm'}
        markers = pd.read_csv(run / 'tables/cluster_markers.csv')
        assert np.isfinite(markers.logfoldchanges).all()
        lr = pd.read_csv(run / 'tables/spatial_lr_all.csv')
        assert lr.padj.between(0, 1).all() and lr['sample'].nunique() == 2
        assert len(list((run / 'notebooks').glob('*.ipynb'))) == 5
        print('All five notebooks completed; saved outputs passed integration assertions.', flush=True)


if __name__ == '__main__':
    main()
