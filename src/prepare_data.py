"""Arrange an existing GEO RAW tar without extracting arbitrary archive paths."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tarfile

import yaml


def prepare_archive(archive, destination, sample_names):
    archive, destination = Path(archive), Path(destination)
    expected = {
        'matrix.mtx': 'filtered_feature_bc_matrix/matrix.mtx',
        'features.tsv': 'filtered_feature_bc_matrix/features.tsv',
        'barcodes.tsv': 'filtered_feature_bc_matrix/barcodes.tsv',
        'tissue_positions.csv': 'spatial/tissue_positions.csv',
        'scalefactors_json.json': 'spatial/scalefactors_json.json',
        'tissue_hires_image.png': 'spatial/tissue_hires_image.png',
        'tissue_lowres_image.png': 'spatial/tissue_lowres_image.png',
    }
    if not sample_names or len(sample_names) != len(set(sample_names)):
        raise ValueError('Sample names must be nonempty and unique')
    if any(Path(name).name != name or name in {'.', '..'} for name in sample_names):
        raise ValueError('Sample names must be plain directory names')
    if destination.exists() and any(destination.iterdir()):
        raise FileExistsError(f'Use an empty destination: {destination}')
    with tarfile.open(archive) as tar:
        members = []
        for sample in sample_names:
            for suffix, relative in expected.items():
                member = tar.getmember(f'{sample}_{suffix}.gz')
                if not member.isfile():
                    raise ValueError(f'Expected a regular file: {member.name}')
                members.append((member, destination / sample / relative))
        for member, target in members:
            target.parent.mkdir(parents=True, exist_ok=True)
            with tar.extractfile(member) as compressed, gzip.GzipFile(fileobj=compressed) as data, target.open('xb') as output:
                shutil.copyfileobj(data, output)
    with archive.open('rb') as handle:
        checksum = hashlib.file_digest(handle, 'sha256').hexdigest()
    (destination / 'source.json').write_text(json.dumps({'accession': 'GSE304669', 'archive_sha256': checksum,
                                                       'samples': sample_names}, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path('data/raw'))
    parser.add_argument('--paths', type=Path, default=Path(__file__).resolve().parents[1] / 'config/paths.yaml')
    args = parser.parse_args()
    settings = yaml.safe_load(args.paths.read_text())
    prepare_archive(args.archive, args.output, [sample['name'] for sample in settings['samples']])
    print(f'Prepared {len(settings["samples"])} samples in {args.output}')


if __name__ == '__main__':
    main()
