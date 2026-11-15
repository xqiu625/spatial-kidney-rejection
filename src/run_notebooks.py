"""Execute the five analysis notebooks with a fresh, explicit run directory."""
import argparse
import os
from pathlib import Path

import nbformat
from nbclient import NotebookClient


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True)
    parser.add_argument('--config-dir', type=Path)
    args = parser.parse_args()
    run = args.run_dir.resolve()
    os.environ['SKR_RUN_DIR'] = str(run)
    if args.config_dir:
        os.environ['SKR_CONFIG_DIR'] = str(args.config_dir.resolve())
    root = Path(__file__).resolve().parents[1]
    notebooks = sorted((root / 'notebooks').glob('[0-9][0-9]_*.ipynb'))
    for path in notebooks:
        print(f'Executing {path.name}', flush=True)
        notebook = nbformat.read(path, as_version=4)
        NotebookClient(notebook, timeout=1800, kernel_name='python3',
                       resources={'metadata': {'path': str(root / 'notebooks')}}).execute()
        output = run / 'notebooks'
        output.mkdir(exist_ok=True)
        nbformat.write(notebook, output / path.name)
        print(f'Completed {path.name}', flush=True)


if __name__ == '__main__':
    main()
