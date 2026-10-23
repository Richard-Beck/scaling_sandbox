"""Keep development PNG figures visible and supporting exports out of their way.

Artifact keys in manifests are relative to the development directory. Manifests
live in supporting/metadata and declare artifact_root='../..'.
"""
from pathlib import Path
import os


def dev_output_path(directory, filename):
    directory = Path(directory)
    filename = Path(filename)
    if filename.name != str(filename):
        raise ValueError('Expected a single development export filename')
    extension = filename.suffix.lower()
    category = {
        '.pdf': 'vectors', '.svg': 'vectors',
        '.csv': 'data', '.npz': 'data',
        '.json': 'metadata', '.txt': 'metadata', '.md': 'metadata',
        '.html': 'viewers',
    }.get(extension)
    if extension == '.png':
        path = directory / filename
    elif category:
        path = directory / 'supporting' / category / filename
    else:
        raise ValueError(f'Unclassified development export: {filename}')
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def dev_artifact_key(directory, path):
    return str(Path(path).relative_to(directory))


def dev_viewer_link(directory, filename):
    return os.path.relpath(dev_output_path(directory, filename),
                           Path(directory) / 'supporting' / 'viewers')
