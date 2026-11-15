"""Temporary, symlink-based staging for file- and directory-oriented QA tools."""

import tempfile
from contextlib import contextmanager
from pathlib import Path

from common.empty_outputs import fasta_paths


@contextmanager
def staged_fasta_directory(source: Path, extension: str | None = None):
    """Stage direct-child FASTAs as symlinks and remove the staging directory."""
    with tempfile.TemporaryDirectory(prefix="meew-annotation-qa-") as temporary:
        staging = Path(temporary)
        used_names = set()
        for fasta in fasta_paths(Path(source)):
            suffix = extension if extension is not None else fasta.suffix
            name = f"{fasta.stem}{suffix}"
            if name in used_names:
                raise ValueError(
                    f"QA input names collide after extension normalization: {name}"
                )
            used_names.add(name)
            (staging / name).symlink_to(fasta.resolve())
        yield staging
