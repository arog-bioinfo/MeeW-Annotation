"""Deterministic native-format outputs for inputs with no FASTA genomes."""

import csv
from pathlib import Path

FASTA_SUFFIXES = {".fa", ".fna", ".fasta"}

CHECKM2_COLUMNS = [
    "Name",
    "Completeness",
    "Contamination",
    "Completeness_Model_Used",
    "Translation_Table_Used",
    "Coding_Density",
    "Contig_N50",
    "Average_Gene_Length",
    "Genome_Size",
    "GC_Content",
    "Total_Coding_Sequences",
    "Total_Contigs",
    "Max_Contig_Length",
    "Additional_Notes",
]

EUKCC_COLUMNS = ["bin", "completeness", "contamination", "ncbi_lng"]


def fasta_paths(source: Path) -> list[Path]:
    """Return deterministic direct-child nucleotide FASTAs from a file/directory."""
    source = Path(source)
    if source.is_file():
        return [source] if source.suffix.lower() in FASTA_SUFFIXES else []
    if not source.exists():
        return []
    return sorted(
        path
        for path in source.iterdir()
        if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES
    )


def has_fasta_inputs(source: Path) -> bool:
    """Report whether a file/directory contains at least one nucleotide FASTA."""
    return bool(fasta_paths(source))


def write_delimited_header(
    path: Path, columns: list[str], *, delimiter: str = "\t"
) -> None:
    """Write a deterministic header-only native table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter=delimiter, lineterminator="\n")
        writer.writerow(columns)


def write_empty_log(path: Path, *, tool: str, domain: str) -> None:
    """Record why a scientific command was not invoked."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"status=empty_domain\ntool={tool}\ndomain={domain}\n"
        "scientific_command_invoked=false\n",
        encoding="utf-8",
    )
