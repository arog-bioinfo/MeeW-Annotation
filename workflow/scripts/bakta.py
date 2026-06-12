"""Snakemake wrapper for Bakta."""

__author__ = "Rodolfo Brandão Dias Ferreira"
__copyright__ = "Copyright 2026, Rodolfo Brandão"
__email__ = "rodolfobrandao88@gmail.com"
__license__ = "MIT"

from pathlib import Path
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

db = snakemake.params.get("db", "")
db_cmd = f"--db {db}" if db else ""

proteins = snakemake.input.get("proteins", "")
proteins_cmd = f"--proteins {proteins}" if proteins else ""

outdir = snakemake.output.outdir
prefix = Path(snakemake.input.fasta).stem

shell(
    "bakta "
    "--output {outdir} "
    "--prefix {prefix} "
    "--force "
    "--threads {snakemake.threads} "
    "{db_cmd} "
    "{proteins_cmd} "
    "{extra} "
    "{snakemake.input.fasta} "
    "{log}"
)
