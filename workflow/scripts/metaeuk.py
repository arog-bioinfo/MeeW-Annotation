"""Snakemake wrapper for MetaEuk."""

__author__ = "Rodolfo Brandão Dias Ferreira"
__copyright__ = "Copyright 2026, Rodolfo Brandão"
__email__ = "rodolfobrandao88@gmail.com"
__license__ = "MIT"

from pathlib import Path
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

first_output = Path(snakemake.output.proteins)
outdir = first_output.parent
prefix = first_output.stem

out_prefix = outdir / prefix
tmp_dir = outdir / "tmp"

tmp_dir.mkdir(parents=True, exist_ok=True)

# 1. Correr o MetaEuk
shell(
    "metaeuk easy-predict "
    "--threads {snakemake.threads} "
    "{extra} "
    "{snakemake.input.fasta} "
    "{snakemake.params.db} "
    "{out_prefix} "
    "{tmp_dir} "
    "{log}"
)

possible_extensions = [".fas", ".fasta", ".faa"]

for ext in possible_extensions:
    generated_file = out_prefix.with_suffix(ext)
    if generated_file.exists():
        if generated_file != first_output:
            generated_file.rename(first_output)
        break
