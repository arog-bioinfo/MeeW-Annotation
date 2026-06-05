"""Snakemake wrapper for UPIMAPI."""

__author__ = "Rodolfo Brandão Dias Ferreira"
__copyright__ = "Copyright 2026, Rodolfo Brandão"
__email__ = "rodolfobrandao88@gmail.com"
__license__ = "MIT"

from pathlib import Path
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

db_input = snakemake.input.get("db", "")
db_param = snakemake.params.get("db", "")

if db_input:
    db_cmd = f"--database {db_input}"
elif db_param:
    db_cmd = f"--database {db_param}"
else:
    db_cmd = "" 

tsv_output = [out for out in snakemake.output if out.endswith(".tsv")][0]
outdir = Path(tsv_output).parent

shell(
    "upimapi "
    "--input {snakemake.input.fasta} "
    "--output {outdir} "
    "{db_cmd} "
    "--threads {snakemake.threads} "
    "{extra} "
    "{log}"
)