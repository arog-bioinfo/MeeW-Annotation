"""Snakemake wrapper for GTDB-Tk classify_wf."""

__author__ = "Rodolfo Brandão Dias Ferreira"
__copyright__ = "Copyright 2026, Rodolfo Brandão"
__email__ = "rodolfobrandao88@gmail.com" 
__license__ = "MIT"

import os
from pathlib import Path
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

db_path = snakemake.params.get("db_path", "")
if db_path:
    os.environ["GTDBTK_DATA_PATH"] = db_path

out_dir = Path(snakemake.output[0]).parent

extension = snakemake.params.get("extension", "fasta")

shell(
    "gtdbtk classify_wf"
    " --genome_dir {snakemake.input[0]:q}"
    " --out_dir {out_dir:q}"
    " --extension {extension}"
    " --cpus {snakemake.threads}"
    " --force"
    " {extra}"
    " {log}"
)