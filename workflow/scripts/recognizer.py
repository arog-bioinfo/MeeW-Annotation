"""Snakemake wrapper for reCOGnizer."""

__author__ = "Rodolfo Brandão Dias Ferreira"
__copyright__ = "Copyright 2026, Rodolfo Brandão"
__email__ = "rodolfobrandao88@gmail.com"
__license__ = "MIT"

from pathlib import Path
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

custom_db_input = snakemake.input.get("custom_db", "")
custom_db_param = snakemake.params.get("custom_db", "")

if custom_db_input:
    custom_db_cmd = f"--custom-databases -dbs {custom_db_input}"
elif custom_db_param:
    custom_db_cmd = f"--custom-databases -dbs {custom_db_param}"
else:
    custom_db_cmd = ""

resources_dir = snakemake.params.get("resources_dir", "")
resources_cmd = f"-rd {resources_dir}" if resources_dir else ""

first_output = Path(snakemake.output[0])
outdir = first_output.parent

shell(
    "recognizer "
    "-f {snakemake.input.fasta} "
    "-o {outdir} "
    "-t {snakemake.threads} "
    "{custom_db_cmd} "
    "{resources_cmd} "
    "{extra} "
    "{log}"
)