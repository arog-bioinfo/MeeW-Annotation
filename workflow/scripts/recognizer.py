"""Snakemake wrapper for reCOGnizer."""

import os
from pathlib import Path

from snakemake.utils import format

from common.native_execution import ensure_lock_coverage, run_native_command
from common.recognizer_database import assess_recognizer_database
from common.recognizer_options import validate_recognizer_extra

extra = validate_recognizer_extra(
    snakemake.params.get("extra", ""),
    config_key=snakemake.params.get("extra_config_key", "config.recognizer_prok.extra"),
)
log_path = Path(snakemake.log[0]).resolve()
log_path.parent.mkdir(parents=True, exist_ok=True)
job_cwd = log_path.parent

resources_dir = snakemake.params.get("resources_dir", "")
resources_dir = Path(resources_dir).resolve() if resources_dir else ""
resources_option = "-rd" if resources_dir else ""
fasta = (Path(snakemake.input.annotations) / str(snakemake.params.fasta_name)).resolve()
outdir = Path(snakemake.output.results).resolve()
databases = str(snakemake.params.databases)
assessment = assess_recognizer_database(
    resources_dir, databases, extra, definitive=False
)


def verify_locked_state():
    current = assess_recognizer_database(
        resources_dir,
        databases,
        extra,
        executable=assessment.native_source,
        definitive=True,
    )
    ensure_lock_coverage(assessment.locks, current.locks)


original_cwd = Path.cwd()
try:
    os.chdir(job_cwd)
    command = format(
        "recognizer "
        "-f {fasta:q} "
        "-o {outdir:q} "
        "-t {snakemake.threads} "
        "-dbs {databases:q} "
        "{resources_option} {resources_dir:q} "
        "{extra} "
        "> {log_path:q} 2>&1"
    )
    run_native_command(
        command,
        assessment.locks,
        prelaunch_check=verify_locked_state,
    )
finally:
    os.chdir(original_cwd)
