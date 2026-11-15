from pathlib import Path

from snakemake.utils import format

from common.native_execution import ensure_lock_coverage, run_native_command
from common.upimapi_database import assess_upimapi_database
from common.upimapi_options import validate_upimapi_extra

extra = validate_upimapi_extra(snakemake.params.get("extra", ""))
Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

db = snakemake.params.get("db", "")
db_custom = snakemake.params.get("db_custom", "")
if db and db_custom:
    raise ValueError(
        "Can choose only one option. Or upimapi options (db) or custom db path (db_custom)"
    )

resources_dir = snakemake.params.get("resources_dir", "")

skip_db_check_if_exists = snakemake.params.get("skip_db_check_if_exists", True)
skip_db = snakemake.params.get("skip_db", False)

db2file = {
    "uniprot": Path(resources_dir) / "uniprot.fasta",
    "swissprot": Path(resources_dir) / "uniprot_sprot.fasta",
    "taxids": Path(resources_dir) / "taxids_database.fasta",
}


def original_skip_db_option():
    db_custom_exists = bool(db_custom) and Path(db_custom).exists()
    db_exists = bool(db) and db in db2file and db2file[db].exists()
    if skip_db or (skip_db_check_if_exists and (db_custom_exists or db_exists)):
        return "--skip-db-check"
    return ""


skip_db_option = original_skip_db_option()

db_option = "--database" if db else ""
db_custom_option = "--database" if db_custom else ""
resources_option = "-rd" if resources_dir else ""

fasta = (Path(snakemake.input.annotations) / str(snakemake.params.fasta_name)).resolve()
outdir = Path(snakemake.output.results).resolve()
outdir.mkdir(parents=True, exist_ok=True)

assessment = assess_upimapi_database(
    resources_dir,
    db,
    db_custom,
    skip_db_check=bool(skip_db_option),
    extra=extra,
    definitive=False,
)
launch_state = {"skip_db_option": skip_db_option}


def verify_locked_state(probe):
    locked_skip_db_option = original_skip_db_option()
    current = assess_upimapi_database(
        resources_dir,
        db,
        db_custom,
        skip_db_check=bool(locked_skip_db_option),
        extra=extra,
        definitive=True,
        diamond_probe=probe,
    )
    ensure_lock_coverage(assessment.locks, current.locks)
    launch_state["skip_db_option"] = locked_skip_db_option


def command():
    skip_db_option = launch_state["skip_db_option"]
    return format(
        "upimapi "
        "--input {fasta:q} "
        "--output {outdir:q} "
        "{db_option} {db:q} "
        "{db_custom_option} {db_custom:q} "
        "{resources_option} {resources_dir:q} "
        "{skip_db_option} "
        "--threads {snakemake.threads} "
        "--max-target-seqs 1 "
        "{extra} "
        "{log}"
    )


run_native_command(
    command,
    assessment.locks,
    prelaunch_probe_check=verify_locked_state,
)
