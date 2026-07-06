from pathlib import Path
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

db = snakemake.params.get("db", "")
db_custom = snakemake.params.get("db_custom", "")
if db and db_custom:
    raise ValueError(
        "Can choose only one option. Or upimapi options (db) or custom db path (db_custom)"
    )

resources_dir = snakemake.params.get("resources_dir")

skip_db_check_if_exists = snakemake.params.get("skip_db_check_if_exists", True)
skip_db = snakemake.params.get("skip_db", False)

db2file = {
    "uniprot": Path(resources_dir) / "uniprot.fasta",
    "swissprot": Path(resources_dir) / "uniprot_sprot.fasta",
    "taxids": Path(resources_dir) / "taxids_database.fasta",
}

db_custom_exists = bool(db_custom) and Path(db_custom).exists()
db_exists = bool(db) and db in db2file and db2file[db].exists()

if skip_db or (skip_db_check_if_exists and (db_custom_exists or db_exists)):
    skip_db = "--skip-db-check"
else:
    skip_db = ""

if db:
    db = f"--database {db}"

if db_custom:
    db_custom = f"--database {db_custom}"

if resources_dir:
    resources_dir = f"-rd {resources_dir}"

tsv_output = [out for out in snakemake.output if out.endswith(".tsv")][0]
outdir = Path(tsv_output).parent

shell(
    "upimapi "
    "--input {snakemake.input.fasta} "
    "--output {outdir} "
    "{db} "
    "{db_custom} "
    "{resources_dir} "
    "{skip_db} "
    "--threads {snakemake.threads} "
    "{extra} "
    "{log}"
)
