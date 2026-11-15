import os
from pathlib import Path

from snakemake.shell import shell

from common.bakta_options import validate_bakta_extra
from common.bakta_database import ensure_bakta_database
from common.identifiers import validate_annotation_identifier

extra = snakemake.params.get("extra", "")
validate_bakta_extra(extra)

db = snakemake.params.get("db", "")
db = Path(db).resolve() if db else ""
db_option = "--db" if db else ""

proteins = snakemake.input.get("proteins", "")
proteins = Path(proteins).resolve() if proteins else ""
proteins_option = "--proteins" if proteins else ""

fasta = Path(snakemake.input.fasta).resolve()
outdir = Path(snakemake.output.results).resolve()
log_path = Path(snakemake.log[0]).resolve()
log_path.parent.mkdir(parents=True, exist_ok=True)
log_path.write_text("", encoding="utf-8")
prefix = validate_annotation_identifier(str(snakemake.params.prefix))
meta_option = "--meta" if snakemake.params.get("meta", False) else ""

if not db:
    raise ValueError("bakta.db must resolve to a persistent database directory.")
ensure_bakta_database(db, snakemake.params.db_type, log_path)
outdir.mkdir(parents=True, exist_ok=True)

original_cwd = Path.cwd()
try:
    os.chdir(outdir)
    shell(
        "bakta "
        "--prefix {prefix:q} "
        "--threads {snakemake.threads} "
        "{meta_option} "
        "{db_option} {db:q} "
        "{proteins_option} {proteins:q} "
        "{extra} "
        "{fasta:q} "
        ">> {log_path:q} 2>&1"
    )
finally:
    os.chdir(original_cwd)
