from pathlib import Path
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

db = snakemake.params.get("db", "")
db_option = "--db" if db else ""

if db:
    db_path = Path(db)
    if not (db_path / "bakta.db").exists():
        db_path.parent.mkdir(parents=True, exist_ok=True)
        shell("bakta_db download --output {db_path.parent:q}")

proteins = snakemake.input.get("proteins", "")
proteins_option = "--proteins" if proteins else ""

outdir = snakemake.output.outdir
prefix = Path(snakemake.input.fasta).stem

shell(
    "bakta "
    "--output {outdir:q} "
    "--prefix {prefix:q} "
    "--force "
    "--threads {snakemake.threads} "
    "{db_option} {db:q} "
    "{proteins_option} {proteins:q} "
    "{extra} "
    "{snakemake.input.fasta:q} "
    "{log}"
)
