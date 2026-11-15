import os
from pathlib import Path

from snakemake.shell import shell

from common.gtdbtk import (
    consolidate_gtdbtk_summaries,
    temporary_gtdbtk_genomes,
)
from common.gtdbtk_database import ensure_gtdbtk_database

extra = snakemake.params.get("extra", "")
log_path = Path(snakemake.log[0]).resolve()
log_path.parent.mkdir(parents=True, exist_ok=True)
log_path.write_text("", encoding="utf-8")

db_path = snakemake.params.get("db_path", "")
if not db_path:
    raise ValueError("gtdbtk.data_dir must resolve to a persistent database directory.")
db_path = Path(db_path).resolve()

out_dir = Path(snakemake.output.classification).resolve().parent
out_dir.mkdir(parents=True, exist_ok=True)

extension = snakemake.params.get("extension", "fasta")

with temporary_gtdbtk_genomes(
    snakemake.params.pairs,
    snakemake.input.fastas,
    snakemake.resources.tmpdir,
) as genome_dir:
    ensure_gtdbtk_database(db_path, log_path)
    os.environ["GTDBTK_DATA_PATH"] = str(db_path)
    shell(
        "gtdbtk classify_wf"
        " --genome_dir {genome_dir:q}"
        " --out_dir {out_dir:q}"
        " --extension {extension:q}"
        " --cpus {snakemake.threads}"
        " --force"
        " {extra}"
        " >> {log_path:q} 2>&1"
    )

summaries = [
    ("bacteria", out_dir / "gtdbtk.bac120.summary.tsv"),
    ("archaea", out_dir / "gtdbtk.ar53.summary.tsv"),
]
consolidate_gtdbtk_summaries(
    summaries,
    snakemake.output.classification,
    snakemake.params.expected_ids,
)
