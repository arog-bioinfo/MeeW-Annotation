"""CheckM2 1.1.0 adapter for one FASTA or a direct-child FASTA directory."""

from pathlib import Path

from snakemake.shell import shell

from common.empty_outputs import (
    CHECKM2_COLUMNS,
    has_fasta_inputs,
    write_delimited_header,
    write_empty_log,
)
from common.qa_inputs import staged_fasta_directory
from common.qa_databases import ensure_checkm2_database

report = Path(snakemake.output.report)
output_directory = report.parent
source = Path(snakemake.input.bins)
log_path = Path(snakemake.log[0])
extension = str(snakemake.params.extension).lstrip(".")

if not has_fasta_inputs(source):
    write_delimited_header(report, CHECKM2_COLUMNS)
    write_empty_log(log_path, tool="checkm2", domain="prok")
else:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("", encoding="utf-8")
    database_file, _ = ensure_checkm2_database(
        snakemake.params.database_path,
        snakemake.params.download_database,
        log_path,
    )

    with staged_fasta_directory(source, extension=f".{extension}") as staged:
        shell(
            "checkm2 predict"
            " --input {staged:q}"
            " --output-directory {output_directory:q}"
            " --database_path {database_file:q}"
            " --threads {snakemake.threads}"
            " --extension {extension:q}"
            " {snakemake.params.extra}"
            " >> {log_path:q} 2>&1"
        )
