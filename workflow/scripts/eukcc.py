"""EukCC 2.1.3 adapter for one FASTA or a direct-child FASTA directory."""

import shutil
import tempfile
from pathlib import Path

from snakemake.shell import shell

from common.empty_outputs import (
    EUKCC_COLUMNS,
    has_fasta_inputs,
    write_delimited_header,
    write_empty_log,
)
from common.qa_inputs import staged_fasta_directory
from common.qa_databases import ensure_eukcc_database

report = Path(snakemake.output.report)
outdir = report.parent
bad_quality = Path(snakemake.output.bad_quality)
missing_markers = Path(snakemake.output.missing_markers)
if outdir.exists():
    shutil.rmtree(outdir)
outdir.mkdir(parents=True, exist_ok=True)
source = Path(snakemake.input.bins)

if not has_fasta_inputs(source):
    write_delimited_header(report, EUKCC_COLUMNS)
    write_delimited_header(bad_quality, EUKCC_COLUMNS)
    write_delimited_header(missing_markers, ["bin"])
    write_empty_log(Path(snakemake.log[0]), tool="eukcc", domain="euk")
else:
    log_path = Path(snakemake.log[0]).resolve()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("", encoding="utf-8")
    resolved_database_path = Path(snakemake.params.database_path).resolve()
    ensure_eukcc_database(
        resolved_database_path, snakemake.params.download_database, log_path
    )
    with staged_fasta_directory(
        source, extension=".fa"
    ) as staged, tempfile.TemporaryDirectory(
        prefix="meew-annotation-eukcc-"
    ) as temporary:
        workspace = Path(temporary) / "output"
        shell(
            "eukcc folder"
            " --out {workspace:q}"
            " --db {resolved_database_path:q}"
            " --threads {snakemake.threads}"
            " {snakemake.params.extra}"
            " --suffix .fa"
            " {staged:q}"
            " >> {log_path:q} 2>&1"
        )
        native_outputs = {
            "eukcc.csv": report,
            "bad_quality.csv": bad_quality,
            "missing_marker_genes.txt": missing_markers,
        }
        for native_name, destination in native_outputs.items():
            native_path = workspace / native_name
            if not native_path.is_file():
                raise FileNotFoundError(
                    f"EukCC completed without expected native output: {native_path}"
                )
            shutil.copyfile(native_path, destination)
