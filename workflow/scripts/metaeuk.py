"""Run MetaEuk directly in its Snakemake-owned result directory."""

import shlex
import subprocess
from pathlib import Path

from common.euk_options import validate_metaeuk_extra
from common.metaeuk_database import ensure_metaeuk_database
from common.versions import version_contains

EXPECTED_TOOL_VERSION = "7.bba0d80"

extra = validate_metaeuk_extra(snakemake.params.get("extra", ""))
database = Path(snakemake.params.db).resolve()
fasta = Path(snakemake.input.fasta).resolve()
outdir = Path(snakemake.output.results).resolve()
log_path = Path(snakemake.log[0]).resolve()
log_path.parent.mkdir(parents=True, exist_ok=True)
log_path.write_text("", encoding="utf-8")

runtime_tool_version = subprocess.run(
    ["metaeuk", "version"], check=True, capture_output=True, text=True
).stdout.strip()
if not version_contains(runtime_tool_version, EXPECTED_TOOL_VERSION):
    raise ValueError(
        f"Expected MetaEuk {EXPECTED_TOOL_VERSION}, found {runtime_tool_version!r}."
    )
ensure_metaeuk_database(
    database,
    snakemake.params.catalog,
    managed=bool(snakemake.params.database_managed),
    threads=snakemake.threads,
    log_path=log_path,
)

outdir.mkdir(parents=True, exist_ok=True)
native_prefix = outdir / "native_predictions"
tmp_dir = outdir / "tmp"
tmp_dir.mkdir(parents=True, exist_ok=True)
argv = [
    "metaeuk",
    "easy-predict",
    str(fasta),
    str(database),
    str(native_prefix),
    str(tmp_dir),
    "--threads",
    str(snakemake.threads),
    "--max-exon-sets",
    str(snakemake.params.max_exon_sets),
    *shlex.split(extra),
]

with log_path.open("a", encoding="utf-8") as log_handle:
    log_handle.write(f"$ {shlex.join(argv)}\n")
    log_handle.flush()
    subprocess.run(argv, check=True, stdout=log_handle, stderr=subprocess.STDOUT)
