"""Validation for workflow-managed UPIMAPI 1.13.3 options."""

from common.managed_options import reject_managed_options

_LONG_OPTIONS = {
    "--bitscore",
    "--blast",
    "--block-size",
    "--columns",
    "--database",
    "--diamond-mode",
    "--evalue",
    "--fasta",
    "--from-db",
    "--full-id",
    "--help",
    "--index-chunks",
    "--input",
    "--local-id-mapping",
    "--max-memory",
    "--max-target-seqs",
    "--max-tries",
    "--mirror",
    "--no-annotation",
    "--output",
    "--output-table",
    "--pident",
    "--resources-directory",
    "--show-available-fields",
    "--skip-db-check",
    "--skip-id-checking",
    "--skip-id-mapping",
    "--sleep",
    "--step",
    "--taxids",
    "--threads",
    "--to-db",
    "--version",
}
_MANAGED_LONG_OPTIONS = {
    "--database",
    "--fasta",
    "--from-db",
    "--input",
    "--max-target-seqs",
    "--no-annotation",
    "--output",
    "--output-table",
    "--resources-directory",
    "--skip-db-check",
    "--skip-id-checking",
    "--skip-id-mapping",
    "--threads",
    "--to-db",
}
_SHORT_OPTIONS = {
    "-i": "--input",
    "-o": "--output",
    "-ot": "--output-table",
    "-rd": "--resources-directory",
    "-db": "--database",
    "-t": "--threads",
    "-mts": "--max-target-seqs",
}


def validate_upimapi_extra(extra):
    """Reserve UPIMAPI core annotation inputs, outputs, and scientific shape."""
    return reject_managed_options(
        extra,
        config_key="config.upimapi.extra",
        tool="UPIMAPI",
        long_options=_LONG_OPTIONS,
        managed_long_options=_MANAGED_LONG_OPTIONS,
        short_options=_SHORT_OPTIONS,
        explanation=(
            "MeeW fixes the input FAA, owned output directory, output table, "
            "resources, database, threads, one target per query, core annotation "
            "and ID-mapping behavior, and manages --skip-db-check itself."
        ),
    )
