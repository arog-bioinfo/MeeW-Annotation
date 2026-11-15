"""Validation for workflow-managed reCOGnizer 1.12.2 options."""

from common.managed_options import reject_managed_options

_LONG_OPTIONS = {
    "--custom-databases",
    "--databases",
    "--debug",
    "--evalue",
    "--file",
    "--help",
    "--keep-intermediates",
    "--keep-spaces",
    "--max-target-seqs",
    "--no-blast-info",
    "--no-output-sequences",
    "--output",
    "--output-rpsbproc-cols",
    "--protein-id-col",
    "--quiet",
    "--resources-directory",
    "--species-taxids",
    "--tax-col",
    "--tax-file",
    "--test-run",
    "--threads",
    "--version",
}
_MANAGED_LONG_OPTIONS = {
    "--custom-databases",
    "--databases",
    "--file",
    "--keep-intermediates",
    "--output",
    "--resources-directory",
    "--threads",
}
_SHORT_OPTIONS = {
    "-f": "--file",
    "-t": "--threads",
    "-o": "--output",
    "-rd": "--resources-directory",
    "-dbs": "--databases",
}


def validate_recognizer_extra(extra, config_key="config.recognizer_prok.extra"):
    """Reserve reCOGnizer I/O, resources, databases, and output variability."""
    return reject_managed_options(
        extra,
        config_key=config_key,
        tool="reCOGnizer",
        long_options=_LONG_OPTIONS,
        managed_long_options=_MANAGED_LONG_OPTIONS,
        short_options=_SHORT_OPTIONS,
        explanation=(
            "MeeW fixes the input FAA, owned output directory, threads, resource "
            "directory, route-specific built-in databases, and disables "
            "--keep-intermediates."
        ),
    )
