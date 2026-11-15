"""Write required native QA contracts for an absent metagenome domain."""

from pathlib import Path

from common.empty_outputs import (
    CHECKM2_COLUMNS,
    EUKCC_COLUMNS,
    write_delimited_header,
    write_empty_log,
)

action = str(snakemake.params.action)
domain = str(snakemake.params.domain)

if action == "checkm2":
    write_delimited_header(Path(snakemake.output.report), CHECKM2_COLUMNS)
elif action == "eukcc":
    write_delimited_header(Path(snakemake.output.report), EUKCC_COLUMNS)
    write_delimited_header(Path(snakemake.output.bad_quality), EUKCC_COLUMNS)
    write_delimited_header(Path(snakemake.output.missing_markers), ["bin"])
else:
    raise ValueError(f"Unsupported absent-domain QA action: {action}")

write_empty_log(Path(snakemake.log[0]), tool=action, domain=domain)
