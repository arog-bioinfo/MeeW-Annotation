"""Write one deterministic annotation summary row per Bakta protein."""

from pathlib import Path

from common.protein_summary import write_protein_summary

bakta = Path(snakemake.input.bakta)
recognizer = Path(snakemake.input.recognizer)
upimapi = Path(snakemake.input.upimapi)
write_protein_summary(
    bakta / "annotation.faa",
    bakta / "annotation.tsv",
    recognizer / "reCOGnizer_results.tsv",
    upimapi / "UPIMAPI_results.tsv",
    snakemake.output.summary,
)

log = Path(snakemake.log[0])
log.parent.mkdir(parents=True, exist_ok=True)
log.write_text(f"Wrote {snakemake.output.summary}.\n", encoding="utf-8")
