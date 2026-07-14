"""Write the standalone Annotation-to-Metabolic-Modeling handoff sheet."""

from pathlib import Path

output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)
with output_path.open("w", newline="", encoding="utf-8") as handle:
    handle.write("mag\tpath\n")
    for mag, protein_fasta in zip(snakemake.params.mags, snakemake.input.proteins):
        handle.write(f"{mag}\t{Path(protein_fasta).resolve()}\n")

Path(snakemake.log[0]).write_text("Stage sheet written.\n", encoding="utf-8")
