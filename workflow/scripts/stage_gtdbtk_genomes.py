"""Symlink prokaryotic genomes into the GTDB-Tk staging directory."""

import shutil
from pathlib import Path

genomes = Path(snakemake.output.genomes)
if genomes.exists():
    shutil.rmtree(genomes)
genomes.mkdir(parents=True)

for sample, source in zip(snakemake.params.samples, snakemake.input):
    target = genomes / f"{sample}.fasta"
    target.symlink_to(Path(source).resolve())

Path(snakemake.log[0]).write_text("GTDB-Tk genomes staged.\n", encoding="utf-8")
