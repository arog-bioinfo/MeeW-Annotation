# Snakemake workflow: MAG Annotation

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

A best-practice Snakemake workflow for the **Annotation of Metagenome-Assembled Genomes (MAGs)**.

## Pipeline Overview

This workflow orchestrates several state-of-the-art bioinformatics tools to perform structural and functional annotation, as well as taxonomic classification. It takes assembled genomes (MAGs or isolates) in `*.fasta` or `*.fna` format as input. FASTA inputs may be produced by upstream `MeeW-Assembly` or `MeeW-Binning` workflows, or supplied from compatible external sources.

1. **Structural Annotation / Gene Calling:** Predict protein-coding genes using `Prodigal` (prokaryotes) and `MetaEuk` (eukaryotes).
2. **Comprehensive Genome Annotation:** Rapid and standard-compliant annotation with `Bakta`.
3. **Functional Domain Annotation:** Identify COG/Pfam/CDD domains using `reCOGnizer`.
4. **Protein Mapping:** Functional annotation via UniProt databases using `UPIMAPI`.
5. **Taxonomic Classification:** Optional prokaryotic-only taxonomic assignment using `GTDB-Tk`.

The workflow also writes an Annotation-owned handoff sheet for downstream Metabolic Modeling at `results/stage_sheets/annotation_to_metabolic_modeling.tsv`. This TSV contains one row per prokaryotic MAG/bin only, with columns `mag` and `path`, where `path` is the absolute path to the Bakta protein FASTA output (`results/bakta/{sample}/{sample}.faa`). Eukaryotic bins are not included in this Metabolic Modeling handoff.

## Configuration & Input Data

Detailed information about input data formats and workflow configuration parameters (such as database paths and tool-specific arguments) can be found in the [`config/README.md`](config/README.md).

By default, the workflow expects:

- A TSV sample sheet (`config/samples.tsv`) containing the paths to your `*.fasta` files.
- A YAML configuration file (`config/config.yaml`) defining tool parameters and database locations.

## Acknowledgements / Authorship

This annotation workflow was originally authored by [@rodolfobrandao8](https://github.com/rodolfobrandao8) as a first-year project component developed under the supervision of [@arog-bioinfo](https://github.com/arog-bioinfo).

This repository is now an independently maintained continuation, maintained since 2026-06-12 by [@arog-bioinfo](https://github.com/arog-bioinfo), as part of second-year master's thesis work.

### 1. Deployment options

To run the workflow from the command line, change to the working directory of the project:

```bash
cd path/to/<repo>
```
