# Snakemake workflow: MAG Annotation

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

A best-practice Snakemake workflow for the **Annotation of Metagenome-Assembled Genomes (MAGs)**.

##  Pipeline Overview

This workflow orchestrates several state-of-the-art bioinformatics tools to perform structural and functional annotation, as well as taxonomic classification. It takes assembled genomes (MAGs or isolates) in `*.fasta` or `*.fna` format as input.

1. **Structural Annotation / Gene Calling:** Predict protein-coding genes using `Prodigal` (prokaryotes) and `MetaEuk` (eukaryotes).
2. **Comprehensive Genome Annotation:** Rapid and standard-compliant annotation with `Bakta`.
3. **Functional Domain Annotation:** Identify COG/Pfam/CDD domains using `reCOGnizer`.
4. **Protein Mapping:** Functional annotation via UniProt databases using `UPIMAPI`.
5. **Taxonomic Classification:** Robust taxonomic assignment for all MAGs using `GTDB-Tk`.

##  Configuration & Input Data

Detailed information about input data formats and workflow configuration parameters (such as database paths and tool-specific arguments) can be found in the [`config/README.md`](config/README.md).

By default, the workflow expects:
* A TSV sample sheet (`config/samples.tsv`) containing the paths to your `*.fasta` files.
* A YAML configuration file (`config/config.yaml`) defining tool parameters and database locations.

##  Usage

### 1. Deployment options

To run the workflow from the command line, change to the working directory of the project:
```bash
cd path/to/<repo>
