# Snakemake pipeline: MeeW-Annotation

[![Snakemake](https://img.shields.io/badge/Snakemake-%E2%89%A59.12.0-039475?logo=snakemake&logoColor=white)](https://snakemake.readthedocs.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/conda/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

<!-- TEMPORARY CATALOG OPT-OUT: later replace every standalone "pipeline" with the standard term formed by "work" + "flow", then remove this comment. -->

A Snakemake pipeline for quality assessment and structural and functional annotation of prokaryotic and eukaryotic genomes.

## Usage

Each invocation processes one isolate or one metagenome project. CheckM2 and EukCC provide optional quality assessment and MAG filtering. Prokaryotic inputs use Bakta, reCOGnizer, UPIMAPI, and optional GTDB-Tk. Eukaryotic isolates use Funannotate2, while eukaryotic MAGs use MetaEuk and reCOGnizer's KOG route.

Install Snakemake `>=9.12.0` and Conda, then select the input mode and sample sheet in [`config/config.yaml`](config/config.yaml). Provision the persistent databases and scientific choices required by that route before execution. In particular, eukaryotic isolate annotation needs a supported species/model route. The default full-UniProt UPIMAPI selection has substantial storage, indexing, memory, and runtime costs.

Do not assume that every annotation route is fully offline. reCOGnizer resource setup remains tool-owned, and UPIMAPI retains its native API behavior; when its selected FASTA is absent, its normal database check may download or update resources. A ready local database therefore does not replace review of the selected tools' network requirements.

Results below `results/{sample}/` include QA selection and reports, route-specific native outputs, prokaryotic protein summaries, optional GTDB-Tk classification, and metagenome annotation inventories. See the [configuration guide](config/README.md) for input schemas, database behavior, algorithms, and exact output contracts. Run the commands below from the repository root.

## Deployment options

Conda is the supported software deployment method. After configuring `config/config.yaml`, inspect the planned jobs and then execute them with:

```bash
snakemake --cores 2 --dry-run
snakemake \
  --cores 2 \
  --software-deployment-method conda
```

Adjust cores and resource overrides for the selected databases and execution environment. This pipeline does not supply a container deployment.

## Execution profiles

The repository supplies [profile guidance](profiles/README.md), but no configured pipeline profile. Add or select a user-maintained profile for executor and site-specific resource settings.

## License

The pipeline is distributed under the [MIT License](LICENSE).

## References

No module-specific DOI is declared. Cite the repository version used and the selected tools and databases, including CheckM2, EukCC, Bakta, reCOGnizer, UPIMAPI, GTDB-Tk, Funannotate2, MetaEuk, MMseqs2, and Snakemake:

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, **2021**. https://doi.org/10.12688/f1000research.29032.2.
