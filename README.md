# MeeW-Annotation

MeeW-Annotation is a Snakemake pipeline for structural and functional annotation of metagenome-assembled genomes (MAGs) and isolates. Its current interfaces support standalone sample-sheet runs and a QA-filtered directory/batch mode imported by main MeeW.

> [!WARNING]
> **Research preview:** MeeW-Annotation is under active development and is not yet a deployment-ready release. Validate results and resource requirements for your use case.

## Pipeline

1. Prokaryotic samples: Prodigal gene calling, Bakta annotation, prokaryotic reCOGnizer, UPIMAPI, and optional GTDB-Tk classification.
2. Eukaryotic MAGs: MetaEuk protein prediction followed by eukaryotic reCOGnizer.
3. Eukaryotic isolates: optional Funannotate2 cleaning, training, prediction, and annotation.
4. In standalone sample-sheet mode, write a prokaryotic protein handoff for Metabolic Modeling.

## Inputs

Standalone mode uses the TSV selected by `sample_sheet`. Required columns are:

- `sample`: unique MAG or isolate identifier
- `path`: input `.fasta`, `.fna`, or `.fa`
- `domain`: `prok` or `euk`

Optional `genome_type` accepts `mag` or `isolate` and defaults to `mag`. Eukaryotic isolates use the Funannotate2 path when enabled; eukaryotic MAGs use MetaEuk and reCOGnizer.

Directory mode instead consumes prokaryotic and eukaryotic Binning directories plus CheckM2 and EukCC reports. Main MeeW supplies this mode. Configured paths expand `~` and environment variables, while relative paths remain relative to the Snakemake working directory. See [`config/README.md`](config/README.md) for both input contracts.

## Configuration and resources

`config/config.yaml` controls optional QA filtering, tool arguments, databases, GTDB-Tk, Funannotate2, and thread presets. Sample-sheet QA filtering can use CheckM2 reports for prokaryotes and EukCC reports for eukaryotes, with configurable completeness, contamination, and missing-sample policy.

Annotation databases are not bundled. Each resource uses its explicit config path, then the relevant tool variable, `MEEW_RESOURCES`, `$XDG_DATA_HOME/meew`, and `~/.local/share/meew`. This applies to Bakta, GTDB-Tk, reCOGnizer, UPIMAPI, MetaEuk, and Funannotate2. Install the resources required by the selected paths before running; configured Bakta/Funannotate2 setup behavior is the only pipeline path that can create its own missing resources. Full lookup names and parameters are in [`config/README.md`](config/README.md).

### Directory/batch mode

When `directory_mode.enabled: true`, the pipeline filters prokaryotic bins with CheckM2 and eukaryotic bins with EukCC using `qa_filter.min_completeness` and `qa_filter.max_contamination`. Passing bins are symlink-staged under `batch/qa/` rather than copied. Batch Bakta, Prodigal, prokaryotic/eukaryotic reCOGnizer, UPIMAPI, and MetaEuk outputs and completion markers are written under `batch/`. Empty domain directories produce empty manifests and completed batch markers.

## Running the pipeline

### Prerequisites

- Snakemake and Conda
- Input genomes and metadata for sample-sheet mode, or Binning directories and QA reports for directory mode
- Databases for each selected annotation tool

Run commands from the repository root.

### Dry-run

The bundled fixture exercises the prokaryotic sample-sheet DAG:

```bash
snakemake \
  --directory .test \
  --configfile .test/config/config.yaml \
  --cores 2 \
  --dry-run
```

### Normal run

After updating `config/config.yaml` and `config/samples.tsv`:

```bash
snakemake --cores 2 --sdm conda
```

### Bundled test run

The existing conda-backed fixture command is:

```bash
snakemake \
  --directory .test \
  --configfile .test/config/config.yaml \
  --cores 2 \
  --sdm conda
```

The fixture still depends on the annotation resources named by its configuration.

## Outputs and stage handoff

Standalone sample-sheet outputs are rooted at `results/{sample}`. Depending on domain and genome type, key outputs include:

- prokaryotic Prodigal, Bakta, reCOGnizer, UPIMAPI, and optional GTDB-Tk results
- eukaryotic MetaEuk/reCOGnizer results for MAGs
- Funannotate2 outputs for enabled eukaryotic isolates

Each prokaryotic sample target also writes:

```text
results/{sample}/stage_sheets/annotation_to_metabolic_modeling.tsv
```

The sheet has `mag` and `path` columns and points to absolute Bakta protein FASTA paths. It contains prokaryotic MAG/bin rows only; eukaryotic samples are excluded because the current Metabolic Modeling stage is prokaryotic-only.

Directory mode writes QA manifests and batch tool directories under `results/{sample}/batch/`. Main MeeW passes the batch Bakta directory to Metabolic Modeling directly rather than using the standalone stage sheet.

## Standalone and integrated use

Standalone Annotation can use MeeW-Binning's `binning_to_annotation.tsv` as its sample sheet, and its prokaryotic handoff can be set as MeeW-Metabolic_Modeling's `input.faa_samples`. Main MeeW imports Annotation from the sibling MeeW-Annotation repository, places outputs under `results/{sample}/annotation`, and enables directory mode with Binning's separated domain directories and QA reports. Passing proteins from the imported batch path then become Metabolic Modeling inputs.

## Current limitations

- Databases are external and must be present before their selected tools run, except for the explicitly configured setup behavior documented for Bakta/Funannotate2.
- GTDB-Tk is optional and prokaryotic-only.
- Funannotate2 targets eukaryotic isolates; eukaryotic MAGs use MetaEuk/reCOGnizer.
- The standalone Metabolic Modeling handoff excludes eukaryotic proteins.
- Bundled fixtures exercise a small prokaryotic input and do not establish full scientific validation across all annotation paths.

## Acknowledgements and authorship

This annotation pipeline was originally authored by [@rodolfobrandao8](https://github.com/rodolfobrandao8) as a first-year project component developed under the supervision of [@arog-bioinfo](https://github.com/arog-bioinfo).

This repository is now an independently maintained continuation, maintained since 2026-06-12 by [@arog-bioinfo](https://github.com/arog-bioinfo), as part of second-year master's thesis work.

## Citation, support, and license

No module-specific DOI is declared. When publishing work that uses this pipeline, cite the [MeeW-Annotation repository](https://github.com/arog-bioinfo/MeeW-Annotation), the tools used by the selected annotation paths, and Snakemake according to its [citation guidance](https://snakemake.readthedocs.io/en/stable/project_info/citations.html).

Use the [MeeW-Annotation issue tracker](https://github.com/arog-bioinfo/MeeW-Annotation/issues) for reproducible bug reports or documentation questions. This repository is licensed under the [MIT License](LICENSE).
