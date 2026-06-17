# Workflow configuration

The workflow processes one or more Metagenome-Assembled Genomes (MAGs) per run. Configure inputs and tool options in `config/config.yaml`.

## General input

- `sample_sheet`: path to a TSV file containing sample names, input FASTA paths, domains, and optional genome types.

## Sample sheet format

The sample sheet is a tab-separated file. Required columns:

- `sample`: unique identifier/name for the MAG or isolate.
- `path`: path to the input genome file in FASTA format (`.fasta`, `.fna`, `.fa`).
- `domain`: annotation path for the sample. Use `prok` for prokaryotic samples and `euk` for eukaryotic samples.

Optional columns:

- `genome_type`: use `mag` for metagenome-assembled genomes or `isolate` for isolate genomes. Missing values default to `mag`. Funannotate2 is only targeted for eukaryotic isolate genomes, so eukaryotic isolates require `domain: euk` and `genome_type: isolate`. Eukaryotic MAGs use the existing MetaEuk/reCOGnizer path.

The workflow dynamically processes all rows defined in this sheet.

Example:

```tsv
sample	path	domain	genome_type
sample_a	data/sample_a.fna	prok	mag
sample_b	data/sample_b.fasta	euk	mag
fungal_isolate	data/fungal_isolate.fasta	euk	isolate
```

## Optional QA filtering

- `qa_filter.enabled`: when `true`, filter samples before annotation targets are expanded.
- `qa_filter.min_completeness`: minimum completeness required for a sample to pass.
- `qa_filter.max_contamination`: maximum contamination allowed for a sample to pass.
- `qa_filter.checkm2_reports`: CheckM2 TSV reports for prokaryotic samples. Reports must include `Name`, `Completeness`, and `Contamination`.
- `qa_filter.eukcc_reports`: EukCC CSV reports for eukaryotic samples. Reports must include `bin`, `completeness`, and `contamination`.
- `qa_filter.missing_sample`: behavior when a sample is absent from QA reports. Supported values are `error`, `keep`, and `drop`.

## Prokaryotic annotation path

- `prodigal.extra`: optional extra options string passed to the Prodigal wrapper.
- `bakta.db`: path to the Bakta database directory.
- `bakta.extra`: optional extra options string passed to the Bakta wrapper.
- `gtdbtk.data_dir`: path to the GTDB-Tk reference database directory.
- `gtdbtk.extra`: optional extra options string passed to the GTDB-Tk wrapper.
- `recognizer_prok.resources_dir`: path to the prokaryotic reCOGnizer resources database directory.
- `recognizer_prok.extra`: optional extra options string passed to the prokaryotic reCOGnizer wrapper.
- `upimapi.db`: UPIMAPI built-in database name to use, for example `swissprot`. Leave empty when using `upimapi.db_custom`.
- `upimapi.db_custom`: path to a custom UPIMAPI database FASTA. Leave empty when using `upimapi.db`.
- `upimapi.resources_dir`: path to the UPIMAPI resources database directory.
- `upimapi.extra`: optional extra options string passed to the UPIMAPI wrapper.
- `upimapi.skip_db_check_if_exists`: when `true`, automatically add `--skip-db-check` only if the selected UPIMAPI database FASTA already exists in `upimapi.resources_dir` or `upimapi.db_custom` exists.

## Eukaryotic annotation path

- `metaeuk.db`: path to the MetaEuk reference database, such as a UniProt database.
- `metaeuk.extra`: optional extra options string passed to the MetaEuk wrapper.
- `recognizer_euk.resources_dir`: path to the eukaryotic reCOGnizer resources database directory.
- `recognizer_euk.custom_db`: path to a KOG/custom database for eukaryotic reCOGnizer. Leave empty to disable a custom eukaryotic database.
- `recognizer_euk.extra`: optional extra options string passed to the eukaryotic reCOGnizer wrapper.
- `funannotate2.enabled`: when `true`, target Funannotate2 for eukaryotic isolate genomes (`domain: euk`, `genome_type: isolate`).
- `funannotate2.db_dir`: Funannotate2 database directory. The default is `/home/argomes/resources/funannotate2_db`.
- `funannotate2.install_db`: when `true`, a real workflow run installs missing Funannotate2 databases into `funannotate2.db_dir`; dry-runs only plan this step.
- `funannotate2.databases`: Funannotate2 database names to install, for example `all`.
- `funannotate2.extra_install`: optional extra CLI options passed to `funannotate2 install`.
- `funannotate2.species`, `funannotate2.strain`, `funannotate2.params`, `funannotate2.pretrained`: optional values passed to Funannotate2 steps when provided.
- `funannotate2.extra_clean`, `funannotate2.extra_train`, `funannotate2.extra_predict`, `funannotate2.extra_annotate`: optional extra CLI options for each Funannotate2 step.

## Thread presets

- `threads`: dictionary containing computational resource presets.
- `threads.high`: thread count for high-resource steps.
- `threads.medium`: thread count for medium-resource steps.
- `threads.low`: thread count for low-resource steps.

## Example config

```yaml
# ====================
# General Input
# ====================
sample_sheet: "config/samples.tsv"

# ====================
# Quality Filtering
# ====================
qa_filter:
  enabled: false
  min_completeness: 50.0
  max_contamination: 10.0
  checkm2_reports: []
  eukcc_reports: []
  missing_sample: "error"

# ====================
# Prokaryotic Annotation
# ====================

# --------------------
# Prodigal
# --------------------
prodigal:
  extra: "-p meta -f gff"

# --------------------
# Bakta
# --------------------
bakta:
  db: "resources/bakta_db/db-light"
  extra: ""

# --------------------
# GTDB-Tk
# --------------------
gtdbtk:
  data_dir: "resources/gtdbtk_db"
  extra: ""

# --------------------
# reCOGnizer Prokaryotic
# --------------------
recognizer_prok:
  resources_dir: "resources/recognizer_db"
  extra: ""

# --------------------
# UPIMAPI
# --------------------
upimapi:
  db: "swissprot"
  db_custom: ""
  resources_dir: "resources/upimapi_db"
  extra: ""
  skip_db_check_if_exists: true

# ====================
# Eukaryotic Annotation
# ====================

# --------------------
# MetaEuk
# --------------------
metaeuk:
  db: "resources/metaeuk_db/uniprot_db"
  extra: ""

# --------------------
# reCOGnizer Eukaryotic
# --------------------
recognizer_euk:
  resources_dir: "resources/recognizer_db"
  custom_db: ""
  extra: ""

# --------------------
# Funannotate2
# --------------------
funannotate2:
  enabled: true
  db_dir: "/home/argomes/resources/funannotate2_db"
  install_db: true
  databases:
    - all
  extra_install: ""
  species: ""
  strain: ""
  params: ""
  pretrained: ""
  extra_clean: ""
  extra_train: ""
  extra_predict: ""
  extra_annotate: ""

# ====================
# Computational Resources
# ====================
threads:
  high: 16
  medium: 8
  low: 1
```
