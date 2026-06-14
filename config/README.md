# Workflow configuration

The workflow processes one or more Metagenome-Assembled Genomes (MAGs) per run.
Set these fields in `config/config.yaml`:

- `sample_sheet`: path to a TSV file containing the sample names and paths.
- `qa_filter`: optional external QA filtering before annotation targets are expanded. When enabled, prokaryotic samples use CheckM2 TSV reports (`Name`, `Completeness`, `Contamination`) and eukaryotic samples use EukCC CSV reports (`bin`, `completeness`, `contamination`). Samples pass when completeness is at least `min_completeness` and contamination is at most `max_contamination`; `missing_sample` can be `error`, `keep`, or `drop`.
- `prodigal.extra`: optional extra options string passed to the Prodigal wrapper (e.g., `-p meta -f gff`).
- `bakta.db`: path to the Bakta database directory.
- `bakta.extra`: optional extra options string passed to the Bakta wrapper.
- `gtdbtk.data_dir`: path to the GTDB-Tk reference database directory.
- `gtdbtk.extra`: optional extra options string passed to the GTDB-Tk wrapper.
- `metaeuk.db`: path to the MetaEuk reference database (UniProt database).
- `metaeuk.extra`: optional extra options string passed to the MetaEuk wrapper.
- `recognizer.resources_dir`: path to the reCOGnizer resources database directory.
- `recognizer.extra`: optional extra options string passed to the reCOGnizer wrapper.
- `recognizer.euk_custom_db`: path to a KOG/custom database for eukaryotic reCOGnizer. Leave empty to disable a custom eukaryotic database.
- `recognizer.euk_extra`: optional extra options string passed to eukaryotic reCOGnizer.
- `upimapi.db`: UPIMAPI built-in database name to use (for example, `swissprot`). Leave empty when using `upimapi.db_custom`.
- `upimapi.db_custom`: path to a custom UPIMAPI database FASTA. Leave empty when using `upimapi.db`.
- `upimapi.resources_dir`: path to the UPIMAPI resources database directory.
- `upimapi.extra`: optional extra options string passed to the UPIMAPI wrapper.
- `upimapi.skip_db_check_if_exists`: when `true` (default), automatically add `--skip-db-check` only if the selected UPIMAPI database FASTA already exists in `upimapi.resources_dir` or `upimapi.db_custom` exists.
- `threads`: dictionary containing computational resource presets (`high`, `medium`, `low`).

## Sample sheet format (TSV)

Required columns:

- `sample`: unique identifier/name for the MAG or isolate.
- `path`: path to the input genome file in FASTA format (`.fasta`, `.fna`, `.fa`).

The workflow will dynamically process all rows defined in this sheet.

## Example files

### `config/config.yaml`:

```yaml
sample_sheet: config/samples.tsv

qa_filter:
  enabled: false
  min_completeness: 50.0
  max_contamination: 10.0
  checkm2_reports: []
  eukcc_reports: []
  missing_sample: "error"

prodigal:
  extra: "-p meta -f gff"

bakta:
  db: "resources/bakta_db"
  extra: ""

gtdbtk:
  data_dir: "resources/gtdbtk_db"
  extra: ""

metaeuk:
  db: "resources/metaeuk_db/uniprot_db"
  extra: "--e 0.0001"

recognizer:
  resources_dir: "resources/recognizer_db"
  extra: "--evalue 0.001"
  euk_custom_db: ""
  euk_extra: ""

upimapi:
  db: "swissprot"
  db_custom: ""
  resources_dir: "resources/upimapi_db"
  extra: "--evalue 1e-5"
  skip_db_check_if_exists: true

threads:
  high: 16
  medium: 8
  low: 1
```
