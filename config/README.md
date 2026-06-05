# Workflow configuration

The workflow processes one or more Metagenome-Assembled Genomes (MAGs) per run.
Set these fields in `config/config.yaml`:

* `sample_sheet`: path to a TSV file containing the sample names and paths.
* `prodigal.extra`: optional extra options string passed to the Prodigal wrapper (e.g., `-p meta -f gff`).
* `bakta.db`: path to the Bakta database directory.
* `bakta.extra`: optional extra options string passed to the Bakta wrapper.
* `gtdbtk.data_dir`: path to the GTDB-Tk reference database directory.
* `gtdbtk.extra`: optional extra options string passed to the GTDB-Tk wrapper.
* `metaeuk.db`: path to the MetaEuk reference database (UniProt database).
* `metaeuk.extra`: optional extra options string passed to the MetaEuk wrapper.
* `recognizer.resources_dir`: path to the reCOGnizer resources database directory.
* `recognizer.extra`: optional extra options string passed to the reCOGnizer wrapper.
* `upimapi.extra`: optional extra options string passed to the UPIMAPI wrapper.
* `threads`: dictionary containing computational resource presets (`high`, `medium`, `low`).

## Sample sheet format (TSV)

Required columns:
* `sample`: unique identifier/name for the MAG or isolate.
* `path`: path to the input genome file in FASTA format (`.fasta`, `.fna`, `.fa`).

The workflow will dynamically process all rows defined in this sheet.

## Example files

### `config/config.yaml`:
```yaml
sample_sheet: config/samples.tsv

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

upimapi:
  extra: "--evalue 1e-5"

threads:
  high: 16
  medium: 8
  low: 1