# Workflow configuration

Each invocation processes exactly one isolate or one metagenomic project.

## General input

```yaml
input:
  mode: isolate # isolate | metagenome
  sample_sheet: config/isolate_samples.tsv
```

`input.mode` selects one of two strict one-record TSV contracts.
`input.sample_sheet` selects that TSV. Exactly one logical header and one
logical data record are required; blank or extra records, column/width errors,
blank required values, unsafe identifiers, and wrong path types fail early.

Isolate exact header and example:

```tsv
sample	domain	path
ecoli1	prok	data/prok/ecoli1.fasta
```

The domain is `prok` or `euk`. The path must be one existing regular
uncompressed nucleotide `.fa`, `.fna`, or `.fasta` file, case-insensitive.

Metagenome exact header and example:

```tsv
sample	prok_bins_dir	euk_bins_dir
study1	results/prok_bins	<blank>
```

Here `<blank>` denotes an empty TSV field; use an actual empty cell in the
metagenome input TSV. At least one domain directory is required. Blank domain
values, supplied empty directories, and absent domains are valid. Supplied
paths must be existing directories. Scanning is sorted, deterministic, and
direct-child only. It rejects duplicate filename-derived IDs across accepted
suffix spellings and domains.

All configured paths expand environment variables and `~`. Relative paths use
the Snakemake execution working directory, not the config or TSV directory.

## QA and routing

- `qa.enabled`: run QA tools and produce native QA reports (default `true`,
  including when omitted from an existing configuration).
- `qa.prok_tools`: the required single-item array `[checkm2]`.
- `qa.euk_tools`: the required single-item array `[eukcc]`.
- `qa.threads`: threads for each QA command.
- `qa.checkm2`: database path/download switch, staging extension, and extra
  options.
- `qa.eukcc`: database path/download switch and extra options.
- `qa_filter.enabled`: apply MAG completeness/contamination thresholds (default
  `true`). Isolates always remain `report_only`.

Exactly three QA/filter combinations are supported:

| `qa.enabled` | `qa_filter.enabled` | Behavior                                                                                     |
| ------------ | ------------------- | -------------------------------------------------------------------------------------------- |
| `true`       | `true`              | Run QA, retain native reports/metrics, and filter MAGs by thresholds.                        |
| `true`       | `false`             | Run QA and retain native reports/metrics, but do not reject MAGs for threshold failures.     |
| `false`      | `false`             | Run no QA tools and produce no native QA reports; write deterministic bypass selection rows. |

`qa.enabled: false` with `qa_filter.enabled: true` is rejected before DAG
construction. In QA-disabled mode the stable `qa/selection.tsv` columns remain:
metric and threshold cells are blank, reason is `qa_disabled`, isolates use
`report_only`, and every direct-child metagenome FASTA uses `pass`.

Prokaryotic isolates route to standard Bakta. Eukaryotic isolates route to
Funannotate2 after selecting one supported model route.
When QA/filtering is enabled, metagenome prokaryotic bins route through CheckM2
filtering to Bakta with exactly one managed `--meta`; eukaryotic bins route
through EukCC filtering to MetaEuk. The other supported modes retain the same
selection/gating interface described above. Absent domains run no scientific
command; empty supplied directories produce header-only contracts. Optional
GTDB-Tk eligibility is `report_only` for a prokaryotic isolate and exactly
`pass` for metagenome prokaryotic bins. Failed bins and every eukaryotic bin are
excluded from GTDB-Tk.

## Annotation settings

- `resources_root`: optional shared root for Annotation database/reference
  defaults. It does not enable any tool or download.
- `bakta.type`, `bakta.db`, `bakta.extra`: database selection, local database
  location, and extra options. `type` is `full` (the normal default) or `light`;
  `db` is an optional explicit local directory and is not a type selector. The
  native `annotation` prefix and `--meta` are workflow-managed. Output-directory
  (`--output`/`-o`) and force-overwrite (`--force`/`-f`) options are rejected in
  `extra`; Bakta runs directly in its declared result directory without either
  option. All argparse-accepted forms of `--skip-cds` and `--skip-plot` are also
  rejected because they suppress declared native outputs.
- `gtdbtk.enabled`, `data_dir`, `extra`: optional prokaryotic isolate or
  passing-MAG classification. A metagenome runs one project-level
  `classify_wf` batch containing all passing prokaryotic MAGs, never one job per
  MAG. Selected FASTAs are declared directly and symlinked by that job in
  unique temporary scratch, which is cleaned after success or failure; no
  persistent staging directory or staging log is produced. The resolved
  `resources.tmpdir` root must not contain ASCII spaces for GTDB-Tk 2.7.1;
  original FASTA source paths with spaces remain supported through validated
  scratch symlink names. If none passes, no GTDB-Tk job or classification
  artifact is requested.
- `recognizer_prok.resources_dir`, `extra`: prokaryotic reCOGnizer settings.
- `upimapi`: built-in/custom database, resources, extra options, and database
  check behavior. Its default `db: uniprot` is full UniProtKB (Swiss-Prot plus
  TrEMBL); `swissprot`, `taxids`, and custom references remain explicit choices.
- `metaeuk`: `database` selects a free-form MMseqs2 catalog name (`UniRef90` by
  default), `db` optionally selects an explicit local database prefix, and
  `cache_dir` selects the persistent managed catalog cache root. It also has
  managed `max_exon_sets` (default `1`) and constrained extras. `--threads`,
  `--max-exon-sets`, and `--exhaustive-search` cannot be supplied through
  `extra`; exhaustive search is not enabled automatically.
- `recognizer_euk`: resources and extra options. The eukaryotic route always
  selects reCOGnizer's built-in `KOG` database and consumes MetaEuk's exact
  stable protein FASTA.
- `funannotate2`: persistent database cache, upstream installation database
  list and `busco_lineage`, optional cleaning
  (`false` by default), species/strain, one reusable model option, and
  phase-specific extras. Supported deterministic routes are: params JSON
  (training skipped; species derived from JSON and checked against any
  configured species), pretrained model plus mandatory species (training
  skipped), or species without either reusable model (genome-specific
  training). Params and pretrained are mutually exclusive. Missing all model
  information is rejected before scientific execution; no generic fallback is
  selected. The shipped generic config deliberately leaves species/model blank
  so users must make that scientific choice for each eukaryotic isolate.
  `busco_lineage` is passed to train and annotate. Funannotate2 owns acquisition
  of that lineage and any lineage it selects internally for predict.
- `threads.high`, `medium`, `low`: tool thread presets.

`funannotate2.databases` keeps its existing nonempty, unique YAML-list
contract. The shipped list uses `- all` (equivalent to `[all]`) to install all
databases. A named subset is accepted as a list and becomes separate argv
elements, for example:

```yaml
databases:
  - uniprot
  - pfam
  - dbCAN
```

This produces `funannotate2 install -d uniprot pfam dbCAN`. Names are
case-sensitive, especially `dbCAN`; allowed values are `all`, `merops`,
`uniprot`, `dbCAN`, `pfam`, `go`, `mibig`, `interpro`, `gene2product`, and
`mito`. Upstream accepts `all` mixed with other names, and MeeW intentionally
does not reject mixed `all` lists or change upstream behavior. Funannotate2
26.6.21 checks whether `all` is present and then replaces the entire requested
list with all named databases, so additional entries are redundant and not
additive.

## Database/resource lookup

Database fields use an explicit tool path, nonblank top-level `resources_root`,
the relevant tool environment variable, `MEEW_RESOURCES`,
`$XDG_DATA_HOME/meew`, then `~/resources`. The first nonblank value wins; an
omitted, empty, or whitespace-only `resources_root` preserves the old chain.
All configured and environment paths expand environment variables and `~`, and
relative paths remain relative to Snakemake's execution working directory. The
shared root supplies only the table's default subpaths: it is never prepended to
an explicit tool path and does not enable a tool or any download switch. With no
configured or environment root, the subpaths remain beneath `~/resources`. This
also applies when `gtdbtk.enabled: true` and `data_dir` is blank.

| Tool          | Variable                | Resource-root subpath |
| ------------- | ----------------------- | --------------------- |
| CheckM2       | `CHECKM2DB`             | `checkm2_db`          |
| EukCC         | `EUKCC2_DB`, `EUKCC_DB` | `eukcc2_db_ver_1.2`   |
| Bakta `full`  | `BAKTA_DB`              | `bakta_db/db`         |
| Bakta `light` | `BAKTA_DB`              | `bakta_db/light/db`   |
| GTDB-Tk       | `GTDBTK_DATA_PATH`      | `gtdbtk_db`           |
| reCOGnizer    | `RECOGNIZER_RESOURCES`  | `recognizer_db`       |
| UPIMAPI       | `UPIMAPI_RESOURCES`     | `upimapi_db`          |
| MetaEuk cache | `METAEUK_DB`            | `metaeuk_db`          |
| Funannotate2  | `FUNANNOTATE2_DB`       | `funannotate2_db`     |

For blank locations, MetaEuk uses `<resolved-resource-root>/metaeuk_db` as its
cache.
Omitted, empty, or whitespace-only `metaeuk.database` uses `UniRef90`; a custom
value is trimmed and otherwise passed unchanged as one `mmseqs databases` argv
element. A blank `metaeuk.db` derives `catalogs/uniref90/db` within the resolved
cache for `UniRef90`. Every other source derives a dedicated
`catalogs/catalog-<encoded-source>/db` prefix, where `encoded-source` is Python
`urllib.parse.quote(trimmed_source, safe="")`; for example,
`UniProtKB/Swiss-Prot` becomes
`catalogs/catalog-UniProtKB%2FSwiss-Prot/db`. Distinct exact trimmed source
strings therefore have distinct cache leaves that cannot escape or collide
with the historical flat prefix layout. Funannotate2 uses
`<resolved-resource-root>/funannotate2_db`, and GTDB-Tk uses
`<resolved-resource-root>/gtdbtk_db`.

Database preparation occurs inside each selected scientific adapter. Native
database material is not a Snakemake output, and separate setup rules/targets no
longer exist. A ready matching cache is reused read-only without installer,
update, or network calls even when annotation reruns for provenance, parameter,
mode, or force reasons. Missing managed material is staged, validated, and
published under a cache-keyed interprocess lock. Existing partial, corrupt,
incompatible, or mismatched material is never automatically overwritten or
repaired. Cache paths that overlap Snakemake-owned analysis results are invalid.

Explicit Bakta, CheckM2, GTDB-Tk, and EukCC database transfers use at most four
attempts, waiting 10, 30, then 60 seconds after identifiable transient transfer
failures. Permanent, local-write, TLS-certificate, validation, checksum,
extraction, configuration, and unclassified failures are not retried, and the
final native exception is preserved. The policy wraps only the listed
transfer/native download operations: native checksum, extraction, installer,
and cache validation behavior remains authoritative. Downloads are not
resumable; in particular, Bakta 1.12.0 reopens its archive from byte zero for
each attempt.

`funannotate2.db_dir` and `metaeuk.cache_dir` select persistent caches.
`metaeuk.db` and `bakta.db` select explicit reused local locations and override
derived locations. When `bakta.db` is blank and `BAKTA_DB` is unset, `full`
uses `<resource-root>/bakta_db/db` while `light` uses
`<resource-root>/bakta_db/light/db`. Bakta requires schema-6 `version.json`, the
selected native `full`/`light` type, and Bakta 1.12.0's readable nonempty file
set. The downloader's `db` (full) and `db-light` (light) directories are handled
explicitly. GTDB-Tk requires `VERSION_DATA=r232` and a usable native reference
layout; fresh preparation uses the official full R232 package from
`https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz`
in isolated staging. Funannotate2 validates requested records in
`funannotate-db-info.json` and their referenced files; it installs only wholly
absent requested entries when they can be merged without replacing existing
material. Pretrained and BUSCO resources remain Funannotate2-owned. No adapter
adds `--force` or `--update`. The bundled environments pin
GTDB-Tk `2.7.1` for R232, Funannotate2 `26.6.21`, MetaEuk `7.bba0d80`, and
MMseqs2 `17.b804f`.

Managed MetaEuk catalog leaves validate nonempty data/index/dbtype files and a
matching identity record. Existing explicit and historical flat prefixes remain
in place. Without a record, a complete legacy prefix is a user assertion whose
catalog identity is reported as unverified; a known recorded mismatch fails.
MeeW never writes a false identity record for a legacy prefix and does not
prepare into a missing explicit/legacy prefix. No taxonomy filtering is added.

CheckM2 and EukCC download only as part of a selected QA job when
`download_database: true`. A direct CheckM2 `.dmnd` path is checked directly,
and EukCC validates native `db_base/db_info.json` references rather than merely
accepting an existing directory.

reCOGnizer receives `-rd <resolved-resource-root>/recognizer_db`; its adapter
has no explicit download command, so provision those resources using supported
external tool instructions and review tool-owned network behavior. UPIMAPI's
default `uniprot` selection uses full UniProtKB (Swiss-Prot plus TrEMBL) at the
native path `<resolved-resource-root>/upimapi_db/uniprot.fasta`, with the
derived index `uniprot.dmnd`. Full UniProtKB needs substantially more disk,
indexing time, memory, and search time than explicit `swissprot`, whose native
FASTA remains `uniprot_sprot.fasta`. UPIMAPI runs only the selected database;
the full-UniProtKB route does not add a separate Swiss-Prot pass. It passes
`--skip-db-check` only when the selected FASTA exists; otherwise its normal
database check runs during annotation and may download or update data. No
unverified standalone installer command is documented for either tool.

Bakta, reCOGnizer, and UPIMAPI produce partially overlapping but complementary
protein evidence. Changing the UPIMAPI database does not merge or replace
source columns: `protein_annotations.tsv` remains one row per exact Bakta
protein ID with separate Bakta, reCOGnizer, and UPIMAPI columns. UPIMAPI ranking
does not automatically prefer reviewed Swiss-Prot hits over better-ranked
TrEMBL hits.

Both adapters guard, rather than replace, their original native invocation.
They retain the configured database labels/custom path, skip decisions,
advanced options, command flags, working directory, logs, and outputs. There is
no MeeW downloader, staging cache, copied native CLI, database substitution, or
mapping/API algorithm change. Consequently the guard does not fix
reCOGnizer 1.12.2's old EggNOG URLs and does not make UPIMAPI 1.13.3's API
retry/poll loops bounded; native setup may still fail on upstream errors.

Verified complete selections use shared locks and can be read concurrently.
Cold or potentially mutating invocations use exclusive locks for the entire
native process-tree lifetime. reCOGnizer readiness includes its native
completion indicator, actual relational data, all members of each selected RPS
index family, and required installation-level COG maps. UPIMAPI readiness
includes both its selected FASTA and derived DIAMOND index, or a structurally
valid directly selected `.dmnd`; `--skip-db-check` alone does not make a missing
index read-only. Custom references coordinate on the canonical native index
target, while the original custom database string is still passed unchanged.
Taxonomic, local-ID, or otherwise unproven advanced modes are conservatively
exclusive because native code may add shared files.

Existing nonempty incomplete or corrupt selected material is diagnosed and
left untouched instead of being considered ready. Preserve it for diagnosis and
use the native tool's supported recovery process before restarting the job;
MeeW does not delete/replace the cache or create a fake readiness timestamp.
After recovery, newly visible database material is a declared optional input and
can cause Snakemake to rerun completed consumers for provenance. Do not suppress
that rerun by modifying result timestamps or metadata.

## Result contracts

All outputs use the concrete root `results/<sample>/`. Isolate QA selection is
one row and annotation follows the domain route. Metagenome selection is one
row per discovered direct-child FASTA. Only passing bins become independent
annotation jobs. The final internal `batch/annotation_results.tsv` remains a
valid header-only table when nothing passes.

Prokaryotic isolate runs retain Bakta's native table at
`bakta/annotation.tsv` and produce `protein_annotations.tsv`. Each passing
prokaryotic bin produces `batch/prok/<bin>/protein_annotations.tsv`; this path
is included in `batch/annotation_results.tsv` alongside the atomic Bakta,
reCOGnizer, and UPIMAPI directory paths. The summary has exactly one row
per exact Bakta FAA record in FAA order and source-prefixed Bakta, selected
reCOGnizer COG, and selected UPIMAPI columns. Full schema and deterministic
selection/hit-count semantics are documented in the main README. Empty Bakta
FAAs are passed to reCOGnizer and UPIMAPI normally; their adapters do not create
synthetic empty outputs or evidence-provenance mappings, and any upstream
failure propagates normally. Rejected bins produce no annotation outputs.

An enabled GTDB-Tk route publishes the consolidated bacterial/archaeal result
at `gtdbtk/classification.tsv` and retains native GTDB-Tk outputs in the same
directory. Consolidation requires the native summary contracts to agree and
requires exactly one classified row for every scratch-staged identifier, with no
missing, duplicate, or unexpected identifier. For a metagenome, the existing
four-column `batch/annotation_results.tsv` schema represents this project-level
output with `bin=<configured project sample>`, `domain=prok`,
`output_type=gtdbtk_classification`, and its absolute `path`. There is one such
row per project, not one per MAG. When no prokaryotic MAG passes, this row and
the classification artifact are absent; the inventory and optional downstream
FAA manifest remain deterministic and may be header-only.

Bakta 1.12.0, reCOGnizer 1.12.2, and UPIMAPI 1.13.3 each own one atomic
`directory()` result. Their native files can vary, so they are not individual
Snakemake outputs or inventory rows. Adapters write directly to those
directories and preserve native filenames, including UPIMAPI's
`unaligned.blast`. Workflow-managed input, output, resources, threads, database
selection, and output-shaping options cannot be supplied through `extra`.

Snakemake determines directory freshness from `.snakemake_timestamp` and
cleans the complete owned directory before a producer reruns. Users must not
manually modify individual internal files; invalidate the whole directory when
a rerun is required. Snakemake `log:` files are outside owned directories so
they survive cleanup and failed jobs.

Funannotate2 is deliberately one stateful job because its phases mutate one
project. It retains the full native project at `funannotate2/project/` and
publishes `annotation.gff3`, `proteins.faa`, `transcripts.fna`, `summary.json`,
and machine-readable `provenance.json`. When cleaning is disabled, the original
input is used directly by training/prediction; when enabled, all later phases
use the cleaned project FASTA.

Each passing eukaryotic MAG retains MetaEuk's direct native outputs under the
atomic `batch/euk/<bin>/metaeuk/` directory:
`native_predictions.fas`, `native_predictions.codon.fas`,
`native_predictions.headersMap.tsv`, `native_predictions.gff`, and tool-owned
temporary/internal files. Native names and identifiers are unchanged; no
synthetic reference-provenance file is published. The corresponding
reCOGnizer directory is also atomic and explicitly uses KOG.
