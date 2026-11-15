# import basic packages
import sys
from pathlib import Path

from snakemake.exceptions import WorkflowError
from snakemake.utils import validate

_common_source_paths = [
    Path(workflow.source_path(f"../scripts/common/{module}.py"))
    for module in (
        "__init__",
        "bakta_options",
        "bakta_database",
        "cache_runtime",
        "database_dependencies",
        "database_selection",
        "download_retries",
        "empty_outputs",
        "euk_options",
        "funannotate2_database",
        "gtdbtk",
        "gtdbtk_database",
        "identifiers",
        "input_config",
        "managed_options",
        "metaeuk_database",
        "native_execution",
        "paths",
        "protein_summary",
        "qa",
        "qa_databases",
        "qa_inputs",
        "recognizer_options",
        "recognizer_database",
        "upimapi_database",
        "upimapi_options",
        "versions",
    )
]
sys.path.insert(0, str(_common_source_paths[0].parent.parent))
from common.bakta_options import validate_bakta_extra
from common.database_dependencies import (
    bakta_database_inputs,
    funannotate2_database_inputs,
    gtdbtk_database_inputs,
    metaeuk_database_inputs,
    recognizer_database_inputs,
    upimapi_database_inputs,
)
from common.database_selection import (
    resolve_bakta_database,
    resolve_metaeuk_database_location,
)
from common.euk_options import (
    validate_funannotate2_extras,
    validate_funannotate2_route,
    validate_metaeuk_extra,
)
from common.identifiers import validate_annotation_identifier
from common.input_config import load_input_record
from common.paths import expand_path, resolve_resource_path
from common.qa import read_selection
from common.recognizer_options import validate_recognizer_extra
from common.upimapi_options import validate_upimapi_extra

ANNOTATION_BIN_WILDCARD_REGEX = r"[A-Za-z0-9][A-Za-z0-9._-]*"
ANNOTATION_RESOURCES_ROOT = config.get("resources_root", "")
GTDBTK_DATA_DIR = resolve_resource_path(
    config.get("gtdbtk", {}).get("data_dir", ""),
    "GTDBTK_DATA_PATH",
    "gtdbtk_db",
    ANNOTATION_RESOURCES_ROOT,
)
FUNANNOTATE2_DB_DIR = resolve_resource_path(
    config.get("funannotate2", {}).get("db_dir", ""),
    "FUNANNOTATE2_DB",
    "funannotate2_db",
    ANNOTATION_RESOURCES_ROOT,
)


QA_ENABLED = config.get("qa", {}).get("enabled", True)
QA_FILTER_ENABLED = config.get("qa_filter", {}).get("enabled", True)
if QA_ENABLED is False and QA_FILTER_ENABLED is True:
    raise WorkflowError(
        "Invalid QA configuration: config.qa.enabled is false while "
        "config.qa_filter.enabled is true. Disable qa_filter when QA is disabled. "
        "Supported modes are QA on/filter on, QA on/filter off, and "
        "QA off/filter off."
    )


validate(config, schema="../schemas/config.schema.yaml")
BAKTA_DB_TYPE, BAKTA_DB_DIR = resolve_bakta_database(
    config.get("bakta", {}).get("type"),
    config.get("bakta", {}).get("db", ""),
    ANNOTATION_RESOURCES_ROOT,
)
METAEUK_CACHE_DIR = resolve_resource_path(
    config.get("metaeuk", {}).get("cache_dir", ""),
    "METAEUK_DB",
    "metaeuk_db",
    ANNOTATION_RESOURCES_ROOT,
)
_metaeuk_location = resolve_metaeuk_database_location(
    config.get("metaeuk", {}).get("database"),
    config.get("metaeuk", {}).get("db", ""),
    METAEUK_CACHE_DIR,
)
METAEUK_DATABASE = _metaeuk_location.catalog
METAEUK_DB_PREFIX = _metaeuk_location.prefix
METAEUK_DB_MANAGED = _metaeuk_location.managed
try:
    validate_bakta_extra(config.get("bakta", {}).get("extra", ""))
    validate_metaeuk_extra(config.get("metaeuk", {}).get("extra", ""))
    validate_recognizer_extra(config.get("recognizer_prok", {}).get("extra", ""))
    validate_recognizer_extra(
        config.get("recognizer_euk", {}).get("extra", ""),
        config_key="config.recognizer_euk.extra",
    )
    validate_upimapi_extra(config.get("upimapi", {}).get("extra", ""))
except ValueError as error:
    raise WorkflowError(str(error)) from error

IMPORT_CONTEXT_ENABLED = bool(config.get("import_context", {}).get("enabled", False))


def _module_input_record():
    """Build the private input context supplied by an importing workflow."""
    context = config.get("import_context", {})
    try:
        sample = context["sample"]
        input_mode = context["mode"]
        if not isinstance(sample, str) or not isinstance(input_mode, str):
            raise TypeError
        sample = validate_annotation_identifier(sample)
    except (KeyError, TypeError, ValueError) as error:
        raise WorkflowError(
            "Invalid imported Annotation context: mode and sample are required, and "
            "sample must satisfy the Annotation identifier contract."
        ) from error

    if input_mode == "metagenome":
        try:
            prok_bins_dir = context["prok_bins_dir"]
            euk_bins_dir = context["euk_bins_dir"]
            if not isinstance(prok_bins_dir, str) or not isinstance(euk_bins_dir, str):
                raise TypeError
        except (KeyError, TypeError, ValueError) as error:
            raise WorkflowError(
                "Imported metagenome Annotation context requires future "
                "prok_bins_dir and euk_bins_dir outputs."
            ) from error
        if not prok_bins_dir.strip() or not euk_bins_dir.strip():
            raise WorkflowError(
                "Imported metagenome Annotation context requires both future "
                "prok_bins_dir and euk_bins_dir outputs."
            )
        return input_mode, {
            "sample": sample,
            "prok_bins_dir": prok_bins_dir,
            "euk_bins_dir": euk_bins_dir,
        }

    if input_mode != "isolate":
        raise WorkflowError(
            "Invalid imported Annotation context mode: expected isolate or metagenome; "
            f"found {input_mode!r}."
        )
    try:
        domain = context["domain"]
        path = context["path"]
        if not isinstance(domain, str) or not isinstance(path, str):
            raise TypeError
    except (KeyError, TypeError, ValueError) as error:
        raise WorkflowError(
            "Imported isolate Annotation context requires a prok/euk domain and a "
            "future Assembly FASTA path."
        ) from error
    if domain not in {"prok", "euk"} or not path.strip():
        raise WorkflowError(
            "Imported isolate Annotation context requires domain=prok or domain=euk "
            "and a non-empty future Assembly FASTA path."
        )
    return input_mode, {"sample": sample, "domain": domain, "path": path}


def _standalone_input_record():
    """Load and validate standalone metadata only for a standalone invocation."""
    input_mode = config["input"]["mode"]
    config["input"]["sample_sheet"] = expand_path(config["input"]["sample_sheet"])
    try:
        record = load_input_record(input_mode, config["input"]["sample_sheet"])
    except ValueError as error:
        raise WorkflowError(str(error)) from error

    sample_schema_source = workflow.current_basedir.join(
        f"../schemas/{input_mode}_samples.schema.yaml"
    )
    try:
        sample_schema = sample_schema_source.get_path_or_uri(secret_free=False)
    except TypeError:
        # Snakemake 9.12 does not yet accept the secret_free keyword.
        sample_schema = sample_schema_source.get_path_or_uri()
    validate(record, schema=sample_schema)
    return input_mode, record


if IMPORT_CONTEXT_ENABLED:
    INPUT_MODE, INPUT_RECORD = _module_input_record()
else:
    INPUT_MODE, INPUT_RECORD = _standalone_input_record()

SAMPLE = INPUT_RECORD["sample"]
RESULTS_ROOT = f"results/{SAMPLE}"
LOGS_ROOT = "<results>/logs"
IS_ISOLATE = INPUT_MODE == "isolate"
IS_METAGENOME = INPUT_MODE == "metagenome"
ISOLATE_DOMAIN = INPUT_RECORD.get("domain", "")
HAS_PROK_INPUT = (
    IS_ISOLATE
    and ISOLATE_DOMAIN == "prok"
    or bool(INPUT_RECORD.get("prok_bins_dir", ""))
)
HAS_EUK_INPUT = (
    IS_ISOLATE and ISOLATE_DOMAIN == "euk" or bool(INPUT_RECORD.get("euk_bins_dir", ""))
)
GTDBTK_ENABLED = config.get("gtdbtk", {}).get("enabled", False)
# Resource fallbacks intentionally follow validation so configured path types
# are checked before expansion while blank values can still select fallbacks.
qa_config = config["qa"]
selected_prok_qa = qa_config["prok_tools"]
selected_euk_qa = qa_config["euk_tools"]
if QA_ENABLED:
    if "checkm2" not in selected_prok_qa:
        raise WorkflowError("config.qa.prok_tools must select checkm2.")
    if "eukcc" not in selected_euk_qa:
        raise WorkflowError("config.qa.euk_tools must select eukcc.")

    checkm2_config = qa_config["checkm2"]
    checkm2_config["database_path"] = resolve_resource_path(
        checkm2_config.get("database_path", ""),
        "CHECKM2DB",
        "checkm2_db",
        ANNOTATION_RESOURCES_ROOT,
    )
    eukcc_config = qa_config["eukcc"]
    eukcc_config["database_path"] = resolve_resource_path(
        eukcc_config.get("database_path", ""),
        ("EUKCC2_DB", "EUKCC_DB"),
        "eukcc2_db_ver_1.2",
        ANNOTATION_RESOURCES_ROOT,
    )
resource_fields = [
    ("recognizer_prok", "resources_dir", "RECOGNIZER_RESOURCES", "recognizer_db"),
    ("upimapi", "resources_dir", "UPIMAPI_RESOURCES", "upimapi_db"),
    ("recognizer_euk", "resources_dir", "RECOGNIZER_RESOURCES", "recognizer_db"),
]
for section, field, tool_environment, subpath in resource_fields:
    section_config = config.setdefault(section, {})
    section_config[field] = resolve_resource_path(
        section_config.get(field, ""),
        tool_environment,
        subpath,
        ANNOTATION_RESOURCES_ROOT,
    )

for section, field in [
    ("upimapi", "db_custom"),
]:
    if field in config.get(section, {}) and config[section][field]:
        config[section][field] = expand_path(config[section][field])

if IS_ISOLATE and ISOLATE_DOMAIN == "euk":
    try:
        validate_funannotate2_route(config.get("funannotate2", {}), inspect_params=True)
        validate_funannotate2_extras(config.get("funannotate2", {}))
    except ValueError as error:
        raise WorkflowError(str(error)) from error

config["bakta"]["db"] = BAKTA_DB_DIR
config["bakta"]["type"] = BAKTA_DB_TYPE
config["gtdbtk"]["data_dir"] = GTDBTK_DATA_DIR
config["funannotate2"]["db_dir"] = FUNANNOTATE2_DB_DIR
config["metaeuk"]["cache_dir"] = METAEUK_CACHE_DIR
config["metaeuk"]["db"] = METAEUK_DB_PREFIX
config["metaeuk"]["database"] = METAEUK_DATABASE


def sample_fasta(wildcards):
    return INPUT_RECORD["path"]


def checked_cache_path(cache_path, results_root, resource_name):
    """Reject cache/output overlap while the DAG is built, before output cleanup."""
    cache = Path(cache_path).resolve()
    results = Path(results_root).resolve()
    try:
        cache.relative_to(results)
        overlaps = True
    except ValueError:
        try:
            results.relative_to(cache)
            overlaps = True
        except ValueError:
            overlaps = False
    if overlaps:
        raise WorkflowError(
            f"{resource_name} persistent cache {cache} overlaps Snakemake-owned "
            f"analysis results {results}. Move the cache outside the analysis output "
            "tree before running; Snakemake may remove directory outputs before the "
            "adapter can inspect them."
        )
    return str(cache)


def funannotate2_model_input(wildcards):
    """Declare reusable params/model files as material DAG inputs when file-backed."""
    route_config = config.get("funannotate2", {})
    params_json = route_config.get("params") or ""
    pretrained = route_config.get("pretrained") or ""
    if params_json:
        return expand_path(params_json)
    if pretrained:
        expanded = expand_path(pretrained)
        if Path(expanded).is_file():
            return expanded
    return []


def qa_domain_source(domain, wildcards):
    if IS_ISOLATE:
        if ISOLATE_DOMAIN != domain:
            raise WorkflowError(
                f"Isolate '{SAMPLE}' has domain '{ISOLATE_DOMAIN}', not '{domain}'."
            )
        return INPUT_RECORD["path"]
    source = INPUT_RECORD[f"{domain}_bins_dir"]
    if not source:
        raise WorkflowError(f"Metagenome '{SAMPLE}' has no {domain} input domain.")
    return source


def prok_qa_source(wildcards):
    return qa_domain_source("prok", wildcards)


def euk_qa_source(wildcards):
    return qa_domain_source("euk", wildcards)


def metagenome_selection_records(wildcards):
    """Read checkpoint-selected bins without retaining runtime state in memory."""
    checkpoint_wildcards = {"sample": SAMPLE} if IMPORT_CONTEXT_ENABLED else {}
    selection = checkpoints.qa_selection.get(**checkpoint_wildcards).output.selection
    rows = read_selection(selection)
    seen = set()
    for row in rows:
        key = (row["sample"], row["domain"])
        if key in seen:
            raise WorkflowError(
                f"Duplicate metagenome QA selection row for {row['domain']} bin "
                f"'{row['sample']}' in {selection}."
            )
        seen.add(key)
    return selection, rows


def metagenome_selection_file(wildcards):
    selection, _ = metagenome_selection_records(wildcards)
    return selection


def metagenome_bin_fasta(domain):
    def selected_fasta(wildcards):
        selection, rows = metagenome_selection_records(wildcards)
        matches = [
            row
            for row in rows
            if row["sample"] == wildcards.bin
            and row["domain"] == domain
            and row["decision"] == "pass"
        ]
        if len(matches) != 1:
            raise WorkflowError(
                f"Expected exactly one passing {domain} QA selection row for bin "
                f"'{wildcards.bin}' in {selection}; found {len(matches)}."
            )
        return matches[0]["path"]

    return selected_fasta


metagenome_prok_bin_fasta = metagenome_bin_fasta("prok")
metagenome_euk_bin_fasta = metagenome_bin_fasta("euk")


def passing_metagenome_bin_ids(wildcards, domain):
    _, rows = metagenome_selection_records(wildcards)
    return sorted(
        row["sample"]
        for row in rows
        if row["domain"] == domain and row["decision"] == "pass"
    )


def format_metagenome_output(path, *, bin_id=None):
    output_wildcards = {}
    if IMPORT_CONTEXT_ENABLED:
        output_wildcards["sample"] = SAMPLE
    if bin_id is not None:
        output_wildcards["bin"] = bin_id
    return str(path).format(**output_wildcards)


def metagenome_annotation_materials(wildcards):
    """Resolve all concrete per-bin outputs after the QA checkpoint."""
    materials = []
    for bin_id in passing_metagenome_bin_ids(wildcards, "prok"):
        for output in (
            rules.bakta_bin.output,
            rules.recognizer_prok_bin.output,
            rules.upimapi_bin.output,
            rules.protein_summary_bin.output,
        ):
            materials.extend(
                format_metagenome_output(path, bin_id=bin_id) for path in output
            )
    for bin_id in passing_metagenome_bin_ids(wildcards, "euk"):
        for output in (rules.metaeuk_bin.output, rules.recognizer_euk_bin.output):
            materials.extend(
                format_metagenome_output(path, bin_id=bin_id) for path in output
            )
    return materials


def isolate_selection_record():
    checkpoint_wildcards = {"sample": SAMPLE} if IMPORT_CONTEXT_ENABLED else {}
    selection = checkpoints.qa_selection.get(**checkpoint_wildcards).output.selection
    rows = read_selection(selection)
    if len(rows) != 1 or rows[0]["sample"] != SAMPLE:
        raise WorkflowError(
            f"Expected exactly one QA selection row for isolate '{SAMPLE}' in {selection}."
        )
    row = rows[0]
    if (
        row["domain"] != ISOLATE_DOMAIN
        or Path(row["path"]).resolve() != Path(INPUT_RECORD["path"]).resolve()
    ):
        raise WorkflowError(
            f"QA selection metadata does not match isolate samplesheet row '{SAMPLE}'."
        )
    return selection, row


def isolate_annotation_gate(wildcards):
    selection, row = isolate_selection_record()
    if row["decision"] == "fail":
        raise WorkflowError(
            f"Annotation is unexpectedly blocked for isolate '{SAMPLE}' by {selection}."
        )
    return selection


def routed_isolate_gate(wildcards, domain):
    selection = isolate_annotation_gate(wildcards)
    if ISOLATE_DOMAIN != domain:
        raise WorkflowError(
            f"Isolate '{SAMPLE}' is not eligible for the {domain} annotation route."
        )
    return selection


def prok_isolate_gate(wildcards):
    return routed_isolate_gate(wildcards, "prok")


def euk_isolate_gate(wildcards):
    return routed_isolate_gate(wildcards, "euk")


def eligible_prokaryotic_gtdbtk_pairs(wildcards):
    if IS_METAGENOME:
        _, rows = metagenome_selection_records(wildcards)
        return sorted(
            (row["sample"], row["path"])
            for row in rows
            if row["domain"] == "prok" and row["decision"] == "pass"
        )
    if ISOLATE_DOMAIN != "prok":
        return []
    _, selection = isolate_selection_record()
    if selection["decision"] == "report_only":
        return [(SAMPLE, INPUT_RECORD["path"])]
    return []


def eligible_prokaryotic_gtdbtk_fastas(wildcards):
    return [path for _, path in eligible_prokaryotic_gtdbtk_pairs(wildcards)]


def gtdbtk_selection_file(wildcards):
    if IS_METAGENOME:
        return metagenome_selection_file(wildcards)
    return prok_isolate_gate(wildcards)


def metagenome_gtdbtk_classification(wildcards):
    if GTDBTK_ENABLED and passing_metagenome_bin_ids(wildcards, "prok"):
        return [format_metagenome_output(rules.gtdbtk.output.classification)]
    return []
