# import basic packages
import sys
from pathlib import Path

import pandas as pd
from snakemake.utils import validate

sys.path.insert(
    0, str(Path(workflow.source_path("../scripts/meew_annotation_paths.py")).parent)
)
sys.path.insert(0, str(Path(workflow.source_path("../scripts/filter_qa.py")).parent))
from filter_qa import filter_samples_by_domain
from meew_annotation_paths import expand_path, resolve_resource_path

if "sample_sheet" in config:
    config["sample_sheet"] = expand_path(config["sample_sheet"])

for field in [
    "prok_bins_dir",
    "euk_bins_dir",
    "checkm2_report",
    "eukcc_report",
]:
    if field in config.get("directory_mode", {}):
        config["directory_mode"][field] = expand_path(config["directory_mode"][field])

qa_filter_config = config.get("qa_filter", {})
for field in ["checkm2_reports", "eukcc_reports"]:
    if field in qa_filter_config:
        qa_filter_config[field] = [
            expand_path(path) for path in qa_filter_config[field]
        ]

resource_fields = [
    ("bakta", "db", "BAKTA_DB", "bakta_db/db"),
    ("gtdbtk", "data_dir", "GTDBTK_DATA_PATH", "gtdbtk_db"),
    (
        "recognizer_prok",
        "resources_dir",
        "RECOGNIZER_RESOURCES",
        "recognizer_db",
    ),
    ("upimapi", "resources_dir", "UPIMAPI_RESOURCES", "upimapi_db"),
    ("metaeuk", "db", "METAEUK_DB", "metaeuk_db/uniprot_db"),
    (
        "recognizer_euk",
        "resources_dir",
        "RECOGNIZER_RESOURCES",
        "recognizer_db",
    ),
    ("funannotate2", "db_dir", "FUNANNOTATE2_DB", "funannotate2_db"),
]
for section, field, tool_environment, subpath in resource_fields:
    section_config = config.setdefault(section, {})
    section_config[field] = resolve_resource_path(
        section_config.get(field, ""), tool_environment, subpath
    )

for section, field in [
    ("upimapi", "db_custom"),
    ("recognizer_euk", "custom_db"),
    ("funannotate2", "params"),
    ("funannotate2", "pretrained"),
]:
    if field in config.get(section, {}) and config[section][field]:
        config[section][field] = expand_path(config[section][field])

SAMPLES_SCHEMA = workflow.current_basedir.join(
    "../schemas/samples.schema.yaml"
).get_path_or_uri(secret_free=False)

_samples_cache = None


def load_samples_from_path(path):
    loaded_samples = (
        pd.read_csv(path, sep="\t", dtype={"sample": str})
        .set_index("sample", drop=False)
        .sort_index()
    )
    if "genome_type" not in loaded_samples.columns:
        loaded_samples["genome_type"] = "mag"
    loaded_samples["genome_type"] = loaded_samples["genome_type"].fillna("mag")
    loaded_samples["path"] = loaded_samples["path"].map(
        lambda value: value if pd.isna(value) else expand_path(value)
    )
    validate(loaded_samples, schema=SAMPLES_SCHEMA)
    return loaded_samples


def load_samples():
    global _samples_cache
    if "sample_sheet" not in config:
        raise ValueError("config.sample_sheet is required to load Annotation samples.")
    if _samples_cache is None:
        _samples_cache = filter_samples_by_domain(
            load_samples_from_path(config["sample_sheet"]), config.get("qa_filter", {})
        )
    return _samples_cache


def filtered_samples_from_path(path):
    return filter_samples_by_domain(
        load_samples_from_path(path), config.get("qa_filter", {})
    )


samples = pd.DataFrame(columns=["sample", "path", "domain", "genome_type"]).set_index(
    "sample", drop=False
)


# validate sample sheet and config file
validate(config, schema="../schemas/config.schema.yaml")


def sample_fasta(wildcards):
    return load_samples().loc[wildcards.sample, "path"]


def prokaryotic_gtdbtk_inputs(wildcards):
    current_samples = load_samples()
    return current_samples.loc[current_samples["domain"] == "prok", "path"].tolist()


def directory_mode_path(key, wildcards):
    value = config.get("directory_mode", {}).get(key, "")
    return str(value).format(sample=wildcards.sample)


def directory_mode_enabled():
    return config.get("directory_mode", {}).get("enabled", False)


def annotation_to_metabolic_modeling_proteins(wildcards):
    return [
        str(rules.bakta.output.faa).format(sample=sample)
        for sample in prokaryotic_samples()
    ]
