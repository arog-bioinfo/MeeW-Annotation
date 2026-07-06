# import basic packages
import sys
from pathlib import Path

import pandas as pd
from snakemake.utils import validate

sys.path.insert(0, str(Path(workflow.source_path("../scripts/filter_qa.py")).parent))
from filter_qa import filter_samples_by_domain

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
