# import basic packages
import sys
from pathlib import Path

import pandas as pd
from snakemake.utils import validate

sys.path.insert(0, str(Path(workflow.source_path("../scripts/filter_qa.py")).parent))
from filter_qa import filter_samples_by_domain

DEFERRED_SAMPLE_SHEET = bool(config.get("deferred_sample_sheet", False))
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
    validate(loaded_samples, schema="../schemas/samples.schema.yaml")
    return loaded_samples


def load_samples():
    global _samples_cache
    if _samples_cache is None:
        _samples_cache = load_samples_from_path(config["sample_sheet"])
    return _samples_cache


def filtered_samples_from_path(path):
    return filter_samples_by_domain(
        load_samples_from_path(path), config.get("qa_filter", {})
    )


if DEFERRED_SAMPLE_SHEET:
    samples = pd.DataFrame(
        columns=["sample", "path", "domain", "genome_type"]
    ).set_index("sample", drop=False)
else:
    samples = load_samples()


# validate sample sheet and config file
validate(config, schema="../schemas/config.schema.yaml")

# filter samples by external QA reports before expanding annotation targets
samples = filter_samples_by_domain(samples, config.get("qa_filter", {}))


def sample_fasta(wildcards):
    return load_samples().loc[wildcards.sample, "path"]
