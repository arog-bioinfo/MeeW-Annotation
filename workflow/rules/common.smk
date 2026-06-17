# import basic packages
import sys
from pathlib import Path

import pandas as pd
from snakemake.utils import validate

sys.path.insert(0, str(Path(workflow.basedir) / "scripts"))
from filter_qa import filter_samples_by_domain

# read sample sheet
samples = (
    pd.read_csv(config["sample_sheet"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)
    .sort_index()
)
if "genome_type" not in samples.columns:
    samples["genome_type"] = "mag"
samples["genome_type"] = samples["genome_type"].fillna("mag")


# validate sample sheet and config file
validate(samples, schema="../schemas/samples.schema.yaml")
validate(config, schema="../schemas/config.schema.yaml")

# filter samples by external QA reports before expanding annotation targets
samples = filter_samples_by_domain(samples, config.get("qa_filter", {}))
