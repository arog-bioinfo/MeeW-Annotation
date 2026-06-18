# ----------------------------------------------------- #
# Annotation workflow targets                           #
# ----------------------------------------------------- #


prok_samples = samples[samples["domain"] == "prok"].index.tolist()
euk_mag_samples = samples[
    (samples["domain"] == "euk") & (samples["genome_type"] == "mag")
].index.tolist()
euk_isolate_samples = samples[
    (samples["domain"] == "euk") & (samples["genome_type"] == "isolate")
].index.tolist()
if not config.get("funannotate2", {}).get("enabled", True):
    euk_isolate_samples = []
gtdbtk_targets = []
if config.get("gtdbtk", {}).get("enabled", False) and prok_samples:
    gtdbtk_targets = ["results/gtdbtk/gtdbtk.done"]


FINAL_TARGETS = (
    expand("results/bakta/{sample}/{sample}.faa", sample=prok_samples)
    + expand(
        "results/recognizer/prok/{sample}/reCOGnizer_results.tsv",
        sample=prok_samples,
    )
    + expand("results/upimapi/{sample}/uniprotinfo.tsv", sample=prok_samples)
    + expand(
        "results/recognizer/euk/{sample}/reCOGnizer_results.tsv",
        sample=euk_mag_samples,
    )
    + expand(
        "results/funannotate2/{sample}/annotate.done",
        sample=euk_isolate_samples,
    )
    + gtdbtk_targets
)
