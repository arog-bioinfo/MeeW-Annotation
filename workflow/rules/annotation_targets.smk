# ----------------------------------------------------- #
# Annotation workflow targets                           #
# ----------------------------------------------------- #


def _format_rule_output(template, **wildcards):
    return str(template).format(**wildcards)


def prokaryotic_samples():
    current_samples = load_samples()
    return current_samples[current_samples["domain"] == "prok"].index.tolist()


def eukaryotic_mag_samples():
    current_samples = load_samples()
    return current_samples[
        (current_samples["domain"] == "euk") & (current_samples["genome_type"] == "mag")
    ].index.tolist()


def eukaryotic_isolate_samples():
    if not config.get("funannotate2", {}).get("enabled", True):
        return []
    current_samples = load_samples()
    return current_samples[
        (current_samples["domain"] == "euk")
        & (current_samples["genome_type"] == "isolate")
    ].index.tolist()


def final_targets(wildcards):
    if directory_mode_enabled():
        batch_sample = config.get("directory_mode", {}).get("sample", "")
        return [
            _format_rule_output(rules.bakta_batch.output.done, sample=batch_sample),
            _format_rule_output(rules.prodigal_batch.output.done, sample=batch_sample),
            _format_rule_output(
                rules.recognizer_prok_batch.output.done, sample=batch_sample
            ),
            _format_rule_output(rules.upimapi_batch.output.done, sample=batch_sample),
            _format_rule_output(rules.metaeuk_batch.output.done, sample=batch_sample),
            _format_rule_output(
                rules.recognizer_euk_batch.output.done, sample=batch_sample
            ),
        ]

    prok_samples = prokaryotic_samples()
    euk_mag_samples = eukaryotic_mag_samples()
    euk_isolate_samples = eukaryotic_isolate_samples()
    gtdbtk_targets = []
    if config.get("gtdbtk", {}).get("enabled", False) and prok_samples:
        gtdbtk_targets = [
            _format_rule_output(rules.gtdbtk.output.done, sample=sample)
            for sample in prok_samples
        ]
    return (
        [
            _format_rule_output(rules.bakta.output.faa, sample=sample)
            for sample in prok_samples
        ]
        + [
            _format_rule_output(rules.recognizer_prok.output.tsv, sample=sample)
            for sample in prok_samples
        ]
        + [
            _format_rule_output(rules.upimapi.output.results, sample=sample)
            for sample in prok_samples
        ]
        + [
            _format_rule_output(rules.recognizer_euk.output.tsv, sample=sample)
            for sample in euk_mag_samples
        ]
        + [
            _format_rule_output(rules.funannotate2_annotate.output.done, sample=sample)
            for sample in euk_isolate_samples
        ]
        + gtdbtk_targets
        + [
            _format_rule_output(
                rules.annotation_to_metabolic_modeling_stage_sheet.output[0],
                sample=sample,
            )
            for sample in prok_samples
        ]
    )


FINAL_TARGETS = final_targets
