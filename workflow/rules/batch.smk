# ----------------------------------------------------- #
# Directory/batch Annotation rules                       #
# ----------------------------------------------------- #


def directory_mode_path(key, wildcards):
    value = config.get("directory_mode", {}).get(key, "")
    return str(value).format(sample=wildcards.sample)


def directory_mode_enabled():
    return config.get("directory_mode", {}).get("enabled", False)


rule filter_passing_prok_bins:
    input:
        bins_dir=lambda wc: directory_mode_path("prok_bins_dir", wc),
        report=lambda wc: directory_mode_path("checkm2_report", wc),
    output:
        passing_dir=directory("<results>/batch/qa/prok/passing_bins"),
        manifest="<results>/batch/qa/prok/passing_bins.tsv",
        done=touch("<results>/batch/qa/prok/passing_bins.done"),
    params:
        action="filter_prok",
        min_completeness=config.get("qa_filter", {}).get("min_completeness", 50.0),
        max_contamination=config.get("qa_filter", {}).get("max_contamination", 10.0),
    message:
        "--- Filtering passing prokaryotic bins from a Binning directory."
    script:
        "../scripts/batch_annotation.py"


rule filter_passing_euk_bins:
    input:
        bins_dir=lambda wc: directory_mode_path("euk_bins_dir", wc),
        report=lambda wc: directory_mode_path("eukcc_report", wc),
    output:
        passing_dir=directory("<results>/batch/qa/euk/passing_bins"),
        manifest="<results>/batch/qa/euk/passing_bins.tsv",
        done=touch("<results>/batch/qa/euk/passing_bins.done"),
    params:
        action="filter_euk",
        min_completeness=config.get("qa_filter", {}).get("min_completeness", 50.0),
        max_contamination=config.get("qa_filter", {}).get("max_contamination", 10.0),
    message:
        "--- Filtering passing eukaryotic bins from a Binning directory."
    script:
        "../scripts/batch_annotation.py"


rule bakta_batch:
    input:
        manifest=rules.filter_passing_prok_bins.output.manifest,
        passing_done=rules.filter_passing_prok_bins.output.done,
    output:
        outdir=directory("<results>/batch/bakta"),
        done=touch("<results>/batch/bakta.done"),
    log:
        "<results>/batch/logs/bakta.log",
    conda:
        "../envs/bakta.yaml"
    threads: config.get("threads", {}).get("high", 16)
    params:
        action="bakta",
        db=config.get("bakta", {}).get("db", ""),
        extra=config.get("bakta", {}).get("extra", ""),
    message:
        "--- Running batch Bakta over passing prokaryotic bins."
    script:
        "../scripts/batch_annotation.py"


rule prodigal_batch:
    input:
        manifest=rules.filter_passing_prok_bins.output.manifest,
        passing_done=rules.filter_passing_prok_bins.output.done,
    output:
        outdir=directory("<results>/batch/prodigal"),
        done=touch("<results>/batch/prodigal.done"),
    log:
        "<results>/batch/logs/prodigal.log",
    conda:
        "../envs/prodigal.yaml"
    threads: config.get("threads", {}).get("low", 1)
    params:
        action="prodigal",
        extra=config.get("prodigal", {}).get("extra", ""),
    message:
        "--- Running batch Prodigal over passing prokaryotic bins."
    script:
        "../scripts/batch_annotation.py"


rule recognizer_prok_batch:
    input:
        manifest=rules.filter_passing_prok_bins.output.manifest,
        proteins_dir=rules.prodigal_batch.output.outdir,
        prodigal_done=rules.prodigal_batch.output.done,
    output:
        outdir=directory("<results>/batch/recognizer/prok"),
        done=touch("<results>/batch/recognizer/prok.done"),
    log:
        "<results>/batch/logs/recognizer_prok.log",
    conda:
        "../envs/recognizer.yaml"
    threads: config.get("threads", {}).get("medium", 8)
    params:
        action="recognizer_prok",
        resources_dir=config.get("recognizer_prok", {}).get("resources_dir", ""),
        extra=config.get("recognizer_prok", {}).get("extra", ""),
    message:
        "--- Running batch prokaryotic reCOGnizer over Prodigal outputs."
    script:
        "../scripts/batch_annotation.py"


rule upimapi_batch:
    input:
        manifest=rules.filter_passing_prok_bins.output.manifest,
        proteins_dir=rules.prodigal_batch.output.outdir,
        prodigal_done=rules.prodigal_batch.output.done,
    output:
        outdir=directory("<results>/batch/upimapi"),
        done=touch("<results>/batch/upimapi.done"),
    log:
        "<results>/batch/logs/upimapi.log",
    conda:
        "../envs/upimapi.yaml"
    threads: config.get("threads", {}).get("medium", 8)
    params:
        action="upimapi",
        db=config.get("upimapi", {}).get("db", "swissprot"),
        db_custom=config.get("upimapi", {}).get("db_custom", ""),
        resources_dir=config.get("upimapi", {}).get("resources_dir", ""),
        extra=config.get("upimapi", {}).get("extra", ""),
        skip_db_check_if_exists=config.get("upimapi", {}).get(
            "skip_db_check_if_exists", True
        ),
    message:
        "--- Running batch UPIMAPI over Prodigal outputs."
    script:
        "../scripts/batch_annotation.py"


rule metaeuk_batch:
    input:
        manifest=rules.filter_passing_euk_bins.output.manifest,
        passing_done=rules.filter_passing_euk_bins.output.done,
    output:
        outdir=directory("<results>/batch/metaeuk"),
        done=touch("<results>/batch/metaeuk.done"),
    log:
        "<results>/batch/logs/metaeuk.log",
    conda:
        "../envs/metaeuk.yaml"
    threads: config.get("threads", {}).get("medium", 8)
    params:
        action="metaeuk",
        db=config.get("metaeuk", {}).get("db", ""),
        extra=config.get("metaeuk", {}).get("extra", ""),
    message:
        "--- Running batch MetaEuk over passing eukaryotic bins."
    script:
        "../scripts/batch_annotation.py"


rule recognizer_euk_batch:
    input:
        manifest=rules.filter_passing_euk_bins.output.manifest,
        proteins_dir=rules.metaeuk_batch.output.outdir,
        metaeuk_done=rules.metaeuk_batch.output.done,
    output:
        outdir=directory("<results>/batch/recognizer/euk"),
        done=touch("<results>/batch/recognizer/euk.done"),
    log:
        "<results>/batch/logs/recognizer_euk.log",
    conda:
        "../envs/recognizer.yaml"
    threads: config.get("threads", {}).get("medium", 8)
    params:
        action="recognizer_euk",
        resources_dir=config.get("recognizer_euk", {}).get("resources_dir", ""),
        custom_db=config.get("recognizer_euk", {}).get("custom_db", ""),
        extra=config.get("recognizer_euk", {}).get("extra", ""),
    message:
        "--- Running batch eukaryotic reCOGnizer over MetaEuk outputs."
    script:
        "../scripts/batch_annotation.py"
