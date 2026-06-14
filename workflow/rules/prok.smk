# ----------------------------------------------------- #
# PROKARYOTIC MAG ANNOTATION RULES                      #
# ----------------------------------------------------- #


# Predict protein-coding genes using Prodigal
# -----------------------------------------------------
rule prodigal:
    input:
        fasta="data/{sample}.fasta",
    output:
        out="results/prodigal/{sample}.gff",
        faa="results/prodigal/{sample}.faa",
        fna="results/prodigal/{sample}.fna",
    log:
        "results/prodigal/{sample}.log",
    conda:
        "../envs/prodigal.yaml"
    threads: config.get("threads", {}).get("low", 1)
    params:
        extra=config.get("prodigal", {}).get("extra", ""),
    message:
        """--- Predicting genes with Prodigal for {wildcards.sample}."""
    script:
        "../scripts/prodigal.py"


# Prokaryotic functional annotation using reCOGnizer
# -----------------------------------------------------
rule recognizer_prok:
    input:
        fasta=rules.prodigal.output.faa,
    output:
        tsv="results/recognizer/prok/{sample}/reCOGnizer_results.tsv",
    log:
        "results/recognizer/prok/{sample}.log",
    conda:
        "../envs/recognizer.yaml"
    threads: config.get("threads", {}).get("medium", 8)
    params:
        resources_dir=config.get("recognizer", {}).get("resources_dir", ""),
        extra=config.get("recognizer", {}).get("extra", ""),
    message:
        """--- Running prokaryotic reCOGnizer domain annotation for {wildcards.sample}."""
    script:
        "../scripts/recognizer.py"


# Functional annotation via UniProt using UPIMAPI
# -----------------------------------------------------
rule upimapi:
    input:
        fasta=rules.prodigal.output.faa,
    output:
        outdir=directory("results/upimapi/{sample}"),
        results="results/upimapi/{sample}/uniprotinfo.tsv",
    log:
        "results/upimapi/{sample}.log",
    conda:
        "../envs/upimapi.yaml"
    threads: config.get("threads", {}).get("medium", 8)
    params:
        db=config.get("upimapi", {}).get("db", "swissprot"),
        db_custom=config.get("upimapi", {}).get("db_custom", ""),
        resources_dir=config.get("upimapi", {}).get(
            "resources_dir", "resources/upimapi_db"
        ),
        extra=config.get("upimapi", {}).get("extra", ""),
        skip_db_check_if_exists=config.get("upimapi", {}).get(
            "skip_db_check_if_exists", True
        ),
    message:
        """--- Running UPIMAPI protein mapping for {wildcards.sample}."""
    script:
        "../scripts/upimapi.py"


# Comprehensive genome annotation using Bakta
# -----------------------------------------------------
rule bakta:
    input:
        fasta="data/{sample}.fasta",
    output:
        outdir=directory("results/bakta/{sample}"),
        gff="results/bakta/{sample}/{sample}.gff3",
        faa="results/bakta/{sample}/{sample}.faa",
    log:
        "results/bakta/{sample}.log",
    conda:
        "../envs/bakta.yaml"
    threads: config.get("threads", {}).get("high", 16)
    params:
        db=config.get("bakta", {}).get("db", ""),
        extra=config.get("bakta", {}).get("extra", ""),
    message:
        """--- Running Bakta comprehensive annotation for {wildcards.sample}."""
    script:
        "../scripts/bakta.py"
