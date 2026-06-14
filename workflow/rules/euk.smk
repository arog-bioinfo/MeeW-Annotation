# ----------------------------------------------------- #
# EUKARYOTIC MAG ANNOTATION RULES                       #
# ----------------------------------------------------- #


# Eukaryotic gene prediction using MetaEuk
# -----------------------------------------------------
rule metaeuk:
    input:
        fasta="data/{sample}.fasta",
    output:
        proteins="results/metaeuk/{sample}.faa",
    log:
        "results/metaeuk/{sample}.log",
    conda:
        "../envs/metaeuk.yaml"
    threads: config.get("threads", {}).get("medium", 8)
    params:
        db=config.get("metaeuk", {}).get("db", ""),
        extra=config.get("metaeuk", {}).get("extra", ""),
    message:
        """--- Running MetaEuk eukaryotic gene prediction for {wildcards.sample}."""
    script:
        "../scripts/metaeuk.py"


# Eukaryotic functional annotation using reCOGnizer
# -----------------------------------------------------
rule recognizer_euk:
    input:
        fasta=rules.metaeuk.output.proteins,
    output:
        tsv="results/recognizer/euk/{sample}/reCOGnizer_results.tsv",
    log:
        "results/recognizer/euk/{sample}.log",
    conda:
        "../envs/recognizer.yaml"
    threads: config.get("threads", {}).get("medium", 8)
    params:
        custom_db=config.get("recognizer_euk", {}).get("custom_db", ""),
        resources_dir=config.get("recognizer_euk", {}).get("resources_dir", ""),
        extra=config.get("recognizer_euk", {}).get("extra", ""),
    message:
        """--- Running eukaryotic reCOGnizer domain annotation for {wildcards.sample}."""
    script:
        "../scripts/recognizer.py"
