# ----------------------------------------------------- #
# MAG ANNOTATION WORKFLOW                               #
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
    threads: config.get("threads", {}).get("low", 1)
    params:
        extra=config.get("prodigal", {}).get("extra", "")
    message:
        """--- Predicting genes with Prodigal for {wildcards.sample}."""
    conda:
        "../envs/prodigal.yaml"
    script:
        "../scripts/prodigal.py"


# Functional annotation using reCOGnizer
# -----------------------------------------------------
rule recognizer:
    input:
        fasta=rules.prodigal.output.faa,
    output:
        tsv="results/recognizer/{sample}/reCOGnizer_results.tsv",
    log:
        "results/recognizer/{sample}.log",
    threads: config.get("threads", {}).get("medium", 8)
    params:
        resources_dir=config.get("recognizer", {}).get("resources_dir", ""),
        extra=config.get("recognizer", {}).get("extra", "")
    message:
        """--- Running reCOGnizer domain annotation for {wildcards.sample}."""
    conda:
        "../envs/recognizer.yaml"
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
    threads: config.get("threads", {}).get("medium", 8)
    params:
        extra=config.get("upimapi", {}).get("extra", "")
    message:
        """--- Running UPIMAPI protein mapping for {wildcards.sample}."""
    conda:
        "../envs/upimapi.yaml"
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
    threads: config.get("threads", {}).get("high", 16)
    params:
        db=config.get("bakta", {}).get("db", ""),
        extra=config.get("bakta", {}).get("extra", "")
    message:
        """--- Running Bakta comprehensive annotation for {wildcards.sample}."""
    conda:
        "../envs/bakta.yaml"
    script:
        "../scripts/bakta.py"


# Taxonomic classification using GTDB-Tk
# -----------------------------------------------------
rule gtdbtk:
    input:
        bins="data/",
    output:
        outdir=directory("results/gtdbtk"),
    log:
        "results/gtdbtk/gtdbtk.log",
    threads: config.get("threads", {}).get("high", 16)
    params:
        db_path=config.get("gtdbtk", {}).get("data_dir", ""), # <-- Mudei de data_dir para db_path aqui
        extra=config.get("gtdbtk", {}).get("extra", "")
    message:
        """--- Running GTDB-Tk taxonomic classification for all MAGs."""
    conda:
        "../envs/gtdb-tk.yaml"
    script:
        "../scripts/gtdb-tk.py"


# Eukaryotic gene prediction using MetaEuk
# -----------------------------------------------------
rule metaeuk:
    input:
        fasta="data/{sample}.fasta",
    output:
        proteins="results/metaeuk/{sample}.faa",
    log:
        "results/metaeuk/{sample}.log",
    threads: config.get("threads", {}).get("medium", 8)
    params:
        db=config.get("metaeuk", {}).get("db", ""),
        extra=config.get("metaeuk", {}).get("extra", "")
    message:
        """--- Running MetaEuk eukaryotic gene prediction for {wildcards.sample}."""
    conda:
        "../envs/metaeuk.yaml"
    script:
        "../scripts/metaeuk.py"

    
    
    
   
