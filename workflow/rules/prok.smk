# ----------------------------------------------------- #
# PROKARYOTIC MAG ANNOTATION RULES                      #
# ----------------------------------------------------- #


def prokaryotic_gtdbtk_inputs(wildcards):
    return samples.loc[samples["domain"] == "prok", "path"].tolist()


# Predict protein-coding genes using Prodigal
# -----------------------------------------------------
rule prodigal:
    input:
        fasta=sample_fasta,
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
        resources_dir=config.get("recognizer_prok", {}).get("resources_dir", ""),
        extra=config.get("recognizer_prok", {}).get("extra", ""),
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
        fasta=sample_fasta,
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


# Stage only prokaryotic genomes for GTDB-Tk batch classification
# -----------------------------------------------------
rule stage_gtdbtk_genomes:
    input:
        prokaryotic_gtdbtk_inputs,
    output:
        genomes=directory("results/gtdbtk_genomes"),
    message:
        """--- Staging prokaryotic genomes for GTDB-Tk."""
    run:
        import shutil
        from pathlib import Path

        genomes = Path(output.genomes)
        if genomes.exists():
            shutil.rmtree(genomes)
        genomes.mkdir(parents=True)
        for _, row in samples[samples["domain"] == "prok"].iterrows():
            source = Path(row["path"])
            target = genomes / f"{row['sample']}.fasta"
            target.symlink_to(source.resolve())


# Optional prokaryotic genome classification using GTDB-Tk
# -----------------------------------------------------
rule gtdbtk:
    input:
        genomes=rules.stage_gtdbtk_genomes.output.genomes,
    output:
        done=touch("results/gtdbtk/gtdbtk.done"),
    log:
        "results/gtdbtk/gtdbtk.log",
    conda:
        "../envs/gtdb-tk.yaml"
    threads: config.get("threads", {}).get("high", 16)
    params:
        db_path=config.get("gtdbtk", {}).get("data_dir", ""),
        extension="fasta",
        extra=config.get("gtdbtk", {}).get("extra", ""),
    message:
        """--- Running GTDB-Tk classification for prokaryotic genomes."""
    script:
        "../scripts/gtdb-tk.py"
