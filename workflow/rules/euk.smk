# ----------------------------------------------------- #
# EUKARYOTIC MAG ANNOTATION RULES                       #
# ----------------------------------------------------- #


# Eukaryotic gene prediction using MetaEuk
# -----------------------------------------------------
rule metaeuk:
    input:
        fasta=sample_fasta,
    output:
        proteins="<results>/metaeuk/{sample}.faa",
    log:
        "<results>/metaeuk/{sample}.log",
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
        tsv="<results>/recognizer/euk/reCOGnizer_results.tsv",
    log:
        "<results>/recognizer/euk/{sample}.log",
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


# Funannotate2 database check/install for eukaryotic isolate annotation
# -----------------------------------------------------
rule funannotate2_db:
    output:
        done=touch("<results>/funannotate2/db/funannotate2_db.done"),
    log:
        "<results>/funannotate2/db/funannotate2_db.log",
    conda:
        "../envs/funannotate2.yaml"
    params:
        step="install_db",
        db_dir=config.get("funannotate2", {}).get(
            "db_dir", "/home/argomes/resources/funannotate2_db"
        ),
        install_db=config.get("funannotate2", {}).get("install_db", False),
        databases=config.get("funannotate2", {}).get("databases", ["all"]),
        extra_install=config.get("funannotate2", {}).get("extra_install", ""),
    message:
        """--- Checking Funannotate2 databases."""
    script:
        "../scripts/funannotate2.py"


# Funannotate2 clean for eukaryotic isolate annotation
# -----------------------------------------------------
rule funannotate2_clean:
    input:
        db=rules.funannotate2_db.output.done,
        fasta=sample_fasta,
    output:
        outdir=directory("<results>/funannotate2/clean"),
        cleaned_fasta="<results>/funannotate2/clean/cleaned.fasta",
        done="<results>/funannotate2/clean.done",
    log:
        "<results>/funannotate2/clean.{sample}.log",
    conda:
        "../envs/funannotate2.yaml"
    threads: config.get("threads", {}).get("medium", 8)
    params:
        step="clean",
        db_dir=config.get("funannotate2", {}).get(
            "db_dir", "/home/argomes/resources/funannotate2_db"
        ),
        extra_clean=config.get("funannotate2", {}).get("extra_clean", ""),
    message:
        """--- Running Funannotate2 clean for {wildcards.sample}."""
    script:
        "../scripts/funannotate2.py"


# Funannotate2 train for eukaryotic isolate annotation
# -----------------------------------------------------
rule funannotate2_train:
    input:
        cleaned_fasta=rules.funannotate2_clean.output.cleaned_fasta,
        clean_done=rules.funannotate2_clean.output.done,
    output:
        done="<results>/funannotate2/train.done",
    log:
        "<results>/funannotate2/train.{sample}.log",
    conda:
        "../envs/funannotate2.yaml"
    threads: config.get("threads", {}).get("medium", 8)
    params:
        step="train",
        db_dir=config.get("funannotate2", {}).get(
            "db_dir", "/home/argomes/resources/funannotate2_db"
        ),
        outdir=lambda wc, output: str(Path(output.done).parent / "run"),
        species=config.get("funannotate2", {}).get("species", ""),
        strain=config.get("funannotate2", {}).get("strain", ""),
        extra_train=config.get("funannotate2", {}).get("extra_train", ""),
    message:
        """--- Running Funannotate2 train for {wildcards.sample}."""
    script:
        "../scripts/funannotate2.py"


# Funannotate2 predict for eukaryotic isolate annotation
# -----------------------------------------------------
rule funannotate2_predict:
    input:
        train_done=rules.funannotate2_train.output.done,
    output:
        done="<results>/funannotate2/predict.done",
    log:
        "<results>/funannotate2/predict.{sample}.log",
    conda:
        "../envs/funannotate2.yaml"
    threads: config.get("threads", {}).get("medium", 8)
    params:
        step="predict",
        db_dir=config.get("funannotate2", {}).get(
            "db_dir", "/home/argomes/resources/funannotate2_db"
        ),
        outdir=lambda wc, output: str(Path(output.done).parent / "run"),
        species=config.get("funannotate2", {}).get("species", ""),
        strain=config.get("funannotate2", {}).get("strain", ""),
        params=config.get("funannotate2", {}).get("params", ""),
        pretrained=config.get("funannotate2", {}).get("pretrained", ""),
        extra_predict=config.get("funannotate2", {}).get("extra_predict", ""),
    message:
        """--- Running Funannotate2 predict for {wildcards.sample}."""
    script:
        "../scripts/funannotate2.py"


# Funannotate2 annotate for eukaryotic isolate annotation
# -----------------------------------------------------
rule funannotate2_annotate:
    input:
        predict_done=rules.funannotate2_predict.output.done,
    output:
        done="<results>/funannotate2/annotate.done",
    log:
        "<results>/funannotate2/annotate.{sample}.log",
    conda:
        "../envs/funannotate2.yaml"
    threads: config.get("threads", {}).get("medium", 8)
    params:
        step="annotate",
        db_dir=config.get("funannotate2", {}).get(
            "db_dir", "/home/argomes/resources/funannotate2_db"
        ),
        outdir=lambda wc, output: str(Path(output.done).parent / "run"),
        species=config.get("funannotate2", {}).get("species", ""),
        strain=config.get("funannotate2", {}).get("strain", ""),
        extra_annotate=config.get("funannotate2", {}).get("extra_annotate", ""),
    message:
        """--- Running Funannotate2 annotate for {wildcards.sample}."""
    script:
        "../scripts/funannotate2.py"
