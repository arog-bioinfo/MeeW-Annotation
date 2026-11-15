# ----------------------------------------------------- #
# Checkpoint-derived metagenome annotation rules         #
# ----------------------------------------------------- #


rule bakta_bin:
    input:
        qa=metagenome_selection_file,
        fasta=metagenome_prok_bin_fasta,
        database=lambda wc: bakta_database_inputs(BAKTA_DB_DIR, BAKTA_DB_TYPE),
    output:
        results=directory("<results>/batch/prok/{bin}/bakta"),
    log:
        "<results>/batch/logs/prok/{bin}/bakta.log",
    conda:
        "../envs/bakta.yaml"
    threads: config.get("threads", {}).get("high", 16)
    params:
        db=lambda wc, output: checked_cache_path(
            BAKTA_DB_DIR, Path(output.results).parents[3], "Bakta"
        ),
        db_type=BAKTA_DB_TYPE,
        prefix="annotation",
        meta=True,
        extra=config.get("bakta", {}).get("extra", ""),
    message:
        "--- Running Bakta for passing prokaryotic bin {wildcards.bin}."
    script:
        "../scripts/bakta.py"


rule recognizer_prok_bin:
    input:
        annotations=rules.bakta_bin.output.results,
        database=lambda wc: recognizer_database_inputs(
            config.get("recognizer_prok", {}).get("resources_dir", ""),
            "NCBI_Curated,Pfam,SMART,KOG,COG,PRK,TIGR",
        ),
    output:
        results=directory("<results>/batch/prok/{bin}/recognizer"),
    log:
        "<results>/batch/logs/prok/{bin}/recognizer.log",
    conda:
        "../envs/recognizer.yaml"
    threads: config.get("threads", {}).get("medium", 8)
    params:
        fasta_name="annotation.faa",
        databases="NCBI_Curated,Pfam,SMART,KOG,COG,PRK,TIGR",
        extra_config_key="config.recognizer_prok.extra",
        resources_dir=config.get("recognizer_prok", {}).get("resources_dir", ""),
        extra=config.get("recognizer_prok", {}).get("extra", ""),
    message:
        "--- Running reCOGnizer for passing prokaryotic bin {wildcards.bin}."
    script:
        "../scripts/recognizer.py"


rule upimapi_bin:
    input:
        annotations=rules.bakta_bin.output.results,
        database=lambda wc: upimapi_database_inputs(
            config.get("upimapi", {}).get("resources_dir", ""),
            config.get("upimapi", {}).get("db", "uniprot"),
            config.get("upimapi", {}).get("db_custom", ""),
        ),
    output:
        results=directory("<results>/batch/prok/{bin}/upimapi"),
    log:
        "<results>/batch/logs/prok/{bin}/upimapi.log",
    conda:
        "../envs/upimapi.yaml"
    threads: config.get("threads", {}).get("medium", 8)
    params:
        fasta_name="annotation.faa",
        db=config.get("upimapi", {}).get("db", "uniprot"),
        db_custom=config.get("upimapi", {}).get("db_custom", ""),
        resources_dir=config.get("upimapi", {}).get("resources_dir", ""),
        extra=config.get("upimapi", {}).get("extra", ""),
        skip_db_check_if_exists=config.get("upimapi", {}).get(
            "skip_db_check_if_exists", True
        ),
    message:
        "--- Running UPIMAPI for passing prokaryotic bin {wildcards.bin}."
    script:
        "../scripts/upimapi.py"


rule protein_summary_bin:
    input:
        bakta=rules.bakta_bin.output.results,
        recognizer=rules.recognizer_prok_bin.output.results,
        upimapi=rules.upimapi_bin.output.results,
    output:
        summary="<results>/batch/prok/{bin}/protein_annotations.tsv",
    log:
        "<results>/batch/logs/prok/{bin}/protein_annotations.log",
    conda:
        "../envs/python.yaml"
    message:
        "--- Summarizing annotations for passing prokaryotic bin {wildcards.bin}."
    script:
        "../scripts/protein_summary.py"


rule metaeuk_bin:
    input:
        qa=metagenome_selection_file,
        fasta=metagenome_euk_bin_fasta,
        database=lambda wc: metaeuk_database_inputs(METAEUK_DB_PREFIX),
    output:
        results=directory("<results>/batch/euk/{bin}/metaeuk"),
    log:
        "<results>/batch/logs/euk/{bin}/metaeuk.log",
    conda:
        "../envs/metaeuk.yaml"
    threads: config.get("threads", {}).get("medium", 8)
    params:
        db=lambda wc, output: checked_cache_path(
            METAEUK_DB_PREFIX, Path(output.results).parents[3], "MetaEuk"
        ),
        catalog=METAEUK_DATABASE,
        database_managed=METAEUK_DB_MANAGED,
        max_exon_sets=config.get("metaeuk", {}).get("max_exon_sets", 1),
        extra=config.get("metaeuk", {}).get("extra", ""),
    message:
        "--- Running MetaEuk for passing eukaryotic bin {wildcards.bin}."
    script:
        "../scripts/metaeuk.py"


rule recognizer_euk_bin:
    input:
        annotations=rules.metaeuk_bin.output.results,
        database=lambda wc: recognizer_database_inputs(
            config.get("recognizer_euk", {}).get("resources_dir", ""),
            "KOG",
        ),
    output:
        results=directory("<results>/batch/euk/{bin}/recognizer"),
    log:
        "<results>/batch/logs/euk/{bin}/recognizer.log",
    conda:
        "../envs/recognizer.yaml"
    threads: config.get("threads", {}).get("medium", 8)
    params:
        fasta_name="native_predictions.fas",
        databases="KOG",
        extra_config_key="config.recognizer_euk.extra",
        resources_dir=config.get("recognizer_euk", {}).get("resources_dir", ""),
        extra=config.get("recognizer_euk", {}).get("extra", ""),
    message:
        "--- Running reCOGnizer for passing eukaryotic bin {wildcards.bin}."
    script:
        "../scripts/recognizer.py"


rule metagenome_annotation_inventory:
    input:
        selection=metagenome_selection_file,
        materials=metagenome_annotation_materials,
        gtdbtk=metagenome_gtdbtk_classification,
    output:
        inventory="<results>/batch/annotation_results.tsv",
    log:
        "<results>/batch/logs/annotation_inventory.log",
    conda:
        "../envs/python.yaml"
    params:
        results_root=lambda wc, output: str(Path(output.inventory).parents[1]),
        project=SAMPLE,
        gtdbtk_enabled=GTDBTK_ENABLED,
    script:
        "../scripts/write_directory_inventories.py"
