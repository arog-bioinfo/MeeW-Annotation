# ----------------------------------------------------- #
# Prokaryotic annotation rules                           #
# ----------------------------------------------------- #


rule bakta:
    input:
        qa=prok_isolate_gate,
        fasta=sample_fasta,
        database=lambda wc: bakta_database_inputs(BAKTA_DB_DIR, BAKTA_DB_TYPE),
    output:
        results=directory("<results>/bakta"),
    log:
        "<results>/logs/bakta.log",
    conda:
        "../envs/bakta.yaml"
    threads: config.get("threads", {}).get("high", 16)
    params:
        db=lambda wc, output: checked_cache_path(
            BAKTA_DB_DIR, Path(output.results).parent, "Bakta"
        ),
        db_type=BAKTA_DB_TYPE,
        prefix="annotation",
        meta=False,
        extra=config.get("bakta", {}).get("extra", ""),
    message:
        f"--- Running Bakta comprehensive annotation for {SAMPLE}."
    script:
        "../scripts/bakta.py"


rule recognizer_prok:
    input:
        annotations=rules.bakta.output.results,
        database=lambda wc: recognizer_database_inputs(
            config.get("recognizer_prok", {}).get("resources_dir", ""),
            "NCBI_Curated,Pfam,SMART,KOG,COG,PRK,TIGR",
        ),
    output:
        results=directory("<results>/recognizer/prok"),
    log:
        "<results>/logs/recognizer_prok.log",
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
        f"--- Running prokaryotic reCOGnizer for {SAMPLE}."
    script:
        "../scripts/recognizer.py"


rule upimapi:
    input:
        annotations=rules.bakta.output.results,
        database=lambda wc: upimapi_database_inputs(
            config.get("upimapi", {}).get("resources_dir", ""),
            config.get("upimapi", {}).get("db", "uniprot"),
            config.get("upimapi", {}).get("db_custom", ""),
        ),
    output:
        results=directory("<results>/upimapi"),
    log:
        "<results>/logs/upimapi.log",
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
        f"--- Running UPIMAPI protein mapping for {SAMPLE}."
    script:
        "../scripts/upimapi.py"


rule protein_summary:
    input:
        bakta=rules.bakta.output.results,
        recognizer=rules.recognizer_prok.output.results,
        upimapi=rules.upimapi.output.results,
    output:
        summary="<results>/protein_annotations.tsv",
    log:
        "<results>/protein_annotations.log",
    conda:
        "../envs/python.yaml"
    message:
        f"--- Summarizing prokaryotic protein annotations for {SAMPLE}."
    script:
        "../scripts/protein_summary.py"
