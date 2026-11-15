# ----------------------------------------------------- #
# Eukaryotic annotation rules                            #
# ----------------------------------------------------- #


rule funannotate2:
    input:
        qa=euk_isolate_gate,
        fasta=sample_fasta,
        model=funannotate2_model_input,
        database=lambda wc: funannotate2_database_inputs(
            FUNANNOTATE2_DB_DIR,
            config.get("funannotate2", {}).get("databases", ["all"]),
        ),
    output:
        project=directory("<results>/funannotate2/project"),
        gff3="<results>/funannotate2/annotation.gff3",
        proteins="<results>/funannotate2/proteins.faa",
        transcripts="<results>/funannotate2/transcripts.fna",
        summary="<results>/funannotate2/summary.json",
        provenance="<results>/funannotate2/provenance.json",
    log:
        "<results>/funannotate2/funannotate2.log",
    conda:
        "../envs/funannotate2.yaml"
    threads: config.get("threads", {}).get("medium", 8)
    params:
        db_dir=lambda wc, output: checked_cache_path(
            FUNANNOTATE2_DB_DIR,
            Path(output.project).parents[1],
            "Funannotate2",
        ),
        databases=config.get("funannotate2", {}).get("databases", ["all"]),
        clean=config.get("funannotate2", {}).get("clean", False),
        species=config.get("funannotate2", {}).get("species", ""),
        strain=config.get("funannotate2", {}).get("strain", ""),
        params_json=config.get("funannotate2", {}).get("params", ""),
        pretrained=config.get("funannotate2", {}).get("pretrained", ""),
        busco_lineage=config.get("funannotate2", {}).get("busco_lineage", "eukaryota"),
        extra_clean=config.get("funannotate2", {}).get("extra_clean", ""),
        extra_train=config.get("funannotate2", {}).get("extra_train", ""),
        extra_predict=config.get("funannotate2", {}).get("extra_predict", ""),
        extra_annotate=config.get("funannotate2", {}).get("extra_annotate", ""),
    message:
        f"--- Running the stateful Funannotate2 isolate workflow for {SAMPLE}."
    script:
        "../scripts/funannotate2.py"
