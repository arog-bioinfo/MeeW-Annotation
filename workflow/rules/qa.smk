# ----------------------------------------------------- #
# Annotation-owned genome quality assessment             #
# ----------------------------------------------------- #


if QA_ENABLED:

    if HAS_PROK_INPUT:

        rule qa_checkm2_prok:
            input:
                bins=prok_qa_source,
            output:
                report="<results>/qa/prok/checkm2/quality_report.tsv",
            log:
                "<results>/logs/qa/checkm2.log",
            conda:
                "../envs/checkm2.yaml"
            threads: config["qa"].get("threads", 8)
            params:
                database_path=lambda wc, output: checked_cache_path(
                    config["qa"]["checkm2"]["database_path"],
                    Path(output.report).parents[3],
                    "CheckM2",
                ),
                download_database=config["qa"]["checkm2"].get(
                    "download_database", False
                ),
                extension=config["qa"]["checkm2"]["extension"],
                extra=config["qa"]["checkm2"].get("extra", ""),
            message:
                f"--- Running Annotation-owned CheckM2 for {SAMPLE}."
            script:
                "../scripts/checkm2.py"

    elif IS_METAGENOME:

        rule qa_checkm2_prok:
            output:
                report="<results>/qa/prok/checkm2/quality_report.tsv",
            log:
                "<results>/logs/qa/checkm2.log",
            conda:
                "../envs/python.yaml"
            params:
                action="checkm2",
                domain="prok",
            message:
                f"--- Writing absent prokaryotic QA contracts for {SAMPLE}."
            script:
                "../scripts/write_empty_qa_outputs.py"

    if HAS_EUK_INPUT:

        rule qa_eukcc_euk:
            input:
                bins=euk_qa_source,
            output:
                report="<results>/qa/euk/eukcc/eukcc.csv",
                bad_quality="<results>/qa/euk/eukcc/bad_quality.csv",
                missing_markers="<results>/qa/euk/eukcc/missing_marker_genes.txt",
            log:
                "<results>/logs/qa/eukcc.log",
            conda:
                "../envs/eukcc.yaml"
            threads: config["qa"].get("threads", 8)
            params:
                database_path=lambda wc, output: checked_cache_path(
                    config["qa"]["eukcc"].get("database_path", ""),
                    Path(output.report).parents[3],
                    "EukCC",
                ),
                download_database=config["qa"]["eukcc"].get("download_database", False),
                extra=config["qa"]["eukcc"].get("extra", ""),
            message:
                f"--- Running Annotation-owned EukCC for {SAMPLE}."
            script:
                "../scripts/eukcc.py"

    elif IS_METAGENOME:

        rule qa_eukcc_euk:
            output:
                report="<results>/qa/euk/eukcc/eukcc.csv",
                bad_quality="<results>/qa/euk/eukcc/bad_quality.csv",
                missing_markers="<results>/qa/euk/eukcc/missing_marker_genes.txt",
            log:
                "<results>/logs/qa/eukcc.log",
            conda:
                "../envs/python.yaml"
            params:
                action="eukcc",
                domain="euk",
            message:
                f"--- Writing absent eukaryotic QA contracts for {SAMPLE}."
            script:
                "../scripts/write_empty_qa_outputs.py"

    if IS_METAGENOME:

        checkpoint qa_selection:
            input:
                prok_bins=lambda wc: (
                    [INPUT_RECORD["prok_bins_dir"]]
                    if INPUT_RECORD.get("prok_bins_dir")
                    else []
                ),
                euk_bins=lambda wc: (
                    [INPUT_RECORD["euk_bins_dir"]]
                    if INPUT_RECORD.get("euk_bins_dir")
                    else []
                ),
                checkm2=rules.qa_checkm2_prok.output.report,
                eukcc=[
                    rules.qa_eukcc_euk.output.report,
                    rules.qa_eukcc_euk.output.bad_quality,
                ],
                eukcc_missing=rules.qa_eukcc_euk.output.missing_markers,
            output:
                selection="<results>/qa/selection.tsv",
            log:
                "<results>/logs/qa/selection.log",
            conda:
                "../envs/python.yaml"
            params:
                action="metagenome",
                qa_filter=config["qa_filter"],
            message:
                "--- Selecting metagenome MAGs after Annotation-owned QA."
            script:
                "../scripts/write_qa_selection.py"

    else:

        checkpoint qa_selection:
            input:
                checkm2=(
                    rules.qa_checkm2_prok.output.report
                    if ISOLATE_DOMAIN == "prok"
                    else []
                ),
                eukcc=(
                    [
                        rules.qa_eukcc_euk.output.report,
                        rules.qa_eukcc_euk.output.bad_quality,
                    ]
                    if ISOLATE_DOMAIN == "euk"
                    else []
                ),
                eukcc_missing=(
                    rules.qa_eukcc_euk.output.missing_markers
                    if ISOLATE_DOMAIN == "euk"
                    else []
                ),
            output:
                selection="<results>/qa/selection.tsv",
            log:
                "<results>/logs/qa/selection.log",
            conda:
                "../envs/python.yaml"
            params:
                action="isolate",
                sample=SAMPLE,
                path=INPUT_RECORD["path"],
                domain=ISOLATE_DOMAIN,
                qa_filter=config["qa_filter"],
            message:
                f"--- Recording the Annotation QA decision for {SAMPLE}."
            script:
                "../scripts/write_qa_selection.py"

else:

    if IS_METAGENOME:

        checkpoint qa_selection:
            input:
                prok_bins=lambda wc: (
                    [INPUT_RECORD["prok_bins_dir"]]
                    if INPUT_RECORD.get("prok_bins_dir")
                    else []
                ),
                euk_bins=lambda wc: (
                    [INPUT_RECORD["euk_bins_dir"]]
                    if INPUT_RECORD.get("euk_bins_dir")
                    else []
                ),
            output:
                selection="<results>/qa/selection.tsv",
            log:
                "<results>/logs/qa/selection.log",
            conda:
                "../envs/python.yaml"
            params:
                action="metagenome_bypass",
            message:
                "--- Selecting every metagenome bin because Annotation QA is disabled."
            script:
                "../scripts/write_qa_selection.py"

    else:

        checkpoint qa_selection:
            input:
                fasta=sample_fasta,
            output:
                selection="<results>/qa/selection.tsv",
            log:
                "<results>/logs/qa/selection.log",
            conda:
                "../envs/python.yaml"
            params:
                action="isolate_bypass",
                sample=SAMPLE,
                path=INPUT_RECORD["path"],
                domain=ISOLATE_DOMAIN,
            message:
                f"--- Recording the QA-disabled Annotation selection for {SAMPLE}."
            script:
                "../scripts/write_qa_selection.py"
