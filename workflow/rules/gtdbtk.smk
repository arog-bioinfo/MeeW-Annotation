# ----------------------------------------------------- #
# Optional shared prokaryotic GTDB-Tk classification     #
# ----------------------------------------------------- #


rule gtdbtk:
    input:
        qa=gtdbtk_selection_file,
        fastas=eligible_prokaryotic_gtdbtk_fastas,
        database=lambda wc: gtdbtk_database_inputs(GTDBTK_DATA_DIR),
    output:
        classification="<results>/gtdbtk/classification.tsv",
    log:
        "<results>/gtdbtk/gtdbtk.log",
    conda:
        "../envs/gtdb-tk.yaml"
    threads: config.get("threads", {}).get("high", 16)
    params:
        db_path=lambda wc, output: checked_cache_path(
            GTDBTK_DATA_DIR,
            Path(output.classification).parents[1],
            "GTDB-Tk",
        ),
        extension="fasta",
        extra=config.get("gtdbtk", {}).get("extra", ""),
        pairs=lambda wc, input: eligible_prokaryotic_gtdbtk_pairs(wc),
        expected_ids=lambda wc, input: [
            sample for sample, _ in eligible_prokaryotic_gtdbtk_pairs(wc)
        ],
    message:
        "--- Running one GTDB-Tk classification batch for eligible prokaryotic genomes."
    script:
        "../scripts/gtdb-tk.py"
