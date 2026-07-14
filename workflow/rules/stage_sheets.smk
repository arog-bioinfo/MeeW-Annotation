# ----------------------------------------------------- #
# Annotation workflow stage sheets                      #
# ----------------------------------------------------- #


rule annotation_to_metabolic_modeling_stage_sheet:
    input:
        proteins=annotation_to_metabolic_modeling_proteins,
    output:
        "<results>/stage_sheets/annotation_to_metabolic_modeling.tsv",
    log:
        "<results>/stage_sheets/annotation_to_metabolic_modeling.log",
    container:
        "docker://python:3.11-slim"
    params:
        mags=lambda wildcards, input: prokaryotic_samples(),
    script:
        "../scripts/write_metabolic_stage_sheet.py"
