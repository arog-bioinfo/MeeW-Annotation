# ----------------------------------------------------- #
# Annotation workflow stage sheets                      #
# ----------------------------------------------------- #


rule annotation_to_metabolic_modeling_stage_sheet:
    input:
        proteins=expand("results/bakta/{sample}/{sample}.faa", sample=prok_samples),
    output:
        "results/stage_sheets/annotation_to_metabolic_modeling.tsv",
    params:
        mags=prok_samples,
    run:
        output_path = Path(output[0])
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", newline="") as handle:
            handle.write("mag\tpath\n")
            for mag, protein_fasta in zip(params.mags, input.proteins):
                handle.write(f"{mag}\t{Path(protein_fasta).resolve()}\n")
