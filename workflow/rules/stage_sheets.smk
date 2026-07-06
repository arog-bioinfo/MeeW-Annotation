# ----------------------------------------------------- #
# Annotation workflow stage sheets                      #
# ----------------------------------------------------- #


def annotation_to_metabolic_modeling_proteins(wildcards):
    return [
        str(rules.bakta.output.faa).format(sample=sample)
        for sample in prokaryotic_samples()
    ]


rule annotation_to_metabolic_modeling_stage_sheet:
    input:
        proteins=annotation_to_metabolic_modeling_proteins,
    output:
        "<results>/stage_sheets/annotation_to_metabolic_modeling.tsv",
    params:
        mags=lambda wildcards, input: prokaryotic_samples(),
    run:
        output_path = Path(output[0])
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", newline="") as handle:
            handle.write("mag\tpath\n")
            for mag, protein_fasta in zip(params.mags, input.proteins):
                handle.write(f"{mag}\t{Path(protein_fasta).resolve()}\n")
