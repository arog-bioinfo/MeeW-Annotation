# ----------------------------------------------------- #
# Annotation workflow stage sheets                      #
# ----------------------------------------------------- #


if DEFERRED_SAMPLE_SHEET:

    checkpoint deferred_annotation_sample_sheet:
        input:
            config["sample_sheet"],
        output:
            "<results>/stage_sheets/deferred_annotation_samples.tsv",
        message:
            "--- Materializing deferred Annotation sample sheet."
        shell:
            "cp {input} {output}"

    def annotation_to_metabolic_modeling_sample_sheet(wildcards):
        return checkpoints.deferred_annotation_sample_sheet.get().output[0]

    def annotation_to_metabolic_modeling_samples(wildcards):
        return filtered_samples_from_path(
            annotation_to_metabolic_modeling_sample_sheet(wildcards)
        )

    def annotation_to_metabolic_modeling_prok_samples(wildcards):
        samples = annotation_to_metabolic_modeling_samples(wildcards)
        return samples[samples["domain"] == "prok"].index.tolist()

    def annotation_to_metabolic_modeling_proteins(wildcards):
        return [
            f"<results>/bakta/{sample}/{sample}.faa"
            for sample in annotation_to_metabolic_modeling_prok_samples(wildcards)
        ]

    rule annotation_to_metabolic_modeling_stage_sheet:
        input:
            sample_sheet=annotation_to_metabolic_modeling_sample_sheet,
            proteins=annotation_to_metabolic_modeling_proteins,
        output:
            "<results>/stage_sheets/annotation_to_metabolic_modeling.tsv",
        params:
            mags=lambda wildcards, input: annotation_to_metabolic_modeling_prok_samples(
                wildcards
            ),
        run:
            output_path = Path(output[0])
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with output_path.open("w", newline="") as handle:
                handle.write("mag\tpath\n")
                for mag, protein_fasta in zip(params.mags, input.proteins):
                    handle.write(f"{mag}\t{Path(protein_fasta).resolve()}\n")

else:

    rule annotation_to_metabolic_modeling_stage_sheet:
        input:
            proteins=expand(
                "<results>/bakta/{sample}/{sample}.faa", sample=prok_samples
            ),
        output:
            "<results>/stage_sheets/annotation_to_metabolic_modeling.tsv",
        params:
            mags=prok_samples,
        run:
            output_path = Path(output[0])
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with output_path.open("w", newline="") as handle:
                handle.write("mag\tpath\n")
                for mag, protein_fasta in zip(params.mags, input.proteins):
                    handle.write(f"{mag}\t{Path(protein_fasta).resolve()}\n")
