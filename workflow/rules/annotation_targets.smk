# ----------------------------------------------------- #
# Annotation workflow targets                           #
# ----------------------------------------------------- #


def isolate_annotation_targets():
    if ISOLATE_DOMAIN == "prok":
        targets = [
            rules.bakta.output.results,
            rules.recognizer_prok.output.results,
            rules.upimapi.output.results,
            rules.protein_summary.output.summary,
        ]
        if GTDBTK_ENABLED:
            targets.extend(rules.gtdbtk.output)
        return targets
    return list(rules.funannotate2.output)


def final_targets(wildcards):
    if IS_METAGENOME:
        return [rules.metagenome_annotation_inventory.output.inventory]

    selection, decision = isolate_selection_record()
    targets = [str(selection)]
    if decision["decision"] == "report_only":
        targets.extend(isolate_annotation_targets())
    return targets


FINAL_TARGETS = final_targets
