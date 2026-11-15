"""Write isolate or metagenome Annotation QA decisions."""

from pathlib import Path

from common.empty_outputs import fasta_paths
from common.identifiers import validate_annotation_identifier
from common.qa import (
    fasta_stem,
    parse_eukcc_missing_markers,
    parse_checkm2_reports,
    parse_eukcc_reports,
    qa_bypass_row,
    selection_row,
    write_selection,
)


def report_for_domain(domain):
    if domain == "prok":
        return parse_checkm2_reports(snakemake.input.checkm2)
    return parse_eukcc_reports(snakemake.input.eukcc)


def missing_markers_for_domain(domain):
    if domain == "prok":
        return set()
    return parse_eukcc_missing_markers(snakemake.input.eukcc_missing)


def metagenome_fastas(action):
    domains = []
    qa_label = "QA-disabled " if action == "metagenome_bypass" else "QA "
    for domain, domain_inputs in (
        ("prok", snakemake.input.prok_bins),
        ("euk", snakemake.input.euk_bins),
    ):
        if not domain_inputs:
            continue
        directory = domain_inputs[0]
        seen = set()
        genomes = []
        for fasta in fasta_paths(Path(directory)):
            sample = validate_annotation_identifier(fasta_stem(fasta))
            if sample in seen:
                raise ValueError(
                    f"Metagenome {qa_label}sample identifier is duplicated in "
                    f"{directory}: {sample}"
                )
            seen.add(sample)
            genomes.append((sample, fasta))
        domains.append((domain, genomes))
    return domains


action = str(snakemake.params.action)
if action not in {"isolate_bypass", "metagenome_bypass", "isolate", "metagenome"}:
    raise ValueError(f"Unsupported QA selection action: {action}")

rows = []
if action in {"isolate_bypass", "isolate"}:
    domain = str(snakemake.params.domain)
    if action == "isolate_bypass":
        row = qa_bypass_row(
            source=snakemake.params.path,
            sample=snakemake.params.sample,
            domain=domain,
            is_isolate=True,
        )
    else:
        row = selection_row(
            source=snakemake.params.path,
            sample=snakemake.params.sample,
            domain=domain,
            is_isolate=True,
            qa=report_for_domain(domain),
            qa_filter=snakemake.params.qa_filter,
            missing_marker_samples=missing_markers_for_domain(domain),
        )
    rows.append(row)
else:
    for domain, genomes in metagenome_fastas(action):
        if action == "metagenome_bypass":
            builder = lambda sample, fasta: qa_bypass_row(
                source=fasta, sample=sample, domain=domain, is_isolate=False
            )
        else:
            qa = report_for_domain(domain)
            missing_marker_samples = missing_markers_for_domain(domain)
            builder = lambda sample, fasta: selection_row(
                source=fasta,
                sample=sample,
                domain=domain,
                is_isolate=False,
                qa=qa,
                qa_filter=snakemake.params.qa_filter,
                missing_marker_samples=missing_marker_samples,
            )
        rows.extend(builder(sample, fasta) for sample, fasta in genomes)

write_selection(rows, snakemake.output.selection)
Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.log[0]).write_text(
    f"selection_rows={len(rows)}\naction={action}\n",
    encoding="utf-8",
)
