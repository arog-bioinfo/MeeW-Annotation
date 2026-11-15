"""Write deterministic material inventories after metagenome QA selection."""

import csv
from pathlib import Path

from common.identifiers import validate_annotation_identifier
from common.qa import read_selection


def path_list(value):
    if isinstance(value, (str, Path)):
        return [value]
    return list(value)


def write_inventory():
    output = Path(snakemake.output.inventory)
    output.parent.mkdir(parents=True, exist_ok=True)
    materials = {
        str(Path(path))
        for path in [
            *path_list(snakemake.input.materials),
            *path_list(snakemake.input.gtdbtk),
        ]
    }
    root = Path(snakemake.params.results_root)
    passing_rows = sorted(
        (
            row
            for row in read_selection(snakemake.input.selection)
            if row["decision"] == "pass"
        ),
        key=lambda row: (row["domain"], row["sample"]),
    )
    expected = []
    for row in passing_rows:
        bin_id = validate_annotation_identifier(row["sample"])
        domain = row["domain"]
        if domain == "prok":
            paths = {
                "bakta_results": root / "batch" / "prok" / bin_id / "bakta",
                "recognizer_results": root / "batch" / "prok" / bin_id / "recognizer",
                "upimapi_results": root / "batch" / "prok" / bin_id / "upimapi",
                "protein_annotations": root
                / "batch"
                / "prok"
                / bin_id
                / "protein_annotations.tsv",
            }
        else:
            paths = {
                "metaeuk_results": root / "batch" / "euk" / bin_id / "metaeuk",
                "recognizer_results": root / "batch" / "euk" / bin_id / "recognizer",
            }
        for output_type, path in paths.items():
            expected.append((bin_id, domain, output_type, path))

    if snakemake.params.gtdbtk_enabled and any(
        row["domain"] == "prok" for row in passing_rows
    ):
        expected.append(
            (
                validate_annotation_identifier(str(snakemake.params.project)),
                "prok",
                "gtdbtk_classification",
                root / "gtdbtk" / "classification.tsv",
            )
        )

    missing_dependencies = [
        str(path) for _, _, _, path in expected if str(path) not in materials
    ]
    missing_paths = [str(path) for _, _, _, path in expected if not path.exists()]
    if missing_dependencies or missing_paths:
        raise FileNotFoundError(
            "Metagenome annotation inventory is missing atomic dependencies/paths. "
            f"Undeclared: {missing_dependencies}; absent: {missing_paths}."
        )

    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["bin", "domain", "output_type", "path"])
        for bin_id, domain, output_type, path in expected:
            writer.writerow([bin_id, domain, output_type, str(path.resolve())])


write_inventory()

Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.log[0]).write_text("Wrote inventory.\n", encoding="utf-8")
