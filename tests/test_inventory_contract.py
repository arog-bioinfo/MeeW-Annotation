import csv
import runpy
from pathlib import Path
from types import SimpleNamespace

import pytest

from common import qa as qa_module
from common.qa import qa_bypass_row, write_selection

REPOSITORY = Path(__file__).resolve().parents[1]
INVENTORY_SCRIPT = (
    REPOSITORY / "workflow" / "scripts" / "write_directory_inventories.py"
)


def selection_row(tmp_path, sample, domain, decision):
    row = qa_bypass_row(
        source=tmp_path / "inputs" / f"{sample}.fa",
        sample=sample,
        domain=domain,
        is_isolate=False,
    )
    row["decision"] = decision
    return row


def run_inventory(
    tmp_path,
    rows,
    materials,
    *,
    gtdbtk=(),
    enabled=False,
):
    selection = tmp_path / "selection.tsv"
    write_selection(rows, selection)
    output = tmp_path / "annotation_results.tsv"
    snakemake = SimpleNamespace(
        input=SimpleNamespace(
            selection=str(selection),
            materials=list(materials),
            gtdbtk=list(gtdbtk),
        ),
        output=SimpleNamespace(inventory=str(output)),
        params=SimpleNamespace(
            results_root=str(tmp_path / "results" / "project"),
            project="project",
            gtdbtk_enabled=enabled,
        ),
        log=[str(tmp_path / "inventory.log")],
    )
    runpy.run_path(str(INVENTORY_SCRIPT), init_globals={"snakemake": snakemake})
    return output


def expected_materials(root, sample, domain):
    if domain == "euk":
        return [
            ("metaeuk_results", root / "batch" / "euk" / sample / "metaeuk"),
            ("recognizer_results", root / "batch" / "euk" / sample / "recognizer"),
        ]
    return [
        ("bakta_results", root / "batch" / "prok" / sample / "bakta"),
        ("recognizer_results", root / "batch" / "prok" / sample / "recognizer"),
        ("upimapi_results", root / "batch" / "prok" / sample / "upimapi"),
        (
            "protein_annotations",
            root / "batch" / "prok" / sample / "protein_annotations.tsv",
        ),
    ]


def create_material(path):
    if path.suffix:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixture\n", encoding="utf-8")
    else:
        path.mkdir(parents=True, exist_ok=True)


def test_annotation_inventory_has_exact_schema_inventory_order_and_single_read(
    tmp_path, monkeypatch
):
    root = tmp_path / "results" / "project"
    rows = [
        selection_row(tmp_path, "z_prok", "prok", "pass"),
        selection_row(tmp_path, "rejected", "prok", "fail"),
        selection_row(tmp_path, "a_euk", "euk", "pass"),
    ]
    expected = [
        ("a_euk", "euk", output_type, path)
        for output_type, path in expected_materials(root, "a_euk", "euk")
    ] + [
        ("z_prok", "prok", output_type, path)
        for output_type, path in expected_materials(root, "z_prok", "prok")
    ]
    for _, _, _, path in expected:
        create_material(path)

    original_read_selection = qa_module.read_selection
    reads = []

    def counted_read_selection(path):
        reads.append(Path(path))
        return original_read_selection(path)

    monkeypatch.setattr(qa_module, "read_selection", counted_read_selection)
    output = run_inventory(
        tmp_path,
        rows,
        [str(path) for _, _, _, path in reversed(expected)],
    )

    with output.open(newline="", encoding="utf-8") as handle:
        observed = list(csv.reader(handle, delimiter="\t"))
    assert observed == [
        ["bin", "domain", "output_type", "path"],
        *[
            [sample, domain, output_type, str(path.resolve())]
            for sample, domain, output_type, path in expected
        ],
    ]
    assert "rejected" not in output.read_text(encoding="utf-8")
    assert reads == [tmp_path / "selection.tsv"]
    assert (tmp_path / "inventory.log").read_text(encoding="utf-8") == (
        "Wrote inventory.\n"
    )


@pytest.mark.parametrize("failure", ["undeclared", "absent"])
def test_annotation_inventory_fails_for_missing_atomic_material(tmp_path, failure):
    root = tmp_path / "results" / "project"
    rows = [selection_row(tmp_path, "mag", "prok", "pass")]
    paths = [path for _, path in expected_materials(root, "mag", "prok")]
    for path in paths:
        create_material(path)
    declared = list(paths)
    if failure == "undeclared":
        declared.pop()
    else:
        missing = paths[-1]
        missing.unlink()

    with pytest.raises(FileNotFoundError, match="Undeclared|absent"):
        run_inventory(tmp_path, rows, declared)


def test_header_only_inventory_is_valid_when_no_bins_pass(tmp_path):
    rows = [selection_row(tmp_path, "rejected", "euk", "fail")]

    output = run_inventory(tmp_path, rows, [])

    assert output.read_text(encoding="utf-8") == "bin\tdomain\toutput_type\tpath\n"
