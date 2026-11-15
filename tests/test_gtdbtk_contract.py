import csv
import importlib
import os
import runpy
from pathlib import Path
from types import SimpleNamespace

import pytest

from common.gtdbtk import (
    consolidate_gtdbtk_summaries,
    temporary_gtdbtk_genomes,
)
from helpers import file_snapshot


def write_summary(path, rows, header=("user_genome", "classification")):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def test_temporary_genomes_use_unique_symlinks_and_preserve_sources(tmp_path):
    source_root = tmp_path / "source genomes"
    source_root.mkdir()
    sources = [source_root / "first unrelated.fa", source_root / "second.fna"]
    for index, source in enumerate(sources, start=1):
        source.write_text(f">contig_{index}\nACGT\n", encoding="utf-8")
        os.utime(source, ns=(1_700_000_000_000_000_000 + index,) * 2)
    source_before = file_snapshot(source_root)

    scratch_root = tmp_path / "scratch-root"
    other_job = scratch_root / "other-job"
    other_job.mkdir(parents=True)
    scratch_sentinel = other_job / "keep.txt"
    scratch_sentinel.write_text("unrelated scratch material\n", encoding="utf-8")
    scratch_before = file_snapshot(scratch_root)
    pairs = [("sample_a", str(sources[0])), ("sample_b", str(sources[1]))]

    with temporary_gtdbtk_genomes(pairs, sources, scratch_root) as first:
        assert first.parent == scratch_root.resolve()
        assert sorted(path.name for path in first.iterdir()) == [
            "sample_a.fasta",
            "sample_b.fasta",
        ]
        assert all(path.is_symlink() for path in first.iterdir())
        assert (first / "sample_a.fasta").resolve() == sources[0].resolve()
        assert (first / "sample_b.fasta").resolve() == sources[1].resolve()

        with temporary_gtdbtk_genomes(pairs, sources, scratch_root) as second:
            assert second.parent == scratch_root.resolve()
            assert second != first
            assert first.exists()

    assert file_snapshot(source_root) == source_before
    assert file_snapshot(scratch_root) == scratch_before
    assert sorted(scratch_root.iterdir()) == [other_job]


def test_temporary_genomes_clean_owned_scratch_after_failure(tmp_path):
    source = tmp_path / "source.fasta"
    source.write_text(">contig\nACGT\n", encoding="utf-8")
    scratch_root = tmp_path / "scratch"
    scratch_root.mkdir()
    sentinel = scratch_root / "keep.txt"
    sentinel.write_text("keep\n", encoding="utf-8")
    before = file_snapshot(tmp_path)

    with pytest.raises(RuntimeError, match="fake GTDB-Tk failure"):
        with temporary_gtdbtk_genomes(
            [("sample", str(source))], [source], scratch_root
        ) as genomes:
            (genomes / "tool-side-effect.txt").write_text("partial\n", encoding="utf-8")
            raise RuntimeError("fake GTDB-Tk failure")

    assert file_snapshot(tmp_path) == before
    assert sorted(scratch_root.iterdir()) == [sentinel]


def test_temporary_genomes_reject_invalid_metadata_before_creating_scratch(tmp_path):
    source = tmp_path / "source.fa"
    source.write_text(">contig\nACGT\n", encoding="utf-8")
    other = tmp_path / "other.fna"
    other.write_text(">other\nACGT\n", encoding="utf-8")
    scratch_root = tmp_path / "unused-scratch"

    invalid_cases = [
        ([], [], "at least one"),
        ([("", str(source))], [source], "Invalid Annotation identifier"),
        (
            [("duplicate", str(source)), ("duplicate", str(other))],
            [source, other],
            "duplicate genome identifiers",
        ),
        (
            [("sample", str(source))],
            [other],
            "do not match declared FASTA inputs",
        ),
        (
            [("sample", str(tmp_path / "missing.fa"))],
            [tmp_path / "missing.fa"],
            "not a usable regular file",
        ),
        (
            [("sample", str(tmp_path))],
            [tmp_path],
            "not a usable regular file",
        ),
    ]
    for pairs, declared, message in invalid_cases:
        with pytest.raises((ValueError, FileNotFoundError), match=message):
            with temporary_gtdbtk_genomes(pairs, declared, scratch_root):
                pass
        assert not scratch_root.exists()


def test_adapter_rejects_spaced_scratch_before_cache_or_classifier(
    tmp_path, monkeypatch
):
    source_root = tmp_path / "source inputs"
    source_root.mkdir()
    source = source_root / "unrelated genome.fna"
    source.write_text(">contig\nACGT\n", encoding="utf-8")
    source_before = file_snapshot(source_root)

    scratch_root = tmp_path / "scratch parent"
    scratch_root.mkdir()
    sentinel = scratch_root / "keep.txt"
    sentinel.write_text("unrelated scratch material\n", encoding="utf-8")
    scratch_before = file_snapshot(scratch_root)
    database = tmp_path / "must-not-prepare-gtdbtk-db"
    database_calls = []
    classifier_calls = []

    database_module = importlib.import_module("common.gtdbtk_database")
    shell_module = importlib.import_module("snakemake.shell")
    monkeypatch.setattr(
        database_module,
        "ensure_gtdbtk_database",
        lambda *args: database_calls.append(args),
    )
    monkeypatch.setattr(
        shell_module,
        "shell",
        lambda *args: classifier_calls.append(args),
    )

    class Params(dict):
        __getattr__ = dict.__getitem__

    snakemake = SimpleNamespace(
        input=SimpleNamespace(fastas=[str(source)]),
        output=SimpleNamespace(classification=str(tmp_path / "classification.tsv")),
        params=Params(
            db_path=str(database),
            extension="fasta",
            extra="",
            pairs=[("sample", str(source))],
            expected_ids=["sample"],
        ),
        resources=SimpleNamespace(tmpdir=str(scratch_root)),
        threads=1,
        log=[str(tmp_path / "gtdbtk.log")],
    )

    with pytest.raises(
        ValueError,
        match=r"resources\.tmpdir resolved to .*scratch parent.*without spaces",
    ):
        runpy.run_path(
            str(
                Path(__file__).resolve().parents[1]
                / "workflow"
                / "scripts"
                / "gtdb-tk.py"
            ),
            init_globals={"snakemake": snakemake},
        )

    assert database_calls == []
    assert classifier_calls == []
    assert not database.exists()
    assert file_snapshot(source_root) == source_before
    assert file_snapshot(scratch_root) == scratch_before
    assert sorted(scratch_root.iterdir()) == [sentinel]


def test_consolidates_bacterial_and_archaeal_summaries(tmp_path):
    bacteria = tmp_path / "gtdbtk.bac120.summary.tsv"
    archaea = tmp_path / "gtdbtk.ar53.summary.tsv"
    output = tmp_path / "classification.tsv"
    write_summary(bacteria, [("mag_b", "d__Bacteria")])
    write_summary(archaea, [("mag_a", "d__Archaea")])

    consolidate_gtdbtk_summaries(
        [("bacteria", bacteria), ("archaea", archaea)],
        output,
        ["mag_a", "mag_b"],
    )

    assert output.read_text(encoding="utf-8") == (
        "domain\tuser_genome\tclassification\n"
        "bacteria\tmag_b\td__Bacteria\n"
        "archaea\tmag_a\td__Archaea\n"
    )


@pytest.mark.parametrize(
    ("bacterial_rows", "archaeal_rows", "expected", "message"),
    [
        ([("mag_a", "d__Bacteria")], [], ["mag_a", "mag_b"], "Missing"),
        (
            [("mag_a", "d__Bacteria"), ("extra", "d__Bacteria")],
            [],
            ["mag_a"],
            "unexpected",
        ),
        (
            [("mag_a", "d__Bacteria")],
            [("mag_a", "d__Archaea")],
            ["mag_a"],
            "duplicated",
        ),
    ],
)
def test_rejects_identifier_mismatches(
    tmp_path, bacterial_rows, archaeal_rows, expected, message
):
    bacteria = tmp_path / "gtdbtk.bac120.summary.tsv"
    archaea = tmp_path / "gtdbtk.ar53.summary.tsv"
    write_summary(bacteria, bacterial_rows)
    write_summary(archaea, archaeal_rows)

    with pytest.raises(ValueError, match=message):
        consolidate_gtdbtk_summaries(
            [("bacteria", bacteria), ("archaea", archaea)],
            tmp_path / "classification.tsv",
            expected,
        )


@pytest.mark.parametrize(
    ("header", "rows", "message"),
    [
        (("classification",), [("d__Bacteria",)], "user_genome"),
        (
            ("user_genome", "user_genome"),
            [("mag_a", "mag_a")],
            "duplicate columns",
        ),
        (("user_genome", "classification"), [("mag_a",)], "malformed row"),
        (("user_genome", "classification"), [("", "d__Bacteria")], "blank"),
    ],
)
def test_rejects_malformed_summary_contract(tmp_path, header, rows, message):
    bacteria = tmp_path / "gtdbtk.bac120.summary.tsv"
    write_summary(bacteria, rows, header=header)

    with pytest.raises(ValueError, match=message):
        consolidate_gtdbtk_summaries(
            [("bacteria", bacteria)],
            tmp_path / "classification.tsv",
            ["mag_a"],
        )


def test_rejects_different_bacterial_and_archaeal_columns(tmp_path):
    bacteria = tmp_path / "gtdbtk.bac120.summary.tsv"
    archaea = tmp_path / "gtdbtk.ar53.summary.tsv"
    write_summary(bacteria, [("mag_b", "d__Bacteria")])
    write_summary(
        archaea,
        [("mag_a", "d__Archaea", "warning")],
        header=("user_genome", "classification", "warnings"),
    )

    with pytest.raises(ValueError, match="columns differ"):
        consolidate_gtdbtk_summaries(
            [("bacteria", bacteria), ("archaea", archaea)],
            tmp_path / "classification.tsv",
            ["mag_a", "mag_b"],
        )


def test_project_classification_inventory_row_keeps_four_column_schema(tmp_path):
    root = tmp_path / "results" / "project"
    selection = root / "qa" / "selection.tsv"
    selection.parent.mkdir(parents=True)
    with selection.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "sample",
                "path",
                "domain",
                "completeness",
                "contamination",
                "min_completeness",
                "max_contamination",
                "decision",
                "reason",
            ]
        )
        writer.writerow(
            [
                "mag_a",
                "/input/mag_a.fasta",
                "prok",
                "95",
                "1",
                "50",
                "10",
                "pass",
                "test",
            ]
        )

    materials = [
        root / "batch" / "prok" / "mag_a" / "bakta",
        root / "batch" / "prok" / "mag_a" / "recognizer",
        root / "batch" / "prok" / "mag_a" / "upimapi",
        root / "batch" / "prok" / "mag_a" / "protein_annotations.tsv",
    ]
    for material in materials:
        if material.suffix:
            material.parent.mkdir(parents=True, exist_ok=True)
            material.write_text("test\n", encoding="utf-8")
        else:
            material.mkdir(parents=True, exist_ok=True)
    classification = root / "gtdbtk" / "classification.tsv"
    classification.parent.mkdir(parents=True)
    classification.write_text("domain\tuser_genome\n", encoding="utf-8")
    output = root / "batch" / "annotation_results.tsv"
    log = root / "batch" / "logs" / "annotation_inventory.log"

    snakemake = SimpleNamespace(
        input=SimpleNamespace(
            selection=str(selection),
            materials=[str(path) for path in materials],
            gtdbtk=[str(classification)],
        ),
        output=SimpleNamespace(inventory=str(output)),
        params=SimpleNamespace(
            results_root=str(root),
            project="project",
            gtdbtk_enabled=True,
        ),
        log=[str(log)],
    )
    runpy.run_path(
        str(
            Path(__file__).resolve().parents[1]
            / "workflow"
            / "scripts"
            / "write_directory_inventories.py"
        ),
        init_globals={"snakemake": snakemake},
    )

    with output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    assert rows[0] == ["bin", "domain", "output_type", "path"]
    assert rows[-1] == [
        "project",
        "prok",
        "gtdbtk_classification",
        str(classification.resolve()),
    ]
