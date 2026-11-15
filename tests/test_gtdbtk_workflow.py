import csv
import json
import os
import time
from pathlib import Path

from helpers import file_snapshot, rule_count


def _records(path):
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]


def _assert_no_persistent_staging(case):
    result_root = case.root / "results" / case.sample
    assert not (result_root / "gtdbtk_genomes").exists()
    assert not (result_root / "gtdbtk_genomes.log").exists()


def _assert_scratch_record(record, scratch_parent, expected_sources):
    args = record["args"]
    assert args[0] == "classify_wf"
    assert "--batchfile" not in args
    assert args[args.index("--extension") + 1] == "fasta"
    assert args[args.index("--cpus") + 1] == "1"
    assert "--force" in args

    genome_dir = Path(record["genome_dir"])
    assert genome_dir.parent == scratch_parent.resolve()
    assert args[args.index("--genome_dir") + 1] == str(genome_dir)
    assert not genome_dir.exists()
    assert record["genomes"] == [
        {
            "name": f"{sample}.fasta",
            "is_symlink": True,
            "target": str(source.resolve()),
            "content": source.read_text(encoding="utf-8"),
        }
        for sample, source in expected_sources
    ]


def test_isolate_gtdbtk_enabled_preserves_route_and_output(
    workflow_case, fake_gtdbtk, offline_guard
):
    case = workflow_case(
        "isolate",
        enabled=True,
        fasta_names={"isolate1": "assembled genome source.fna"},
    )
    output = case.run(dry_run=True, extra_arguments=("--rerun-triggers", "mtime"))

    assert "stage_gtdbtk_genomes" not in output
    assert rule_count(output, "gtdbtk") == 1
    assert "results/isolate1/gtdbtk/classification.tsv" in output
    assert "gtdbtk_genomes" not in output

    source_before = file_snapshot(case.root / "input")
    scratch_parent = case.root / "scratch-parent"
    scratch_parent.mkdir()
    scratch_sentinel = scratch_parent / "keep.txt"
    scratch_sentinel.write_text("unrelated scratch material\n", encoding="utf-8")
    scratch_before = file_snapshot(scratch_parent)

    env, calls, records = fake_gtdbtk
    scratch_arguments = (
        "--default-resources",
        f"tmpdir={str(scratch_parent)!r}",
    )
    case.run(
        dry_run=False,
        targets=("results/isolate1/gtdbtk/classification.tsv",),
        env=env,
        extra_arguments=scratch_arguments,
    )
    assert calls.read_text(encoding="utf-8").splitlines() == ["classify_wf"]
    _assert_scratch_record(
        _records(records)[0],
        scratch_parent,
        [("isolate1", case.fastas["isolate1"])],
    )
    with case.classification.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows == [
        {
            "domain": "bacteria",
            "user_genome": "isolate1",
            "classification": "d__Bacteria",
        }
    ]
    assert file_snapshot(case.root / "input") == source_before
    assert file_snapshot(scratch_parent) == scratch_before
    _assert_no_persistent_staging(case)

    database = case.root / "gtdbtk_db"
    metadata = database / "metadata" / "metadata.txt"
    metadata.write_text("VERSION_DATA=r232\nfixture-change=true\n", encoding="utf-8")
    current = metadata.stat()
    metadata_time = time.time_ns()
    os.utime(
        metadata,
        ns=(current.st_atime_ns, metadata_time),
    )
    changed_cache = file_snapshot(database)
    case.run(
        dry_run=False,
        targets=("results/isolate1/gtdbtk/classification.tsv",),
        env=env,
        extra_arguments=scratch_arguments,
    )
    assert calls.read_text(encoding="utf-8").splitlines() == [
        "classify_wf",
        "classify_wf",
    ]
    assert file_snapshot(database) == changed_cache
    assert file_snapshot(case.root / "input") == source_before
    assert file_snapshot(scratch_parent) == scratch_before
    assert len(_records(records)) == 2
    for record in _records(records):
        _assert_scratch_record(
            record,
            scratch_parent,
            [("isolate1", case.fastas["isolate1"])],
        )
    _assert_no_persistent_staging(case)
    assert not offline_guard["attempt_log"].exists()


def test_isolate_gtdbtk_disabled_has_no_annotation_jobs(workflow_case):
    case = workflow_case("isolate", enabled=False)
    output = case.run(dry_run=True)

    assert "stage_gtdbtk_genomes" not in output
    assert rule_count(output, "gtdbtk") == 0
    assert "gtdbtk_genomes" not in output
    _assert_no_persistent_staging(case)


def test_metagenome_one_batch_stages_only_passing_prokaryotes(
    workflow_case, fake_gtdbtk
):
    case = workflow_case(
        "metagenome",
        enabled=True,
        rows=[
            ("pass_prok_a", "prok", "pass"),
            ("pass_prok_b", "prok", "pass"),
            ("fail_prok", "prok", "fail"),
            ("pass_euk", "euk", "pass"),
        ],
        fasta_names={
            "pass_prok_a": "pass_prok_a.fa",
            "pass_prok_b": "pass_prok_b.fna",
            "fail_prok": "failed.fasta",
            "pass_euk": "pass_euk.fa",
        },
    )
    dry_run = case.run(dry_run=True)
    assert "stage_gtdbtk_genomes" not in dry_run
    assert rule_count(dry_run, "gtdbtk") == 1
    assert "gtdbtk_genomes" not in dry_run

    source_before = file_snapshot(case.root / "input")
    scratch_parent = case.root / "scratch-parent"
    scratch_parent.mkdir()
    scratch_sentinel = scratch_parent / "keep.txt"
    scratch_sentinel.write_text("unrelated scratch material\n", encoding="utf-8")
    scratch_before = file_snapshot(scratch_parent)

    env, calls, records = fake_gtdbtk
    case.run(
        dry_run=False,
        targets=("results/project/gtdbtk/classification.tsv",),
        env=env,
        extra_arguments=(
            "--default-resources",
            f"tmpdir={str(scratch_parent)!r}",
        ),
    )

    assert calls.read_text(encoding="utf-8").splitlines() == ["classify_wf"]
    _assert_scratch_record(
        _records(records)[0],
        scratch_parent,
        [
            ("pass_prok_a", case.fastas["pass_prok_a"]),
            ("pass_prok_b", case.fastas["pass_prok_b"]),
        ],
    )
    with case.classification.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows == [
        {
            "domain": "bacteria",
            "user_genome": "pass_prok_a",
            "classification": "d__Bacteria",
        },
        {
            "domain": "bacteria",
            "user_genome": "pass_prok_b",
            "classification": "d__Bacteria",
        },
    ]
    assert file_snapshot(case.root / "input") == source_before
    assert file_snapshot(scratch_parent) == scratch_before
    _assert_no_persistent_staging(case)


def test_metagenome_gtdbtk_disabled_has_no_annotation_jobs(workflow_case):
    case = workflow_case(
        "metagenome",
        enabled=False,
        rows=[("pass_prok", "prok", "pass")],
    )
    output = case.run(dry_run=True)

    assert "stage_gtdbtk_genomes" not in output
    assert rule_count(output, "gtdbtk") == 0
    assert "gtdbtk_genomes" not in output
    _assert_no_persistent_staging(case)


def test_no_passing_prokaryotes_omit_gtdbtk_and_complete_final_targets(
    workflow_case,
):
    case = workflow_case(
        "metagenome",
        enabled=True,
        rows=[
            ("fail_prok", "prok", "fail"),
            ("fail_euk", "euk", "fail"),
        ],
    )
    dry_run = case.run(dry_run=True)
    assert "stage_gtdbtk_genomes" not in dry_run
    assert rule_count(dry_run, "gtdbtk") == 0
    assert "gtdbtk_genomes" not in dry_run

    case.run(dry_run=False)

    inventory = case.root / "results" / "project" / "batch" / "annotation_results.tsv"
    assert inventory.read_text(encoding="utf-8") == ("bin\tdomain\toutput_type\tpath\n")
    assert not case.classification.exists()
    _assert_no_persistent_staging(case)


def test_gtdbtk_tool_failure_cleans_only_owned_scratch(workflow_case, fake_gtdbtk):
    case = workflow_case(
        "metagenome",
        enabled=True,
        rows=[("pass_prok", "prok", "pass")],
        fasta_names={"pass_prok": "pass_prok.fasta"},
    )
    source_before = file_snapshot(case.root / "input")
    scratch_parent = case.root / "scratch-parent"
    other_job = scratch_parent / "other-job"
    other_job.mkdir(parents=True)
    scratch_sentinel = other_job / "keep.txt"
    scratch_sentinel.write_text("unrelated scratch material\n", encoding="utf-8")
    scratch_before = file_snapshot(scratch_parent)

    env, calls, records = fake_gtdbtk
    env["FAKE_GTDBTK_FAIL"] = "1"
    result = case.run_result(
        dry_run=False,
        targets=("results/project/gtdbtk/classification.tsv",),
        env=env,
        extra_arguments=(
            "--default-resources",
            f"tmpdir={str(scratch_parent)!r}",
        ),
    )

    assert result.returncode != 0
    assert calls.read_text(encoding="utf-8").splitlines() == ["classify_wf"]
    _assert_scratch_record(
        _records(records)[0],
        scratch_parent,
        [("pass_prok", case.fastas["pass_prok"])],
    )
    assert file_snapshot(case.root / "input") == source_before
    assert file_snapshot(scratch_parent) == scratch_before
    assert sorted(scratch_parent.iterdir()) == [other_job]
    assert not case.classification.exists()
    _assert_no_persistent_staging(case)
