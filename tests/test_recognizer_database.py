import json
import multiprocessing
import os
import runpy
import stat
import subprocess
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

from common.cache_runtime import CacheIntegrityError
from common.recognizer_database import (
    RECOGNIZER_COG2EC_SOURCES,
    RECOGNIZER_COG2KO_SOURCES,
    RECOGNIZER_COG_MAP_SOURCES,
    RECOGNIZER_COMMON_RELATIONAL_FILES,
    RECOGNIZER_COMPLETION_INDICATOR,
    RECOGNIZER_DATABASE_RELATIONAL_FILES,
    RECOGNIZER_INDEX_EXTENSIONS,
    RECOGNIZER_PREFIXES,
    assess_recognizer_database,
)
from helpers import material_snapshot

ROOT = Path(__file__).resolve().parents[1]
WRAPPER = ROOT / "workflow" / "scripts" / "recognizer.py"


class Values(dict):
    __getattr__ = dict.__getitem__


class FakeSnakemake(SimpleNamespace):
    pass


def seed_native(bin_dir):
    bin_dir.mkdir(parents=True)
    executable = bin_dir / "recognizer"
    executable.write_text(
        "#!/usr/bin/env python3\n"
        "import json, os, sys\n"
        "from pathlib import Path\n"
        "record = {'argv': sys.argv[1:], 'cwd': os.getcwd()}\n"
        "with Path(os.environ['FAKE_RECOGNIZER_RECORD']).open('a') as handle:\n"
        "    handle.write(json.dumps(record) + '\\n')\n"
        "print('fake recognizer output')\n"
        "if os.environ.get('FAKE_RECOGNIZER_PARTIAL') == '1':\n"
        "    root = Path(sys.argv[sys.argv.index('-rd') + 1])\n"
        "    root.mkdir(parents=True, exist_ok=True)\n"
        "    (root / 'partial.native').write_text('partial\\n')\n"
        "if os.environ.get('FAKE_RECOGNIZER_FAIL'):\n"
        "    raise SystemExit(int(os.environ['FAKE_RECOGNIZER_FAIL']))\n"
        "out = Path(sys.argv[sys.argv.index('-o') + 1])\n"
        "out.mkdir(parents=True, exist_ok=True)\n"
        "(out / 'reCOGnizer_results.tsv').write_text('qseqid\\n')\n",
        encoding="utf-8",
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    (bin_dir / "fun.tsv").write_text("fixture\n", encoding="utf-8")
    for name in ("cog2ec_recognizer.tsv", "cog2ko_recognizer.tsv"):
        (bin_dir / name).write_text("fixture\n", encoding="utf-8")
    return executable


def seed_ready(root, native, databases):
    root.mkdir(parents=True)
    (root / RECOGNIZER_COMPLETION_INDICATOR).write_text(
        "2026-01-01 00:00:00", encoding="utf-8"
    )
    for name in RECOGNIZER_COMMON_RELATIONAL_FILES:
        (root / name).write_text("fixture\n", encoding="utf-8")
    for database in databases.split(","):
        for name in RECOGNIZER_DATABASE_RELATIONAL_FILES.get(database, ()):
            (root / name).write_text("fixture\n", encoding="utf-8")
        prefix, smp_prefixes = RECOGNIZER_PREFIXES[database]
        smps = root / "smps"
        smps.mkdir(exist_ok=True)
        (smps / f"{smp_prefixes[0]}fixture.smp").write_text(
            "fixture\n", encoding="utf-8"
        )
        indexes = root / "dbs"
        indexes.mkdir(exist_ok=True)
        for extension in RECOGNIZER_INDEX_EXTENSIONS:
            (indexes / f"{prefix}.{extension}").write_text(
                "fixture\n", encoding="utf-8"
            )
    if "COG" in databases.split(","):
        for name in RECOGNIZER_COG_MAP_SOURCES:
            (root / name).write_text("fixture\n", encoding="utf-8")
        for name in ("cog2ec_recognizer.tsv", "cog2ko_recognizer.tsv"):
            (native.parent / name).write_text("fixture\n", encoding="utf-8")


def make_snakemake(tmp_path, resources, extra="", *, databases=None, label=""):
    annotations = tmp_path / "bakta"
    annotations.mkdir(exist_ok=True)
    (annotations / "annotation.faa").write_text(">p\nM\n", encoding="utf-8")
    output_name = "recognizer" if not label else f"recognizer-{label}"
    output = tmp_path / "results" / output_name
    log = tmp_path / "results" / "logs" / f"{output_name}.log"
    return FakeSnakemake(
        input=Values(annotations=str(annotations)),
        output=Values(results=str(output)),
        params=Values(
            fasta_name="annotation.faa",
            databases=databases or "NCBI_Curated,Pfam,SMART,KOG,COG,PRK,TIGR",
            extra_config_key="config.recognizer_prok.extra",
            resources_dir=str(resources),
            extra=extra,
        ),
        log=[str(log)],
        threads=3,
    )


def run_wrapper_worker(snakemake, queue):
    try:
        runpy.run_path(str(WRAPPER), init_globals={"snakemake": snakemake})
    except BaseException as error:
        queue.put(("error", type(error).__name__, str(error)))
    else:
        queue.put(("ok",))


def seed_initializing_native(bin_dir):
    executable = seed_native(bin_dir)
    executable.write_text(
        "#!/usr/bin/env python3\n"
        "import json, os, sys, time\n"
        "from pathlib import Path\n"
        f"common = {RECOGNIZER_COMMON_RELATIONAL_FILES!r}\n"
        f"extensions = {RECOGNIZER_INDEX_EXTENSIONS!r}\n"
        "record = {'argv': sys.argv[1:], 'cwd': os.getcwd()}\n"
        "with Path(os.environ['FAKE_RECOGNIZER_RECORD']).open('a') as handle:\n"
        "    handle.write(json.dumps(record) + '\\n')\n"
        "root = Path(sys.argv[sys.argv.index('-rd') + 1])\n"
        "timestamp = root / 'recognizer_dwnl.timestamp'\n"
        "index = root / 'dbs' / 'Kog.aux'\n"
        "if not timestamp.exists() or not index.exists():\n"
        "    root.mkdir(parents=True, exist_ok=True)\n"
        "    (root / 'dbs').mkdir(exist_ok=True)\n"
        "    index.write_text('partial\\n')\n"
        "    Path(os.environ['FAKE_RECOGNIZER_INIT_MARKER']).touch()\n"
        "    time.sleep(0.35)\n"
        "    if os.environ.get('FAKE_RECOGNIZER_INIT_FAIL') == '1':\n"
        "        raise SystemExit(31)\n"
        "    for name in common:\n"
        "        (root / name).write_text('fixture\\n')\n"
        "    (root / 'kog').write_text('fixture\\n')\n"
        "    (root / 'smps').mkdir(exist_ok=True)\n"
        "    (root / 'smps' / 'KOGfixture.smp').write_text('fixture\\n')\n"
        "    for extension in extensions:\n"
        "        (root / 'dbs' / ('Kog.' + extension)).write_text('fixture\\n')\n"
        "    timestamp.write_text('2026-01-01 00:00:00\\n')\n"
        "out = Path(sys.argv[sys.argv.index('-o') + 1])\n"
        "out.mkdir(parents=True, exist_ok=True)\n"
        "(out / 'reCOGnizer_results.tsv').write_text('qseqid\\n')\n",
        encoding="utf-8",
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    return executable


def start_wrapper(context, queue, snakemake):
    process = context.Process(target=run_wrapper_worker, args=(snakemake, queue))
    process.start()
    return process


def test_ready_selection_is_shared_and_immutable(tmp_path):
    native = seed_native(tmp_path / "bin")
    resources = tmp_path / "resources"
    seed_ready(resources, native, "KOG,COG")
    before_resources = material_snapshot(resources)
    before_native = material_snapshot(native.parent)

    assessment = assess_recognizer_database(resources, "KOG,COG", executable=native)

    assert assessment.ready
    assert not assessment.may_mutate
    assert all(not request.exclusive for request in assessment.locks)
    assert material_snapshot(resources) == before_resources
    assert material_snapshot(native.parent) == before_native


def test_ready_selection_accepts_empty_cdtrack(tmp_path):
    native = seed_native(tmp_path / "bin")
    resources = tmp_path / "resources"
    seed_ready(resources, native, "NCBI_Curated")
    (resources / "cdtrack.txt").write_bytes(b"")

    assessment = assess_recognizer_database(
        resources, "NCBI_Curated", executable=native
    )

    assert assessment.ready
    assert not assessment.may_mutate


def test_ready_selection_rejects_missing_cdtrack(tmp_path):
    native = seed_native(tmp_path / "bin")
    resources = tmp_path / "resources"
    seed_ready(resources, native, "NCBI_Curated")
    (resources / "cdtrack.txt").unlink()

    with pytest.raises(CacheIntegrityError, match=r"relational resource cdtrack\.txt"):
        assess_recognizer_database(resources, "NCBI_Curated", executable=native)


def test_ready_selection_rejects_other_empty_common_relational_file(tmp_path):
    native = seed_native(tmp_path / "bin")
    resources = tmp_path / "resources"
    seed_ready(resources, native, "NCBI_Curated")
    (resources / "cddannot.dat").write_bytes(b"")

    with pytest.raises(CacheIntegrityError, match=r"relational resource cddannot\.dat"):
        assess_recognizer_database(resources, "NCBI_Curated", executable=native)


def test_missing_cache_and_lazy_material_require_exclusive_native_execution(tmp_path):
    native = seed_native(tmp_path / "bin")
    absent = assess_recognizer_database(tmp_path / "absent", "KOG", executable=native)
    assert absent.may_mutate
    assert absent.locks[0].exclusive

    resources = tmp_path / "ready"
    seed_ready(resources, native, "KOG,COG")
    (resources / "dbs" / "Kog.aux").unlink()
    with pytest.raises(CacheIntegrityError, match="index family"):
        assess_recognizer_database(resources, "KOG,COG", executable=native)
    preliminary = assess_recognizer_database(
        resources,
        "KOG,COG",
        executable=native,
        definitive=False,
    )
    assert preliminary.may_mutate
    assert preliminary.locks[0].exclusive

    # A wholly absent family is native 1.12.2's supported lazy-build case.
    for extension in RECOGNIZER_INDEX_EXTENSIONS:
        candidate = resources / "dbs" / f"Kog.{extension}"
        if candidate.exists():
            candidate.unlink()
    lazy = assess_recognizer_database(resources, "KOG,COG", executable=native)
    assert lazy.may_mutate
    assert lazy.locks[0].exclusive


def test_timestamp_only_and_interrupted_initialization_are_not_ready(tmp_path):
    native = seed_native(tmp_path / "bin")
    timestamp_only = tmp_path / "timestamp-only"
    timestamp_only.mkdir()
    (timestamp_only / RECOGNIZER_COMPLETION_INDICATOR).write_text("date\n")
    with pytest.raises(CacheIntegrityError, match="relational resource"):
        assess_recognizer_database(timestamp_only, "KOG", executable=native)

    interrupted = tmp_path / "interrupted"
    interrupted.mkdir()
    (interrupted / "cdd.tar.gz").write_bytes(b"partial")
    with pytest.raises(CacheIntegrityError, match="incomplete native initialization"):
        assess_recognizer_database(interrupted, "KOG", executable=native)


@pytest.mark.parametrize(
    ("missing_map_name", "required_sources", "unneeded_sources"),
    [
        (
            "cog2ec_recognizer.tsv",
            RECOGNIZER_COG2EC_SOURCES,
            RECOGNIZER_COG2KO_SOURCES,
        ),
        (
            "cog2ko_recognizer.tsv",
            RECOGNIZER_COG2KO_SOURCES,
            RECOGNIZER_COG2EC_SOURCES,
        ),
    ],
)
def test_each_missing_map_requires_only_its_own_sources(
    tmp_path, missing_map_name, required_sources, unneeded_sources
):
    native = seed_native(tmp_path / "bin")
    resources = tmp_path / "resources"
    seed_ready(resources, native, "COG")
    (native.parent / missing_map_name).unlink()
    for source in unneeded_sources:
        (resources / source).unlink()
    assessment = assess_recognizer_database(resources, "COG", executable=native)
    assert all(request.exclusive for request in assessment.locks)
    assert all((resources / source).is_file() for source in required_sources)


@pytest.mark.parametrize(
    ("missing_map_name", "missing_source"),
    [
        ("cog2ec_recognizer.tsv", RECOGNIZER_COG2EC_SOURCES[0]),
        ("cog2ko_recognizer.tsv", RECOGNIZER_COG2KO_SOURCES[0]),
    ],
)
def test_missing_map_still_requires_its_actual_source_pair(
    tmp_path, missing_map_name, missing_source
):
    native = seed_native(tmp_path / "bin")
    resources = tmp_path / "resources"
    seed_ready(resources, native, "COG")
    (native.parent / missing_map_name).unlink()
    (resources / missing_source).unlink()
    with pytest.raises(CacheIntegrityError, match=missing_source):
        assess_recognizer_database(resources, "COG", executable=native)


def test_advanced_modes_are_exclusive(tmp_path):
    native = seed_native(tmp_path / "bin")
    resources = tmp_path / "resources"
    seed_ready(resources, native, "COG")

    advanced = assess_recognizer_database(
        resources, "COG", "--tax-file taxonomy.tsv", executable=native
    )
    assert all(request.exclusive for request in advanced.locks)


def test_wrapper_preserves_original_command_cwd_extras_logging_and_failure(
    tmp_path, monkeypatch
):
    native = seed_native(tmp_path / "bin")
    resources = tmp_path / "resource root"
    seed_ready(resources, native, "NCBI_Curated,Pfam,SMART,KOG,COG,PRK,TIGR")
    snakemake = make_snakemake(tmp_path, resources, "--quiet --evalue 1e-5")
    records = tmp_path / "records.jsonl"
    monkeypatch.setenv("PATH", f"{native.parent}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_RECOGNIZER_RECORD", str(records))
    original_cwd = Path.cwd()
    resources_before = material_snapshot(resources)
    native_before = material_snapshot(native.parent)

    runpy.run_path(str(WRAPPER), init_globals={"snakemake": snakemake})

    assert Path.cwd() == original_cwd
    record = json.loads(records.read_text(encoding="utf-8").strip())
    fasta = str((Path(snakemake.input.annotations) / "annotation.faa").resolve())
    outdir = str(Path(snakemake.output.results).resolve())
    assert record["argv"] == [
        "-f",
        fasta,
        "-o",
        outdir,
        "-t",
        "3",
        "-dbs",
        "NCBI_Curated,Pfam,SMART,KOG,COG,PRK,TIGR",
        "-rd",
        str(resources),
        "--quiet",
        "--evalue",
        "1e-5",
    ]
    assert record["cwd"] == str(Path(snakemake.log[0]).resolve().parent)
    assert "fake recognizer output" in Path(snakemake.log[0]).read_text()
    assert material_snapshot(resources) == resources_before
    native_after = material_snapshot(native.parent)
    assert all(native_after[path] == value for path, value in native_before.items())
    assert all(
        path in native_before or path.name.startswith(".meew-recognizer-native-")
        for path in native_after
    )

    monkeypatch.setenv("FAKE_RECOGNIZER_FAIL", "23")
    with pytest.raises(subprocess.CalledProcessError) as failure:
        runpy.run_path(str(WRAPPER), init_globals={"snakemake": snakemake})
    assert failure.value.returncode == 23


def test_failed_native_initialization_leaves_evidence_and_blocks_later_job(
    tmp_path, monkeypatch
):
    native = seed_native(tmp_path / "bin")
    resources = tmp_path / "cold"
    snakemake = make_snakemake(tmp_path, resources)
    records = tmp_path / "records.jsonl"
    monkeypatch.setenv("PATH", f"{native.parent}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_RECOGNIZER_RECORD", str(records))
    monkeypatch.setenv("FAKE_RECOGNIZER_PARTIAL", "1")
    monkeypatch.setenv("FAKE_RECOGNIZER_FAIL", "19")

    with pytest.raises(subprocess.CalledProcessError):
        runpy.run_path(str(WRAPPER), init_globals={"snakemake": snakemake})
    assert (resources / "partial.native").is_file()
    assert not (resources / RECOGNIZER_COMPLETION_INDICATOR).exists()
    calls = records.read_text(encoding="utf-8").splitlines()

    with pytest.raises(CacheIntegrityError, match="incomplete native initialization"):
        runpy.run_path(str(WRAPPER), init_globals={"snakemake": snakemake})
    assert records.read_text(encoding="utf-8").splitlines() == calls


@pytest.mark.parametrize("initial_state", ["cold", "missing-index"])
def test_concurrent_full_wrappers_wait_for_native_initialization(
    tmp_path, monkeypatch, initial_state
):
    context = multiprocessing.get_context("fork")
    queue = context.Queue()
    native = seed_initializing_native(tmp_path / "bin")
    resources = tmp_path / "resources"
    if initial_state == "missing-index":
        seed_ready(resources, native, "KOG")
        for extension in RECOGNIZER_INDEX_EXTENSIONS:
            (resources / "dbs" / f"Kog.{extension}").unlink()
    marker = tmp_path / "initializer-partial"
    records = tmp_path / "records.jsonl"
    monkeypatch.setenv("PATH", f"{native.parent}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_RECOGNIZER_RECORD", str(records))
    monkeypatch.setenv("FAKE_RECOGNIZER_INIT_MARKER", str(marker))

    first = start_wrapper(
        context,
        queue,
        make_snakemake(tmp_path, resources, databases="KOG", label="first"),
    )
    deadline = time.monotonic() + 3
    while not marker.exists() and time.monotonic() < deadline:
        time.sleep(0.02)
    assert marker.exists()
    second = start_wrapper(
        context,
        queue,
        make_snakemake(tmp_path, resources, databases="KOG", label="second"),
    )
    for process in (first, second):
        process.join(5)
        assert process.exitcode == 0

    assert [queue.get(timeout=1)[0] for _ in range(2)] == ["ok", "ok"]
    assert len(records.read_text(encoding="utf-8").splitlines()) == 2
    assert (resources / RECOGNIZER_COMPLETION_INDICATOR).is_file()
    assert all(
        (resources / "dbs" / f"Kog.{extension}").is_file()
        for extension in RECOGNIZER_INDEX_EXTENSIONS
    )


def test_waiting_full_wrapper_blocks_after_failed_native_initialization(
    tmp_path, monkeypatch
):
    context = multiprocessing.get_context("fork")
    queue = context.Queue()
    native = seed_initializing_native(tmp_path / "bin")
    resources = tmp_path / "resources"
    marker = tmp_path / "initializer-partial"
    records = tmp_path / "records.jsonl"
    monkeypatch.setenv("PATH", f"{native.parent}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_RECOGNIZER_RECORD", str(records))
    monkeypatch.setenv("FAKE_RECOGNIZER_INIT_MARKER", str(marker))
    monkeypatch.setenv("FAKE_RECOGNIZER_INIT_FAIL", "1")

    first = start_wrapper(
        context,
        queue,
        make_snakemake(tmp_path, resources, databases="KOG", label="first"),
    )
    deadline = time.monotonic() + 3
    while not marker.exists() and time.monotonic() < deadline:
        time.sleep(0.02)
    assert marker.exists()
    second = start_wrapper(
        context,
        queue,
        make_snakemake(tmp_path, resources, databases="KOG", label="second"),
    )
    for process in (first, second):
        process.join(5)
        assert process.exitcode == 0

    results = [queue.get(timeout=1) for _ in range(2)]
    assert {result[1] for result in results} == {
        "CalledProcessError",
        "CacheIntegrityError",
    }
    assert len(records.read_text(encoding="utf-8").splitlines()) == 1
    assert not (resources / RECOGNIZER_COMPLETION_INDICATOR).exists()
