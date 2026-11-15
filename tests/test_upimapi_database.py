import base64
import hashlib
import json
import multiprocessing
import os
import runpy
import shlex
import shutil
import stat
import struct
import subprocess
import time
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pytest

from common.cache_runtime import CacheIntegrityError
from common.native_execution import run_native_command
from common.upimapi_database import (
    DIAMOND_MAGIC,
    assess_upimapi_database,
    selected_upimapi_reference,
)
from helpers import material_snapshot

ROOT = Path(__file__).resolve().parents[1]
WRAPPER = ROOT / "workflow" / "scripts" / "upimapi.py"
DIAMOND_FIXTURES = ROOT / "tests" / "fixtures" / "diamond"
PINNED_DIAMOND = Path(
    "/home/argomes/data/meta/EST6/.snakemake/conda/"
    "269ab901aa2ce5806b168821fc68d8ec_/bin/diamond"
)


class Values(dict):
    __getattr__ = dict.__getitem__


class FakeSnakemake(SimpleNamespace):
    def log_fmt_shell(self, stdout, stderr):
        assert stdout and stderr
        return f"> {shlex.quote(str(self.log[0]))} 2>&1"


def diamond_executable():
    configured = os.environ.get("MEEW_TEST_DIAMOND")
    candidates = [
        Path(configured) if configured else None,
        PINNED_DIAMOND,
        Path(shutil.which("diamond")) if shutil.which("diamond") else None,
    ]
    for candidate in candidates:
        if candidate is not None and candidate.is_file():
            return candidate.resolve()
    pytest.skip("A DIAMOND executable is required for native index validation tests")


def fixture_index(name="small"):
    encoded = (DIAMOND_FIXTURES / f"{name}.dmnd.b64").read_text(encoding="ascii")
    return base64.b64decode(encoded)


def write_index(path, name="small"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(fixture_index(name))


def seed_pair(fasta):
    fasta.parent.mkdir(parents=True, exist_ok=True)
    fasta.write_text(">ref\nM\n", encoding="utf-8")
    index = fasta.with_suffix(".dmnd")
    write_index(index)
    return index


def seed_native(bin_dir):
    bin_dir.mkdir(parents=True)
    (bin_dir / "diamond").symlink_to(diamond_executable())
    executable = bin_dir / "upimapi"
    executable.write_text(
        "#!/usr/bin/env python3\n"
        "import json, os, sys\n"
        "from pathlib import Path\n"
        "with Path(os.environ['FAKE_UPIMAPI_RECORD']).open('a') as handle:\n"
        "    handle.write(json.dumps({'argv': sys.argv[1:], 'cwd': os.getcwd()}) + '\\n')\n"
        "print('fake upimapi output')\n"
        "if os.environ.get('FAKE_UPIMAPI_FAIL'):\n"
        "    raise SystemExit(int(os.environ['FAKE_UPIMAPI_FAIL']))\n"
        "out = Path(sys.argv[sys.argv.index('--output') + 1])\n"
        "out.mkdir(parents=True, exist_ok=True)\n"
        "(out / 'UPIMAPI_results.tsv').write_text('qseqid\\n')\n",
        encoding="utf-8",
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    return executable


def without_guard_locks(snapshot):
    return {
        path: value
        for path, value in snapshot.items()
        if not path.name.startswith(".meew-upimapi-material-")
    }


def make_snakemake(
    tmp_path,
    resources,
    *,
    db="swissprot",
    db_custom="",
    skip_if_exists=True,
    skip_db=False,
    extra="",
    label="",
):
    annotations = tmp_path / "bakta"
    annotations.mkdir(exist_ok=True)
    (annotations / "annotation.faa").write_text(">p\nM\n", encoding="utf-8")
    output_name = "upimapi" if not label else f"upimapi-{label}"
    output = tmp_path / "results" / output_name
    log = tmp_path / "results" / "logs" / f"{output_name}.log"
    return FakeSnakemake(
        input=Values(annotations=str(annotations)),
        output=Values(results=str(output)),
        params=Values(
            fasta_name="annotation.faa",
            db=db,
            db_custom=str(db_custom),
            resources_dir=str(resources),
            extra=extra,
            skip_db_check_if_exists=skip_if_exists,
            skip_db=skip_db,
        ),
        log=[str(log)],
        threads=4,
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
    complete_index = fixture_index()
    executable.write_text(
        "#!/usr/bin/env python3\n"
        "import json, os, sys, time\n"
        "from pathlib import Path\n"
        f"magic = bytes.fromhex('{DIAMOND_MAGIC.hex()}')\n"
        f"complete_index = {complete_index!r}\n"
        "with Path(os.environ['FAKE_UPIMAPI_RECORD']).open('a') as handle:\n"
        "    handle.write(json.dumps({'argv': sys.argv[1:], 'cwd': os.getcwd()}) + '\\n')\n"
        "root = Path(sys.argv[sys.argv.index('-rd') + 1])\n"
        "fasta = root / 'uniprot_sprot.fasta'\n"
        "index = root / 'uniprot_sprot.dmnd'\n"
        "if not index.exists() or index.read_bytes() != complete_index:\n"
        "    Path(os.environ['FAKE_UPIMAPI_EARLY_MARKER']).touch()\n"
        "    time.sleep(0.2)\n"
        "    root.mkdir(parents=True, exist_ok=True)\n"
        "    if not fasta.exists():\n"
        "        fasta.write_text('>ref\\nM\\n')\n"
        "    index.write_bytes(magic)\n"
        "    Path(os.environ['FAKE_UPIMAPI_PARTIAL_MARKER']).touch()\n"
        "    time.sleep(0.3)\n"
        "    if os.environ.get('FAKE_UPIMAPI_INIT_FAIL') == '1':\n"
        "        raise SystemExit(37)\n"
        "    index.write_bytes(complete_index)\n"
        "out = Path(sys.argv[sys.argv.index('--output') + 1])\n"
        "out.mkdir(parents=True, exist_ok=True)\n"
        "(out / 'UPIMAPI_results.tsv').write_text('qseqid\\n')\n",
        encoding="utf-8",
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    return executable


def start_wrapper(context, queue, snakemake):
    process = context.Process(target=run_wrapper_worker, args=(snakemake, queue))
    process.start()
    return process


def assess_locked(resources_dir, db, db_custom, **kwargs):
    preliminary = assess_upimapi_database(
        resources_dir,
        db,
        db_custom,
        definitive=False,
        **kwargs,
    )
    result = {}

    def verify(probe):
        result["assessment"] = assess_upimapi_database(
            resources_dir,
            db,
            db_custom,
            definitive=True,
            diamond_probe=probe,
            **kwargs,
        )

    path = f"{diamond_executable().parent}{os.pathsep}{os.environ['PATH']}"
    with mock.patch.dict(os.environ, {"PATH": path}):
        run_native_command(
            ":",
            preliminary.locks,
            prelaunch_probe_check=verify,
        )
    return result["assessment"]


def test_ready_builtin_is_shared_and_missing_index_is_exclusive(tmp_path):
    resources = tmp_path / "resources"
    fasta = resources / "uniprot_sprot.fasta"
    seed_pair(fasta)
    before = material_snapshot(resources)

    preliminary = assess_upimapi_database(
        resources,
        "swissprot",
        "",
        skip_db_check=True,
        definitive=False,
    )
    assert not preliminary.ready
    assert not preliminary.may_mutate
    assert all(not request.exclusive for request in preliminary.locks)

    ready = assess_locked(resources, "swissprot", "", skip_db_check=True)
    assert ready.ready
    assert all(not request.exclusive for request in ready.locks)
    assert without_guard_locks(material_snapshot(resources)) == before

    fasta.with_suffix(".dmnd").unlink()
    missing = assess_upimapi_database(
        resources,
        "swissprot",
        "",
        skip_db_check=True,
        definitive=False,
    )
    assert missing.may_mutate
    assert missing.locks[0].exclusive


def test_corrupt_index_and_incomplete_pair_are_never_ready(tmp_path):
    resources = tmp_path / "resources"
    fasta = resources / "uniprot_sprot.fasta"
    seed_pair(fasta)
    fasta.with_suffix(".dmnd").write_bytes(b"not diamond")
    with pytest.raises(CacheIntegrityError, match="incomplete or invalid"):
        assess_locked(resources, "swissprot", "", skip_db_check=True)

    fasta.unlink()
    write_index(fasta.with_suffix(".dmnd"))
    with pytest.raises(CacheIntegrityError, match="without its selected FASTA"):
        assess_locked(resources, "swissprot", "", skip_db_check=True)


@pytest.mark.parametrize(
    ("fixture_name", "size", "sha256", "sequences", "letters"),
    [
        (
            "small",
            139,
            "ced6b1eea34766b32e1eea20ade2f1868f31fed8403d800171d46f595a68e848",
            "Sequences  1",
            "Letters  1",
        ),
        (
            "typical",
            444,
            "1073bc45a66e51ee799b0abb8e0a95cfcbc19a1d7dfe7708bb1f5480af2b2f2c",
            "Sequences  4",
            "Letters  165",
        ),
    ],
)
def test_native_derived_small_and_typical_indexes_are_accepted_under_locks(
    tmp_path, fixture_name, size, sha256, sequences, letters
):
    resources = tmp_path / fixture_name
    fasta = resources / "uniprot_sprot.fasta"
    fasta.parent.mkdir(parents=True)
    fasta.write_bytes((DIAMOND_FIXTURES / f"{fixture_name}.faa").read_bytes())
    index = fasta.with_suffix(".dmnd")
    write_index(index, fixture_name)

    assert index.stat().st_size == size
    assert hashlib.sha256(index.read_bytes()).hexdigest() == sha256
    native = subprocess.run(
        [diamond_executable(), "dbinfo", "--db", index],
        check=True,
        capture_output=True,
        text=True,
    )
    assert sequences in native.stdout
    assert letters in native.stdout

    assessment = assess_locked(resources, "swissprot", "", skip_db_check=True)
    assert assessment.ready
    assert not assessment.may_mutate
    assert all(not request.exclusive for request in assessment.locks)


@pytest.mark.parametrize(
    "case",
    [
        "magic-only",
        "truncated-header",
        "truncated-body",
        "padded-fake",
        "invalid-offset",
        "invalid-content",
    ],
)
def test_truncated_or_structurally_invalid_index_is_not_ready(tmp_path, case):
    resources = tmp_path / "resources"
    fasta = resources / "uniprot_sprot.fasta"
    seed_pair(fasta)
    native = bytearray(fixture_index())
    payloads = {
        "magic-only": bytes(native[:8]),
        "truncated-header": bytes(native[:31]),
        "truncated-body": bytes(native[:106]),
        "padded-fake": bytes(native[:32]) + b"\x01" + b"\x00" * (4096 - 33),
        "invalid-offset": bytes(native[:32])
        + struct.pack("<Q", 10_000)
        + bytes(native[40:]),
        "invalid-content": bytes(native[:98]) + b"\x00" + bytes(native[99:]),
    }
    payload = payloads[case]
    fasta.with_suffix(".dmnd").write_bytes(payload)

    preliminary = assess_upimapi_database(
        resources,
        "swissprot",
        "",
        skip_db_check=True,
        definitive=False,
    )
    assert preliminary.may_mutate
    assert any(request.exclusive for request in preliminary.locks)
    with pytest.raises(CacheIntegrityError, match="incomplete or invalid"):
        assess_locked(resources, "swissprot", "", skip_db_check=True)


def test_structural_guard_rejects_padded_header_that_dbinfo_accepts(tmp_path):
    resources = tmp_path / "resources"
    fasta = resources / "uniprot_sprot.fasta"
    seed_pair(fasta)
    native = fixture_index()
    fake = fasta.with_suffix(".dmnd")
    fake.write_bytes(native[:32] + b"\x01" + b"\x00" * (4096 - 33))
    extracted = tmp_path / "native-read.faa"

    dbinfo = subprocess.run(
        [diamond_executable(), "dbinfo", "--db", fake],
        check=False,
        capture_output=True,
        text=True,
    )
    getseq = subprocess.run(
        [
            diamond_executable(),
            "getseq",
            "--db",
            fake,
            "--seq",
            "1",
            "--out",
            extracted,
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert dbinfo.returncode == 0
    assert getseq.returncode == 0
    assert extracted.read_bytes().startswith(b">\nAAAA")
    with pytest.raises(CacheIntegrityError, match="incomplete or invalid"):
        assess_locked(resources, "swissprot", "", skip_db_check=True)


@pytest.mark.parametrize(
    "payload",
    [DIAMOND_MAGIC, fixture_index()[:98] + b"\x00" + fixture_index()[99:]],
    ids=["magic-only", "invalid-body"],
)
def test_wrapper_rejects_truncated_index_before_native_command(
    tmp_path, monkeypatch, payload
):
    native = seed_native(tmp_path / "bin")
    resources = tmp_path / "resources"
    fasta = resources / "uniprot_sprot.fasta"
    seed_pair(fasta)
    fasta.with_suffix(".dmnd").write_bytes(payload)
    records = tmp_path / "records.jsonl"
    monkeypatch.setenv("PATH", f"{native.parent}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_UPIMAPI_RECORD", str(records))

    with pytest.raises(CacheIntegrityError, match="incomplete or invalid"):
        runpy.run_path(
            str(WRAPPER),
            init_globals={"snakemake": make_snakemake(tmp_path, resources)},
        )
    assert not records.exists()


def test_custom_reference_aliases_share_index_identity_and_advanced_modes_are_exclusive(
    tmp_path,
):
    resources = tmp_path / "resources"
    resources.mkdir()
    custom = tmp_path / "references" / "custom.fasta"
    seed_pair(custom)
    alias = tmp_path / "custom-alias.fasta"
    alias.symlink_to(custom)

    direct = assess_locked(resources, "", custom, skip_db_check=True)
    through_alias = assess_locked(resources, "", alias, skip_db_check=True)
    assert direct.ready
    # Native derives a different next-to-the-CLI-string alias.dmnd, so this
    # alias needs an exclusive index build while still coordinating on the
    # canonical source inode.
    assert not through_alias.ready
    assert direct.locks[1].lock_path == through_alias.locks[1].lock_path
    assert any(
        request.exclusive and request.namespace == "upimapi-material"
        for request in through_alias.locks
    )

    local = assess_locked(
        resources,
        "",
        custom,
        skip_db_check=True,
        extra="--local-id-mapping",
    )
    assert all(request.exclusive for request in local.locks)


def _material_lock(assessment, payload):
    payload = Path(payload).resolve()
    return next(
        request
        for request in assessment.locks
        if request.namespace == "upimapi-material" and request.payload == payload
    )


def test_builtin_and_custom_selections_share_canonical_material_locks(tmp_path):
    builtin_root = tmp_path / "builtin"
    fasta = builtin_root / "uniprot_sprot.fasta"
    index = seed_pair(fasta)
    other_root = tmp_path / "other-resources"
    other_root.mkdir()

    builtin_reader = assess_upimapi_database(
        builtin_root, "swissprot", "", skip_db_check=True, definitive=False
    )
    custom_fasta_reader = assess_upimapi_database(
        other_root, "", fasta, skip_db_check=True, definitive=False
    )
    custom_index_reader = assess_upimapi_database(
        other_root, "", index, skip_db_check=True, definitive=False
    )
    builtin_writer = assess_upimapi_database(
        builtin_root, "swissprot", "", skip_db_check=False, definitive=False
    )

    assert (
        _material_lock(builtin_reader, fasta).lock_path
        == _material_lock(custom_fasta_reader, fasta).lock_path
    )
    builtin_index_lock = _material_lock(builtin_reader, index)
    assert (
        builtin_index_lock.lock_path
        == _material_lock(custom_fasta_reader, index).lock_path
    )
    assert (
        builtin_index_lock.lock_path
        == _material_lock(custom_index_reader, index).lock_path
    )
    assert _material_lock(builtin_writer, fasta).exclusive
    assert _material_lock(builtin_writer, index).exclusive

    fasta_alias = tmp_path / "reference-alias.fasta"
    fasta_alias.symlink_to(fasta)
    index_alias = tmp_path / "reference-alias.dmnd"
    index_alias.symlink_to(index)
    alias_reader = assess_upimapi_database(
        other_root, "", fasta_alias, skip_db_check=True, definitive=False
    )
    direct_index_alias_reader = assess_upimapi_database(
        other_root, "", index_alias, skip_db_check=True, definitive=False
    )
    assert (
        _material_lock(alias_reader, fasta).lock_path
        == _material_lock(builtin_reader, fasta).lock_path
    )
    assert _material_lock(alias_reader, index).lock_path == builtin_index_lock.lock_path
    assert (
        _material_lock(direct_index_alias_reader, index).lock_path
        == builtin_index_lock.lock_path
    )


@pytest.mark.parametrize(
    ("skip_if_exists", "skip_db", "expected_skip"),
    [(True, False, True), (False, False, False), (False, True, True)],
)
def test_wrapper_preserves_builtin_skip_branches_command_cwd_and_logging(
    tmp_path,
    monkeypatch,
    skip_if_exists,
    skip_db,
    expected_skip,
):
    native = seed_native(tmp_path / "bin")
    resources = tmp_path / "resource root"
    seed_pair(resources / "uniprot_sprot.fasta")
    snakemake = make_snakemake(
        tmp_path,
        resources,
        skip_if_exists=skip_if_exists,
        skip_db=skip_db,
        extra="--sleep 0 --max-tries 2",
    )
    records = tmp_path / "records.jsonl"
    monkeypatch.setenv("PATH", f"{native.parent}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_UPIMAPI_RECORD", str(records))
    invocation_cwd = Path.cwd()
    resources_before = material_snapshot(resources)

    runpy.run_path(str(WRAPPER), init_globals={"snakemake": snakemake})

    record = json.loads(records.read_text(encoding="utf-8").strip())
    assert record["cwd"] == str(invocation_cwd)
    assert record["argv"] == [
        "--input",
        str((Path(snakemake.input.annotations) / "annotation.faa").resolve()),
        "--output",
        str(Path(snakemake.output.results).resolve()),
        "--database",
        "swissprot",
        "-rd",
        str(resources),
        *(["--skip-db-check"] if expected_skip else []),
        "--threads",
        "4",
        "--max-target-seqs",
        "1",
        "--sleep",
        "0",
        "--max-tries",
        "2",
    ]
    assert "fake upimapi output" in Path(snakemake.log[0]).read_text()
    assert without_guard_locks(material_snapshot(resources)) == resources_before


def test_uniprot_default_reference_and_wrapper_invocation(tmp_path, monkeypatch):
    native = seed_native(tmp_path / "bin")
    resources = tmp_path / "resource root"
    reference = resources / "uniprot.fasta"
    index = seed_pair(reference)
    snakemake = make_snakemake(tmp_path, resources, db="uniprot")
    records = tmp_path / "records.jsonl"
    monkeypatch.setenv("PATH", f"{native.parent}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_UPIMAPI_RECORD", str(records))

    selected_reference, selected_index = selected_upimapi_reference(
        resources, "uniprot", ""
    )
    assert selected_reference == reference.resolve()
    assert selected_index == index.resolve()

    runpy.run_path(str(WRAPPER), init_globals={"snakemake": snakemake})

    args = json.loads(records.read_text(encoding="utf-8").strip())["argv"]
    database_at = args.index("--database")
    assert args[database_at + 1] == "uniprot"
    assert "--skip-db-check" in args
    assert args[args.index("--max-target-seqs") + 1] == "1"


@pytest.mark.parametrize("index_ready", [True, False], ids=["ready", "missing-index"])
def test_wrapper_preserves_custom_database_argument(tmp_path, monkeypatch, index_ready):
    native = seed_native(tmp_path / "bin")
    resources = tmp_path / "resources"
    resources.mkdir()
    custom = tmp_path / "custom reference.fasta"
    index = seed_pair(custom)
    if not index_ready:
        index.unlink()
    snakemake = make_snakemake(tmp_path, resources, db="", db_custom=custom)
    records = tmp_path / "records.jsonl"
    monkeypatch.setenv("PATH", f"{native.parent}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_UPIMAPI_RECORD", str(records))

    runpy.run_path(str(WRAPPER), init_globals={"snakemake": snakemake})

    args = json.loads(records.read_text(encoding="utf-8").strip())["argv"]
    database_at = args.index("--database")
    assert args[database_at + 1] == str(custom)
    assert "swissprot" not in args
    assert ("--skip-db-check" in args) is True
    assert not index.exists() if not index_ready else index.exists()


def test_wrapper_preserves_direct_dmnd_selection(tmp_path, monkeypatch):
    native = seed_native(tmp_path / "bin")
    resources = tmp_path / "resources"
    resources.mkdir()
    direct = tmp_path / "direct.dmnd"
    write_index(direct)
    snakemake = make_snakemake(tmp_path, resources, db="", db_custom=direct)
    records = tmp_path / "records.jsonl"
    monkeypatch.setenv("PATH", f"{native.parent}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_UPIMAPI_RECORD", str(records))

    assessment = assess_locked(resources, "", direct, skip_db_check=True)
    assert assessment.ready
    assert assessment.reference == assessment.diamond_index == direct.resolve()

    runpy.run_path(str(WRAPPER), init_globals={"snakemake": snakemake})

    args = json.loads(records.read_text(encoding="utf-8").strip())["argv"]
    database_at = args.index("--database")
    assert args[database_at + 1] == str(direct)
    assert "--skip-db-check" in args


def test_wrapper_keeps_db_custom_exclusion_and_native_failure(tmp_path, monkeypatch):
    native = seed_native(tmp_path / "bin")
    resources = tmp_path / "resources"
    seed_pair(resources / "uniprot_sprot.fasta")
    records = tmp_path / "records.jsonl"
    monkeypatch.setenv("PATH", f"{native.parent}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_UPIMAPI_RECORD", str(records))

    invalid = make_snakemake(
        tmp_path, resources, db="swissprot", db_custom=resources / "custom.fasta"
    )
    with pytest.raises(ValueError, match="Can choose only one option"):
        runpy.run_path(str(WRAPPER), init_globals={"snakemake": invalid})
    assert not records.exists()

    failing = make_snakemake(tmp_path, resources)
    monkeypatch.setenv("FAKE_UPIMAPI_FAIL", "29")
    with pytest.raises(subprocess.CalledProcessError) as failure:
        runpy.run_path(str(WRAPPER), init_globals={"snakemake": failing})
    assert failure.value.returncode == 29


def test_native_index_path_derivation_matches_upimapi_1133(tmp_path):
    resources = tmp_path / "resources"
    reference, index = selected_upimapi_reference(resources, "swissprot", "")
    assert reference == (resources / "uniprot_sprot.fasta").resolve()
    assert index == (resources / "uniprot_sprot.dmnd").resolve()

    no_suffix = tmp_path / "custom"
    _, native_oddity = selected_upimapi_reference(resources, "", no_suffix)
    assert native_oddity == (Path.cwd() / ".dmnd").resolve()


def test_concurrent_cold_wrappers_wait_and_recompute_original_skip_branch(
    tmp_path, monkeypatch
):
    context = multiprocessing.get_context("fork")
    queue = context.Queue()
    native = seed_initializing_native(tmp_path / "bin")
    resources = tmp_path / "resources"
    records = tmp_path / "records.jsonl"
    early_marker = tmp_path / "initializer-locked"
    partial_marker = tmp_path / "initializer-partial"
    monkeypatch.setenv("PATH", f"{native.parent}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_UPIMAPI_RECORD", str(records))
    monkeypatch.setenv("FAKE_UPIMAPI_EARLY_MARKER", str(early_marker))
    monkeypatch.setenv("FAKE_UPIMAPI_PARTIAL_MARKER", str(partial_marker))

    first = start_wrapper(
        context, queue, make_snakemake(tmp_path, resources, label="initializer")
    )
    deadline = time.monotonic() + 3
    while not early_marker.exists() and time.monotonic() < deadline:
        time.sleep(0.02)
    assert early_marker.exists()
    stale_waiter = start_wrapper(
        context, queue, make_snakemake(tmp_path, resources, label="stale-waiter")
    )
    deadline = time.monotonic() + 3
    while not partial_marker.exists() and time.monotonic() < deadline:
        time.sleep(0.02)
    assert partial_marker.exists()
    partial_waiter = start_wrapper(
        context, queue, make_snakemake(tmp_path, resources, label="partial-waiter")
    )

    for process in (first, stale_waiter, partial_waiter):
        process.join(6)
        assert process.exitcode == 0
    assert [queue.get(timeout=1)[0] for _ in range(3)] == ["ok", "ok", "ok"]

    invocations = [
        json.loads(line)["argv"] for line in records.read_text().splitlines()
    ]
    assert len(invocations) == 3
    assert "--skip-db-check" not in invocations[0]
    assert all("--skip-db-check" in args for args in invocations[1:])


def test_concurrent_wrappers_wait_for_missing_diamond_index_build(
    tmp_path, monkeypatch
):
    context = multiprocessing.get_context("fork")
    queue = context.Queue()
    native = seed_initializing_native(tmp_path / "bin")
    resources = tmp_path / "resources"
    fasta = resources / "uniprot_sprot.fasta"
    seed_pair(fasta).unlink()
    records = tmp_path / "records.jsonl"
    early_marker = tmp_path / "initializer-locked"
    partial_marker = tmp_path / "initializer-partial"
    monkeypatch.setenv("PATH", f"{native.parent}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_UPIMAPI_RECORD", str(records))
    monkeypatch.setenv("FAKE_UPIMAPI_EARLY_MARKER", str(early_marker))
    monkeypatch.setenv("FAKE_UPIMAPI_PARTIAL_MARKER", str(partial_marker))

    first = start_wrapper(
        context, queue, make_snakemake(tmp_path, resources, label="initializer")
    )
    deadline = time.monotonic() + 3
    while not partial_marker.exists() and time.monotonic() < deadline:
        time.sleep(0.02)
    assert partial_marker.exists()
    second = start_wrapper(
        context, queue, make_snakemake(tmp_path, resources, label="waiter")
    )
    for process in (first, second):
        process.join(6)
        assert process.exitcode == 0
    assert [queue.get(timeout=1)[0] for _ in range(2)] == ["ok", "ok"]
    assert len(records.read_text().splitlines()) == 2


def test_waiting_wrapper_blocks_after_failed_diamond_initialization(
    tmp_path, monkeypatch
):
    context = multiprocessing.get_context("fork")
    queue = context.Queue()
    native = seed_initializing_native(tmp_path / "bin")
    resources = tmp_path / "resources"
    records = tmp_path / "records.jsonl"
    early_marker = tmp_path / "initializer-locked"
    partial_marker = tmp_path / "initializer-partial"
    monkeypatch.setenv("PATH", f"{native.parent}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_UPIMAPI_RECORD", str(records))
    monkeypatch.setenv("FAKE_UPIMAPI_EARLY_MARKER", str(early_marker))
    monkeypatch.setenv("FAKE_UPIMAPI_PARTIAL_MARKER", str(partial_marker))
    monkeypatch.setenv("FAKE_UPIMAPI_INIT_FAIL", "1")

    first = start_wrapper(
        context, queue, make_snakemake(tmp_path, resources, label="initializer")
    )
    deadline = time.monotonic() + 3
    while not partial_marker.exists() and time.monotonic() < deadline:
        time.sleep(0.02)
    assert partial_marker.exists()
    second = start_wrapper(
        context, queue, make_snakemake(tmp_path, resources, label="waiter")
    )
    for process in (first, second):
        process.join(6)
        assert process.exitcode == 0

    results = [queue.get(timeout=1) for _ in range(2)]
    assert {result[1] for result in results} == {
        "CalledProcessError",
        "CacheIntegrityError",
    }
    assert len(records.read_text().splitlines()) == 1
