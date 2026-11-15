import json
import multiprocessing
import os
import shutil
import ssl
import subprocess
import tarfile
import time
import urllib.error
from pathlib import Path

import pytest

from common import bakta_database as bakta_module
from common import download_retries as retry_module
from common import funannotate2_database as funannotate2_module
from common import gtdbtk_database as gtdbtk_module
from common import metaeuk_database as metaeuk_module
from common import qa_databases as qa_module
from common.bakta_database import (
    BAKTA_COMMON_FILES,
    ensure_bakta_database,
    validate_bakta_database,
)
from common.cache_runtime import CacheIntegrityError
from common.database_dependencies import (
    bakta_database_inputs,
    funannotate2_database_inputs,
    gtdbtk_database_inputs,
    metaeuk_database_inputs,
)
from common.funannotate2_database import (
    FUNANNOTATE2_ALL_DATABASES,
    ensure_funannotate2_database,
    expand_funannotate2_databases,
)
from common.gtdbtk_database import (
    GTDBTK_REQUIRED_DIRECTORIES,
    GTDBTK_REQUIRED_FILES,
    ensure_gtdbtk_database,
    validate_gtdbtk_database,
)
from common.metaeuk_database import (
    METAEUK_IDENTITY_NAME,
    ensure_metaeuk_database,
    validate_metaeuk_native_database,
)
from common.qa_databases import (
    ensure_checkm2_database,
    ensure_eukcc_database,
    selected_checkm2_database,
    validate_eukcc_database,
)
from helpers import file_snapshot as cache_snapshot


def seed_bakta(path, database_type="full"):
    path.mkdir(parents=True)
    (path / "version.json").write_text(
        json.dumps(
            {"date": "2026-01-01", "major": 6, "minor": 1, "type": database_type}
        ),
        encoding="utf-8",
    )
    for name in (
        *BAKTA_COMMON_FILES,
        "psc.dmnd" if database_type == "full" else "pscc.dmnd",
    ):
        (path / name).write_bytes(b"fixture\n")
    return path


def seed_metaeuk(prefix, catalog=None):
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for path in (prefix, Path(f"{prefix}.index"), Path(f"{prefix}.dbtype")):
        path.write_bytes(b"fixture\n")
    if catalog is not None:
        (prefix.parent / METAEUK_IDENTITY_NAME).write_text(
            json.dumps(
                {
                    "schema": 1,
                    "catalog": catalog,
                    "database_prefix": prefix.name,
                    "mmseqs_version": "17.b804f",
                    "managed_by": "MeeW-Annotation",
                }
            ),
            encoding="utf-8",
        )
    return prefix


def seed_funannotate_entry(root, database, filename=None):
    root.mkdir(parents=True, exist_ok=True)
    filename = filename or f"{database}.db"
    material = root / filename
    material.write_bytes(b"fixture\n")
    return {
        "db": str(material.resolve()),
        "type": "fixture",
        "version": "1",
        "date": "2026-01-01",
        "n_records": 1,
        "md5": "fixture",
    }


def seed_gtdbtk(root, release="r232"):
    for relative in GTDBTK_REQUIRED_FILES:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        content = f"VERSION_DATA={release}\n" if path.name == "metadata.txt" else "x\n"
        path.write_text(content, encoding="utf-8")
    for relative in GTDBTK_REQUIRED_DIRECTORIES:
        path = root / relative
        path.mkdir(parents=True, exist_ok=True)
        (path / "fixture").write_text("x\n", encoding="utf-8")
    return root


def seed_eukcc(root):
    base = root / "db_base"
    base.mkdir(parents=True)
    (base / "material.hmm").write_text("fixture\n", encoding="utf-8")
    (base / "placement").mkdir()
    (base / "placement" / "marker.hmm").write_text("fixture\n", encoding="utf-8")
    (base / "db_info.json").write_text(
        json.dumps(
            {
                "files": {
                    "hmm_db": "material.hmm",
                    "hmm_placement_place": "placement",
                }
            }
        ),
        encoding="utf-8",
    )
    return root


def append_attempt_output(log_path, message):
    with Path(log_path).open("a", encoding="utf-8") as handle:
        handle.write(message + "\n")


def test_read_only_database_dependencies_are_optional_and_select_native_material(
    tmp_path,
):
    assert bakta_database_inputs(tmp_path / "missing-bakta", "full") == []
    assert gtdbtk_database_inputs(tmp_path / "missing-gtdbtk") == []
    assert (
        funannotate2_database_inputs(tmp_path / "missing-funannotate", ["merops"]) == []
    )
    assert metaeuk_database_inputs(tmp_path / "missing-metaeuk" / "db") == []

    bakta = seed_bakta(tmp_path / "bakta")
    assert str((bakta / "version.json").resolve()) in bakta_database_inputs(
        bakta, "full"
    )
    assert str((bakta / "bakta.db").resolve()) in bakta_database_inputs(bakta, "full")

    gtdbtk = seed_gtdbtk(tmp_path / "gtdbtk")
    assert set(gtdbtk_database_inputs(gtdbtk)) == {
        str((gtdbtk / relative).resolve()) for relative in GTDBTK_REQUIRED_FILES
    }

    funannotate = tmp_path / "funannotate"
    merops = seed_funannotate_entry(funannotate, "merops")
    fun_info = funannotate / "funannotate-db-info.json"
    fun_info.write_text(json.dumps({"merops": merops}), encoding="utf-8")
    assert set(funannotate2_database_inputs(funannotate, ["merops"])) == {
        str(fun_info.resolve()),
        str((funannotate / "merops.db").resolve()),
    }

    metaeuk = seed_metaeuk(
        tmp_path / "metaeuk" / "catalogs" / "uniref90" / "db", "UniRef90"
    )
    assert set(metaeuk_database_inputs(metaeuk)) == {
        str(metaeuk.resolve()),
        str(Path(f"{metaeuk}.index").resolve()),
        str(Path(f"{metaeuk}.dbtype").resolve()),
        str((metaeuk.parent / METAEUK_IDENTITY_NAME).resolve()),
    }


def _ensure_bakta_worker(database, log, queue):
    try:
        result = ensure_bakta_database(database, "light", log)
        queue.put(result.status)
    except BaseException as error:  # pragma: no cover - reported in parent
        queue.put(f"error:{error!r}")


def test_bakta_ready_reuse_is_read_only_and_never_calls_installer(
    tmp_path, monkeypatch
):
    database = seed_bakta(tmp_path / "db", "full")
    before = cache_snapshot(database)
    for path in database.rglob("*"):
        path.chmod(0o555 if path.is_dir() else 0o444)
    database.chmod(0o555)

    def unexpected(*args, **kwargs):
        raise AssertionError("ready cache must not invoke bakta_db")

    monkeypatch.setattr(bakta_module, "run_logged", unexpected)
    result = ensure_bakta_database(database, "full", tmp_path / "annotation.log")

    assert result.status == "reused"
    assert cache_snapshot(database) == before
    database.chmod(0o755)
    for path in database.rglob("*"):
        path.chmod(0o755 if path.is_dir() else 0o644)


def test_bakta_rejects_wrong_type_schema_and_partial_cache(tmp_path):
    light = seed_bakta(tmp_path / "light", "light")
    with pytest.raises(CacheIntegrityError, match="type"):
        validate_bakta_database(light, "full")

    wrong_schema = seed_bakta(tmp_path / "schema", "full")
    metadata = json.loads((wrong_schema / "version.json").read_text())
    metadata["major"] = 5
    (wrong_schema / "version.json").write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(CacheIntegrityError, match="schema 6"):
        validate_bakta_database(wrong_schema, "full")

    partial = seed_bakta(tmp_path / "partial", "full")
    (partial / "pfam.h3f").unlink()
    with pytest.raises(CacheIntegrityError, match="pfam.h3f"):
        ensure_bakta_database(partial, "full", tmp_path / "partial.log")


def test_bakta_light_install_uses_native_db_light_layout(tmp_path, monkeypatch):
    database = tmp_path / "cache" / "db"
    calls = []

    def fake_run(argv, log_path, **kwargs):
        calls.append([str(value) for value in argv])
        output = Path(argv[argv.index("--output") + 1])
        seed_bakta(output / "db-light", "light")

    monkeypatch.setattr(bakta_module, "run_logged", fake_run)
    result = ensure_bakta_database(database, "light", tmp_path / "bakta.log")

    assert result.status == "prepared"
    assert calls == [
        ["bakta_db", "download", "--type", "light", "--output", calls[0][-1]]
    ]
    assert "--force" not in calls[0]
    assert "--update" not in calls[0]
    assert validate_bakta_database(database, "light")["type"] == "light"


@pytest.mark.parametrize("database_type", ["full", "light"])
def test_bakta_retries_observed_download_marker_and_preserves_layout(
    tmp_path, monkeypatch, database_type
):
    calls = []
    sleeps = []

    def fake_run(argv, log_path, **kwargs):
        calls.append([str(value) for value in argv])
        if len(calls) == 1:
            append_attempt_output(
                log_path,
                f"{bakta_module.BAKTA_DOWNLOAD_FAILURE}\n"
                "url=https://zenodo.org/record/403123/files/db.tar.xz, "
                "path=/cache/project404/db.tar.xz",
            )
            raise subprocess.CalledProcessError(1, argv)
        output = Path(argv[argv.index("--output") + 1])
        seed_bakta(
            output / ("db" if database_type == "full" else "db-light"),
            database_type,
        )

    monkeypatch.setattr(bakta_module, "run_logged", fake_run)
    monkeypatch.setattr(retry_module.time, "sleep", sleeps.append)
    database = tmp_path / f"bakta-{database_type}"

    result = ensure_bakta_database(database, database_type, tmp_path / "bakta.log")

    assert result.status == "prepared"
    assert len(calls) == 2
    assert sleeps == [10]
    assert validate_bakta_database(database, database_type)["type"] == database_type


@pytest.mark.parametrize(
    ("previous_output", "current_output"),
    [
        (bakta_module.BAKTA_DOWNLOAD_FAILURE, "unclassified native failure"),
        (
            "",
            f"{bakta_module.BAKTA_DOWNLOAD_FAILURE}\nchecksum mismatch",
        ),
        (
            "",
            f"{bakta_module.BAKTA_DOWNLOAD_FAILURE}\nHTTP Error 404: Not Found",
        ),
    ],
    ids=[
        "earlier-marker-isolation",
        "fatal-marker-precedence",
        "permanent-http-status",
    ],
)
def test_bakta_does_not_retry_unclassified_or_fatal_current_failure(
    tmp_path, monkeypatch, previous_output, current_output
):
    log_path = tmp_path / "bakta.log"
    log_path.write_text(previous_output + "\n", encoding="utf-8")
    calls = []
    sleeps = []
    error = subprocess.CalledProcessError(1, ["bakta_db"])

    def fail(argv, selected_log, **kwargs):
        calls.append(argv)
        append_attempt_output(selected_log, current_output)
        raise error

    monkeypatch.setattr(bakta_module, "run_logged", fail)
    monkeypatch.setattr(retry_module.time, "sleep", sleeps.append)
    with pytest.raises(subprocess.CalledProcessError) as raised:
        ensure_bakta_database(tmp_path / "bakta", "full", log_path)

    assert raised.value is error
    assert len(calls) == 1
    assert sleeps == []


def test_bakta_missing_executable_is_not_retried(tmp_path, monkeypatch):
    calls = []
    sleeps = []
    error = FileNotFoundError("bakta_db")

    def fail(*args, **kwargs):
        calls.append(1)
        raise error

    monkeypatch.setattr(bakta_module, "run_logged", fail)
    monkeypatch.setattr(retry_module.time, "sleep", sleeps.append)
    with pytest.raises(FileNotFoundError) as raised:
        ensure_bakta_database(tmp_path / "bakta", "full", tmp_path / "bakta.log")

    assert raised.value is error
    assert calls == [1]
    assert sleeps == []


def test_bakta_failed_install_does_not_publish_payload(tmp_path, monkeypatch):
    database = tmp_path / "db"

    def fail(*args, **kwargs):
        raise subprocess.CalledProcessError(7, ["bakta_db"])

    monkeypatch.setattr(bakta_module, "run_logged", fail)
    with pytest.raises(subprocess.CalledProcessError):
        ensure_bakta_database(database, "full", tmp_path / "bakta.log")
    assert not database.exists()


def test_concurrent_bakta_first_install_runs_downloader_once(tmp_path, monkeypatch):
    database = tmp_path / "db"
    calls = tmp_path / "calls"

    def fake_run(argv, log_path, **kwargs):
        with calls.open("a", encoding="utf-8") as handle:
            handle.write("download\n")
            handle.flush()
            os.fsync(handle.fileno())
        time.sleep(0.2)
        output = Path(argv[argv.index("--output") + 1])
        seed_bakta(output / "db-light", "light")

    monkeypatch.setattr(bakta_module, "run_logged", fake_run)
    context = multiprocessing.get_context("fork")
    queue = context.Queue()
    processes = [
        context.Process(
            target=_ensure_bakta_worker,
            args=(database, tmp_path / f"job-{index}.log", queue),
        )
        for index in range(2)
    ]
    for process in processes:
        process.start()
    for process in processes:
        process.join(10)
        assert process.exitcode == 0

    assert sorted(queue.get(timeout=1) for _ in processes) == ["prepared", "reused"]
    assert calls.read_text(encoding="utf-8").splitlines() == ["download"]


def test_gtdbtk_requires_r232_and_representative_layout(tmp_path):
    ready = seed_gtdbtk(tmp_path / "ready")
    assert validate_gtdbtk_database(ready) == {"release": "r232"}

    wrong = seed_gtdbtk(tmp_path / "wrong", "r220")
    with pytest.raises(CacheIntegrityError, match="requires VERSION_DATA=r232"):
        validate_gtdbtk_database(wrong)

    (ready / "markers" / "tigrfam" / "tigrfam.hmm").unlink()
    with pytest.raises(CacheIntegrityError, match="tigrfam.hmm"):
        validate_gtdbtk_database(ready)


def test_funannotate_all_semantics_ready_reuse_and_partial_failure(
    tmp_path, monkeypatch
):
    assert (
        expand_funannotate2_databases(["merops", "all"]) == FUNANNOTATE2_ALL_DATABASES
    )
    root = tmp_path / "funannotate"
    info = {
        database: seed_funannotate_entry(root, database)
        for database in FUNANNOTATE2_ALL_DATABASES
    }
    (root / "funannotate-db-info.json").write_text(json.dumps(info), encoding="utf-8")
    before = cache_snapshot(root)

    monkeypatch.setattr(
        "common.funannotate2_database.run_logged",
        lambda *args, **kwargs: pytest.fail("ready selection invoked installer"),
    )
    result = ensure_funannotate2_database(root, ["merops", "all"], tmp_path / "log")
    assert result.status == "reused"
    assert cache_snapshot(root) == before

    Path(info["merops"]["db"]).unlink()
    with pytest.raises(CacheIntegrityError, match="merops.*missing"):
        ensure_funannotate2_database(root, ["merops"], tmp_path / "partial.log")


def test_funannotate_first_and_additive_installs_rebase_paths_and_preserve_entries(
    tmp_path, monkeypatch
):
    calls = []

    def fake_run(argv, log_path, *, env=None, **kwargs):
        calls.append([str(value) for value in argv])
        staging = Path(env["FUNANNOTATE2_DB"])
        info = {}
        for database in argv[argv.index("-d") + 1 :]:
            info[database] = seed_funannotate_entry(staging, database)
        (staging / "funannotate-db-info.json").write_text(
            json.dumps(info), encoding="utf-8"
        )

    monkeypatch.setattr(funannotate2_module, "run_logged", fake_run)
    root = tmp_path / "funannotate"
    result = ensure_funannotate2_database(root, ["merops"], tmp_path / "fresh.log")
    assert result.status == "prepared"
    installed = json.loads((root / "funannotate-db-info.json").read_text())
    assert installed["merops"]["db"] == str((root / "merops.db").resolve())

    pretrained = root / "pretrained" / "model.json"
    pretrained.parent.mkdir()
    pretrained.write_text("model\n", encoding="utf-8")
    before_merops = (root / "merops.db").read_bytes()
    result = ensure_funannotate2_database(
        root, ["merops", "uniprot"], tmp_path / "additive.log"
    )
    assert result.prepared == ("uniprot",)
    installed = json.loads((root / "funannotate-db-info.json").read_text())
    assert set(installed) == {"merops", "uniprot"}
    assert (root / "merops.db").read_bytes() == before_merops
    assert pretrained.read_text(encoding="utf-8") == "model\n"
    assert calls[1][calls[1].index("-d") + 1 :] == ["uniprot"]


def test_funannotate_additive_install_preserves_dangling_symlink_collision(
    tmp_path, monkeypatch
):
    root = tmp_path / "funannotate"
    merops = seed_funannotate_entry(root, "merops")
    info_path = root / "funannotate-db-info.json"
    info_path.write_text(json.dumps({"merops": merops}), encoding="utf-8")
    dangling = root / "uniprot.db"
    dangling.symlink_to("missing-user-material")
    before = {
        "merops": (root / "merops.db").read_bytes(),
        "metadata": (info_path.read_bytes(), info_path.stat().st_mtime_ns),
        "link": os.readlink(dangling),
    }
    calls = []

    def fake_run(argv, log_path, *, env=None, **kwargs):
        calls.append([str(value) for value in argv])
        staging = Path(env["FUNANNOTATE2_DB"])
        record = seed_funannotate_entry(staging, "uniprot")
        (staging / "funannotate-db-info.json").write_text(
            json.dumps({"uniprot": record}), encoding="utf-8"
        )

    monkeypatch.setattr(funannotate2_module, "run_logged", fake_run)
    with pytest.raises(CacheIntegrityError, match="collides.*uniprot.db"):
        ensure_funannotate2_database(
            root, ["merops", "uniprot"], tmp_path / "additive.log"
        )

    assert len(calls) == 1
    assert dangling.is_symlink()
    assert os.readlink(dangling) == before["link"]
    assert (root / "merops.db").read_bytes() == before["merops"]
    assert (info_path.read_bytes(), info_path.stat().st_mtime_ns) == before["metadata"]


def test_funannotate_failed_install_does_not_publish(tmp_path, monkeypatch):
    root = tmp_path / "funannotate"
    monkeypatch.setattr(
        funannotate2_module,
        "run_logged",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            subprocess.CalledProcessError(8, ["funannotate2"])
        ),
    )
    with pytest.raises(subprocess.CalledProcessError):
        ensure_funannotate2_database(root, ["merops"], tmp_path / "failed.log")
    assert not root.exists()


def test_gtdbtk_missing_cache_prepares_once_and_ready_reuse_skips_download(
    tmp_path, monkeypatch
):
    calls = []

    def fake_run(argv, log_path, **kwargs):
        calls.append([str(value) for value in argv])
        if argv[0] == "wget":
            Path(argv[argv.index("-O") + 1]).write_bytes(b"archive")
        else:
            seed_gtdbtk(Path(argv[argv.index("-C") + 1]))

    monkeypatch.setattr(gtdbtk_module, "run_logged", fake_run)
    root = tmp_path / "gtdbtk"
    result = ensure_gtdbtk_database(root, tmp_path / "gtdbtk.log")
    assert result.status == "prepared"
    assert calls[0][1:3] == ["--tries=1", gtdbtk_module.GTDBTK_DATABASE_URL]
    ensure_gtdbtk_database(root, tmp_path / "reuse.log")
    assert len(calls) == 2


def test_gtdbtk_retries_only_wget_and_extracts_once(tmp_path, monkeypatch):
    wget_calls = []
    tar_calls = []
    sleeps = []

    def fake_run(argv, log_path, **kwargs):
        if argv[0] == "wget":
            wget_calls.append([str(value) for value in argv])
            if len(wget_calls) == 1:
                append_attempt_output(log_path, "Temporary failure in name resolution")
                raise subprocess.CalledProcessError(4, argv)
            Path(argv[argv.index("-O") + 1]).write_bytes(b"archive")
        else:
            tar_calls.append([str(value) for value in argv])
            seed_gtdbtk(Path(argv[argv.index("-C") + 1]))

    monkeypatch.setattr(gtdbtk_module, "run_logged", fake_run)
    monkeypatch.setattr(retry_module.time, "sleep", sleeps.append)
    result = ensure_gtdbtk_database(tmp_path / "gtdbtk", tmp_path / "gtdbtk.log")

    assert result.status == "prepared"
    assert len(wget_calls) == 2
    assert all("--tries=1" in call for call in wget_calls)
    assert len(tar_calls) == 1
    assert sleeps == [10]


def test_gtdbtk_retry_exhaustion_never_runs_extraction(tmp_path, monkeypatch):
    errors = [subprocess.CalledProcessError(8, ["wget"], output=None) for _ in range(4)]
    calls = []
    tar_calls = []
    sleeps = []

    def fail(argv, log_path, **kwargs):
        if argv[0] != "wget":
            tar_calls.append(argv)
            return
        error = errors[len(calls)]
        calls.append(error)
        append_attempt_output(log_path, "ERROR 503: Service Unavailable")
        raise error

    monkeypatch.setattr(gtdbtk_module, "run_logged", fail)
    monkeypatch.setattr(retry_module.time, "sleep", sleeps.append)
    with pytest.raises(subprocess.CalledProcessError) as raised:
        ensure_gtdbtk_database(tmp_path / "gtdbtk", tmp_path / "gtdbtk.log")

    assert raised.value is errors[-1]
    assert calls == errors
    assert tar_calls == []
    assert sleeps == [10, 30, 60]


@pytest.mark.parametrize(
    ("returncode", "output"),
    [
        (8, "ERROR 404: Not Found"),
        (4, "certificate verification failed"),
        (8, "unclassified server response"),
    ],
    ids=["permanent-http", "fatal-precedence", "unclassified"],
)
def test_gtdbtk_permanent_and_unclassified_failures_are_immediate(
    tmp_path, monkeypatch, returncode, output
):
    calls = []
    sleeps = []
    error = subprocess.CalledProcessError(returncode, ["wget"])

    def fail(argv, log_path, **kwargs):
        calls.append(argv)
        append_attempt_output(log_path, output)
        raise error

    monkeypatch.setattr(gtdbtk_module, "run_logged", fail)
    monkeypatch.setattr(retry_module.time, "sleep", sleeps.append)
    with pytest.raises(subprocess.CalledProcessError) as raised:
        ensure_gtdbtk_database(tmp_path / "gtdbtk", tmp_path / "gtdbtk.log")

    assert raised.value is error
    assert len(calls) == 1
    assert sleeps == []


def test_metaeuk_managed_siblings_legacy_reuse_and_mismatch(tmp_path, monkeypatch):
    calls = []

    def fake_run(argv, log_path, **kwargs):
        calls.append([str(value) for value in argv])
        seed_metaeuk(Path(argv[3]))

    monkeypatch.setattr(metaeuk_module, "run_logged", fake_run)
    first = tmp_path / "catalogs" / "one" / "db"
    second = tmp_path / "catalogs" / "two" / "db"
    ensure_metaeuk_database(
        first, "Catalog One", managed=True, threads=2, log_path=tmp_path / "one.log"
    )
    first_before = cache_snapshot(first.parent)
    ensure_metaeuk_database(
        second, "Catalog Two", managed=True, threads=2, log_path=tmp_path / "two.log"
    )
    assert cache_snapshot(first.parent) == first_before
    assert [call[2] for call in calls] == ["Catalog One", "Catalog Two"]

    ensure_metaeuk_database(
        first,
        "Catalog One",
        managed=True,
        threads=2,
        log_path=tmp_path / "reuse.log",
    )
    assert len(calls) == 2

    legacy = seed_metaeuk(tmp_path / "legacy")
    legacy_before = cache_snapshot(tmp_path)
    result = ensure_metaeuk_database(
        legacy,
        "User Asserted",
        managed=False,
        threads=1,
        log_path=tmp_path / "legacy.log",
    )
    assert result.identity == "asserted_unverified"
    assert (
        cache_snapshot(legacy.parent)[Path("legacy")] == legacy_before[Path("legacy")]
    )
    assert not (legacy.parent / METAEUK_IDENTITY_NAME).exists()

    with pytest.raises(CacheIntegrityError, match="identity mismatch"):
        ensure_metaeuk_database(
            first,
            "Wrong Catalog",
            managed=True,
            threads=1,
            log_path=tmp_path / "wrong.log",
        )


def test_metaeuk_rejects_lone_dbtype_and_missing_external_prefix(tmp_path):
    prefix = tmp_path / "catalog"
    Path(f"{prefix}.dbtype").write_text("fixture\n", encoding="utf-8")
    with pytest.raises(CacheIntegrityError, match="native family"):
        validate_metaeuk_native_database(prefix)
    with pytest.raises(FileNotFoundError, match="will not prepare"):
        ensure_metaeuk_database(
            tmp_path / "missing",
            "Catalog",
            managed=False,
            threads=1,
            log_path=tmp_path / "log",
        )


def test_checkm2_direct_file_guard_and_download_switch(tmp_path, monkeypatch):
    direct = tmp_path / "direct.dmnd"
    direct.write_bytes(b"database\n")
    selected, root = selected_checkm2_database(direct)
    assert selected == direct.resolve()
    assert root is None
    monkeypatch.setattr(
        "common.qa_databases.run_logged",
        lambda *args, **kwargs: pytest.fail("valid direct file triggered download"),
    )
    database, result = ensure_checkm2_database(direct, True, tmp_path / "log")
    assert database == direct.resolve()
    assert result.status == "reused"

    missing_direct = tmp_path / "missing.dmnd"
    with pytest.raises(FileNotFoundError, match="will not create a directory"):
        ensure_checkm2_database(missing_direct, True, tmp_path / "missing.log")
    assert not missing_direct.exists()

    with pytest.raises(FileNotFoundError, match="download_database"):
        ensure_checkm2_database(tmp_path / "root", False, tmp_path / "disabled.log")

    monkeypatch.setattr(
        qa_module,
        "run_logged",
        lambda argv, log_path, **kwargs: (
            Path(argv[argv.index("--path") + 1])
            / "CheckM2_database"
            / "uniref100.KO.1.dmnd"
        ).parent.mkdir(parents=True)
        or (
            Path(argv[argv.index("--path") + 1])
            / "CheckM2_database"
            / "uniref100.KO.1.dmnd"
        ).write_bytes(b"database\n"),
    )
    prepared, result = ensure_checkm2_database(
        tmp_path / "downloaded", True, tmp_path / "download.log"
    )
    assert result.status == "prepared"
    assert prepared.read_bytes() == b"database\n"


@pytest.mark.parametrize(
    "failure_output",
    [
        "Error during download: (Connection broken: IncompleteRead(1048576 bytes "
        "read, 123456 more expected), IncompleteRead(1048576 bytes read, 123456 "
        "more expected))\ncheckm2.zenodo_backpack.ZenodoConnectionException",
        "requests.exceptions.ChunkedEncodingError: ('Connection broken: "
        "IncompleteRead(4096 bytes read, 2048 more expected)', "
        "IncompleteRead(4096 bytes read, 2048 more expected))",
    ],
    ids=["native-incomplete-read", "requests-chunked-encoding"],
)
def test_checkm2_retries_identifiable_network_failure(
    tmp_path, monkeypatch, failure_output
):
    calls = []
    sleeps = []

    def fake_run(argv, log_path, **kwargs):
        calls.append([str(value) for value in argv])
        if len(calls) == 1:
            append_attempt_output(log_path, failure_output)
            raise subprocess.CalledProcessError(1, argv)
        database = (
            Path(argv[argv.index("--path") + 1])
            / "CheckM2_database"
            / "uniref100.KO.1.dmnd"
        )
        database.parent.mkdir(parents=True)
        database.write_bytes(b"database\n")

    monkeypatch.setattr(qa_module, "run_logged", fake_run)
    monkeypatch.setattr(retry_module.time, "sleep", sleeps.append)
    database, result = ensure_checkm2_database(
        tmp_path / "checkm2", True, tmp_path / "checkm2.log"
    )

    assert result.status == "prepared"
    assert database.read_bytes() == b"database\n"
    assert len(calls) == 2
    assert sleeps == [10]


@pytest.mark.parametrize(
    "output",
    [
        "Connection error: timed out\nChecksum is incorrect",
        "invalid native command usage",
    ],
    ids=["fatal-marker-precedence", "unclassified"],
)
def test_checkm2_fatal_and_unclassified_failures_are_immediate(
    tmp_path, monkeypatch, output
):
    calls = []
    sleeps = []
    error = subprocess.CalledProcessError(1, ["checkm2"])

    def fail(argv, log_path, **kwargs):
        calls.append(argv)
        append_attempt_output(log_path, output)
        raise error

    monkeypatch.setattr(qa_module, "run_logged", fail)
    monkeypatch.setattr(retry_module.time, "sleep", sleeps.append)
    with pytest.raises(subprocess.CalledProcessError) as raised:
        ensure_checkm2_database(tmp_path / "checkm2", True, tmp_path / "checkm2.log")

    assert raised.value is error
    assert len(calls) == 1
    assert sleeps == []


def test_eukcc_native_validation_and_download_switch(tmp_path, monkeypatch):
    ready = seed_eukcc(tmp_path / "eukcc")
    assert "files" in validate_eukcc_database(ready)
    monkeypatch.setattr(
        qa_module.urllib.request,
        "urlretrieve",
        lambda *args, **kwargs: pytest.fail("ready EukCC cache triggered transfer"),
    )
    assert ensure_eukcc_database(ready, True, tmp_path / "ready.log").status == "reused"

    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(CacheIntegrityError, match="db_base/db_info.json"):
        ensure_eukcc_database(empty, False, tmp_path / "empty.log")

    with pytest.raises(FileNotFoundError, match="download_database"):
        ensure_eukcc_database(tmp_path / "missing", False, tmp_path / "missing.log")

    archive_source = tmp_path / "eukcc.tar.gz"
    archive_root = tmp_path / "archive-source" / "eukcc2_db_ver_1.2"
    seed_eukcc(archive_root)
    with tarfile.open(archive_source, "w:gz") as archive:
        archive.add(archive_root, arcname="eukcc2_db_ver_1.2")

    monkeypatch.setattr(
        qa_module.urllib.request,
        "urlretrieve",
        lambda url, destination: shutil.copy2(archive_source, destination),
    )
    downloaded = tmp_path / "downloaded-eukcc"
    result = ensure_eukcc_database(downloaded, True, tmp_path / "download.log")
    assert result.status == "prepared"
    validate_eukcc_database(downloaded)


def test_eukcc_retries_transfer_without_repeating_extraction(tmp_path, monkeypatch):
    transfers = []
    extractions = []
    sleeps = []

    def fake_retrieve(url, destination):
        transfers.append((url, destination))
        if len(transfers) == 1:
            raise urllib.error.ContentTooShortError(
                "incomplete database transfer", b"partial"
            )
        Path(destination).write_bytes(b"archive")

    def fake_extract(archive, destination):
        extractions.append((archive, destination))
        seed_eukcc(Path(destination) / qa_module.EUKCC_ARCHIVE_ROOT)

    monkeypatch.setattr(qa_module.urllib.request, "urlretrieve", fake_retrieve)
    monkeypatch.setattr(qa_module, "_safe_extract", fake_extract)
    monkeypatch.setattr(retry_module.time, "sleep", sleeps.append)
    result = ensure_eukcc_database(tmp_path / "eukcc", True, tmp_path / "eukcc.log")

    assert result.status == "prepared"
    assert len(transfers) == 2
    assert len(extractions) == 1
    assert sleeps == [10]


@pytest.mark.parametrize(
    "error",
    [
        urllib.error.HTTPError("https://example.invalid", 404, "Not Found", None, None),
        urllib.error.URLError(
            ssl.SSLCertVerificationError(1, "certificate verify failed")
        ),
        PermissionError("local destination is not writable"),
        ValueError("unclassified failure"),
    ],
    ids=["permanent-http", "tls-certificate", "local-write", "unclassified"],
)
def test_eukcc_permanent_local_and_unclassified_failures_are_immediate(
    tmp_path, monkeypatch, error
):
    calls = []
    sleeps = []

    def fail(url, destination):
        calls.append((url, destination))
        raise error

    monkeypatch.setattr(qa_module.urllib.request, "urlretrieve", fail)
    monkeypatch.setattr(retry_module.time, "sleep", sleeps.append)
    with pytest.raises(type(error)) as raised:
        ensure_eukcc_database(tmp_path / "eukcc", True, tmp_path / "eukcc.log")

    assert raised.value is error
    assert len(calls) == 1
    assert sleeps == []
