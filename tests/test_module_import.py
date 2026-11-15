import json
import os
import stat
import subprocess
import time
from pathlib import Path

import pytest
import yaml

from common.bakta_database import BAKTA_COMMON_FILES
from helpers import material_snapshot

ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "workflow" / "Snakefile"


def seed_bakta(database):
    database.mkdir(parents=True)
    (database / "version.json").write_text(
        json.dumps({"date": "2026-01-01", "major": 6, "minor": 0, "type": "full"}),
        encoding="utf-8",
    )
    for name in (*BAKTA_COMMON_FILES, "psc.dmnd"):
        (database / name).write_text("fixture\n", encoding="utf-8")


@pytest.mark.parametrize(
    ("context", "expected_rules", "excluded_rules"),
    [
        (
            {
                "mode": "metagenome",
                "sample": "imported_project",
                "prok_bins_dir": "future/prok_bins",
                "euk_bins_dir": "future/euk_bins",
            },
            {"imported_qa_selection", "imported_metagenome_annotation_inventory"},
            {"imported_bakta", "imported_funannotate2"},
        ),
        (
            {
                "mode": "isolate",
                "sample": "imported_prok",
                "domain": "prok",
                "path": "future/assembly-prok.fasta",
            },
            {"imported_qa_selection", "imported_bakta", "imported_protein_summary"},
            {"imported_funannotate2", "imported_metagenome_annotation_inventory"},
        ),
        (
            {
                "mode": "isolate",
                "sample": "imported_euk",
                "domain": "euk",
                "path": "future/assembly-euk.fasta",
            },
            {"imported_qa_selection", "imported_funannotate2"},
            {"imported_bakta", "imported_metagenome_annotation_inventory"},
        ),
    ],
    ids=["metagenome", "prok", "euk"],
)
def test_module_import_uses_future_inputs_without_standalone_sample_sheet(
    tmp_path, context, expected_rules, excluded_rules
):
    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    config.pop("input")
    config["funannotate2"]["species"] = "Imported fixture"
    config["import_context"] = {"enabled": True, **context}
    harness = tmp_path / "Snakefile"
    harness.write_text(
        "MODULE_CONFIG = " + repr(config) + "\n\n"
        "module annotation:\n"
        f"    snakefile: {json.dumps(str(SNAKEFILE))}\n"
        "    config: MODULE_CONFIG\n"
        "    skip_validation: True\n\n"
        "use rule * from annotation exclude all as imported_*\n",
        encoding="utf-8",
    )

    result = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(harness),
            "--directory",
            str(tmp_path),
            "--list-rules",
            "--nolock",
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    rules = set(result.stdout.splitlines())
    assert expected_rules <= rules
    assert rules.isdisjoint(excluded_rules)
    assert not (tmp_path / "config" / "isolate_samples.tsv").exists()
    assert not (tmp_path / "config" / "metagenome_samples.tsv").exists()
    for path in context.values():
        if isinstance(path, str) and path.startswith("future/"):
            assert not (tmp_path / path).exists()


def test_imported_adapter_loads_cache_and_retry_helpers_from_module_source(
    tmp_path, monkeypatch, offline_guard
):
    fasta = tmp_path / "future" / "assembly.fasta"
    fasta.parent.mkdir()
    fasta.write_text(">contig\nACGT\n", encoding="utf-8")
    database = tmp_path / "cache" / "bakta"
    seed_bakta(database)

    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    config.pop("input")
    config["qa"]["enabled"] = False
    config["qa_filter"]["enabled"] = False
    config["bakta"]["db"] = str(database)
    config["import_context"] = {
        "enabled": True,
        "mode": "isolate",
        "sample": "imported",
        "domain": "prok",
        "path": str(fasta),
    }
    harness = tmp_path / "Snakefile"
    harness.write_text(
        "MODULE_CONFIG = " + repr(config) + "\n\n"
        "module annotation:\n"
        "    pathvars:\n"
        "        results='results/imported',\n"
        f"    snakefile: {json.dumps(str(SNAKEFILE))}\n"
        "    config: MODULE_CONFIG\n"
        "    skip_validation: True\n\n"
        "use rule * from annotation exclude all as imported_*\n",
        encoding="utf-8",
    )
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    calls = tmp_path / "bakta-calls.txt"
    fake_bakta = fake_bin / "bakta"
    fake_bakta.write_text(
        "#!/usr/bin/env python3\n"
        "import os\n"
        "from pathlib import Path\n"
        "with Path(os.environ['FAKE_BAKTA_CALLS']).open('a') as handle:\n"
        "    handle.write('annotation\\n')\n"
        "Path('annotation.faa').write_text('>protein\\nM\\n')\n",
        encoding="utf-8",
    )
    fake_bakta.chmod(fake_bakta.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", f"{fake_bin}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_BAKTA_CALLS", str(calls))

    command = [
        "snakemake",
        "--snakefile",
        str(harness),
        "--directory",
        str(tmp_path),
        "--cores",
        "1",
        "--nolock",
        "results/imported/bakta",
    ]
    result = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        env=os.environ.copy(),
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert calls.read_text(encoding="utf-8").splitlines() == ["annotation"]
    assert (tmp_path / "results" / "imported" / "bakta" / "annotation.faa").is_file()

    material = database / "bakta.db"
    material.write_bytes(b"changed-imported-compatible-cache\n")
    current = material.stat()
    changed_mtime = time.time_ns()
    os.utime(
        material,
        ns=(current.st_atime_ns, changed_mtime),
    )
    changed_cache = material_snapshot(database)
    rerun = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        env=os.environ.copy(),
    )

    assert rerun.returncode == 0, rerun.stdout + rerun.stderr
    assert str(material.resolve()) in rerun.stdout + rerun.stderr
    assert calls.read_text(encoding="utf-8").splitlines() == [
        "annotation",
        "annotation",
    ]
    assert material_snapshot(database) == changed_cache
    assert not offline_guard["attempt_log"].exists()


def test_imported_overlap_rejection_precedes_forced_output_cleanup(tmp_path):
    fasta = tmp_path / "future" / "assembly.fasta"
    fasta.parent.mkdir()
    fasta.write_text(">contig\nACGT\n", encoding="utf-8")
    output = tmp_path / "results" / "imported" / "bakta"
    database = output / "persistent-cache"
    seed_bakta(database)
    (database / "sentinel.bin").write_bytes(b"preserve-imported-cache\n")
    (database / "dangling-native-link").symlink_to("missing-user-target")
    before = material_snapshot(database)

    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    config.pop("input")
    config["qa"]["enabled"] = False
    config["qa_filter"]["enabled"] = False
    config["bakta"]["db"] = str(database)
    config["import_context"] = {
        "enabled": True,
        "mode": "isolate",
        "sample": "imported",
        "domain": "prok",
        "path": str(fasta),
    }
    harness = tmp_path / "overlap-import.smk"
    harness.write_text(
        "MODULE_CONFIG = " + repr(config) + "\n\n"
        "module annotation:\n"
        "    pathvars:\n"
        "        results='results/imported',\n"
        f"    snakefile: {json.dumps(str(SNAKEFILE))}\n"
        "    config: MODULE_CONFIG\n"
        "    skip_validation: True\n\n"
        "use rule * from annotation exclude all as imported_*\n",
        encoding="utf-8",
    )
    calls = tmp_path / "unexpected-imported-calls.txt"
    fake_bin = tmp_path / "overlap-fake-bin"
    fake_bin.mkdir()
    for name in ("bakta", "bakta_db"):
        executable = fake_bin / name
        executable.write_text(
            "#!/usr/bin/env python3\n"
            "import os\n"
            "from pathlib import Path\n"
            "with Path(os.environ['OVERLAP_CALLS']).open('a') as handle:\n"
            "    handle.write(Path(__file__).name + '\\n')\n"
            "raise SystemExit(98)\n",
            encoding="utf-8",
        )
        executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    environment = os.environ.copy()
    environment["PATH"] = f"{fake_bin}{os.pathsep}{environment['PATH']}"
    environment["OVERLAP_CALLS"] = str(calls)

    result = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(harness),
            "--directory",
            str(tmp_path),
            "--cores",
            "1",
            "--nolock",
            "--forceall",
            "results/imported/bakta",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )

    assert result.returncode != 0
    assert "overlaps Snakemake-owned analysis results" in result.stdout + result.stderr
    assert not calls.exists()
    assert database.is_dir()
    assert material_snapshot(database) == before
