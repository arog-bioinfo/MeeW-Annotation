import json
import os
import shutil
import stat
import subprocess
import time
from pathlib import Path

import yaml

from common.bakta_database import BAKTA_COMMON_FILES
from helpers import file_snapshot as snapshot

ROOT = Path(__file__).resolve().parents[1]


def seed_bakta(path):
    path.mkdir(parents=True)
    (path / "version.json").write_text(
        json.dumps({"date": "2026-01-01", "major": 6, "minor": 1, "type": "full"}),
        encoding="utf-8",
    )
    for name in (*BAKTA_COMMON_FILES, "psc.dmnd"):
        (path / name).write_bytes(b"fixture\n")


def run_snakemake(workdir, environment, *arguments):
    result = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(workdir / "workflow" / "Snakefile"),
            "--directory",
            str(workdir),
            "--cores",
            "1",
            "--nolock",
            "--printshellcmds",
            "--software-deployment-method",
            "env-modules",
            *arguments,
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    return result, result.stdout + result.stderr


def test_default_provenance_and_forced_reruns_never_reinstall_ready_cache(tmp_path):
    workflow = tmp_path / "workflow"
    scripts = workflow / "scripts"
    scripts.mkdir(parents=True)
    shutil.copy2(ROOT / "workflow" / "scripts" / "bakta.py", scripts / "bakta.py")
    shutil.copytree(ROOT / "workflow" / "scripts" / "common", scripts / "common")
    (workflow / "Snakefile").write_text(
        "import sys\n"
        f"sys.path.insert(0, {str(scripts)!r})\n"
        "from common.database_dependencies import bakta_database_inputs\n\n"
        "configfile: 'config.yaml'\n\n"
        "rule annotate:\n"
        "    input:\n"
        "        fasta='input.fa',\n"
        "        database=lambda wc: bakta_database_inputs(config['db'], 'full'),\n"
        "    output:\n"
        "        results=directory('results/bakta'),\n"
        "    log:\n"
        "        'results/bakta.log',\n"
        "    envmodules:\n"
        "        config['environment_module']\n"
        "    threads: 1\n"
        "    params:\n"
        "        db=config['db'],\n"
        "        db_type='full',\n"
        "        prefix='annotation',\n"
        "        meta=config['meta'],\n"
        "        extra=config['extra'],\n"
        "    script:\n"
        "        'scripts/bakta.py'\n",
        encoding="utf-8",
    )
    (tmp_path / "input.fa").write_text(">contig\nACGT\n", encoding="utf-8")
    database = tmp_path / "cache" / "bakta"
    seed_bakta(database)
    config = {
        "db": str(database),
        "meta": False,
        "extra": "",
        "environment_module": "bakta/1.12.0-a",
    }
    (tmp_path / "config.yaml").write_text(
        yaml.safe_dump(config, sort_keys=False), encoding="utf-8"
    )

    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_module = fake_bin / "module"
    fake_module.write_text("#!/usr/bin/env sh\nexit 0\n", encoding="utf-8")
    fake_module.chmod(fake_module.stat().st_mode | stat.S_IXUSR)
    fake_bakta = fake_bin / "bakta"
    fake_bakta.write_text(
        "#!/usr/bin/env python3\n"
        "import os\n"
        "from pathlib import Path\n"
        "with Path(os.environ['FAKE_BAKTA_CALLS']).open('a') as handle:\n"
        "    handle.write('annotation\\n')\n"
        "if os.environ.get('FAKE_BAKTA_FAIL') == '1':\n"
        "    raise SystemExit(9)\n"
        "Path('annotation.faa').write_text('>protein\\nM\\n')\n",
        encoding="utf-8",
    )
    fake_bakta.chmod(fake_bakta.stat().st_mode | stat.S_IXUSR)
    fake_installer = fake_bin / "bakta_db"
    fake_installer.write_text(
        "#!/usr/bin/env python3\n"
        "import os\n"
        "from pathlib import Path\n"
        "with Path(os.environ['FAKE_INSTALLER_CALLS']).open('a') as handle:\n"
        "    handle.write('unexpected install\\n')\n"
        "raise SystemExit(97)\n",
        encoding="utf-8",
    )
    fake_installer.chmod(fake_installer.stat().st_mode | stat.S_IXUSR)
    calls = tmp_path / "bakta-calls.txt"
    installer_calls = tmp_path / "installer-calls.txt"
    environment = os.environ.copy()
    environment["PATH"] = f"{fake_bin}{os.pathsep}{environment['PATH']}"
    environment["FAKE_BAKTA_CALLS"] = str(calls)
    environment["FAKE_INSTALLER_CALLS"] = str(installer_calls)
    cache_before = snapshot(database)

    first, first_output = run_snakemake(tmp_path, environment)
    assert first.returncode == 0, first_output
    assert calls.read_text(encoding="utf-8").splitlines() == ["annotation"]

    unchanged, unchanged_output = run_snakemake(tmp_path, environment)
    assert unchanged.returncode == 0, unchanged_output
    assert calls.read_text(encoding="utf-8").splitlines() == ["annotation"]
    assert not installer_calls.exists()
    assert snapshot(database) == cache_before

    material = database / "bakta.db"
    material_stat = material.stat()
    material.write_bytes(b"changed-compatible-fixture\n")
    changed_mtime = time.time_ns()
    os.utime(
        material,
        ns=(material_stat.st_atime_ns, changed_mtime),
    )
    changed_cache = snapshot(database)
    database_change, database_change_output = run_snakemake(tmp_path, environment)
    assert database_change.returncode == 0, database_change_output
    assert str(material) in database_change_output
    assert calls.read_text(encoding="utf-8").splitlines() == [
        "annotation",
        "annotation",
    ]
    assert not installer_calls.exists()
    assert snapshot(database) == changed_cache

    config["meta"] = True
    (tmp_path / "config.yaml").write_text(
        yaml.safe_dump(config, sort_keys=False), encoding="utf-8"
    )
    mode, mode_output = run_snakemake(tmp_path, environment)
    assert mode.returncode == 0, mode_output
    assert "Params have changed" in mode_output

    config["extra"] = "--verbose"
    (tmp_path / "config.yaml").write_text(
        yaml.safe_dump(config, sort_keys=False), encoding="utf-8"
    )
    params, params_output = run_snakemake(tmp_path, environment)
    assert params.returncode == 0, params_output
    assert "Params have changed" in params_output

    with (scripts / "bakta.py").open("a", encoding="utf-8") as handle:
        handle.write("\n# lifecycle code-trigger fixture\n")
    code, code_output = run_snakemake(tmp_path, environment)
    assert code.returncode == 0, code_output
    assert "Code has changed" in code_output

    config["environment_module"] = "bakta/1.12.0-b"
    (tmp_path / "config.yaml").write_text(
        yaml.safe_dump(config, sort_keys=False), encoding="utf-8"
    )
    software, software_output = run_snakemake(tmp_path, environment)
    assert software.returncode == 0, software_output
    assert "Software environment definition has changed" in software_output

    forced, forced_output = run_snakemake(tmp_path, environment, "--forceall")
    assert forced.returncode == 0, forced_output
    assert len(calls.read_text(encoding="utf-8").splitlines()) == 7
    assert not installer_calls.exists()
    assert snapshot(database) == changed_cache

    failing_environment = dict(environment)
    failing_environment["FAKE_BAKTA_FAIL"] = "1"
    failed, failed_output = run_snakemake(tmp_path, failing_environment, "--forceall")
    assert failed.returncode != 0, failed_output
    assert snapshot(database) == changed_cache
    assert not installer_calls.exists()
