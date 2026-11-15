import copy
import os
import subprocess
from pathlib import Path

import pytest
import yaml
from jsonschema import Draft202012Validator, ValidationError

from common.bakta_database import BAKTA_COMMON_FILES
from common.database_dependencies import (
    bakta_database_inputs,
    funannotate2_database_inputs,
    metaeuk_database_inputs,
)
from common.database_selection import (
    metaeuk_database_child,
    resolve_bakta_database,
    resolve_metaeuk_database_location,
)
from helpers import material_snapshot, rule_count

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = ROOT / "workflow" / "schemas" / "config.schema.yaml"
SNAKEFILE = ROOT / "workflow" / "Snakefile"


@pytest.mark.parametrize("value", [None, "", "   ", "\t\n"])
def test_metaeuk_omission_and_blank_values_use_uniref90(value, tmp_path):
    location = resolve_metaeuk_database_location(value, "", tmp_path)

    assert location.catalog == "UniRef90"
    assert location.prefix == str(tmp_path / "catalogs" / "uniref90" / "db")
    assert location.managed is True


def test_metaeuk_custom_source_is_trimmed_and_safely_encoded(tmp_path):
    location = resolve_metaeuk_database_location(
        "  UniProtKB/Swiss-Prot  ", "", tmp_path
    )

    assert location.catalog == "UniProtKB/Swiss-Prot"
    assert location.prefix == str(
        tmp_path / "catalogs" / "catalog-UniProtKB%2FSwiss-Prot" / "db"
    )
    assert location.managed is True
    assert Path(location.prefix).parents[2] == tmp_path
    assert len(Path(location.prefix).relative_to(tmp_path).parts) == 3


def test_distinct_metaeuk_sources_have_distinct_cache_children():
    sources = ["catalog/name", "catalog%2Fname", "catalog-name", "../catalog"]

    children = [metaeuk_database_child(source) for source in sources]

    assert len(set(children)) == len(sources)
    assert all(len(Path(child).parts) == 1 for child in children)


def test_explicit_metaeuk_prefix_overrides_source_derivation(tmp_path):
    explicit = tmp_path / "reused" / "database"

    location = resolve_metaeuk_database_location(
        "Custom/Catalog", str(explicit), tmp_path / "cache"
    )

    assert location.catalog == "Custom/Catalog"
    assert location.prefix == str(explicit)
    assert location.managed is False


def test_metaeuk_preserves_legacy_flat_prefix_and_prefers_managed_leaf(tmp_path):
    legacy = tmp_path / "uniref90"
    Path(f"{legacy}.dbtype").write_text("legacy\n", encoding="utf-8")

    legacy_location = resolve_metaeuk_database_location("UniRef90", "", tmp_path)
    assert legacy_location.catalog == "UniRef90"
    assert legacy_location.prefix == str(legacy)
    assert legacy_location.managed is False

    managed = tmp_path / "catalogs" / "uniref90" / "db"
    managed.parent.mkdir(parents=True)
    managed_location = resolve_metaeuk_database_location("UniRef90", "", tmp_path)
    assert managed_location.catalog == "UniRef90"
    assert managed_location.prefix == str(managed)
    assert managed_location.managed is True


def test_bakta_default_and_light_paths_are_distinct(tmp_path, monkeypatch):
    monkeypatch.delenv("BAKTA_DB", raising=False)
    monkeypatch.setenv("MEEW_RESOURCES", str(tmp_path))

    assert resolve_bakta_database(None, "") == (
        "full",
        str(tmp_path / "bakta_db" / "db"),
    )
    assert resolve_bakta_database("light", "") == (
        "light",
        str(tmp_path / "bakta_db" / "light" / "db"),
    )


def test_explicit_bakta_directory_overrides_type_and_environment(tmp_path, monkeypatch):
    explicit = tmp_path / "explicit"
    environment = tmp_path / "environment"
    monkeypatch.setenv("BAKTA_DB", str(environment))

    assert resolve_bakta_database("light", str(explicit)) == (
        "light",
        str(explicit),
    )
    assert resolve_bakta_database("light", "") == ("light", str(environment))


def test_database_selection_schema_contracts():
    schema = yaml.safe_load(SCHEMA.read_text(encoding="utf-8"))
    validator = Draft202012Validator(schema)
    shipped = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    validator.validate(shipped)
    assert shipped["bakta"]["type"] == "full"
    assert shipped["metaeuk"]["database"] == "UniRef90"
    assert shipped["upimapi"]["db"] == "uniprot"
    assert schema["properties"]["upimapi"]["properties"]["db"]["default"] == ("uniprot")

    for rule_file in ("prok.smk", "batch.smk"):
        rule_source = (ROOT / "workflow" / "rules" / rule_file).read_text(
            encoding="utf-8"
        )
        assert rule_source.count('.get("db", "uniprot")') == 2
        assert '.get("db", "swissprot")' not in rule_source

    omitted = copy.deepcopy(shipped)
    omitted["bakta"].pop("type")
    omitted["metaeuk"].pop("database")
    validator.validate(omitted)

    blank_metaeuk = copy.deepcopy(shipped)
    blank_metaeuk["metaeuk"]["database"] = "   "
    validator.validate(blank_metaeuk)

    invalid_bakta = copy.deepcopy(shipped)
    invalid_bakta["bakta"]["type"] = "compact"
    with pytest.raises(ValidationError):
        validator.validate(invalid_bakta)

    removed_stage_sheets = copy.deepcopy(shipped)
    removed_stage_sheets["stage_sheets"] = {"enabled": False}
    with pytest.raises(ValidationError):
        validator.validate(removed_stage_sheets)


def test_standalone_invalid_bakta_type_fails_workflow_schema_validation(workflow_case):
    case = workflow_case("isolate", enabled=False)
    case.config["bakta"]["type"] = "compact"
    case.config_path.write_text(
        yaml.safe_dump(case.config, sort_keys=False), encoding="utf-8"
    )

    result = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(SNAKEFILE),
            "--directory",
            str(case.root),
            "--configfile",
            str(case.config_path),
            "--list-rules",
            "--nolock",
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "ValidationError" in result.stderr
    assert "'compact' is not one of ['full', 'light']" in result.stderr


def test_database_cache_overlap_fails_before_forced_job_cleanup(
    workflow_case, tmp_path
):
    case = workflow_case("isolate", enabled=False)
    output_dir = case.root / "results" / case.sample / "bakta"
    database = output_dir / "persistent-bakta-cache"
    database.mkdir(parents=True)
    (database / "version.json").write_text(
        '{"date":"2026-01-01","major":6,"minor":1,"type":"full"}\n',
        encoding="utf-8",
    )
    for name in (*BAKTA_COMMON_FILES, "psc.dmnd"):
        (database / name).write_bytes(f"native:{name}\n".encode())
    (database / "sentinel.bin").write_bytes(b"preserve-me\n")
    (database / "dangling-native-link").symlink_to("missing-user-target")
    before = material_snapshot(database)

    case.config["bakta"]["db"] = str(database)
    case.config_path.write_text(
        yaml.safe_dump(case.config, sort_keys=False), encoding="utf-8"
    )

    fake_bin = tmp_path / "overlap-fake-bin"
    fake_bin.mkdir()
    calls = tmp_path / "overlap-calls.txt"
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
        executable.chmod(executable.stat().st_mode | 0o100)
    environment = os.environ.copy()
    environment["PATH"] = f"{fake_bin}{os.pathsep}{environment['PATH']}"
    environment["OVERLAP_CALLS"] = str(calls)

    result = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(SNAKEFILE),
            "--directory",
            str(case.root),
            "--configfile",
            str(case.config_path),
            "--cores",
            "1",
            "--nolock",
            "--forceall",
            "results/isolate1/bakta",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )

    assert result.returncode != 0
    assert "overlaps Snakemake-owned analysis results" in (
        result.stdout + result.stderr
    )
    assert not calls.exists()
    assert database.is_dir()
    assert material_snapshot(database) == before


def test_funannotate2_ready_native_metadata_is_a_read_only_standalone_input(
    workflow_case,
):
    case = workflow_case(
        "isolate", enabled=False, isolate_domain="euk", qa_enabled=False
    )
    database = case.root / "funannotate2_db"
    database.mkdir()
    material = database / "merops.dmnd"
    material.write_text("fixture\n", encoding="utf-8")
    info = database / "funannotate-db-info.json"
    info.write_text(
        '{"merops":{"db":"merops.dmnd","type":"fixture","version":"1",'
        '"date":"2026-01-01","n_records":1,"md5":"fixture"}}\n',
        encoding="utf-8",
    )
    case.config["funannotate2"]["databases"] = ["merops"]
    case.config_path.write_text(
        yaml.safe_dump(case.config, sort_keys=False), encoding="utf-8"
    )

    output = case.run(extra_arguments=("--rerun-triggers", "mtime"))
    assert rule_count(output, "funannotate2") == 1
    assert funannotate2_database_inputs(database, ["merops"]) == [
        str(info.resolve()),
        str(material.resolve()),
    ]


@pytest.mark.parametrize(
    ("database", "database_child", "bakta_type", "bakta_subpath"),
    [
        (None, "uniref90", None, "bakta_db/db"),
        (
            "  UniProtKB/Swiss-Prot  ",
            "catalog-UniProtKB%2FSwiss-Prot",
            "light",
            "bakta_db/light/db",
        ),
    ],
    ids=["normal-defaults", "alternate-catalog-and-light"],
)
def test_standalone_dag_resolves_default_and_alternate_database_resources(
    workflow_case,
    database,
    database_child,
    bakta_type,
    bakta_subpath,
):
    case = workflow_case(
        "metagenome",
        enabled=False,
        rows=[("prok_pass", "prok", "pass"), ("euk_pass", "euk", "pass")],
    )
    case.config["bakta"]["db"] = ""
    case.config["metaeuk"]["db"] = ""
    case.config["metaeuk"]["cache_dir"] = ""
    if bakta_type is None:
        case.config["bakta"].pop("type")
    else:
        case.config["bakta"]["type"] = bakta_type
    if database is None:
        case.config["metaeuk"].pop("database")
    else:
        case.config["metaeuk"]["database"] = database
    case.config_path.write_text(
        yaml.safe_dump(case.config, sort_keys=False), encoding="utf-8"
    )

    resources = case.root / "selected_resources"
    env = os.environ.copy()
    env.pop("BAKTA_DB", None)
    env.pop("METAEUK_DB", None)
    env["MEEW_RESOURCES"] = str(resources)
    bakta_database = resources / bakta_subpath
    bakta_database.mkdir(parents=True)
    (bakta_database / "version.json").write_text("{}\n", encoding="utf-8")
    (bakta_database / "bakta.db").write_text("fixture\n", encoding="utf-8")
    metaeuk_prefix = resources / "metaeuk_db" / "catalogs" / database_child / "db"
    metaeuk_prefix.parent.mkdir(parents=True)
    for material in (
        metaeuk_prefix,
        Path(f"{metaeuk_prefix}.index"),
        Path(f"{metaeuk_prefix}.dbtype"),
    ):
        material.write_text("fixture\n", encoding="utf-8")
    (metaeuk_prefix.parent / "catalog-identity.json").write_text(
        "{}\n", encoding="utf-8"
    )
    output = case.run(env=env, extra_arguments=("--rerun-triggers", "mtime"))

    assert rule_count(output, "bakta_bin") == 1
    assert rule_count(output, "metaeuk_bin") == 1
    bakta_material = str((bakta_database / "bakta.db").resolve())
    metaeuk_material = str(Path(f"{metaeuk_prefix}.dbtype").resolve())
    assert bakta_material in bakta_database_inputs(bakta_database, bakta_type or "full")
    assert metaeuk_material in metaeuk_database_inputs(metaeuk_prefix)


@pytest.mark.parametrize(
    ("source", "rule", "selector"),
    [
        ("prok.smk", "bakta", "bakta_database_inputs"),
        ("batch.smk", "bakta_bin", "bakta_database_inputs"),
        ("gtdbtk.smk", "gtdbtk", "gtdbtk_database_inputs"),
        ("euk.smk", "funannotate2", "funannotate2_database_inputs"),
        ("batch.smk", "metaeuk_bin", "metaeuk_database_inputs"),
    ],
)
def test_database_material_is_declared_only_as_consumer_input(source, rule, selector):
    rule_source = (ROOT / "workflow" / "rules" / source).read_text(encoding="utf-8")
    block = rule_source.split(f"rule {rule}:", 1)[1].split("\n\nrule ", 1)[0]
    inputs, outputs = block.split("    output:", 1)

    assert selector in inputs
    assert selector not in outputs
