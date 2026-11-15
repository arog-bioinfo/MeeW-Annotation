import copy
import json
import os
import subprocess
from pathlib import Path

import pytest
import yaml
from jsonschema import Draft202012Validator, ValidationError

from common.paths import resolve_resource_path

ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "workflow" / "Snakefile"
SCHEMA = ROOT / "workflow" / "schemas" / "config.schema.yaml"
RESOURCE_ENVIRONMENT = (
    "BAKTA_DB",
    "CHECKM2DB",
    "EUKCC2_DB",
    "EUKCC_DB",
    "FUNANNOTATE2_DB",
    "GTDBTK_DATA_PATH",
    "METAEUK_DB",
    "RECOGNIZER_RESOURCES",
    "UPIMAPI_RESOURCES",
    "MEEW_RESOURCES",
    "XDG_DATA_HOME",
)


def _clear_resource_environment(monkeypatch):
    for name in RESOURCE_ENVIRONMENT:
        monkeypatch.delenv(name, raising=False)


@pytest.mark.parametrize(
    ("selected", "expected_leaf"),
    [
        ("explicit", "explicit"),
        ("configured", "configured/checkm2_db"),
        ("tool", "tool"),
        ("meew", "meew/checkm2_db"),
        ("xdg", "xdg/meew/checkm2_db"),
        ("home", "home/resources/checkm2_db"),
    ],
)
def test_resource_path_precedence(selected, expected_leaf, tmp_path, monkeypatch):
    _clear_resource_environment(monkeypatch)
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    explicit = ""
    configured = ""
    if selected == "explicit":
        explicit = str(tmp_path / "explicit")
        configured = str(tmp_path / "configured")
        monkeypatch.setenv("CHECKM2DB", str(tmp_path / "tool"))
        monkeypatch.setenv("MEEW_RESOURCES", str(tmp_path / "meew"))
        monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    elif selected == "configured":
        configured = str(tmp_path / "configured")
        monkeypatch.setenv("CHECKM2DB", str(tmp_path / "tool"))
        monkeypatch.setenv("MEEW_RESOURCES", str(tmp_path / "meew"))
        monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    elif selected == "tool":
        monkeypatch.setenv("CHECKM2DB", str(tmp_path / "tool"))
        monkeypatch.setenv("MEEW_RESOURCES", str(tmp_path / "meew"))
        monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    elif selected == "meew":
        monkeypatch.setenv("MEEW_RESOURCES", str(tmp_path / "meew"))
        monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    elif selected == "xdg":
        monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))

    resolved = resolve_resource_path(explicit, "CHECKM2DB", "checkm2_db", configured)

    assert resolved == str(tmp_path / expected_leaf)
    assert not (tmp_path / expected_leaf).exists()


@pytest.mark.parametrize("configured", [None, "", "   ", "\t\n"])
def test_omitted_or_blank_root_preserves_tool_environment(
    configured, tmp_path, monkeypatch
):
    _clear_resource_environment(monkeypatch)
    tool_path = tmp_path / "tool-path"
    monkeypatch.setenv("CHECKM2DB", str(tool_path))

    assert resolve_resource_path("", "CHECKM2DB", "checkm2_db", configured) == str(
        tool_path
    )
    assert not tool_path.exists()


def test_configured_root_and_explicit_paths_expand_without_absolutizing(
    tmp_path, monkeypatch
):
    _clear_resource_environment(monkeypatch)
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("RESOURCE_SEGMENT", "shared resources")
    monkeypatch.setenv("EXPLICIT_SEGMENT", "explicit database")

    assert resolve_resource_path(
        "", "CHECKM2DB", "checkm2_db", "~/$RESOURCE_SEGMENT"
    ) == str(tmp_path / "home" / "shared resources" / "checkm2_db")
    assert resolve_resource_path(
        "~/$EXPLICIT_SEGMENT",
        "CHECKM2DB",
        "checkm2_db",
        str(tmp_path / "configured"),
    ) == str(tmp_path / "home" / "explicit database")
    assert (
        resolve_resource_path("", "CHECKM2DB", "checkm2_db", "relative/resources")
        == "relative/resources/checkm2_db"
    )


def test_resources_root_schema_is_optional_string():
    schema = yaml.safe_load(SCHEMA.read_text(encoding="utf-8"))
    validator = Draft202012Validator(schema)
    shipped = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    validator.validate(shipped)
    assert shipped["resources_root"] == ""

    omitted = copy.deepcopy(shipped)
    omitted.pop("resources_root")
    validator.validate(omitted)

    configured = copy.deepcopy(shipped)
    configured["resources_root"] = "~/shared/$RESOURCE_SET"
    validator.validate(configured)

    invalid = copy.deepcopy(shipped)
    invalid["resources_root"] = 42
    with pytest.raises(ValidationError):
        validator.validate(invalid)


@pytest.mark.parametrize(
    "explicit_override", [False, True], ids=["defaults", "override"]
)
def test_workflow_resolves_every_consumer_under_root_without_side_effects(
    workflow_case, tmp_path, monkeypatch, explicit_override
):
    case = workflow_case(
        "metagenome",
        enabled=False,
        rows=[("prok_pass", "prok", "pass"), ("euk_pass", "euk", "pass")],
        qa_enabled=True,
        filter_enabled=True,
    )
    resources = case.root / "shared-resources"
    explicit_recognizer = case.root / "explicit-recognizer"
    custom_upimapi = case.root / "custom-upimapi.fasta"
    case.config["resources_root"] = str(resources)
    case.config["qa"]["checkm2"]["database_path"] = ""
    case.config["qa"]["eukcc"]["database_path"] = ""
    for section, field in (
        ("bakta", "db"),
        ("gtdbtk", "data_dir"),
        ("funannotate2", "db_dir"),
        ("metaeuk", "cache_dir"),
        ("metaeuk", "db"),
        ("recognizer_prok", "resources_dir"),
        ("recognizer_euk", "resources_dir"),
        ("upimapi", "resources_dir"),
    ):
        case.config[section][field] = ""
    if explicit_override:
        case.config["recognizer_prok"]["resources_dir"] = str(explicit_recognizer)
    case.config["upimapi"]["db"] = ""
    case.config["upimapi"]["db_custom"] = "$CUSTOM_UPIMAPI"
    case.config_path.write_text(
        yaml.safe_dump(case.config, sort_keys=False), encoding="utf-8"
    )
    _clear_resource_environment(monkeypatch)
    monkeypatch.setenv("CUSTOM_UPIMAPI", str(custom_upimapi))

    harness = case.root / "inspect-resource-root.smk"
    harness.write_text(
        "import json\n\n"
        f"include: {json.dumps(str(SNAKEFILE))}\n\n"
        "print('RESOURCE_PATHS=' + json.dumps({\n"
        "    'bakta': BAKTA_DB_DIR,\n"
        "    'checkm2': config['qa']['checkm2']['database_path'],\n"
        "    'eukcc': config['qa']['eukcc']['database_path'],\n"
        "    'funannotate2': FUNANNOTATE2_DB_DIR,\n"
        "    'gtdbtk': GTDBTK_DATA_DIR,\n"
        "    'metaeuk_cache': METAEUK_CACHE_DIR,\n"
        "    'metaeuk_prefix': METAEUK_DB_PREFIX,\n"
        "    'recognizer_prok': config['recognizer_prok']['resources_dir'],\n"
        "    'recognizer_euk': config['recognizer_euk']['resources_dir'],\n"
        "    'upimapi': config['upimapi']['resources_dir'],\n"
        "    'upimapi_custom': config['upimapi']['db_custom'],\n"
        "}, sort_keys=True))\n",
        encoding="utf-8",
    )
    result = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(harness),
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
        env=os.environ.copy(),
    )

    assert result.returncode == 0, result.stdout + result.stderr
    payload = result.stdout.split("RESOURCE_PATHS=", 1)[1].splitlines()[0]
    paths = json.loads(payload)
    assert paths == {
        "bakta": str(resources / "bakta_db" / "db"),
        "checkm2": str(resources / "checkm2_db"),
        "eukcc": str(resources / "eukcc2_db_ver_1.2"),
        "funannotate2": str(resources / "funannotate2_db"),
        "gtdbtk": str(resources / "gtdbtk_db"),
        "metaeuk_cache": str(resources / "metaeuk_db"),
        "metaeuk_prefix": str(
            resources / "metaeuk_db" / "catalogs" / "uniref90" / "db"
        ),
        "recognizer_prok": str(
            explicit_recognizer if explicit_override else resources / "recognizer_db"
        ),
        "recognizer_euk": str(resources / "recognizer_db"),
        "upimapi": str(resources / "upimapi_db"),
        "upimapi_custom": str(custom_upimapi),
    }
    assert not resources.exists()
    assert not explicit_recognizer.exists()
    assert not custom_upimapi.exists()
