import csv
from copy import deepcopy
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github" / "workflows" / "main.yaml"
BASE_CONFIG = ROOT / ".test" / "config" / "config.yaml"

EXPECTED_ROUTES = [
    {
        "name": "prokaryotic isolate",
        "id": "prokaryotic-isolate",
        "config": ".test/config/config.prokaryotic-isolate.yaml",
    },
    {
        "name": "eukaryotic isolate",
        "id": "eukaryotic-isolate",
        "config": ".test/config/config.eukaryotic-isolate.yaml",
    },
    {
        "name": "mixed prokaryotic/eukaryotic metagenome",
        "id": "mixed-metagenome",
        "config": ".test/config/config.mixed-metagenome.yaml",
    },
]


def _load_yaml(path):
    return yaml.load(path.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)


def _named_step(job, name):
    return next(step for step in job["steps"] if step.get("name") == name)


def _merge_config(base, override):
    merged = deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge_config(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged


def _read_one_row_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    return rows[0]


def _assert_tiny_fasta(path):
    assert path.is_file()
    lines = path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 2
    assert lines[0].startswith(">") and len(lines[0]) > 1
    assert lines[1] and set(lines[1].upper()) <= set("ACGTN")


def test_template_plus_ci_contract():
    workflow = _load_yaml(WORKFLOW)

    triggers = workflow["on"]
    assert set(triggers) == {"push", "pull_request", "workflow_dispatch"}
    assert triggers["push"]["branches"] == ["main"]
    assert triggers["pull_request"]["branches"] == ["main"]
    assert workflow["permissions"] == {"contents": "read"}
    assert workflow["concurrency"] == {
        "group": "${{ github.workflow }}-${{ github.ref }}",
        "cancel-in-progress": "true",
    }

    jobs = workflow["jobs"]
    assert set(jobs) == {"formatting", "linting", "python-tests", "testing"}
    assert "needs" not in jobs["formatting"]
    assert "needs" not in jobs["linting"]
    assert "needs" not in jobs["python-tests"]
    assert jobs["testing"]["needs"] == [
        "formatting",
        "linting",
        "python-tests",
    ]
    assert {name: job["timeout-minutes"] for name, job in jobs.items()} == {
        "formatting": "15",
        "linting": "15",
        "python-tests": "20",
        "testing": "20",
    }
    for job in jobs.values():
        checkout_refs = [
            step["uses"]
            for step in job["steps"]
            if step.get("uses", "").startswith("actions/checkout@")
        ]
        assert checkout_refs == ["actions/checkout@v6"]

    formatting = jobs["formatting"]
    assert formatting["steps"][0] == {
        "uses": "actions/checkout@v6",
        "with": {"fetch-depth": "0"},
    }
    formatter = _named_step(formatting, "Check formatting")
    assert formatter["uses"] == "super-linter/super-linter@v8"
    assert {
        key: value
        for key, value in formatter["env"].items()
        if key.startswith("VALIDATE_")
    } == {
        "VALIDATE_ALL_CODEBASE": "true",
        "VALIDATE_PYTHON_BLACK": "true",
        "VALIDATE_SNAKEMAKE_SNAKEFMT": "true",
        "VALIDATE_YAML_PRETTIER": "true",
    }

    lint = _named_step(jobs["linting"], "Lint workflow")
    assert lint["uses"] == "snakemake/snakemake-github-action@v2"
    assert lint["with"] == {
        "directory": ".test",
        "snakefile": "workflow/Snakefile",
        "args": "--lint --configfile .test/config/config.prokaryotic-isolate.yaml",
        "snakemake-version": "9.19.0",
    }

    python_job = jobs["python-tests"]
    setup_python = next(
        step
        for step in python_job["steps"]
        if step.get("uses") == "actions/setup-python@v6"
    )
    assert setup_python["with"]["python-version"] == "3.13"
    assert setup_python["with"]["cache"] == "pip"
    assert setup_python["with"]["cache-dependency-path"].splitlines() == [
        "tests/requirements.txt",
        "tests/requirements-workflow.txt",
    ]
    assert _named_step(python_job, "Install test dependencies")["run"] == (
        "python -m pip install -r tests/requirements-workflow.txt"
    )
    assert _named_step(python_job, "Run Python tests")["run"] == (
        "python -m pytest -q -ra tests"
    )

    testing = jobs["testing"]
    assert testing["strategy"]["matrix"]["route"] == EXPECTED_ROUTES
    assert testing["steps"][0] == {"uses": "actions/checkout@v6"}
    assert len(testing["steps"]) == 2

    dry_run = _named_step(testing, "Dry-run workflow")
    assert dry_run["uses"] == "snakemake/snakemake-github-action@v2"
    assert dry_run["with"]["directory"] == ".test"
    assert dry_run["with"]["snakefile"] == "workflow/Snakefile"
    assert dry_run["with"]["snakemake-version"] == "9.19.0"
    dag_args = " ".join(dry_run["with"]["args"].split())
    assert dag_args == (
        "--configfile ${{ matrix.route.config }} "
        "--dry-run --printshellcmds --cores 2"
    )
    for forbidden in (
        "--software-deployment-method",
        "--sdm",
        "conda",
        "--create-envs-only",
        "--report",
        "--containerize",
    ):
        assert forbidden not in dag_args.lower()
    assert all(
        "upload-artifact" not in step.get("uses", "") for step in testing["steps"]
    )


def test_ci_dependencies_are_pinned():
    requirements = (ROOT / "tests" / "requirements.txt").read_text(encoding="utf-8")
    workflow_requirements = (ROOT / "tests" / "requirements-workflow.txt").read_text(
        encoding="utf-8"
    )

    assert requirements.splitlines() == [
        "pytest==9.0.2",
        "PyYAML==6.0.3",
        "jsonschema==4.26.0",
        "pandas==2.3.3",
    ]
    assert workflow_requirements.splitlines() == [
        "-r requirements.txt",
        "snakemake==9.19.0",
    ]


def test_ci_route_overrides_merge_with_valid_tiny_inputs():
    base = _load_yaml(BASE_CONFIG)
    test_root = ROOT / ".test"
    expected_inputs = {
        "prokaryotic-isolate": {
            "mode": "isolate",
            "sample": "ci_prokaryotic_isolate",
            "domain": "prok",
        },
        "eukaryotic-isolate": {
            "mode": "isolate",
            "sample": "ci_eukaryotic_isolate",
            "domain": "euk",
        },
        "mixed-metagenome": {
            "mode": "metagenome",
            "sample": "ci_mixed_metagenome",
        },
    }

    for route in EXPECTED_ROUTES:
        override = _load_yaml(ROOT / route["config"])
        expected_override_sections = {"input", "qa", "qa_filter"}
        if route["id"] == "eukaryotic-isolate":
            expected_override_sections.add("funannotate2")
            assert override["funannotate2"] == {"species": "CI fixture species"}
        assert set(override) == expected_override_sections
        assert override["qa"] == {"enabled": "true"}
        assert override["qa_filter"] == {"enabled": "true"}

        merged = _merge_config(base, override)
        expected = expected_inputs[route["id"]]
        assert merged["input"]["mode"] == expected["mode"]
        assert {
            "input",
            "qa",
            "qa_filter",
            "bakta",
            "gtdbtk",
            "metaeuk",
            "recognizer_prok",
            "recognizer_euk",
            "funannotate2",
            "upimapi",
            "threads",
        } <= set(merged)
        assert {"prok_tools", "euk_tools", "checkm2", "eukcc"} <= set(merged["qa"])

        sheet = test_root / merged["input"]["sample_sheet"]
        row = _read_one_row_tsv(sheet)
        assert row["sample"] == expected["sample"]
        if expected["mode"] == "isolate":
            assert set(row) == {"sample", "domain", "path"}
            assert row["domain"] == expected["domain"]
            _assert_tiny_fasta(test_root / row["path"])
        else:
            assert set(row) == {"sample", "prok_bins_dir", "euk_bins_dir"}
            for domain in ("prok", "euk"):
                bin_dir = test_root / row[f"{domain}_bins_dir"]
                fastas = sorted(bin_dir.glob("*.fasta"))
                assert len(fastas) == 1
                _assert_tiny_fasta(fastas[0])
