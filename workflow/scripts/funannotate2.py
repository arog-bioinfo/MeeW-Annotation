"""Run one stateful Funannotate2 isolate project and publish stable outputs."""

import json
import os
import shlex
import shutil
import subprocess
from pathlib import Path

from common.euk_options import select_funannotate2_route, validate_funannotate2_extras
from common.funannotate2_database import ensure_funannotate2_database
from common.versions import version_contains

EXPECTED_TOOL_VERSION = "26.6.21"


def value(mapping, key, default=""):
    result = mapping.get(key, default)
    return default if result is None else result


def nonempty_text(raw):
    return str(raw).strip() if raw else ""


def command_argv(program, phase, *arguments, extra=""):
    return [
        program,
        phase,
        *[str(argument) for argument in arguments],
        *shlex.split(extra),
    ]


def run_command(argv, phase, log_handle, commands):
    commands.append({"phase": phase, "argv": argv})
    log_handle.write(f"$ {shlex.join(argv)}\n")
    log_handle.flush()
    subprocess.run(argv, check=True, stdout=log_handle, stderr=subprocess.STDOUT)


def runtime_version(program):
    completed = subprocess.run(
        [program, "--version"], check=True, capture_output=True, text=True
    )
    return (completed.stdout or completed.stderr).strip()


def unique_native(project, pattern, producer):
    matches = sorted(project.glob(pattern))
    if len(matches) != 1 or not matches[0].is_file():
        raise FileNotFoundError(
            f"{producer} expected exactly one native file matching {project / pattern}; "
            f"found {len(matches)}."
        )
    return matches[0]


def publish(source, destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


configured_species = value(snakemake.params, "species")
validate_funannotate2_extras(snakemake.params)
params_json = value(snakemake.params, "params_json")
pretrained = value(snakemake.params, "pretrained")
model_input = snakemake.input.get("model", "")
if params_json:
    params_json = model_input
elif pretrained and model_input:
    pretrained = model_input
route, species, reusable_model = select_funannotate2_route(
    configured_species, params_json, pretrained
)

project = Path(snakemake.output.project)
fasta = Path(snakemake.input.fasta)
db_dir = nonempty_text(value(snakemake.params, "db_dir"))
clean_enabled = bool(value(snakemake.params, "clean", False))
strain = nonempty_text(value(snakemake.params, "strain"))
configured_busco_lineage = nonempty_text(value(snakemake.params, "busco_lineage"))
log_path = Path(snakemake.log[0]).resolve()
log_path.parent.mkdir(parents=True, exist_ok=True)
log_path.write_text("", encoding="utf-8")

if not db_dir:
    raise ValueError("funannotate2.db_dir must resolve to a persistent resource cache.")
db_dir = str(Path(db_dir).resolve())
os.environ["FUNANNOTATE2_DB"] = db_dir
runtime_tool_version = runtime_version("funannotate2")
if not version_contains(runtime_tool_version, EXPECTED_TOOL_VERSION):
    raise ValueError(
        f"Expected Funannotate2 {EXPECTED_TOOL_VERSION}, found {runtime_tool_version!r}."
    )
database_status = ensure_funannotate2_database(
    db_dir, value(snakemake.params, "databases", ["all"]), log_path
)
db_info = Path(db_dir) / "funannotate-db-info.json"
project.mkdir(parents=True, exist_ok=True)
commands = []
phases = []
with log_path.open("a", encoding="utf-8") as log_handle:
    phases.append(
        {
            "phase": "database_prepare",
            **database_status.as_dict(),
        }
    )

    scientific_fasta = fasta
    if clean_enabled:
        scientific_fasta = project / "cleaned.fasta"
        argv = command_argv(
            "funannotate2",
            "clean",
            "-f",
            fasta,
            "-o",
            scientific_fasta,
            "--cpus",
            snakemake.threads,
            extra=value(snakemake.params, "extra_clean"),
        )
        run_command(argv, "clean", log_handle, commands)
        if not scientific_fasta.is_file():
            raise FileNotFoundError(
                f"Funannotate2 clean completed without {scientific_fasta}."
            )
        phases.append({"phase": "clean", "status": "completed"})
    else:
        phases.append(
            {
                "phase": "clean",
                "status": "disabled",
                "downstream_input": str(fasta.resolve()),
            }
        )

    model_argument = reusable_model
    if route == "genome_specific_training":
        argv = command_argv(
            "funannotate2",
            "train",
            "-f",
            scientific_fasta,
            "-s",
            species,
            "-o",
            project,
            "--cpus",
            snakemake.threads,
            *(["--strain", strain] if strain else []),
            "--busco-lineage",
            configured_busco_lineage,
            extra=value(snakemake.params, "extra_train"),
        )
        run_command(argv, "train", log_handle, commands)
        model_argument = str(
            unique_native(project, "train_results/*.params.json", "Funannotate2 train")
        )
        phases.append({"phase": "model", "status": "trained"})
    else:
        phases.append({"phase": "model", "status": "reused", "route": route})

    predict_args = [
        "-f",
        scientific_fasta,
        "-o",
        project,
        "-p",
        model_argument,
        "-s",
        species,
        "--cpus",
        snakemake.threads,
    ]
    if strain:
        predict_args.extend(["--strain", strain])
    argv = command_argv(
        "funannotate2",
        "predict",
        *predict_args,
        extra=value(snakemake.params, "extra_predict"),
    )
    run_command(argv, "predict", log_handle, commands)
    phases.append({"phase": "predict", "status": "completed"})

    annotate_args = [
        "-i",
        project,
        "--cpus",
        snakemake.threads,
        "-s",
        species,
        "--busco-lineage",
        configured_busco_lineage,
    ]
    if strain:
        annotate_args.extend(["-st", strain])
    argv = command_argv(
        "funannotate2",
        "annotate",
        *annotate_args,
        extra=value(snakemake.params, "extra_annotate"),
    )
    run_command(argv, "annotate", log_handle, commands)
    phases.append({"phase": "annotate", "status": "completed"})

native = {
    "gff3": unique_native(project, "annotate_results/*.gff3", "Funannotate2 annotate"),
    "proteins": unique_native(
        project, "annotate_results/*.proteins.fa", "Funannotate2 annotate"
    ),
    "transcripts": unique_native(
        project, "annotate_results/*.transcripts.fa", "Funannotate2 annotate"
    ),
    "summary": unique_native(
        project, "annotate_results/*.summary.json", "Funannotate2 annotate"
    ),
}
for key, source in native.items():
    publish(source, snakemake.output[key])

database = {
    "resolved_path": str(Path(db_dir).resolve()),
    "info_file": str(db_info),
    "requested_databases": sorted(value(snakemake.params, "databases", ["all"])),
    "cache_status": database_status.as_dict(),
}

provenance = {
    "tool": {"name": "Funannotate2", "version": runtime_tool_version},
    "route": route,
    "species": species,
    "strain": strain or None,
    "input_fasta": str(fasta.resolve()),
    "scientific_fasta": str(scientific_fasta.resolve()),
    "cleaning_enabled": clean_enabled,
    "model": {
        "kind": route,
        "source": (
            str(Path(model_argument).resolve())
            if route != "pretrained"
            else model_argument
        ),
    },
    "database": database,
    "busco_lineage": configured_busco_lineage,
    "phases": phases,
    "commands": commands,
    "native_outputs": {key: str(path.resolve()) for key, path in native.items()},
    "stable_outputs": {
        key: str(Path(snakemake.output[key]).resolve()) for key in native
    },
    "project_state": str(project.resolve()),
}
Path(snakemake.output.provenance).write_text(
    json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
)
