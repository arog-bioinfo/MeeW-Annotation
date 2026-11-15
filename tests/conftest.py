import csv
import hashlib
import os
import socket
import stat
import subprocess
import sys
import urllib.request
from pathlib import Path

import pytest
import yaml

REPOSITORY = Path(__file__).resolve().parents[1]
SNAKEFILE = REPOSITORY / "workflow" / "Snakefile"
sys.path.insert(0, str(REPOSITORY / "workflow" / "scripts"))

from common.gtdbtk_database import (  # noqa: E402
    GTDBTK_REQUIRED_DIRECTORIES,
    GTDBTK_REQUIRED_FILES,
)

SELECTION_COLUMNS = [
    "sample",
    "path",
    "domain",
    "completeness",
    "contamination",
    "min_completeness",
    "max_contamination",
    "decision",
    "reason",
]

OFFLINE_EXECUTABLES = (
    "aria2c",
    "bakta_db",
    "checkm2",
    "conda",
    "curl",
    "funannotate2",
    "mamba",
    "micromamba",
    "mmseqs",
    "wget",
)


@pytest.fixture(scope="session")
def offline_guard_assets(tmp_path_factory):
    root = tmp_path_factory.mktemp("annotation-offline-guard")
    executable_dir = root / "bin"
    executable_dir.mkdir()
    wrapper = (
        "#!/usr/bin/env python3\n"
        "import os\n"
        "import sys\n"
        "from pathlib import Path\n"
        "log = os.environ.get('MEEW_ANNOTATION_OFFLINE_ATTEMPT_LOG')\n"
        "if log:\n"
        "    with Path(log).open('a', encoding='utf-8') as handle:\n"
        "        handle.write(Path(sys.argv[0]).name + ' ' + ' '.join(sys.argv[1:]) + '\\n')\n"
        "sys.stderr.write('MeeW Annotation test offline guard blocked executable: ' + Path(sys.argv[0]).name + '\\n')\n"
        "raise SystemExit(97)\n"
    )
    for name in OFFLINE_EXECUTABLES:
        executable = executable_dir / name
        executable.write_text(wrapper, encoding="utf-8")
        executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    return {"root": root, "bin": executable_dir}


@pytest.fixture(autouse=True)
def offline_guard(monkeypatch, offline_guard_assets, request):
    """Block network/download commands unless a test injects an explicit fake."""

    def blocked(operation):
        def reject(*args, **kwargs):
            raise RuntimeError(
                f"MeeW Annotation test offline guard blocked {operation}"
            )

        return reject

    monkeypatch.setattr(
        socket, "create_connection", blocked("socket.create_connection")
    )
    monkeypatch.setattr(socket.socket, "connect", blocked("socket.connect"))
    monkeypatch.setattr(urllib.request, "urlopen", blocked("urllib.request.urlopen"))
    monkeypatch.setattr(
        urllib.request, "urlretrieve", blocked("urllib.request.urlretrieve")
    )

    attempt_name = hashlib.sha256(request.node.nodeid.encode()).hexdigest()[:20]
    attempt_log = offline_guard_assets["root"] / f"attempts-{attempt_name}.log"
    sitecustomize = REPOSITORY / "tests" / "offline_sitecustomize"
    existing_pythonpath = os.environ.get("PYTHONPATH", "")
    pythonpath = str(sitecustomize)
    if existing_pythonpath:
        pythonpath = f"{pythonpath}{os.pathsep}{existing_pythonpath}"
    monkeypatch.setenv("MEEW_ANNOTATION_TEST_OFFLINE", "1")
    monkeypatch.setenv("MEEW_ANNOTATION_OFFLINE_ATTEMPT_LOG", str(attempt_log))
    monkeypatch.setenv("PYTHONPATH", pythonpath)
    monkeypatch.setenv(
        "PATH", f"{offline_guard_assets['bin']}{os.pathsep}{os.environ['PATH']}"
    )
    return {**offline_guard_assets, "attempt_log": attempt_log}


def _write_fasta(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(">contig\nACGTACGT\n", encoding="utf-8")


def _selection_row(sample, path, domain, decision):
    return {
        "sample": sample,
        "path": str(path.resolve()),
        "domain": domain,
        "completeness": "95",
        "contamination": "1",
        "min_completeness": "50",
        "max_contamination": "10",
        "decision": decision,
        "reason": "test_fixture",
    }


class WorkflowCase:
    def __init__(
        self,
        root,
        mode,
        enabled,
        rows,
        *,
        isolate_domain="prok",
        qa_enabled=False,
        filter_enabled=False,
        supplied_domains=None,
        fasta_names=None,
    ):
        self.root = root
        self.sample = "project" if mode == "metagenome" else "isolate1"
        self.mode = mode
        self.rows = rows
        fasta_names = {} if fasta_names is None else fasta_names
        self.config = yaml.safe_load(
            (REPOSITORY / "config" / "config.yaml").read_text(encoding="utf-8")
        )
        self.config["input"]["mode"] = mode
        self.config["qa"]["enabled"] = qa_enabled
        self.config["qa_filter"]["enabled"] = filter_enabled
        self.config["gtdbtk"]["enabled"] = enabled
        self.config["gtdbtk"]["data_dir"] = str(self.root / "gtdbtk_db")
        self.config["bakta"]["db"] = str(self.root / "bakta_db")
        self.config["metaeuk"]["cache_dir"] = str(self.root / "metaeuk_db")
        self.config["metaeuk"]["db"] = ""
        self.config["funannotate2"]["db_dir"] = str(self.root / "funannotate2_db")
        self.config["funannotate2"]["species"] = "Fixture species"
        self.config["recognizer_prok"]["resources_dir"] = str(
            self.root / "recognizer_db"
        )
        self.config["recognizer_euk"]["resources_dir"] = str(
            self.root / "recognizer_db"
        )
        self.config["upimapi"]["resources_dir"] = str(self.root / "upimapi_db")

        if mode == "isolate":
            fasta = self.root / "input" / fasta_names.get("isolate1", "isolate1.fasta")
            _write_fasta(fasta)
            self.fastas = {"isolate1": fasta}
            sheet = self.root / "isolate_samples.tsv"
            sheet.write_text(
                f"sample\tdomain\tpath\n{self.sample}\t{isolate_domain}\t{fasta}\n",
                encoding="utf-8",
            )
            selection_rows = [
                _selection_row("isolate1", fasta, isolate_domain, "report_only")
            ]
        else:
            supplied_domains = (
                {"prok", "euk"} if supplied_domains is None else set(supplied_domains)
            )
            prok_dir = self.root / "input" / "prok"
            euk_dir = self.root / "input" / "euk"
            if "prok" in supplied_domains:
                prok_dir.mkdir(parents=True)
            if "euk" in supplied_domains:
                euk_dir.mkdir(parents=True)
            self.fastas = {}
            selection_rows = []
            for sample, domain, decision in rows:
                directory = prok_dir if domain == "prok" else euk_dir
                if domain not in supplied_domains:
                    raise ValueError(f"Rows supplied for absent {domain} domain")
                fasta = directory / fasta_names.get(sample, f"{sample}.fasta")
                _write_fasta(fasta)
                self.fastas[sample] = fasta
                selection_rows.append(_selection_row(sample, fasta, domain, decision))
            sheet = self.root / "metagenome_samples.tsv"
            sheet.write_text(
                "sample\tprok_bins_dir\teuk_bins_dir\n"
                f"{self.sample}\t"
                f"{prok_dir if 'prok' in supplied_domains else ''}\t"
                f"{euk_dir if 'euk' in supplied_domains else ''}\n",
                encoding="utf-8",
            )

        self.config["input"]["sample_sheet"] = str(sheet)
        self.config_path = self.root / "config.yaml"
        self.config_path.write_text(
            yaml.safe_dump(self.config, sort_keys=False), encoding="utf-8"
        )

        if enabled:
            # Two integration cases execute the GTDB-Tk adapter with a fake
            # classifier, so seed a complete tiny native R232 cache for reuse.
            for relative in GTDBTK_REQUIRED_FILES:
                material = self.root / "gtdbtk_db" / relative
                material.parent.mkdir(parents=True, exist_ok=True)
                material.write_text(
                    (
                        "VERSION_DATA=r232\n"
                        if material.name == "metadata.txt"
                        else "fixture\n"
                    ),
                    encoding="utf-8",
                )
            for relative in GTDBTK_REQUIRED_DIRECTORIES:
                directory = self.root / "gtdbtk_db" / relative
                directory.mkdir(parents=True, exist_ok=True)
                (directory / "fixture").write_text("fixture\n", encoding="utf-8")

        if qa_enabled:
            qa_outputs = {
                self.root
                / "results"
                / self.sample
                / "qa"
                / "prok"
                / "checkm2"
                / "quality_report.tsv": "Name\tCompleteness\tContamination\n",
                self.root
                / "results"
                / self.sample
                / "qa"
                / "euk"
                / "eukcc"
                / "eukcc.csv": "bin\tcompleteness\tcontamination\tncbi_lng\n",
                self.root
                / "results"
                / self.sample
                / "qa"
                / "euk"
                / "eukcc"
                / "bad_quality.csv": "bin\tcompleteness\tcontamination\tncbi_lng\n",
                self.root
                / "results"
                / self.sample
                / "qa"
                / "euk"
                / "eukcc"
                / "missing_marker_genes.txt": "bin\n",
            }
            for report, content in qa_outputs.items():
                report.parent.mkdir(parents=True, exist_ok=True)
                report.write_text(content, encoding="utf-8")

        selection = self.root / "results" / self.sample / "qa" / "selection.tsv"
        selection.parent.mkdir(parents=True, exist_ok=True)
        with selection.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=SELECTION_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(selection_rows)
        self.selection = selection

    @property
    def classification(self):
        return self.root / "results" / self.sample / "gtdbtk" / "classification.tsv"

    def run_result(self, *, dry_run=True, targets=(), env=None, extra_arguments=()):
        command = [
            "snakemake",
            "--snakefile",
            str(SNAKEFILE),
            "--directory",
            str(self.root),
            "--configfile",
            str(self.config_path),
            "--cores",
            "1",
            "--nolock",
            "--printshellcmds",
        ]
        if dry_run:
            command.append("--dry-run")
        command.extend(targets)
        command.extend(extra_arguments)
        result = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            env=env,
        )
        return result

    def run(self, *, dry_run=True, targets=(), env=None, extra_arguments=()):
        result = self.run_result(
            dry_run=dry_run,
            targets=targets,
            env=env,
            extra_arguments=extra_arguments,
        )
        output = result.stdout + result.stderr
        assert result.returncode == 0, output
        return output


@pytest.fixture
def workflow_case(tmp_path):
    def make(mode, *, enabled, rows=(), **kwargs):
        return WorkflowCase(tmp_path, mode, enabled, rows, **kwargs)

    return make


@pytest.fixture
def fake_gtdbtk(tmp_path):
    executable_dir = tmp_path / "fake_bin"
    executable_dir.mkdir()
    executable = executable_dir / "gtdbtk"
    executable.write_text(
        """#!/usr/bin/env python3
import json
import os
import sys
from pathlib import Path

args = sys.argv[1:]
if not args or args[0] != "classify_wf":
    raise SystemExit("expected classify_wf")

def option(name):
    return Path(args[args.index(name) + 1])

genome_dir = option("--genome_dir")
out_dir = option("--out_dir")
if " " in str(genome_dir):
    raise SystemExit("GTDB-Tk 2.7.1 rejects genome paths containing spaces")
record = {
    "args": args,
    "genome_dir": str(genome_dir.resolve()),
    "genomes": [
        {
            "name": path.name,
            "is_symlink": path.is_symlink(),
            "target": str(path.resolve()),
            "content": path.read_text(encoding="utf-8"),
        }
        for path in sorted(genome_dir.glob("*.fasta"))
    ],
}
with Path(os.environ["FAKE_GTDBTK_RECORDS"]).open("a", encoding="utf-8") as handle:
    handle.write(json.dumps(record) + "\\n")
with Path(os.environ["FAKE_GTDBTK_CALLS"]).open("a", encoding="utf-8") as handle:
    handle.write("classify_wf\\n")
if os.environ.get("FAKE_GTDBTK_FAIL") == "1":
    raise SystemExit(23)

out_dir.mkdir(parents=True, exist_ok=True)
genomes = sorted(path.stem for path in genome_dir.glob("*.fasta"))
header = "user_genome\\tclassification\\n"
(out_dir / "gtdbtk.bac120.summary.tsv").write_text(
    header + "".join(f"{genome}\\td__Bacteria\\n" for genome in genomes),
    encoding="utf-8",
)
(out_dir / "gtdbtk.ar53.summary.tsv").write_text(header, encoding="utf-8")
""",
        encoding="utf-8",
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    calls = tmp_path / "gtdbtk_calls.txt"
    records = tmp_path / "gtdbtk_records.jsonl"
    env = os.environ.copy()
    env["PATH"] = f"{executable_dir}{os.pathsep}{env['PATH']}"
    env["FAKE_GTDBTK_CALLS"] = str(calls)
    env["FAKE_GTDBTK_RECORDS"] = str(records)
    return env, calls, records
