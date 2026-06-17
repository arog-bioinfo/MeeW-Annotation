"""Run Funannotate2 steps for eukaryotic isolate genome annotation."""

import os

from pathlib import Path

from snakemake.shell import shell


def param_value(key, default=""):
    """Return a Snakemake parameter, normalizing YAML null to a string default."""
    value = snakemake.params.get(key, default)
    return default if value is None else value


def bool_param(key, default=False):
    """Return a Snakemake boolean parameter with YAML null normalized."""
    value = snakemake.params.get(key, default)
    return default if value is None else bool(value)


def list_param(key, default=None):
    """Return a Snakemake list parameter with YAML null/empty normalized."""
    if default is None:
        default = []
    value = snakemake.params.get(key, default)
    return default if value is None or value == "" else value


step = param_value("step")

species = param_value("species")
strain = param_value("strain")
params = param_value("params")
pretrained = param_value("pretrained")
db_dir = param_value("db_dir")
install_db = bool_param("install_db", False)
databases = list_param("databases", ["all"])

# These extra CLI options are supplied by workflow configuration.
extra_install = param_value("extra_install")
extra_clean = param_value("extra_clean")
extra_train = param_value("extra_train")
extra_predict = param_value("extra_predict")
extra_annotate = param_value("extra_annotate")

if db_dir:
    os.environ["FUNANNOTATE2_DB"] = str(db_dir)


def named_value(values, key, default=""):
    """Return a Snakemake named-list value without requiring every key."""
    value = values.get(key, default)
    return value if value else default


def first_value(values, default=""):
    """Return the first Snakemake list value when present."""
    return values[0] if len(values) else default


def output_path(key, default):
    """Prefer named outputs, then output.outdir, then a predictable default."""
    value = named_value(snakemake.output, key)
    if value:
        return Path(value)

    outdir = named_value(snakemake.output, "outdir")
    if outdir:
        return Path(outdir) / default

    return Path(default)


def output_outdir(default="funannotate2"):
    """Return the named output directory when available."""
    return Path(
        named_value(snakemake.output, "outdir") or param_value("outdir") or default
    )


def input_fasta(default=""):
    """Prefer snakemake.input.fasta, falling back to first input/default."""
    value = named_value(snakemake.input, "fasta")
    return Path(value) if value else Path(first_value(snakemake.input, default))


def log_file():
    """Return a log target suitable for append redirection."""
    value = first_value(snakemake.log, "/dev/null")
    if value != "/dev/null":
        Path(value).parent.mkdir(parents=True, exist_ok=True)
    return value


def touch_done():
    """Touch a named done marker when the rule provides one."""
    done = named_value(snakemake.output, "done")
    if done:
        Path(done).parent.mkdir(parents=True, exist_ok=True)
        Path(done).touch()


outdir = output_outdir()
cleaned_fasta = output_path("cleaned_fasta", Path("preprocess") / "cleaned.fasta")
log = log_file()


def run_install_db():
    """Check for Funannotate2 databases and optionally install missing DBs."""
    if not db_dir:
        raise ValueError("funannotate2.db_dir is required for the install_db step.")

    db_path = Path(db_dir)
    db_info = db_path / "funannotate-db-info.json"
    os.environ["FUNANNOTATE2_DB"] = str(db_path)

    if db_info.exists():
        return

    if not install_db:
        raise FileNotFoundError(
            f"Funannotate2 database info file is missing at {db_info}; "
            "set funannotate2.install_db: true or install databases manually."
        )

    db_path.mkdir(parents=True, exist_ok=True)
    for database in databases:
        shell("funannotate2 install -d {database:q} {extra_install} >> {log:q} 2>&1")

    if not db_info.exists():
        raise FileNotFoundError(
            f"Funannotate2 database install completed but {db_info} is missing."
        )


def run_clean():
    """Run Funannotate2 clean on the input FASTA."""
    fasta = input_fasta()
    cleaned_fasta.parent.mkdir(parents=True, exist_ok=True)
    shell(
        "funannotate2 clean "
        "-f {fasta:q} "
        "-o {cleaned_fasta:q} "
        "--cpus {snakemake.threads} "
        "{extra_clean} "
        ">> {log:q} 2>&1"
    )


def run_train():
    """Run Funannotate2 training on the cleaned FASTA."""
    if not species:
        raise ValueError("funannotate2 species is required for the train step.")

    fasta = input_fasta(cleaned_fasta)
    outdir.mkdir(parents=True, exist_ok=True)
    command = (
        "funannotate2 train "
        "-f {fasta:q} "
        "-s {species:q} "
        "-o {outdir:q} "
        "--cpus {snakemake.threads} "
    )
    if strain:
        command += "--strain {strain:q} "
    command += "{extra_train} >> {log:q} 2>&1"
    shell(command)


def run_predict():
    """Run Funannotate2 prediction in an existing Funannotate2 outdir."""
    # Funannotate2 documents -p/--params/--pretrained for supplying training
    # parameters. Prefer the explicit params setting when both are supplied.
    command = "funannotate2 predict -i {outdir:q} --cpus {snakemake.threads} "
    if species:
        command += "-s {species:q} "
    if strain:
        command += "--strain {strain:q} "
    if params:
        command += "-p {params:q} "
    elif pretrained:
        command += "--pretrained {pretrained:q} "
    command += "{extra_predict} >> {log:q} 2>&1"
    shell(command)


def run_annotate():
    """Run Funannotate2 annotation in an existing Funannotate2 outdir."""
    command = "funannotate2 annotate -i {outdir:q} --cpus {snakemake.threads} "
    if species:
        command += "-s {species:q} "
    if strain:
        command += "--strain {strain:q} "
    command += "{extra_annotate} >> {log:q} 2>&1"
    shell(command)


steps = {
    "install_db": run_install_db,
    "clean": run_clean,
    "train": run_train,
    "predict": run_predict,
    "annotate": run_annotate,
}

if step not in steps:
    raise ValueError(
        f"Invalid Funannotate2 step {step!r}; expected one of: "
        f"{', '.join(sorted(steps))}"
    )

steps[step]()
touch_done()
