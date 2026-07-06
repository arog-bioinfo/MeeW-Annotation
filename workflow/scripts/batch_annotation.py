from pathlib import Path
import shutil

import pandas as pd
from snakemake.shell import shell

from filter_qa import fasta_stem, parse_checkm2_reports, parse_eukcc_reports

FASTA_SUFFIXES = (".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz")


def fasta_paths(directory):
    directory = Path(directory)
    if not directory.exists():
        return []
    return sorted(
        path
        for path in directory.iterdir()
        if path.is_file()
        and any(path.name.endswith(suffix) for suffix in FASTA_SUFFIXES)
    )


def reset_directory(path):
    path = Path(path)
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def qa_row_for_bin(bin_path, qa):
    names = [bin_path.name, fasta_stem(bin_path), Path(bin_path).stem]
    for name in dict.fromkeys(names):
        if name in qa.index:
            return qa.loc[name]
    return None


def write_manifest(rows, manifest):
    manifest = Path(manifest)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", encoding="utf-8") as handle:
        handle.write("sample\tpath\n")
        for sample, path in rows:
            handle.write(f"{sample}\t{path}\n")


def read_manifest(manifest):
    manifest = Path(manifest)
    if not manifest.exists():
        return []
    rows = pd.read_csv(manifest, sep="\t")
    if rows.empty:
        return []
    return [(row["sample"], row["path"]) for _, row in rows.iterrows()]


def passing_bins(domain):
    if domain == "prok":
        qa = parse_checkm2_reports(snakemake.input.report)
    else:
        qa = parse_eukcc_reports(snakemake.input.report)

    min_completeness = float(snakemake.params.min_completeness)
    max_contamination = float(snakemake.params.max_contamination)
    outdir = reset_directory(snakemake.output.passing_dir)
    rows = []

    for bin_path in fasta_paths(snakemake.input.bins_dir):
        qa_row = qa_row_for_bin(bin_path, qa)
        if qa_row is None:
            continue
        if (
            qa_row["completeness"] >= min_completeness
            and qa_row["contamination"] <= max_contamination
        ):
            # Symlink staging gives directory-oriented batch tools a filtered bin
            # directory without copying accepted bin FASTAs or duplicating data.
            target = outdir / bin_path.name
            target.symlink_to(bin_path.resolve())
            rows.append((fasta_stem(bin_path), target.resolve()))

    write_manifest(rows, snakemake.output.manifest)
    Path(snakemake.output.done).touch()


def log_redirect():
    log = Path(snakemake.log[0])
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text("", encoding="utf-8")
    return f">> {log} 2>&1"


def run_bakta(rows):
    outdir = reset_directory(snakemake.output.outdir)
    db = snakemake.params.get("db", "")
    db_cmd = f"--db {db}" if db else ""
    extra = snakemake.params.get("extra", "")
    redir = log_redirect()
    for sample, fasta in rows:
        sample_outdir = outdir / sample
        shell(
            "bakta --output {sample_outdir} --prefix {sample} --force "
            "--threads {snakemake.threads} {db_cmd} {extra} {fasta} {redir}"
        )


def run_prodigal(rows):
    outdir = reset_directory(snakemake.output.outdir)
    extra = snakemake.params.get("extra", "")
    redir = log_redirect()
    for sample, fasta in rows:
        shell(
            "prodigal -i {fasta} -o {outdir}/{sample}.gff "
            "-a {outdir}/{sample}.faa -d {outdir}/{sample}.fna {extra} {redir}"
        )


def discover_protein_path(proteins_dir, sample):
    proteins_dir = Path(proteins_dir)
    for protein_path in [
        proteins_dir / f"{sample}.faa",
        proteins_dir / sample / f"{sample}.faa",
    ]:
        if protein_path.exists():
            return protein_path
    return None


def protein_rows_from_manifest():
    proteins_dir = Path(snakemake.input.proteins_dir)
    rows = []
    for sample, _ in read_manifest(snakemake.input.manifest):
        protein_path = discover_protein_path(proteins_dir, sample)
        if protein_path is not None:
            rows.append((sample, protein_path))
    return rows


def run_recognizer(rows, euk=False):
    outdir = reset_directory(snakemake.output.outdir)
    resources_dir = snakemake.params.get("resources_dir", "")
    resources_cmd = f"-rd {resources_dir}" if resources_dir else ""
    custom_db = snakemake.params.get("custom_db", "") if euk else ""
    custom_db_cmd = f"--custom-databases -dbs {custom_db}" if custom_db else ""
    extra = snakemake.params.get("extra", "")
    redir = log_redirect()
    for sample, fasta in rows:
        sample_outdir = outdir / sample
        sample_outdir.mkdir(parents=True, exist_ok=True)
        shell(
            "recognizer -f {fasta} -o {sample_outdir} -t {snakemake.threads} "
            "{custom_db_cmd} {resources_cmd} {extra} {redir}"
        )


def run_upimapi(rows):
    outdir = reset_directory(snakemake.output.outdir)
    db = snakemake.params.get("db", "")
    db_custom = snakemake.params.get("db_custom", "")
    if db and db_custom:
        raise ValueError("Only one UPIMAPI database mode can be configured.")
    db_cmd = f"--database {db}" if db else ""
    db_custom_cmd = f"--database {db_custom}" if db_custom else ""
    resources_dir = snakemake.params.get("resources_dir", "")
    resources_cmd = f"-rd {resources_dir}" if resources_dir else ""
    skip_db_check_if_exists = snakemake.params.get("skip_db_check_if_exists", True)
    db2file = {
        "uniprot": Path(resources_dir) / "uniprot.fasta",
        "swissprot": Path(resources_dir) / "uniprot_sprot.fasta",
        "taxids": Path(resources_dir) / "taxids_database.fasta",
    }
    db_custom_exists = bool(db_custom) and Path(db_custom).exists()
    db_exists = bool(db) and db in db2file and db2file[db].exists()
    skip_db = (
        "--skip-db-check"
        if skip_db_check_if_exists and (db_custom_exists or db_exists)
        else ""
    )
    extra = snakemake.params.get("extra", "")
    redir = log_redirect()
    for sample, fasta in rows:
        sample_outdir = outdir / sample
        sample_outdir.mkdir(parents=True, exist_ok=True)
        shell(
            "upimapi --input {fasta} --output {sample_outdir} {db_cmd} "
            "{db_custom_cmd} {resources_cmd} --threads {snakemake.threads} "
            "{skip_db} {extra} {redir}"
        )


def run_metaeuk(rows):
    outdir = reset_directory(snakemake.output.outdir)
    db = snakemake.params.db
    extra = snakemake.params.get("extra", "")
    redir = log_redirect()
    for sample, fasta in rows:
        sample_outdir = outdir / sample
        tmp_dir = sample_outdir / "tmp"
        tmp_dir.mkdir(parents=True, exist_ok=True)
        prefix = sample_outdir / sample
        shell(
            "metaeuk easy-predict --threads {snakemake.threads} {extra} "
            "{fasta} {db} {prefix} {tmp_dir} {redir}"
        )
        expected = sample_outdir / f"{sample}.faa"
        for suffix in [".fas", ".fasta", ".faa"]:
            generated = prefix.with_suffix(suffix)
            if generated.exists():
                if generated != expected:
                    generated.rename(expected)
                break


action = snakemake.params.action

if action == "filter_prok":
    passing_bins("prok")
elif action == "filter_euk":
    passing_bins("euk")
else:
    manifest_rows = read_manifest(snakemake.input.manifest)
    if action == "bakta":
        run_bakta(manifest_rows)
    elif action == "prodigal":
        run_prodigal(manifest_rows)
    elif action == "recognizer_prok":
        run_recognizer(protein_rows_from_manifest())
    elif action == "upimapi":
        run_upimapi(protein_rows_from_manifest())
    elif action == "metaeuk":
        run_metaeuk(manifest_rows)
    elif action == "recognizer_euk":
        run_recognizer(protein_rows_from_manifest(), euk=True)
    else:
        raise ValueError(f"Unsupported batch action: {action}")

if "done" in snakemake.output:
    Path(snakemake.output.done).parent.mkdir(parents=True, exist_ok=True)
    Path(snakemake.output.done).touch()
