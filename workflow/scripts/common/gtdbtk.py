"""Scratch-stage inputs and consolidate native GTDB-Tk summaries."""

import csv
import tempfile
from contextlib import contextmanager
from pathlib import Path

from common.identifiers import validate_annotation_identifier


def _validated_staging_pairs(pairs, declared_fastas):
    pairs = [
        (validate_annotation_identifier(str(sample)), str(source))
        for sample, source in pairs
    ]
    if not pairs:
        raise ValueError(
            "GTDB-Tk staging requires at least one eligible prokaryotic genome."
        )

    sample_ids = [sample for sample, _ in pairs]
    duplicates = sorted(
        {sample for sample in sample_ids if sample_ids.count(sample) > 1}
    )
    if duplicates:
        raise ValueError(
            "GTDB-Tk staging received duplicate genome identifiers: "
            + ", ".join(duplicates)
        )

    paired_fastas = [source for _, source in pairs]
    if paired_fastas != [str(path) for path in declared_fastas]:
        raise ValueError(
            "GTDB-Tk eligible sample/path pairs do not match declared FASTA inputs."
        )

    validated = []
    for sample, source in pairs:
        source_path = Path(source)
        try:
            resolved = source_path.resolve(strict=True)
        except OSError as error:
            raise FileNotFoundError(
                f"GTDB-Tk source FASTA is not a usable regular file: {source_path}."
            ) from error
        if not resolved.is_file():
            raise FileNotFoundError(
                f"GTDB-Tk source FASTA is not a usable regular file: {source_path}."
            )
        try:
            with resolved.open("rb"):
                pass
        except OSError as error:
            raise OSError(
                f"GTDB-Tk source FASTA is not readable: {source_path}."
            ) from error
        validated.append((sample, resolved))
    return validated


@contextmanager
def temporary_gtdbtk_genomes(pairs, declared_fastas, tmpdir):
    """Symlink exactly the declared genomes in unique, automatically cleaned scratch."""
    validated = _validated_staging_pairs(pairs, declared_fastas)
    scratch_root = Path(tmpdir).resolve()
    if " " in str(scratch_root):
        raise ValueError(
            "GTDB-Tk 2.7.1 cannot use genome paths containing an ASCII space: "
            f"resources.tmpdir resolved to {scratch_root}. Select a resources.tmpdir "
            "path without spaces."
        )
    scratch_root.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="gtdbtk-", dir=scratch_root) as workdir:
        genomes = Path(workdir)
        for sample, source in validated:
            (genomes / f"{sample}.fasta").symlink_to(source)
        yield genomes


def _read_summary(domain, path):
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        fieldnames = reader.fieldnames or []
        if not fieldnames:
            raise ValueError(f"GTDB-Tk summary has no header: {path}.")
        if any(not field for field in fieldnames) or len(fieldnames) != len(
            set(fieldnames)
        ):
            raise ValueError(
                f"GTDB-Tk summary has invalid or duplicate columns: {path}."
            )
        if "user_genome" not in fieldnames:
            raise ValueError(
                f"GTDB-Tk summary is missing required user_genome column: {path}."
            )

        rows = []
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(
                    f"GTDB-Tk summary has a malformed row at {path}:{line_number}."
                )
            genome_id = row["user_genome"].strip()
            if not genome_id:
                raise ValueError(
                    f"GTDB-Tk summary has a blank user_genome at "
                    f"{path}:{line_number}."
                )
            rows.append({"domain": domain, **row})
    return fieldnames, rows


def consolidate_gtdbtk_summaries(summaries, output, expected_ids):
    """Publish one validated table for exactly the staged genome identifiers."""
    expected_ids = [str(genome_id) for genome_id in expected_ids]
    duplicate_expected = sorted(
        {genome_id for genome_id in expected_ids if expected_ids.count(genome_id) > 1}
    )
    if duplicate_expected:
        raise ValueError(
            "GTDB-Tk staged identifiers are duplicated: "
            + ", ".join(duplicate_expected)
        )
    if not expected_ids:
        raise ValueError("GTDB-Tk classification requires at least one staged genome.")

    available = [
        (str(domain), Path(path)) for domain, path in summaries if Path(path).is_file()
    ]
    if not available:
        summary_dir = Path(summaries[0][1]).parent if summaries else Path(".")
        raise FileNotFoundError(
            "GTDB-Tk completed without a bacterial or archaeal summary in "
            f"{summary_dir}."
        )

    fieldnames = None
    rows = []
    for domain, path in available:
        current_fields, current_rows = _read_summary(domain, path)
        if fieldnames is None:
            fieldnames = current_fields
        elif current_fields != fieldnames:
            raise ValueError(
                f"GTDB-Tk bacterial/archaeal summary columns differ: {path}."
            )
        rows.extend(current_rows)

    classified_ids = [row["user_genome"] for row in rows]
    duplicate_classified = sorted(
        {
            genome_id
            for genome_id in classified_ids
            if classified_ids.count(genome_id) > 1
        }
    )
    if duplicate_classified:
        raise ValueError(
            "GTDB-Tk classified genome identifiers are duplicated: "
            + ", ".join(duplicate_classified)
        )

    expected = set(expected_ids)
    classified = set(classified_ids)
    missing = sorted(expected - classified)
    unexpected = sorted(classified - expected)
    if missing or unexpected:
        raise ValueError(
            "GTDB-Tk classified identifiers do not match staged genomes. "
            f"Missing: {missing}; unexpected: {unexpected}."
        )

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["domain"] + fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
