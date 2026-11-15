"""Strict one-record input contracts for MeeW-Annotation."""

import csv
from pathlib import Path

from common.empty_outputs import FASTA_SUFFIXES
from common.identifiers import validate_annotation_identifier
from common.paths import expand_path

INPUT_HEADERS = {
    "isolate": ["sample", "domain", "path"],
    "metagenome": ["sample", "prok_bins_dir", "euk_bins_dir"],
}
NUCLEOTIDE_CHARACTERS = set("ACGTURYSWKMBDHVNacgturyswkmbdhvn-.")


def _expected_header(mode):
    return "\t".join(INPUT_HEADERS[mode])


def _validate_header(header, mode, sheet):
    expected = INPUT_HEADERS[mode]
    if not header or all(value == "" for value in header):
        raise ValueError(f"Input samplesheet '{sheet}' has a blank header record.")

    duplicates = sorted({value for value in header if header.count(value) > 1})
    if duplicates:
        raise ValueError(
            f"Input samplesheet '{sheet}' contains duplicate column(s): "
            + ", ".join(duplicates)
            + f". Expected exact {mode} header: {_expected_header(mode)}"
        )

    unknown = sorted(set(header) - set(expected))
    missing = sorted(set(expected) - set(header))
    details = []
    if unknown:
        details.append("unknown column(s): " + ", ".join(unknown))
    if missing:
        details.append("missing column(s): " + ", ".join(missing))
    if details:
        raise ValueError(
            f"Input samplesheet '{sheet}' has "
            + "; ".join(details)
            + f". Expected exact {mode} header: {_expected_header(mode)}"
        )
    if header != expected:
        raise ValueError(
            f"Input samplesheet '{sheet}' columns are in the wrong order. "
            f"Expected exact {mode} header: {_expected_header(mode)}"
        )


def _read_exact_record(sheet, mode):
    path = Path(sheet)
    if not path.exists():
        raise ValueError(f"Input samplesheet does not exist: {path}")
    if not path.is_file():
        raise ValueError(f"Input samplesheet is not a regular file: {path}")

    try:
        with path.open(newline="", encoding="utf-8") as handle:
            records = list(csv.reader(handle, delimiter="\t", strict=True))
    except csv.Error as error:
        raise ValueError(
            f"Malformed TSV input samplesheet '{path}': {error}"
        ) from error

    if any(not record or all(value == "" for value in record) for record in records):
        raise ValueError(
            f"Input samplesheet '{path}' contains a blank logical record; "
            "exactly one header and one data record are required."
        )
    if not records:
        raise ValueError(
            f"Input samplesheet '{path}' is empty; exactly one header and one data "
            "record are required."
        )

    if len(records) != 2:
        found = max(0, len(records) - 1)
        raise ValueError(
            f"Input samplesheet '{path}' must contain exactly one data record; "
            f"found {found}."
        )

    _validate_header(records[0], mode, path)
    row = records[1]
    expected_width = len(records[0])
    if len(row) > expected_width:
        raise ValueError(
            f"Input samplesheet '{path}' data record has too many fields: "
            f"expected {expected_width}, found {len(row)}."
        )
    if len(row) < expected_width:
        raise ValueError(
            f"Input samplesheet '{path}' data record has too few fields: "
            f"expected {expected_width}, found {len(row)}."
        )
    return dict(zip(records[0], row))


def _validated_sample(value, sheet):
    try:
        return validate_annotation_identifier(value)
    except ValueError as error:
        raise ValueError(
            f"Input samplesheet '{sheet}' has invalid sample identifier {value!r}; "
            "use an ASCII identifier accepted by the Annotation identifier contract."
        ) from error


def _validate_isolate(record, sheet):
    normalized = {
        "sample": record["sample"],
        "domain": record["domain"].strip(),
        "path": expand_path(record["path"].strip()),
    }
    blank = [
        field
        for field, value in normalized.items()
        if not value or not str(value).strip()
    ]
    if blank:
        raise ValueError(
            f"Isolate samplesheet '{sheet}' has blank required field(s): "
            + ", ".join(blank)
            + "."
        )
    normalized["sample"] = _validated_sample(normalized["sample"], sheet)
    if normalized["domain"] not in {"prok", "euk"}:
        raise ValueError(
            f"Isolate samplesheet '{sheet}' has invalid domain "
            f"{normalized['domain']!r}; expected prok or euk."
        )

    fasta = Path(normalized["path"])
    if not fasta.exists():
        raise ValueError(f"Isolate FASTA does not exist: {fasta}")
    if not fasta.is_file():
        raise ValueError(f"Isolate FASTA is not a regular file: {fasta}")
    if fasta.suffix.lower() not in FASTA_SUFFIXES:
        raise ValueError(
            f"Isolate FASTA has unsupported suffix: {fasta}. Expected an "
            "uncompressed .fa, .fna, or .fasta file."
        )
    _validate_nucleotide_fasta(fasta)
    return normalized


def _validate_nucleotide_fasta(path):
    has_header = False
    has_sequence = False
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            value = line.strip()
            if not value:
                continue
            if value.startswith(">"):
                if len(value) == 1:
                    raise ValueError(
                        f"Isolate FASTA has a blank header at {path}:{line_number}."
                    )
                has_header = True
                continue
            if not has_header:
                raise ValueError(
                    f"Isolate FASTA sequence precedes its first header at "
                    f"{path}:{line_number}."
                )
            invalid = sorted(set(value) - NUCLEOTIDE_CHARACTERS)
            if invalid:
                raise ValueError(
                    f"Isolate FASTA contains non-nucleotide character(s) at "
                    f"{path}:{line_number}: {''.join(invalid)!r}."
                )
            has_sequence = True
    if not has_header or not has_sequence:
        raise ValueError(
            f"Isolate FASTA must contain at least one header and nucleotide "
            f"sequence: {path}."
        )


def _bin_ids(directory, domain):
    identifiers = {}
    for fasta in sorted(
        path
        for path in directory.iterdir()
        if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES
    ):
        identifier = _validated_sample(fasta.stem, fasta)
        identifiers.setdefault(identifier, []).append((domain, fasta))
    return identifiers


def _validate_metagenome(record, sheet):
    normalized = {
        "sample": record["sample"],
        "prok_bins_dir": expand_path(record["prok_bins_dir"].strip()),
        "euk_bins_dir": expand_path(record["euk_bins_dir"].strip()),
    }
    if not normalized["sample"] or not normalized["sample"].strip():
        raise ValueError(f"Metagenome samplesheet '{sheet}' has blank sample field.")
    normalized["sample"] = _validated_sample(normalized["sample"], sheet)
    if not normalized["prok_bins_dir"] and not normalized["euk_bins_dir"]:
        raise ValueError(
            f"Metagenome samplesheet '{sheet}' must provide at least one nonblank "
            "prok_bins_dir or euk_bins_dir."
        )

    all_identifiers = {}
    for domain in ("prok", "euk"):
        field = f"{domain}_bins_dir"
        value = normalized[field]
        if not value:
            continue
        directory = Path(value)
        if not directory.exists():
            raise ValueError(f"Metagenome {field} does not exist: {directory}")
        if not directory.is_dir():
            raise ValueError(f"Metagenome {field} is not a directory: {directory}")
        for identifier, sources in _bin_ids(directory, domain).items():
            all_identifiers.setdefault(identifier, []).extend(sources)

    duplicates = {
        identifier: sources
        for identifier, sources in all_identifiers.items()
        if len(sources) > 1
    }
    if duplicates:
        details = "; ".join(
            f"{identifier}: " + ", ".join(str(path) for _, path in sources)
            for identifier, sources in sorted(duplicates.items())
        )
        raise ValueError(
            "Metagenome input contains duplicate filename-derived bin identifier(s), "
            f"including accepted-suffix collisions: {details}."
        )
    return normalized


def load_input_record(mode, samplesheet_path):
    """Load and manually validate one isolate or metagenome TSV record."""
    if mode not in INPUT_HEADERS:
        raise ValueError(
            f"Invalid input.mode {mode!r}; expected isolate or metagenome."
        )
    sheet = expand_path(samplesheet_path)
    record = _read_exact_record(sheet, mode)
    if mode == "isolate":
        return _validate_isolate(record, sheet)
    return _validate_metagenome(record, sheet)
