"""Deterministic Bakta-anchored prokaryotic protein annotation summaries."""

import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

BAKTA_HEADER = [
    "#Sequence Id",
    "Type",
    "Start",
    "Stop",
    "Strand",
    "Locus Tag",
    "Gene",
    "Product",
    "DbXrefs",
]
RECOGNIZER_REQUIRED = [
    "qseqid",
    "DB ID",
]
UPIMAPI_REQUIRED = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "Entry",
    "Entry Name",
]

SUMMARY_COLUMNS = [
    "protein_id",
    "sequence_length_aa",
    "bakta_contig",
    "bakta_feature_type",
    "bakta_start",
    "bakta_stop",
    "bakta_strand",
    "bakta_locus_tag",
    "bakta_gene",
    "bakta_product",
    "bakta_dbxrefs",
    "recognizer_cog_id",
    "recognizer_cog_general_categories",
    "recognizer_cog_categories",
    "recognizer_cog_description",
    "recognizer_cog_ec",
    "recognizer_cog_ko",
    "recognizer_subject_id",
    "recognizer_evalue",
    "recognizer_bitscore",
    "recognizer_percent_identity",
    "recognizer_alignment_length",
    "recognizer_hit_count",
    "upimapi_subject_id",
    "upimapi_uniprot_entry",
    "upimapi_uniprot_entry_name",
    "upimapi_gene_names",
    "upimapi_protein_names",
    "upimapi_ec",
    "upimapi_function",
    "upimapi_pathway",
    "upimapi_go",
    "upimapi_ko",
    "upimapi_evalue",
    "upimapi_bitscore",
    "upimapi_percent_identity",
    "upimapi_alignment_length",
    "upimapi_hit_count",
]

_RECOGNIZER_AGGREGATES = {
    "recognizer_cog_general_categories": "General functional category",
    "recognizer_cog_categories": "Functional category",
    "recognizer_cog_description": "Protein description",
    "recognizer_cog_ec": "EC number",
    "recognizer_cog_ko": "KO",
}
_RECOGNIZER_COG_ANNOTATION_COLUMNS = tuple(_RECOGNIZER_AGGREGATES.values())
_RECOGNIZER_COG_REQUIRED_COLUMNS = (
    "CDD ID",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    *_RECOGNIZER_COG_ANNOTATION_COLUMNS,
)
_ALIGNMENT_NUMERIC_COLUMNS = (
    "evalue",
    "bitscore",
    "pident",
    "length",
    "qstart",
    "qend",
    "sstart",
    "send",
    "mismatch",
    "gapopen",
)
_RECOGNIZER_ALIGNMENT_COLUMNS = (
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
)
_RECOGNIZER_HIT_IDENTITY_COLUMNS = (
    "DB ID",
    "CDD ID",
    "qstart",
    "qend",
    "sstart",
    "send",
    "mismatch",
    "gapopen",
    "pident",
    "length",
    "evalue",
    "bitscore",
)
_UPIMAPI_ALIASES = {
    "upimapi_gene_names": ("Gene Names", "Gene names"),
    "upimapi_protein_names": ("Protein names", "Protein Names"),
    "upimapi_ec": ("EC number", "EC Number"),
    "upimapi_function": ("Function [CC]",),
    "upimapi_pathway": ("Pathway",),
    "upimapi_go": ("Gene Ontology (GO)", "Gene ontology (GO)"),
    "upimapi_ko": ("Cross-reference (KO)", "KO"),
}


def _read_fasta_records_with_headers(path):
    """Read ordered first-token IDs, lengths, and full protein FASTA headers."""
    records = []
    identifier = None
    header = None
    length = 0
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if line.startswith(">"):
                if identifier is not None:
                    records.append((identifier, length, header))
                header = line[1:].strip()
                identifier = header.split(maxsplit=1)[0] if header else ""
                if not identifier:
                    raise ValueError(
                        f"Blank FASTA identifier in {path} on line {line_number}."
                    )
                length = 0
            elif line:
                if identifier is None:
                    raise ValueError(
                        f"Protein sequence before first FASTA header in {path} on line "
                        f"{line_number}."
                    )
                length += len("".join(line.split()))
    if identifier is not None:
        records.append((identifier, length, header))
    duplicates = sorted(
        key
        for key, count in Counter(identifier for identifier, _, _ in records).items()
        if count > 1
    )
    if duplicates:
        raise ValueError(
            f"Duplicate protein identifier(s) in {path}: {', '.join(duplicates)}."
        )
    return records


def _recognizer_query_aliases(fasta_records):
    """Map native reCOGnizer full-header qseqids to canonical Bakta IDs."""
    candidates = defaultdict(set)
    for identifier, _, header in fasta_records:
        candidates[identifier].add(identifier)
        candidates[header.replace(" ", "_")].add(identifier)

    aliases = {}
    ambiguous = {}
    for query, identifiers in candidates.items():
        if len(identifiers) == 1:
            aliases[query] = next(iter(identifiers))
        else:
            ambiguous[query] = sorted(identifiers)
    return aliases, ambiguous


def _read_table(path, *, producer, required, commented_header=False):
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as handle:
        if commented_header:
            rows = [line for line in handle if line.startswith("#Sequence Id\t")]
            if len(rows) != 1:
                raise ValueError(
                    f"{producer} table {path} must contain exactly one native "
                    "'#Sequence Id' header."
                )
            header = rows[0].rstrip("\r\n").split("\t")
            if header != required:
                raise ValueError(
                    f"{producer} table {path} has malformed header: expected "
                    f"{required}, found {header}."
                )
            handle.seek(0)
            data_lines = (
                line for line in handle if line.strip() and not line.startswith("#")
            )
            reader = csv.DictReader(
                data_lines, fieldnames=header, delimiter="\t", strict=True
            )
            parsed = list(reader)
            if any(None in row.values() or None in row for row in parsed):
                raise ValueError(f"Malformed {producer} table row in {path}.")
            return header, parsed

        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        header = reader.fieldnames or []
        if not header:
            raise ValueError(f"{producer} table {path} has no header.")
        duplicates = sorted(key for key, count in Counter(header).items() if count > 1)
        if duplicates:
            raise ValueError(
                f"{producer} table {path} has duplicate header(s): "
                f"{', '.join(duplicates)}."
            )
        missing = [column for column in required if column not in header]
        if missing:
            raise ValueError(
                f"{producer} table {path} is missing required column(s): "
                f"{', '.join(missing)}."
            )
        parsed = list(reader)
        if any(None in row.values() or None in row for row in parsed):
            raise ValueError(f"Malformed {producer} table row in {path}.")
        return header, parsed


def read_bakta_proteins(path):
    """Read unique protein-producing Bakta rows indexed by native locus tag."""
    _, rows = _read_table(
        path, producer="Bakta", required=BAKTA_HEADER, commented_header=True
    )
    proteins = {}
    for row_number, row in enumerate(rows, start=1):
        if row["Type"] not in {"cds", "sorf"}:
            continue
        locus = row["Locus Tag"]
        if not locus:
            raise ValueError(
                f"Bakta protein row {row_number} in {path} has a blank Locus Tag."
            )
        if locus in proteins:
            raise ValueError(f"Duplicate Bakta protein locus tag '{locus}' in {path}.")
        proteins[locus] = row
    return proteins


def _number(value, *, producer, column, path, row_number):
    value = value.strip()
    if not value:
        return None
    try:
        number = float(value)
    except ValueError as error:
        raise ValueError(
            f"{producer} table {path} row {row_number} has malformed numeric "
            f"{column} value {value!r}."
        ) from error
    if not math.isfinite(number):
        raise ValueError(
            f"{producer} table {path} row {row_number} has non-finite numeric "
            f"{column} value {value!r}."
        )
    return number


def _canonical(row, header):
    return "\t".join(row.get(column, "") or "" for column in header)


def _missing_last_ascending(value):
    return (value is None, value if value is not None else 0.0)


def _missing_last_descending(value):
    return (value is None, -value if value is not None else 0.0)


def _prepare_rows(
    path,
    producer,
    header,
    rows,
    valid_ids,
    *,
    eligible=None,
    eligible_header_required=(),
    eligible_required=(),
    numeric_columns=_ALIGNMENT_NUMERIC_COLUMNS,
    query_aliases=None,
    ambiguous_query_aliases=None,
):
    prepared = defaultdict(list)
    unknown = set()
    ambiguous = set()
    if query_aliases is None:
        query_aliases = {identifier: identifier for identifier in valid_ids}
    if ambiguous_query_aliases is None:
        ambiguous_query_aliases = {}
    for row in rows:
        query = row["qseqid"]
        if not query.strip():
            unknown.add("<blank>")
        elif query in ambiguous_query_aliases:
            ambiguous.add(query)
        elif query not in query_aliases:
            unknown.add(query)
    if ambiguous:
        details = "; ".join(
            f"{query} -> {', '.join(ambiguous_query_aliases[query])}"
            for query in sorted(ambiguous)
        )
        raise ValueError(
            f"{producer} reported ambiguous qseqid value(s) for the Bakta FAA: "
            f"{details}. Table: {path}."
        )
    if unknown:
        raise ValueError(
            f"{producer} reported qseqid value(s) absent from the Bakta FAA: "
            f"{', '.join(sorted(unknown))}. Table: {path}."
        )

    has_eligible_rows = eligible is None or any(eligible(row) for row in rows)
    if has_eligible_rows:
        missing_eligible_headers = [
            column for column in eligible_header_required if column not in header
        ]
        if missing_eligible_headers:
            raise ValueError(
                f"{producer} table {path} is missing required column(s) for "
                f"eligible evidence: {', '.join(missing_eligible_headers)}."
            )

    for row_number, row in enumerate(rows, start=2):
        query = row["qseqid"]
        is_eligible = eligible is None or eligible(row)
        if not is_eligible:
            continue
        blank_required = [
            column for column in eligible_required if not row.get(column, "").strip()
        ]
        if blank_required:
            raise ValueError(
                f"{producer} table {path} row {row_number} has blank core hit "
                f"identifier(s): {', '.join(blank_required)}."
            )
        missing_numeric_headers = [
            column for column in numeric_columns if column not in header
        ]
        if missing_numeric_headers:
            raise ValueError(
                f"{producer} table {path} is missing alignment column(s) required "
                f"for eligible evidence on row {row_number}: "
                f"{', '.join(missing_numeric_headers)}."
            )
        numeric = {
            column: _number(
                row.get(column, ""),
                producer=producer,
                column=column,
                path=path,
                row_number=row_number,
            )
            for column in numeric_columns
        }
        prepared[query_aliases[query]].append((row, numeric, _canonical(row, header)))
    return prepared


def _distinct_hit_count(rows, columns):
    """Count distinct reported alternatives by their source identity columns."""
    return len({tuple(row.get(column, "") for column in columns) for row, _, _ in rows})


def _recognizer_summary(rows):
    if not rows:
        return {}
    selected, numbers, _ = min(
        rows,
        key=lambda item: (
            _missing_last_ascending(item[1]["evalue"]),
            _missing_last_descending(item[1]["bitscore"]),
            _missing_last_descending(item[1]["pident"]),
            _missing_last_descending(item[1]["length"]),
            item[0]["DB ID"],
            item[0]["CDD ID"],
            item[2],
        ),
    )
    alignment_key = tuple(
        selected[column] for column in _RECOGNIZER_HIT_IDENTITY_COLUMNS
    )
    duplicate_alignment_rows = [
        row
        for row, _, _ in rows
        if tuple(row[column] for column in _RECOGNIZER_HIT_IDENTITY_COLUMNS)
        == alignment_key
    ]
    summary = {
        "recognizer_cog_id": selected["DB ID"],
        "recognizer_subject_id": selected["CDD ID"],
        "recognizer_evalue": selected["evalue"],
        "recognizer_bitscore": selected["bitscore"],
        "recognizer_percent_identity": selected["pident"],
        "recognizer_alignment_length": selected["length"],
        "recognizer_hit_count": str(
            _distinct_hit_count(rows, _RECOGNIZER_HIT_IDENTITY_COLUMNS)
        ),
    }
    for output_column, source_column in _RECOGNIZER_AGGREGATES.items():
        summary[output_column] = "; ".join(
            sorted(
                {
                    row.get(source_column, "").strip()
                    for row in duplicate_alignment_rows
                    if row.get(source_column, "").strip()
                }
            )
        )
    return summary


def _first_alias(row, aliases):
    return next((row[column] for column in aliases if row.get(column, "")), "")


def _upimapi_summary(rows):
    if not rows:
        return {}
    selected, _, _ = min(
        rows,
        key=lambda item: (
            _missing_last_descending(item[1]["bitscore"]),
            _missing_last_ascending(item[1]["evalue"]),
            _missing_last_descending(item[1]["pident"]),
            _missing_last_descending(item[1]["length"]),
            item[0]["sseqid"],
            item[0]["Entry"],
            item[2],
        ),
    )
    summary = {
        "upimapi_subject_id": selected["sseqid"],
        "upimapi_uniprot_entry": selected["Entry"],
        "upimapi_uniprot_entry_name": selected["Entry Name"],
        "upimapi_evalue": selected["evalue"],
        "upimapi_bitscore": selected["bitscore"],
        "upimapi_percent_identity": selected["pident"],
        "upimapi_alignment_length": selected["length"],
        "upimapi_hit_count": str(_distinct_hit_count(rows, ("sseqid", "Entry"))),
    }
    for output_column, aliases in _UPIMAPI_ALIASES.items():
        summary[output_column] = _first_alias(selected, aliases)
    return summary


def build_protein_summary(faa, bakta, recognizer, upimapi):
    """Build summary rows, anchored to exact Bakta FAA order and coverage."""
    fasta_records_with_headers = _read_fasta_records_with_headers(faa)
    fasta_records = [
        (identifier, length) for identifier, length, _ in fasta_records_with_headers
    ]
    protein_ids = {identifier for identifier, _ in fasta_records}
    recognizer_aliases, ambiguous_recognizer_aliases = _recognizer_query_aliases(
        fasta_records_with_headers
    )
    bakta_rows = read_bakta_proteins(bakta)
    missing = [
        identifier for identifier, _ in fasta_records if identifier not in bakta_rows
    ]
    foreign = sorted(set(bakta_rows) - protein_ids)
    if missing or foreign:
        raise ValueError(
            "Bakta protein TSV/FAA coverage mismatch. "
            f"Missing Bakta protein records: {missing}; "
            f"Bakta protein records absent from FAA: {foreign}. TSV: {bakta}; FAA: {faa}."
        )

    recognizer_header, recognizer_table_rows = _read_table(
        recognizer, producer="reCOGnizer", required=RECOGNIZER_REQUIRED
    )
    upimapi_header, upimapi_table_rows = _read_table(
        upimapi, producer="UPIMAPI", required=UPIMAPI_REQUIRED
    )

    recognizer_rows = _prepare_rows(
        recognizer,
        "reCOGnizer",
        recognizer_header,
        recognizer_table_rows,
        protein_ids,
        eligible=lambda row: row["DB ID"].startswith("COG"),
        eligible_header_required=_RECOGNIZER_COG_REQUIRED_COLUMNS,
        eligible_required=("CDD ID",),
        numeric_columns=_RECOGNIZER_ALIGNMENT_COLUMNS,
        query_aliases=recognizer_aliases,
        ambiguous_query_aliases=ambiguous_recognizer_aliases,
    )
    upimapi_rows = _prepare_rows(
        upimapi,
        "UPIMAPI",
        upimapi_header,
        upimapi_table_rows,
        protein_ids,
        eligible_required=("sseqid", "Entry"),
    )

    output_rows = []
    for identifier, length in fasta_records:
        bakta_row = bakta_rows[identifier]
        row = {column: "" for column in SUMMARY_COLUMNS}
        row.update(
            {
                "protein_id": identifier,
                "sequence_length_aa": str(length),
                "bakta_contig": bakta_row["#Sequence Id"],
                "bakta_feature_type": bakta_row["Type"],
                "bakta_start": bakta_row["Start"],
                "bakta_stop": bakta_row["Stop"],
                "bakta_strand": bakta_row["Strand"],
                "bakta_locus_tag": bakta_row["Locus Tag"],
                "bakta_gene": bakta_row["Gene"],
                "bakta_product": bakta_row["Product"],
                "bakta_dbxrefs": bakta_row["DbXrefs"],
            }
        )
        row.update(_recognizer_summary(recognizer_rows[identifier]))
        row.update(_upimapi_summary(upimapi_rows[identifier]))
        output_rows.append(row)
    return output_rows


def write_protein_summary(faa, bakta, recognizer, upimapi, output):
    rows = build_protein_summary(faa, bakta, recognizer, upimapi)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    return rows
