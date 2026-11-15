import csv
from pathlib import Path

import pytest

from common import protein_summary as protein_summary_module
from common.protein_summary import (
    BAKTA_HEADER,
    RECOGNIZER_REQUIRED,
    SUMMARY_COLUMNS,
    UPIMAPI_REQUIRED,
    build_protein_summary,
    write_protein_summary,
)

RECOGNIZER_HEADER = RECOGNIZER_REQUIRED + [
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
    "General functional category",
    "Functional category",
    "Protein description",
    "EC number",
    "KO",
]
UPIMAPI_HEADER = UPIMAPI_REQUIRED + [
    "Gene Names",
    "Protein names",
    "EC number",
    "Function [CC]",
    "Pathway",
    "Gene Ontology (GO)",
    "Cross-reference (KO)",
]


def write_table(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=header, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def summary_fixture(tmp_path, faa=">P2 second protein\nMAA\n>P1 first protein\nM\n"):
    files = {
        "faa": tmp_path / "annotation.faa",
        "bakta": tmp_path / "annotation.tsv",
        "recognizer": tmp_path / "recognizer.tsv",
        "upimapi": tmp_path / "upimapi.tsv",
    }
    files["faa"].write_text(faa, encoding="utf-8")
    files["bakta"].write_text(
        "# Annotated with Bakta\n"
        + "\t".join(BAKTA_HEADER)
        + "\n"
        + "contig1\tcds\t1\t3\t+\tP1\tgene1\tproduct 1\tUniRef:1\n"
        + "contig2\tsorf\t10\t18\t-\tP2\tgene2\tproduct 2\tCOG:2\n",
        encoding="utf-8",
    )
    write_table(files["recognizer"], RECOGNIZER_HEADER, [])
    write_table(files["upimapi"], UPIMAPI_HEADER, [])
    return files


def build(files):
    return build_protein_summary(
        files["faa"], files["bakta"], files["recognizer"], files["upimapi"]
    )


def recognizer_hit(**updates):
    row = {
        "qseqid": "P1",
        "DB ID": "COG0001",
        "CDD ID": "CDD:1",
        "pident": "95",
        "length": "100",
        "mismatch": "0",
        "gapopen": "0",
        "qstart": "1",
        "qend": "100",
        "sstart": "1",
        "send": "100",
        "evalue": "1e-20",
        "bitscore": "80",
        "General functional category": "Information",
        "Functional category": "Translation",
        "Protein description": "chosen hit",
        "EC number": "1.1.1.1",
        "KO": "K00001",
    }
    row.update(updates)
    return row


def upimapi_hit(**updates):
    row = {
        "qseqid": "P1",
        "sseqid": "sp|P1|ENTRY",
        "pident": "90",
        "length": "90",
        "mismatch": "0",
        "gapopen": "0",
        "qstart": "1",
        "qend": "90",
        "sstart": "1",
        "send": "90",
        "evalue": "1e-10",
        "bitscore": "100",
        "Entry": "P1",
        "Entry Name": "ENTRY",
        "Gene Names": "geneA",
        "Protein names": "protein A",
        "EC number": "2.2.2.2",
        "Function [CC]": "function A",
        "Pathway": "pathway A",
        "Gene Ontology (GO)": "GO:1",
        "Cross-reference (KO)": "K00002",
    }
    row.update(updates)
    return row


def test_summary_uses_exact_faa_order_lengths_and_blank_missing_evidence(tmp_path):
    files = summary_fixture(tmp_path)

    rows = build(files)

    assert [row["protein_id"] for row in rows] == ["P2", "P1"]
    assert [row["sequence_length_aa"] for row in rows] == ["3", "1"]
    assert rows[0]["bakta_contig"] == "contig2"
    assert all(
        not value
        for key, value in rows[0].items()
        if key.startswith(("recognizer_", "upimapi_"))
    )


def test_evidence_ranking_aggregation_and_output_are_row_order_independent(tmp_path):
    files = summary_fixture(tmp_path)
    recognizer_rows = [
        recognizer_hit(
            **{
                "DB ID": "COG0002",
                "CDD ID": "CDD:2",
                "evalue": "1e-10",
                "bitscore": "999",
            }
        ),
        recognizer_hit(**{"General functional category": "Zeta"}),
        recognizer_hit(**{"General functional category": "Alpha"}),
    ]
    upimapi_rows = [
        upimapi_hit(sseqid="sp|LOW|ENTRY", Entry="LOW", bitscore="90"),
        upimapi_hit(),
    ]
    write_table(files["recognizer"], RECOGNIZER_HEADER, recognizer_rows)
    write_table(files["upimapi"], UPIMAPI_HEADER, upimapi_rows)
    first = tmp_path / "first.tsv"
    write_protein_summary(
        files["faa"],
        files["bakta"],
        files["recognizer"],
        files["upimapi"],
        first,
    )

    with first.open(newline="", encoding="utf-8") as handle:
        selected = list(csv.DictReader(handle, delimiter="\t"))[1]
    assert selected == {
        "protein_id": "P1",
        "sequence_length_aa": "1",
        "bakta_contig": "contig1",
        "bakta_feature_type": "cds",
        "bakta_start": "1",
        "bakta_stop": "3",
        "bakta_strand": "+",
        "bakta_locus_tag": "P1",
        "bakta_gene": "gene1",
        "bakta_product": "product 1",
        "bakta_dbxrefs": "UniRef:1",
        "recognizer_cog_id": "COG0001",
        "recognizer_cog_general_categories": "Alpha; Zeta",
        "recognizer_cog_categories": "Translation",
        "recognizer_cog_description": "chosen hit",
        "recognizer_cog_ec": "1.1.1.1",
        "recognizer_cog_ko": "K00001",
        "recognizer_subject_id": "CDD:1",
        "recognizer_evalue": "1e-20",
        "recognizer_bitscore": "80",
        "recognizer_percent_identity": "95",
        "recognizer_alignment_length": "100",
        "recognizer_hit_count": "2",
        "upimapi_subject_id": "sp|P1|ENTRY",
        "upimapi_uniprot_entry": "P1",
        "upimapi_uniprot_entry_name": "ENTRY",
        "upimapi_gene_names": "geneA",
        "upimapi_protein_names": "protein A",
        "upimapi_ec": "2.2.2.2",
        "upimapi_function": "function A",
        "upimapi_pathway": "pathway A",
        "upimapi_go": "GO:1",
        "upimapi_ko": "K00002",
        "upimapi_evalue": "1e-10",
        "upimapi_bitscore": "100",
        "upimapi_percent_identity": "90",
        "upimapi_alignment_length": "90",
        "upimapi_hit_count": "2",
    }

    write_table(files["recognizer"], RECOGNIZER_HEADER, reversed(recognizer_rows))
    write_table(files["upimapi"], UPIMAPI_HEADER, reversed(upimapi_rows))
    second = tmp_path / "second.tsv"
    write_protein_summary(
        files["faa"],
        files["bakta"],
        files["recognizer"],
        files["upimapi"],
        second,
    )
    assert first.read_bytes() == second.read_bytes()


def test_source_columns_remain_separate_without_reviewed_hit_preference(tmp_path):
    files = summary_fixture(tmp_path)
    write_table(files["recognizer"], RECOGNIZER_HEADER, [recognizer_hit()])
    write_table(
        files["upimapi"],
        UPIMAPI_HEADER,
        [
            upimapi_hit(
                sseqid="sp|REVIEWED|ENTRY",
                Entry="REVIEWED",
                bitscore="90",
            ),
            upimapi_hit(
                sseqid="tr|TREMBL|ENTRY",
                Entry="TREMBL",
                bitscore="110",
            ),
        ],
    )

    rows = build(files)

    assert len(rows) == 2
    selected = rows[1]
    assert selected["protein_id"] == "P1"
    assert selected["bakta_product"] == "product 1"
    assert selected["recognizer_cog_id"] == "COG0001"
    assert selected["upimapi_subject_id"] == "tr|TREMBL|ENTRY"
    assert selected["upimapi_uniprot_entry"] == "TREMBL"


def test_summary_reads_each_evidence_table_once(tmp_path, monkeypatch):
    files = summary_fixture(tmp_path)
    write_table(files["recognizer"], RECOGNIZER_HEADER, [recognizer_hit()])
    write_table(files["upimapi"], UPIMAPI_HEADER, [upimapi_hit()])
    original_read_table = protein_summary_module._read_table
    reads = []

    def counted_read_table(path, **kwargs):
        reads.append(Path(path))
        return original_read_table(path, **kwargs)

    monkeypatch.setattr(protein_summary_module, "_read_table", counted_read_table)

    build(files)

    assert reads.count(files["recognizer"]) == 1
    assert reads.count(files["upimapi"]) == 1


def test_summary_validates_both_evidence_headers_before_eligible_rows(tmp_path):
    files = summary_fixture(tmp_path)
    write_table(
        files["recognizer"],
        RECOGNIZER_HEADER,
        [recognizer_hit(bitscore="not-a-number")],
    )
    write_table(files["upimapi"], ["qseqid"], [])

    with pytest.raises(ValueError, match="UPIMAPI.*missing required column"):
        build(files)


@pytest.mark.parametrize(
    ("producer", "header"),
    [
        ("recognizer", ["qseqid"]),
        ("upimapi", ["qseqid", "sseqid"]),
    ],
)
def test_summary_rejects_malformed_evidence_headers(tmp_path, producer, header):
    files = summary_fixture(tmp_path)
    write_table(files[producer], header, [])

    with pytest.raises(ValueError, match="missing required column"):
        build(files)


@pytest.mark.parametrize("producer", ["recognizer", "upimapi"])
def test_summary_rejects_evidence_for_unknown_proteins(tmp_path, producer):
    files = summary_fixture(tmp_path)
    if producer == "recognizer":
        write_table(
            files[producer],
            RECOGNIZER_HEADER,
            [recognizer_hit(qseqid="UNKNOWN")],
        )
    else:
        write_table(files[producer], UPIMAPI_HEADER, [upimapi_hit(qseqid="UNKNOWN")])

    with pytest.raises(ValueError, match="absent from the Bakta FAA"):
        build(files)


@pytest.mark.parametrize(
    ("producer", "row", "message"),
    [
        ("recognizer", recognizer_hit(**{"CDD ID": ""}), "blank core hit"),
        ("recognizer", recognizer_hit(bitscore="nan"), "non-finite numeric"),
        ("upimapi", upimapi_hit(Entry=""), "blank core hit"),
        ("upimapi", upimapi_hit(evalue="malformed"), "malformed numeric"),
    ],
)
def test_summary_rejects_malformed_eligible_evidence(tmp_path, producer, row, message):
    files = summary_fixture(tmp_path)
    header = RECOGNIZER_HEADER if producer == "recognizer" else UPIMAPI_HEADER
    write_table(files[producer], header, [row])

    with pytest.raises(ValueError, match=message):
        build(files)


def test_summary_rejects_ambiguous_recognizer_full_header_alias(tmp_path):
    files = summary_fixture(
        tmp_path,
        faa=">P1 P2\nM\n>P1_P2\nMM\n",
    )
    files["bakta"].write_text(
        "# Annotated with Bakta\n"
        + "\t".join(BAKTA_HEADER)
        + "\n"
        + "contig\tcds\t1\t3\t+\tP1\t\tproduct\t\n"
        + "contig\tcds\t4\t9\t+\tP1_P2\t\tproduct\t\n",
        encoding="utf-8",
    )
    write_table(
        files["recognizer"],
        RECOGNIZER_HEADER,
        [recognizer_hit(qseqid="P1_P2")],
    )

    with pytest.raises(ValueError, match="ambiguous qseqid"):
        build(files)


@pytest.mark.parametrize(
    ("faa", "message"),
    [
        (">P1\nM\n>P1\nM\n", "Duplicate protein identifier"),
        (">MISSING\nM\n", "coverage mismatch"),
    ],
)
def test_summary_rejects_duplicate_or_mismatched_faa_material(tmp_path, faa, message):
    files = summary_fixture(tmp_path, faa=faa)

    with pytest.raises(ValueError, match=message):
        build(files)


def test_empty_faa_writes_valid_header_only_summary(tmp_path):
    files = summary_fixture(tmp_path, faa="")
    files["bakta"].write_text(
        "# Annotated with Bakta\n" + "\t".join(BAKTA_HEADER) + "\n",
        encoding="utf-8",
    )
    output = tmp_path / "summary.tsv"

    assert (
        write_protein_summary(
            files["faa"],
            files["bakta"],
            files["recognizer"],
            files["upimapi"],
            output,
        )
        == []
    )
    assert output.read_text(encoding="utf-8") == "\t".join(SUMMARY_COLUMNS) + "\n"
