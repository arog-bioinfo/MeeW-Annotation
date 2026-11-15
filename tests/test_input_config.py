from pathlib import Path

import pytest

from common.empty_outputs import fasta_paths
from common.input_config import load_input_record


def write_sheet(path, header, row):
    path.write_text("\t".join(header) + "\n" + "\t".join(row) + "\n", encoding="utf-8")
    return path


def isolate_sheet(tmp_path, fasta, *, sample="isolate", domain="prok"):
    return write_sheet(
        tmp_path / "isolate.tsv",
        ["sample", "domain", "path"],
        [sample, domain, str(fasta)],
    )


@pytest.mark.parametrize("suffix", [".fa", ".fna", ".fasta", ".FA", ".FNA", ".FASTA"])
def test_isolate_accepts_supported_uncompressed_nucleotide_fastas(tmp_path, suffix):
    fasta = tmp_path / f"genome{suffix}"
    fasta.write_text(">contig one\nACGTURYSWKMBDHVN-.\n", encoding="utf-8")

    assert load_input_record("isolate", isolate_sheet(tmp_path, fasta)) == {
        "sample": "isolate",
        "domain": "prok",
        "path": str(fasta),
    }


@pytest.mark.parametrize(
    ("header", "message"),
    [
        (["sample", "domain", "domain"], "duplicate column"),
        (["sample", "path"], "missing column"),
        (["sample", "domain", "path", "extra"], "unknown column"),
        (["domain", "sample", "path"], "wrong order"),
    ],
)
def test_isolate_requires_exact_header_and_order(tmp_path, header, message):
    sheet = write_sheet(tmp_path / "bad.tsv", header, ["sample", "prok", "input.fa"])

    with pytest.raises(ValueError, match=message):
        load_input_record("isolate", sheet)


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("sample\tdomain\tpath\n", "exactly one data record"),
        (
            "sample\tdomain\tpath\na\tprok\ta.fa\nb\tprok\tb.fa\n",
            "exactly one data record",
        ),
        ("sample\tdomain\tpath\na\tprok\ta.fa\n\n", "blank logical record"),
        ("sample\tdomain\tpath\na\tprok\n", "too few fields"),
        ("sample\tdomain\tpath\na\tprok\ta.fa\textra\n", "too many fields"),
    ],
)
def test_input_sheet_requires_exactly_one_well_formed_record(tmp_path, text, message):
    sheet = tmp_path / "bad.tsv"
    sheet.write_text(text, encoding="utf-8")

    with pytest.raises(ValueError, match=message):
        load_input_record("isolate", sheet)


@pytest.mark.parametrize("sample", ["sample_1", "A.b-c", "A" + "b" * 127])
def test_input_preserves_valid_annotation_identifiers(tmp_path, sample):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">contig\nACGT\n", encoding="utf-8")

    assert (
        load_input_record("isolate", isolate_sheet(tmp_path, fasta, sample=sample))[
            "sample"
        ]
        == sample
    )


@pytest.mark.parametrize("sample", ["", ".sample", "bad/id", "mág"])
def test_input_rejects_invalid_annotation_identifiers(tmp_path, sample):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">contig\nACGT\n", encoding="utf-8")

    with pytest.raises(ValueError, match="blank required|invalid sample identifier"):
        load_input_record("isolate", isolate_sheet(tmp_path, fasta, sample=sample))


@pytest.mark.parametrize(
    ("content", "message"),
    [
        ("", "at least one header"),
        ("ACGT\n", "precedes its first header"),
        (">\nACGT\n", "blank header"),
        (">protein\nMPEPTIDE\n", "non-nucleotide"),
        (">contig\n", "at least one header"),
    ],
)
def test_isolate_rejects_malformed_or_non_nucleotide_fasta(tmp_path, content, message):
    fasta = tmp_path / "genome.fna"
    fasta.write_text(content, encoding="utf-8")

    with pytest.raises(ValueError, match=message):
        load_input_record("isolate", isolate_sheet(tmp_path, fasta))


def test_isolate_requires_existing_regular_file_with_supported_suffix(tmp_path):
    missing = tmp_path / "missing.fa"
    with pytest.raises(ValueError, match="does not exist"):
        load_input_record("isolate", isolate_sheet(tmp_path, missing))

    directory = tmp_path / "directory.fa"
    directory.mkdir()
    with pytest.raises(ValueError, match="not a regular file"):
        load_input_record("isolate", isolate_sheet(tmp_path, directory))

    compressed = tmp_path / "genome.fna.gz"
    compressed.write_text(">contig\nACGT\n", encoding="utf-8")
    with pytest.raises(ValueError, match="unsupported suffix"):
        load_input_record("isolate", isolate_sheet(tmp_path, compressed))


def metagenome_sheet(tmp_path, prok="", euk="", *, sample="project"):
    return write_sheet(
        tmp_path / "metagenome.tsv",
        ["sample", "prok_bins_dir", "euk_bins_dir"],
        [sample, str(prok), str(euk)],
    )


def test_metagenome_allows_one_blank_domain_and_empty_supplied_directory(tmp_path):
    prok = tmp_path / "empty-prok"
    prok.mkdir()

    assert load_input_record("metagenome", metagenome_sheet(tmp_path, prok=prok)) == {
        "sample": "project",
        "prok_bins_dir": str(prok),
        "euk_bins_dir": "",
    }


@pytest.mark.parametrize(
    ("header", "message"),
    [
        (["sample", "prok_bins_dir", "prok_bins_dir"], "duplicate column"),
        (["sample", "prok_bins_dir"], "missing column"),
        (
            ["sample", "prok_bins_dir", "euk_bins_dir", "extra"],
            "unknown column",
        ),
        (["sample", "euk_bins_dir", "prok_bins_dir"], "wrong order"),
    ],
)
def test_metagenome_requires_exact_header_and_order(tmp_path, header, message):
    bins = tmp_path / "bins"
    bins.mkdir()
    row = ["project", str(bins), "", "unused"][: len(header)]
    sheet = write_sheet(tmp_path / "bad.tsv", header, row)

    with pytest.raises(ValueError, match=message):
        load_input_record("metagenome", sheet)


def test_metagenome_rejects_both_domains_blank_or_invalid_path_types(tmp_path):
    with pytest.raises(ValueError, match="at least one nonblank"):
        load_input_record("metagenome", metagenome_sheet(tmp_path))

    missing = tmp_path / "missing"
    with pytest.raises(ValueError, match="does not exist"):
        load_input_record("metagenome", metagenome_sheet(tmp_path, prok=missing))

    fasta = tmp_path / "bins.fa"
    fasta.write_text(">contig\nACGT\n", encoding="utf-8")
    with pytest.raises(ValueError, match="not a directory"):
        load_input_record("metagenome", metagenome_sheet(tmp_path, prok=fasta))


def test_fasta_discovery_is_sorted_direct_child_only(tmp_path):
    bins = tmp_path / "bins"
    nested = bins / "nested"
    nested.mkdir(parents=True)
    for name in ("z.fna", "a.FASTA", "m.fa", "ignored.txt"):
        (bins / name).write_text(">contig\nACGT\n", encoding="utf-8")
    (nested / "nested.fa").write_text(">contig\nACGT\n", encoding="utf-8")

    assert [path.name for path in fasta_paths(bins)] == ["a.FASTA", "m.fa", "z.fna"]


def test_metagenome_rejects_duplicate_filename_ids_within_and_across_domains(tmp_path):
    prok = tmp_path / "prok"
    euk = tmp_path / "euk"
    prok.mkdir()
    euk.mkdir()
    (prok / "same.fa").write_text(">contig\nACGT\n", encoding="utf-8")
    (prok / "same.fna").write_text(">contig\nACGT\n", encoding="utf-8")

    with pytest.raises(ValueError, match="duplicate filename-derived"):
        load_input_record("metagenome", metagenome_sheet(tmp_path, prok, euk))

    (prok / "same.fna").unlink()
    (euk / "same.FASTA").write_text(">contig\nACGT\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate filename-derived"):
        load_input_record("metagenome", metagenome_sheet(tmp_path, prok, euk))


def test_metagenome_rejects_invalid_filename_derived_bin_id(tmp_path):
    bins = tmp_path / "bins"
    bins.mkdir()
    (bins / "bad id.fa").write_text(">contig\nACGT\n", encoding="utf-8")

    with pytest.raises(ValueError, match="invalid sample identifier"):
        load_input_record("metagenome", metagenome_sheet(tmp_path, prok=bins))


def test_sample_sheet_path_must_be_an_existing_regular_file(tmp_path):
    with pytest.raises(ValueError, match="does not exist"):
        load_input_record("isolate", tmp_path / "missing.tsv")

    with pytest.raises(ValueError, match="not a regular file"):
        load_input_record("isolate", tmp_path)


def test_invalid_input_mode_is_rejected_before_loading(tmp_path):
    with pytest.raises(ValueError, match="Invalid input.mode"):
        load_input_record("reads", Path(tmp_path / "unused.tsv"))
