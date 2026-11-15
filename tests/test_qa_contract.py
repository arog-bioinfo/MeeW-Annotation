import csv
import runpy
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from common import qa as qa_module
from common.qa import (
    SELECTION_COLUMNS,
    parse_checkm2_reports,
    parse_eukcc_missing_markers,
    parse_eukcc_reports,
    qa_bypass_row,
    read_selection,
    selection_row,
    write_selection,
)

ROOT = Path(__file__).resolve().parents[1]
SELECTION_WRITER = ROOT / "workflow" / "scripts" / "write_qa_selection.py"

FILTER = {
    "enabled": True,
    "min_completeness": 50,
    "max_contamination": 10,
}


def report(path, header, rows, delimiter="\t"):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter=delimiter, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return path


def qa_frame(sample="sample", completeness=90, contamination=1):
    return pd.DataFrame(
        {
            "completeness": [completeness],
            "contamination": [contamination],
        },
        index=[sample],
    )


def decide(tmp_path, completeness, contamination, **kwargs):
    return selection_row(
        source=tmp_path / "sample.fna",
        sample="sample",
        domain=kwargs.pop("domain", "prok"),
        is_isolate=kwargs.pop("is_isolate", False),
        qa=qa_frame("sample", completeness, contamination),
        qa_filter=kwargs.pop("qa_filter", FILTER),
        **kwargs,
    )


def test_checkm2_parser_combines_reports_and_converts_metrics(tmp_path):
    first = report(
        tmp_path / "first.tsv",
        ["Name", "Completeness", "Contamination", "Other"],
        [["b", "91.5", "1", "ignored"]],
    )
    second = report(
        tmp_path / "second.tsv",
        ["Name", "Completeness", "Contamination"],
        [["a", "80", "2.5"]],
    )

    parsed = parse_checkm2_reports([first, second])

    assert list(parsed.index) == ["b", "a"]
    assert parsed.loc["b"].to_dict() == {
        "completeness": 91.5,
        "contamination": 1.0,
    }


@pytest.mark.parametrize("delimiter", ["\t", ","])
def test_eukcc_parser_accepts_native_tsv_and_csv(tmp_path, delimiter):
    path = report(
        tmp_path / "eukcc.csv",
        ["bin", "completeness", "contamination", "ncbi_lng"],
        [["euk.fa", "75", "3", "unknown"]],
        delimiter,
    )

    parsed = parse_eukcc_reports(path)

    assert parsed.loc["euk.fa"].to_dict() == {
        "completeness": 75.0,
        "contamination": 3.0,
    }


@pytest.mark.parametrize("parser", [parse_checkm2_reports, parse_eukcc_reports])
def test_qa_parsers_reject_missing_required_columns(tmp_path, parser):
    path = report(tmp_path / "bad.tsv", ["sample", "score"], [["a", "90"]])

    with pytest.raises(ValueError, match="missing required columns"):
        parser(path)


def test_qa_parser_rejects_duplicate_sample_ids_across_reports(tmp_path):
    paths = [
        report(
            tmp_path / f"report-{index}.tsv",
            ["Name", "Completeness", "Contamination"],
            [["duplicate", "90", "1"]],
        )
        for index in range(2)
    ]

    with pytest.raises(ValueError, match="duplicate sample ids: duplicate"):
        parse_checkm2_reports(paths)


def test_qa_aliases_match_fasta_names_and_reject_conflicting_aliases(tmp_path):
    source = tmp_path / "sample.fna"
    one_alias = qa_frame("sample.fna", 90, 1)
    assert (
        selection_row(
            source=source,
            sample="sample",
            domain="prok",
            is_isolate=False,
            qa=one_alias,
            qa_filter=FILTER,
        )["decision"]
        == "pass"
    )

    conflicting = pd.DataFrame(
        {
            "completeness": [90, 20],
            "contamination": [1, 30],
        },
        index=["sample", "sample.fna"],
    )
    with pytest.raises(ValueError, match="conflicting aliases"):
        selection_row(
            source=source,
            sample="sample",
            domain="prok",
            is_isolate=False,
            qa=conflicting,
            qa_filter=FILTER,
        )


@pytest.mark.parametrize(
    ("completeness", "contamination", "decision", "reason"),
    [
        (50, 10, "pass", "meets_thresholds"),
        (49.9, 10, "fail", "below_min_completeness"),
        (50, 10.1, "fail", "above_max_contamination"),
        (
            49.9,
            10.1,
            "fail",
            "below_min_completeness_and_above_max_contamination",
        ),
    ],
)
def test_current_threshold_boundaries_are_minimum_and_inclusive_maximum(
    tmp_path, completeness, contamination, decision, reason
):
    row = decide(tmp_path, completeness, contamination)

    assert (row["decision"], row["reason"]) == (decision, reason)


def test_disabled_filter_keeps_low_quality_mag_with_metrics(tmp_path):
    row = decide(
        tmp_path,
        1,
        99,
        qa_filter={
            "enabled": False,
            "min_completeness": 50,
            "max_contamination": 10,
        },
    )

    assert (row["decision"], row["reason"]) == ("pass", "filtering_disabled")
    assert (row["completeness"], row["contamination"]) == ("1", "99")


@pytest.mark.parametrize("value", ["NA", "bad", "nan", "inf", "-inf", None])
def test_mag_requires_finite_completeness_and_contamination(tmp_path, value):
    with pytest.raises(ValueError, match="finite completeness and contamination"):
        decide(tmp_path, value, 1)

    with pytest.raises(ValueError, match="finite completeness and contamination"):
        decide(tmp_path, 90, value)


def test_missing_mag_metric_row_is_actionable(tmp_path):
    empty = pd.DataFrame(columns=["completeness", "contamination"])

    with pytest.raises(ValueError, match="missing qa entry.*sample 'sample'"):
        selection_row(
            source=tmp_path / "sample.fna",
            sample="sample",
            domain="prok",
            is_isolate=False,
            qa=empty,
            qa_filter=FILTER,
        )


def test_isolate_routing_is_report_only_for_pass_fail_missing_and_invalid(tmp_path):
    complete = decide(tmp_path, 90, 1, is_isolate=True)
    failed = decide(tmp_path, 1, 99, is_isolate=True)
    invalid = decide(tmp_path, "NA", 1, is_isolate=True)
    missing = selection_row(
        source=tmp_path / "sample.fna",
        sample="sample",
        domain="prok",
        is_isolate=True,
        qa=pd.DataFrame(columns=["completeness", "contamination"]),
        qa_filter=FILTER,
    )

    assert [row["decision"] for row in (complete, failed, invalid, missing)] == [
        "report_only",
        "report_only",
        "report_only",
        "report_only",
    ]
    assert invalid["reason"] == "isolate_report_only_invalid_qa_metrics"
    assert missing["reason"] == "isolate_report_only_missing_qa_entry"


def test_eukcc_missing_marker_parser_and_routing_use_fasta_aliases(tmp_path):
    markers = report(tmp_path / "missing.txt", ["bin"], [["markerless.fa"]])
    missing = parse_eukcc_missing_markers(markers)

    mag = selection_row(
        source=tmp_path / "markerless.fna",
        sample="markerless",
        domain="euk",
        is_isolate=False,
        qa=pd.DataFrame(columns=["completeness", "contamination"]),
        qa_filter=FILTER,
        missing_marker_samples=missing,
    )
    isolate = selection_row(
        source=tmp_path / "markerless.fna",
        sample="markerless",
        domain="euk",
        is_isolate=True,
        qa=pd.DataFrame(columns=["completeness", "contamination"]),
        qa_filter=FILTER,
        missing_marker_samples=missing,
    )

    assert (mag["decision"], mag["reason"]) == ("fail", "missing_marker_genes")
    assert (isolate["decision"], isolate["reason"]) == (
        "report_only",
        "missing_marker_genes",
    )


@pytest.mark.parametrize(
    ("content", "message"),
    [
        ("sample\nvalue\n", "invalid header"),
        ("bin\n \n", "blank bin id"),
        ("bin\nduplicate\nduplicate\n", "duplicate sample ids"),
    ],
)
def test_eukcc_missing_marker_reports_are_strict(tmp_path, content, message):
    path = tmp_path / "missing.txt"
    path.write_text(content, encoding="utf-8")

    with pytest.raises(ValueError, match=message):
        parse_eukcc_missing_markers(path)


def test_qa_bypass_rows_have_blank_metrics_and_mode_specific_decisions(tmp_path):
    isolate = qa_bypass_row(
        source=tmp_path / "isolate.fa",
        sample="isolate",
        domain="prok",
        is_isolate=True,
    )
    mag = qa_bypass_row(
        source=tmp_path / "mag.fa",
        sample="mag",
        domain="euk",
        is_isolate=False,
    )

    assert (isolate["decision"], mag["decision"]) == ("report_only", "pass")
    for row in (isolate, mag):
        assert row["reason"] == "qa_disabled"
        assert all(
            row[column] == ""
            for column in (
                "completeness",
                "contamination",
                "min_completeness",
                "max_contamination",
            )
        )


def test_selection_round_trip_has_exact_header_and_deterministic_order(tmp_path):
    rows = [
        qa_bypass_row(
            source=tmp_path / "z.fa",
            sample="z",
            domain="prok",
            is_isolate=False,
        ),
        qa_bypass_row(
            source=tmp_path / "a.fa",
            sample="a",
            domain="euk",
            is_isolate=False,
        ),
    ]
    output = tmp_path / "selection.tsv"

    write_selection(rows, output)

    assert output.read_text(encoding="utf-8").splitlines()[0] == "\t".join(
        SELECTION_COLUMNS
    )
    assert [row["sample"] for row in read_selection(output)] == ["a", "z"]


@pytest.mark.parametrize(
    ("column", "value", "message"),
    [
        ("decision", "maybe", "invalid decision"),
        ("domain", "viral", "invalid domain"),
    ],
)
def test_selection_reader_rejects_invalid_routing_values(
    tmp_path, column, value, message
):
    row = qa_bypass_row(
        source=tmp_path / "sample.fa",
        sample="sample",
        domain="prok",
        is_isolate=False,
    )
    row[column] = value
    path = tmp_path / "selection.tsv"
    write_selection([row], path)

    with pytest.raises(ValueError, match=message):
        read_selection(path)


def test_selection_reader_rejects_noncontract_header(tmp_path):
    path = tmp_path / "selection.tsv"
    path.write_text("sample\tdomain\n", encoding="utf-8")

    with pytest.raises(ValueError, match="invalid header"):
        read_selection(path)


def selection_writer_fixture(tmp_path, action):
    prok = tmp_path / "inputs" / "prok"
    euk = tmp_path / "inputs" / "euk"
    prok.mkdir(parents=True)
    euk.mkdir(parents=True)
    for path in (prok / "z.fasta", prok / "a.fa", euk / "markerless.fna"):
        path.write_text(">contig\nACGT\n", encoding="utf-8")
    nested = prok / "nested"
    nested.mkdir()
    (nested / "ignored.fa").write_text(">contig\nACGT\n", encoding="utf-8")

    checkm2 = report(
        tmp_path / "checkm2.tsv",
        ["Name", "Completeness", "Contamination"],
        [["a", "90", "1"], ["z.fasta", "40", "11"], ["isolate", "80", "2"]],
    )
    eukcc = report(
        tmp_path / "eukcc.csv",
        ["bin", "completeness", "contamination", "ncbi_lng"],
        [["markerless.fna", "99", "0", "unknown"]],
        delimiter=",",
    )
    eukcc_bad = report(
        tmp_path / "eukcc-bad.csv",
        ["bin", "completeness", "contamination", "ncbi_lng"],
        [],
        delimiter=",",
    )
    missing = report(
        tmp_path / "missing_marker_genes.txt", ["bin"], [["markerless.fna"]]
    )
    isolate = tmp_path / "inputs" / "isolate.fa"
    isolate.write_text(">contig\nACGT\n", encoding="utf-8")
    output = tmp_path / "selection.tsv"
    log = tmp_path / "selection.log"
    snakemake = SimpleNamespace(
        input=SimpleNamespace(
            prok_bins=[str(prok)],
            euk_bins=[str(euk)],
            checkm2=[str(checkm2)],
            eukcc=[str(eukcc), str(eukcc_bad)],
            eukcc_missing=[str(missing)],
        ),
        output=SimpleNamespace(selection=str(output)),
        params=SimpleNamespace(
            action=action,
            path=str(isolate),
            sample="isolate",
            domain="prok",
            qa_filter=FILTER,
        ),
        log=[str(log)],
    )
    return snakemake, output, log, prok, euk, isolate


@pytest.mark.parametrize(
    ("action", "expected"),
    [
        (
            "isolate_bypass",
            [("isolate", "prok", "", "", "report_only", "qa_disabled")],
        ),
        (
            "isolate",
            [("isolate", "prok", "80", "2", "report_only", "isolate_report_only")],
        ),
        (
            "metagenome_bypass",
            [
                ("a", "prok", "", "", "pass", "qa_disabled"),
                ("markerless", "euk", "", "", "pass", "qa_disabled"),
                ("z", "prok", "", "", "pass", "qa_disabled"),
            ],
        ),
        (
            "metagenome",
            [
                ("a", "prok", "90", "1", "pass", "meets_thresholds"),
                ("markerless", "euk", "99", "0", "fail", "missing_marker_genes"),
                (
                    "z",
                    "prok",
                    "40",
                    "11",
                    "fail",
                    "below_min_completeness_and_above_max_contamination",
                ),
            ],
        ),
    ],
)
def test_selection_writer_preserves_all_actions_and_parses_each_report_once(
    tmp_path, monkeypatch, action, expected
):
    snakemake, output, log, _, _, _ = selection_writer_fixture(tmp_path, action)
    calls = {"checkm2": 0, "eukcc": 0, "missing": 0}

    def counted(name, parser):
        def wrapped(paths):
            calls[name] += 1
            return parser(paths)

        return wrapped

    monkeypatch.setattr(
        qa_module,
        "parse_checkm2_reports",
        counted("checkm2", qa_module.parse_checkm2_reports),
    )
    monkeypatch.setattr(
        qa_module,
        "parse_eukcc_reports",
        counted("eukcc", qa_module.parse_eukcc_reports),
    )
    monkeypatch.setattr(
        qa_module,
        "parse_eukcc_missing_markers",
        counted("missing", qa_module.parse_eukcc_missing_markers),
    )

    runpy.run_path(str(SELECTION_WRITER), init_globals={"snakemake": snakemake})

    rows = read_selection(output)
    assert [
        (
            row["sample"],
            row["domain"],
            row["completeness"],
            row["contamination"],
            row["decision"],
            row["reason"],
        )
        for row in rows
    ] == expected
    expected_calls = {
        "isolate_bypass": {"checkm2": 0, "eukcc": 0, "missing": 0},
        "isolate": {"checkm2": 1, "eukcc": 0, "missing": 0},
        "metagenome_bypass": {"checkm2": 0, "eukcc": 0, "missing": 0},
        "metagenome": {"checkm2": 1, "eukcc": 1, "missing": 1},
    }
    assert calls == expected_calls[action]
    assert log.read_text(encoding="utf-8") == (
        f"selection_rows={len(expected)}\naction={action}\n"
    )


@pytest.mark.parametrize("action", ["metagenome_bypass", "metagenome"])
def test_selection_writer_rejects_duplicate_direct_child_identifiers(tmp_path, action):
    snakemake, _, _, prok, _, _ = selection_writer_fixture(tmp_path, action)
    (prok / "a.fna").write_text(">contig\nACGT\n", encoding="utf-8")

    with pytest.raises(ValueError, match="sample identifier is duplicated"):
        runpy.run_path(str(SELECTION_WRITER), init_globals={"snakemake": snakemake})


def test_selection_writer_rejects_unknown_action_before_dispatch(tmp_path):
    snakemake = SimpleNamespace(
        input=SimpleNamespace(),
        output=SimpleNamespace(selection=str(tmp_path / "selection.tsv")),
        params=SimpleNamespace(action="unknown"),
        log=[str(tmp_path / "selection.log")],
    )

    with pytest.raises(ValueError, match="Unsupported QA selection action: unknown"):
        runpy.run_path(str(SELECTION_WRITER), init_globals={"snakemake": snakemake})
