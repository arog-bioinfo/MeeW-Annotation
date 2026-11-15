"""Parse native QA reports and build the Annotation selection contract."""

import csv
import math
from pathlib import Path

import pandas as pd

CHECKM2_COLUMNS = {
    "Name": "qa_sample",
    "Completeness": "completeness",
    "Contamination": "contamination",
}
EUKCC_COLUMNS = {
    "bin": "qa_sample",
    "completeness": "completeness",
    "contamination": "contamination",
}
DEFAULT_MIN_COMPLETENESS = 50.0
DEFAULT_MAX_CONTAMINATION = 10.0
FASTA_SUFFIXES = (".fasta", ".fa", ".fna")
COMPRESSION_SUFFIXES = (".gz", ".bz2", ".xz")
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


def _as_list(paths):
    if paths is None:
        return []
    if isinstance(paths, (str, Path)):
        return [paths]
    return list(paths)


def _read_report(path, sep):
    if sep is None:
        return pd.read_csv(path, sep=None, engine="python", dtype=str)
    return pd.read_csv(path, sep=sep, dtype=str)


def _read_reports(paths, sep, columns, report_name):
    reports = []
    for path in _as_list(paths):
        report = None
        if isinstance(sep, (list, tuple)):
            for candidate_sep in sep:
                candidate_report = _read_report(path, candidate_sep)
                if all(column in candidate_report.columns for column in columns):
                    report = candidate_report
                    break
                report = candidate_report
        else:
            report = _read_report(path, sep)
        missing = [column for column in columns if column not in report.columns]
        if missing:
            raise ValueError(
                f"{report_name} report '{path}' is missing required columns: "
                f"{', '.join(missing)}"
            )
        reports.append(report[list(columns)].rename(columns=columns))

    if not reports:
        return pd.DataFrame(columns=columns.values())

    qa = pd.concat(reports, ignore_index=True)
    qa["qa_sample"] = qa["qa_sample"].astype(str)
    qa["completeness"] = pd.to_numeric(qa["completeness"], errors="coerce")
    qa["contamination"] = pd.to_numeric(qa["contamination"], errors="coerce")
    duplicates = qa.loc[qa["qa_sample"].duplicated(), "qa_sample"].unique()
    if len(duplicates):
        raise ValueError(
            f"{report_name} reports contain duplicate sample ids: "
            f"{', '.join(sorted(duplicates))}"
        )
    return qa.set_index("qa_sample")


def parse_checkm2_reports(paths):
    """Parse one or more CheckM2 TSV reports."""
    return _read_reports(paths, "\t", CHECKM2_COLUMNS, "CheckM2")


def parse_eukcc_reports(paths):
    """Parse one or more EukCC CSV/TSV reports."""
    return _read_reports(paths, ("\t", ",", None), EUKCC_COLUMNS, "EukCC")


def parse_eukcc_missing_markers(paths):
    """Parse EukCC's native missing_marker_genes.txt outcome channel."""
    missing = []
    for path in _as_list(paths):
        with Path(path).open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t", strict=True)
            if reader.fieldnames != ["bin"]:
                raise ValueError(
                    f"EukCC missing-marker report '{path}' has an invalid header; "
                    "expected: bin"
                )
            for row in reader:
                value = row["bin"].strip()
                if not value:
                    raise ValueError(
                        f"EukCC missing-marker report '{path}' contains a blank bin id."
                    )
                missing.append(value)
    duplicates = sorted({value for value in missing if missing.count(value) > 1})
    if duplicates:
        raise ValueError(
            "EukCC missing-marker reports contain duplicate sample ids: "
            + ", ".join(duplicates)
        )
    return set(missing)


def fasta_stem(path):
    """Return a FASTA filename stem, ignoring common compression suffixes."""
    name = Path(path).name
    for suffix in COMPRESSION_SUFFIXES:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    for suffix in FASTA_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def canonical_qa_identifier(value):
    """Normalize only recognized FASTA suffixes from a native QA identifier."""
    name = Path(str(value)).name
    for compression in COMPRESSION_SUFFIXES:
        if name.endswith(compression):
            name = name[: -len(compression)]
            break
    for suffix in FASTA_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def qa_row_for_source(source, sample, qa):
    """Find exactly one report row using stable sample/file aliases."""
    aliases = list(
        dict.fromkeys(
            [str(sample), Path(source).name, fasta_stem(source), Path(source).stem]
        )
    )
    canonical_aliases = {canonical_qa_identifier(alias) for alias in aliases}
    matches = [
        qa_id
        for qa_id in qa.index
        if qa_id in aliases or canonical_qa_identifier(qa_id) in canonical_aliases
    ]
    if not matches:
        return None
    if len(set(matches)) > 1:
        rows = qa.loc[matches][["completeness", "contamination"]].drop_duplicates()
        if len(rows) > 1:
            raise ValueError(
                f"QA report has conflicting aliases for sample '{sample}': "
                + ", ".join(matches)
            )
    return qa.loc[matches[0]]


def source_has_missing_markers(source, sample, missing_marker_samples):
    aliases = {
        canonical_qa_identifier(value)
        for value in (sample, Path(source).name, fasta_stem(source), Path(source).stem)
    }
    return any(
        canonical_qa_identifier(value) in aliases for value in missing_marker_samples
    )


def _metric(value):
    try:
        metric = float(value)
    except (TypeError, ValueError):
        return None
    return metric if math.isfinite(metric) else None


def _format_metric(value):
    return "" if value is None else format(value, ".15g")


def selection_row(
    *,
    source,
    sample,
    domain,
    is_isolate,
    qa,
    qa_filter,
    missing_marker_samples=None,
):
    """Create one deterministic selection row from CheckM2 or EukCC metrics."""
    qa_filter = qa_filter or {}
    minimum = float(qa_filter.get("min_completeness", DEFAULT_MIN_COMPLETENESS))
    maximum = float(qa_filter.get("max_contamination", DEFAULT_MAX_CONTAMINATION))
    missing_marker_samples = missing_marker_samples or set()

    report_row = qa_row_for_source(source, sample, qa)
    completeness = (
        _metric(report_row["completeness"]) if report_row is not None else None
    )
    contamination = (
        _metric(report_row["contamination"]) if report_row is not None else None
    )
    problem = None
    if report_row is None:
        problem = "missing_qa_entry"
    elif completeness is None or contamination is None:
        problem = "invalid_qa_metrics"

    missing_markers = source_has_missing_markers(source, sample, missing_marker_samples)

    if missing_markers:
        decision = "report_only" if is_isolate else "fail"
        reason = "missing_marker_genes"
    elif is_isolate:
        decision = "report_only"
        reason = (
            "isolate_report_only" if not problem else f"isolate_report_only_{problem}"
        )
    elif problem:
        tool = "CheckM2" if domain == "prok" else "EukCC"
        raise ValueError(
            f"{tool} {problem.replace('_', ' ')} for sample '{sample}' "
            f"({source}). The generated primary QA output must contain finite "
            "completeness and contamination metrics for every MAG."
        )
    elif not qa_filter.get("enabled", True):
        decision = "pass"
        reason = "filtering_disabled"
    else:
        low_completeness = completeness < minimum
        high_contamination = contamination > maximum
        if not low_completeness and not high_contamination:
            decision = "pass"
            reason = "meets_thresholds"
        elif low_completeness and high_contamination:
            decision = "fail"
            reason = "below_min_completeness_and_above_max_contamination"
        elif low_completeness:
            decision = "fail"
            reason = "below_min_completeness"
        else:
            decision = "fail"
            reason = "above_max_contamination"

    return {
        "sample": str(sample),
        "path": str(Path(source).resolve()),
        "domain": domain,
        "completeness": _format_metric(completeness),
        "contamination": _format_metric(contamination),
        "min_completeness": _format_metric(minimum),
        "max_contamination": _format_metric(maximum),
        "decision": decision,
        "reason": reason,
    }


def qa_bypass_row(*, source, sample, domain, is_isolate):
    """Create one deterministic selection row without QA metrics or thresholds."""
    return {
        "sample": str(sample),
        "path": str(Path(source).resolve()),
        "domain": domain,
        "completeness": "",
        "contamination": "",
        "min_completeness": "",
        "max_contamination": "",
        "decision": "report_only" if is_isolate else "pass",
        "reason": "qa_disabled",
    }


def write_selection(rows, path):
    """Write rows with a fixed header and deterministic ordering."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(rows, key=lambda row: (row["sample"], row["domain"], row["path"]))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=SELECTION_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(ordered)


def read_selection(path):
    """Read and strictly validate an Annotation selection TSV."""
    path = Path(path)
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        if reader.fieldnames != SELECTION_COLUMNS:
            raise ValueError(
                f"Selection file '{path}' has an invalid header; expected: "
                + "\t".join(SELECTION_COLUMNS)
            )
        rows = list(reader)
    for row in rows:
        if row["decision"] not in {"pass", "fail", "report_only"}:
            raise ValueError(
                f"Selection file '{path}' has invalid decision '{row['decision']}'."
            )
        if row["domain"] not in {"prok", "euk"}:
            raise ValueError(f"Selection file '{path}' has invalid domain.")
    return rows
