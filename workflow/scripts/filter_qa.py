"""Helpers for filtering samples by external QA reports."""

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
DEFAULT_MISSING_SAMPLE = "error"
FASTA_SUFFIXES = (".fasta", ".fa", ".fna")
COMPRESSION_SUFFIXES = (".gz", ".bz2", ".xz")


def _as_list(paths):
    if paths is None:
        return []
    if isinstance(paths, (str, Path)):
        return [paths]
    return list(paths)


def _read_report(path, sep):
    if sep is None:
        return pd.read_csv(path, sep=None, engine="python")
    return pd.read_csv(path, sep=sep)


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


def _sample_ids(row):
    ids = []
    if "sample" in row.index:
        ids.append(str(row["sample"]))
    ids.append(str(row.name))
    if "path" in row.index:
        ids.append(fasta_stem(row["path"]))
    return list(dict.fromkeys(ids))


def _qa_for_sample(row, qa):
    for sample_id in _sample_ids(row):
        if sample_id in qa.index:
            return qa.loc[sample_id]
    return None


def filter_samples_by_domain(
    samples,
    qa_config,
    min_completeness=None,
    max_contamination=None,
    missing_sample=None,
):
    """Filter a samples DataFrame using CheckM2 for prok and EukCC for euk."""
    qa_config = qa_config or {}
    if not qa_config.get("enabled", False):
        return samples

    min_completeness = qa_config.get(
        "min_completeness", min_completeness or DEFAULT_MIN_COMPLETENESS
    )
    max_contamination = qa_config.get(
        "max_contamination", max_contamination or DEFAULT_MAX_CONTAMINATION
    )
    missing_sample = qa_config.get(
        "missing_sample", missing_sample or DEFAULT_MISSING_SAMPLE
    )
    if missing_sample not in {"error", "keep", "drop"}:
        raise ValueError("qa_filter.missing_sample must be one of: error, keep, drop")

    qa_by_domain = {
        "prok": parse_checkm2_reports(qa_config.get("checkm2_reports", [])),
        "euk": parse_eukcc_reports(qa_config.get("eukcc_reports", [])),
    }

    keep = []
    missing = []
    for _, row in samples.iterrows():
        domain = row.get("domain")
        if domain not in qa_by_domain:
            keep.append(True)
            continue

        qa = _qa_for_sample(row, qa_by_domain[domain])
        if qa is None:
            missing.append(str(row.get("sample", row.name)))
            keep.append(missing_sample == "keep")
            continue

        keep.append(
            qa["completeness"] >= float(min_completeness)
            and qa["contamination"] <= float(max_contamination)
        )

    if missing and missing_sample == "error":
        raise ValueError(
            "QA report entries are missing for samples: " + ", ".join(sorted(missing))
        )

    return samples.loc[keep].copy()
