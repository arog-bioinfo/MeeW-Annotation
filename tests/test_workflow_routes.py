import pytest

from helpers import rule_count

SUPPORTED_QA_MODES = [
    pytest.param(True, True, id="qa-filter"),
    pytest.param(True, False, id="qa-report-no-filter"),
    pytest.param(False, False, id="qa-bypass"),
]


@pytest.mark.parametrize(("qa_enabled", "filter_enabled"), SUPPORTED_QA_MODES)
@pytest.mark.parametrize(
    ("domain", "expected_rules", "excluded_rules"),
    [
        (
            "prok",
            ["bakta", "recognizer_prok", "upimapi", "protein_summary"],
            ["funannotate2", "metaeuk_bin"],
        ),
        (
            "euk",
            ["funannotate2"],
            ["bakta", "recognizer_prok", "upimapi", "protein_summary"],
        ),
    ],
)
def test_isolate_dag_routes_report_only_selection_across_supported_qa_modes(
    workflow_case,
    qa_enabled,
    filter_enabled,
    domain,
    expected_rules,
    excluded_rules,
):
    case = workflow_case(
        "isolate",
        enabled=False,
        isolate_domain=domain,
        qa_enabled=qa_enabled,
        filter_enabled=filter_enabled,
    )

    dry_run = case.run(dry_run=True)

    for rule in expected_rules:
        assert rule_count(dry_run, rule) == 1
    for rule in excluded_rules:
        assert rule_count(dry_run, rule) == 0


@pytest.mark.parametrize(("qa_enabled", "filter_enabled"), SUPPORTED_QA_MODES)
def test_mixed_metagenome_dag_annotates_only_passing_bins(
    workflow_case, qa_enabled, filter_enabled
):
    low_quality_decision = "fail" if filter_enabled else "pass"
    case = workflow_case(
        "metagenome",
        enabled=False,
        qa_enabled=qa_enabled,
        filter_enabled=filter_enabled,
        rows=[
            ("pass_prok", "prok", "pass"),
            ("low_quality_prok", "prok", low_quality_decision),
            ("pass_euk", "euk", "pass"),
            ("low_quality_euk", "euk", low_quality_decision),
        ],
    )

    dry_run = case.run(dry_run=True)

    expected_per_domain = 1 if filter_enabled else 2
    assert rule_count(dry_run, "bakta_bin") == expected_per_domain
    assert rule_count(dry_run, "recognizer_prok_bin") == expected_per_domain
    assert rule_count(dry_run, "upimapi_bin") == expected_per_domain
    assert rule_count(dry_run, "protein_summary_bin") == expected_per_domain
    assert rule_count(dry_run, "metaeuk_bin") == expected_per_domain
    assert rule_count(dry_run, "recognizer_euk_bin") == expected_per_domain
    assert "pass_prok" in dry_run
    assert "pass_euk" in dry_run
    assert ("batch/prok/low_quality_prok" in dry_run) is not filter_enabled
    assert ("batch/euk/low_quality_euk" in dry_run) is not filter_enabled


@pytest.mark.parametrize("domain", ["prok", "euk"])
def test_metagenome_absent_domain_has_no_annotation_jobs(workflow_case, domain):
    case = workflow_case(
        "metagenome",
        enabled=False,
        rows=[(f"pass_{domain}", domain, "pass")],
        supplied_domains={domain},
    )

    dry_run = case.run(dry_run=True)

    if domain == "prok":
        assert rule_count(dry_run, "bakta_bin") == 1
        assert rule_count(dry_run, "metaeuk_bin") == 0
        assert rule_count(dry_run, "recognizer_euk_bin") == 0
    else:
        assert rule_count(dry_run, "metaeuk_bin") == 1
        assert rule_count(dry_run, "bakta_bin") == 0
        assert rule_count(dry_run, "recognizer_prok_bin") == 0


def test_empty_supplied_domains_complete_without_annotation_jobs(workflow_case):
    case = workflow_case(
        "metagenome",
        enabled=False,
        rows=[],
        supplied_domains={"prok", "euk"},
    )

    dry_run = case.run(dry_run=True)

    for rule in (
        "bakta_bin",
        "recognizer_prok_bin",
        "upimapi_bin",
        "protein_summary_bin",
        "metaeuk_bin",
        "recognizer_euk_bin",
    ):
        assert rule_count(dry_run, rule) == 0
    assert rule_count(dry_run, "metagenome_annotation_inventory") == 1


def test_all_failed_metagenome_bins_create_no_annotation_jobs(workflow_case):
    case = workflow_case(
        "metagenome",
        enabled=False,
        qa_enabled=True,
        filter_enabled=True,
        rows=[
            ("failed_prok", "prok", "fail"),
            ("failed_euk", "euk", "fail"),
        ],
    )

    dry_run = case.run(dry_run=True)

    for rule in (
        "bakta_bin",
        "recognizer_prok_bin",
        "upimapi_bin",
        "protein_summary_bin",
        "metaeuk_bin",
        "recognizer_euk_bin",
    ):
        assert rule_count(dry_run, rule) == 0
    assert "batch/prok/failed_prok" not in dry_run
    assert "batch/euk/failed_euk" not in dry_run


def test_unsupported_qa_off_filter_on_combination_is_rejected(workflow_case):
    case = workflow_case(
        "metagenome",
        enabled=False,
        qa_enabled=False,
        filter_enabled=True,
        rows=[("pass_prok", "prok", "pass")],
    )

    with pytest.raises(AssertionError, match="Invalid QA configuration"):
        case.run(dry_run=True)
