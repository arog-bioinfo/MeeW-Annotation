import json

import pytest

from common.bakta_options import validate_bakta_extra
from common.euk_options import (
    params_json_context,
    validate_funannotate2_extras,
    validate_funannotate2_route,
    validate_metaeuk_extra,
)
from common.recognizer_options import validate_recognizer_extra
from common.upimapi_options import validate_upimapi_extra


@pytest.mark.parametrize(
    "extra",
    [
        "--keep-contig-headers",
        "--description 'text containing --prefix=value and -pvalue'",
        "--translation-table 11",
    ],
)
def test_bakta_accepts_unmanaged_options_without_substring_false_positives(extra):
    assert validate_bakta_extra(extra) == extra


@pytest.mark.parametrize(
    "extra",
    [
        "--prefix value",
        "--prefix=value",
        "--pref=value",
        "-pvalue",
        "-p=value",
        "--meta",
        "--output out",
        "-oout",
        "--force",
        "-f",
        "--skip-cd",
        "--skip-plo",
    ],
)
def test_bakta_rejects_managed_argparse_forms(extra):
    with pytest.raises(ValueError, match="managed by MeeW"):
        validate_bakta_extra(extra)


@pytest.mark.parametrize("extra", ["--min-length 20", "--min-intron 5", ""])
def test_metaeuk_accepts_unmanaged_options(extra):
    assert validate_metaeuk_extra(extra) == extra


@pytest.mark.parametrize(
    "extra",
    [
        "--threads 2",
        "--threads=2",
        "--max-exon-sets 4",
        "--max-exon-sets=4",
        "--exhaustive-search",
    ],
)
def test_metaeuk_rejects_managed_long_option_forms(extra):
    with pytest.raises(ValueError, match="managed by MeeW"):
        validate_metaeuk_extra(extra)


def test_funannotate2_accepts_unmanaged_phase_options():
    validate_funannotate2_extras(
        {
            "extra_clean": "--minlen 1000",
            "extra_train": "--optimize-augustus",
            "extra_predict": "--organism fungi",
            "extra_annotate": "--sbt template.sbt",
        }
    )


def test_funannotate2_params_validation_expands_supported_configured_paths(
    tmp_path, monkeypatch
):
    home = tmp_path / "home"
    environment_root = tmp_path / "environment"
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("MODEL_ROOT", str(environment_root))
    monkeypatch.chdir(tmp_path)
    configured_paths = {
        "~/models/tilde.json": home / "models" / "tilde.json",
        "$MODEL_ROOT/environment.json": environment_root / "environment.json",
        "models/relative.json": tmp_path / "models" / "relative.json",
        str(tmp_path / "normal.json"): tmp_path / "normal.json",
    }

    for configured, model in configured_paths.items():
        model.parent.mkdir(parents=True, exist_ok=True)
        model.write_text(
            json.dumps(
                {"species": "Fixture species", "taxonomy": {"kingdom": "Fungi"}}
            ),
            encoding="utf-8",
        )
        assert params_json_context(configured) == (
            "Fixture species",
            {"kingdom": "Fungi"},
        )
        validate_funannotate2_route(
            {"params": configured, "species": "Fixture species"},
            inspect_params=True,
        )

    with pytest.raises(ValueError, match="disagrees with params JSON species"):
        validate_funannotate2_route(
            {"params": "~/models/tilde.json", "species": "Other species"},
            inspect_params=True,
        )

    with pytest.raises(ValueError) as error:
        params_json_context("$MODEL_ROOT/missing.json")
    assert str(environment_root / "missing.json") in str(error.value)


@pytest.mark.parametrize(
    ("field", "extra"),
    [
        ("extra_clean", "-finput.fa"),
        ("extra_clean", "--out=clean.fa"),
        ("extra_train", "-sSpecies"),
        ("extra_train", "--busco-lineage=fungi"),
        ("extra_predict", "-pparams.json"),
        ("extra_predict", "-ststrain"),
        ("extra_predict", "--pretrained=model"),
        ("extra_annotate", "-iresults"),
        ("extra_annotate", "--cpus=4"),
    ],
)
def test_funannotate2_rejects_managed_phase_argparse_forms(field, extra):
    with pytest.raises(ValueError, match=rf"config.funannotate2.{field}.*managed"):
        validate_funannotate2_extras({field: extra})


@pytest.mark.parametrize("extra", ["--evalue 1e-10", "--quiet", "--keep-spaces"])
def test_recognizer_accepts_unmanaged_options(extra):
    assert validate_recognizer_extra(extra) == extra


@pytest.mark.parametrize(
    "extra",
    [
        "--threads 2",
        "--threads=2",
        "--thr=2",
        "-t2",
        "-t=2",
        "--file=input.faa",
        "-finput.faa",
        "--output out",
        "-oout",
        "--resources-directory=/db",
        "-rd=/db",
        "--databases COG",
        "-dbs=COG",
        "--custom-databases custom",
        "--keep-intermediates",
    ],
)
def test_recognizer_rejects_managed_argparse_forms(extra):
    with pytest.raises(ValueError, match="managed by MeeW"):
        validate_recognizer_extra(extra)


@pytest.mark.parametrize("extra", ["--evalue 1e-10", "--pident=80", "--mirror 1"])
def test_upimapi_accepts_unmanaged_options(extra):
    assert validate_upimapi_extra(extra) == extra


@pytest.mark.parametrize(
    "extra",
    [
        "--input input.faa",
        "--input=input.faa",
        "--inp=input.faa",
        "-iinput.faa",
        "--output out",
        "-o=out",
        "--output-table=result.tsv",
        "-ot=result.tsv",
        "--resources-directory /db",
        "-rd=/db",
        "--database=swissprot",
        "-db=swissprot",
        "--threads=2",
        "-t2",
        "--max-target-seqs=1",
        "-mts=1",
        "--skip-db-check",
        "--skip-id-checking",
        "--skip-id-mapping",
        "--no-annotation",
        "--from-db source",
        "--to-db target",
    ],
)
def test_upimapi_rejects_managed_argparse_forms(extra):
    with pytest.raises(ValueError, match="managed by MeeW"):
        validate_upimapi_extra(extra)


@pytest.mark.parametrize(
    ("validator", "config_key"),
    [
        (validate_bakta_extra, "config.bakta.extra"),
        (validate_metaeuk_extra, "config.metaeuk.extra"),
        (validate_recognizer_extra, "config.recognizer_prok.extra"),
        (validate_upimapi_extra, "config.upimapi.extra"),
    ],
)
def test_managed_option_helpers_report_malformed_shell_syntax(validator, config_key):
    with pytest.raises(ValueError, match=rf"{config_key} is not valid shell syntax"):
        validator("'unterminated")
