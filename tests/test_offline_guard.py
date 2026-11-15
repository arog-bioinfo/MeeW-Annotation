import os
import subprocess
import sys
import urllib.request
from pathlib import Path
from urllib.request import urlretrieve as imported_urlretrieve

import pytest

BLOCKED_URL = "https://example.com/meew-test-database.tar.gz"


def test_offline_guard_blocks_parent_alias_direct_urllib_and_child_downloaders(
    offline_guard, tmp_path
):
    parent_direct = tmp_path / "parent-direct.tar.gz"
    with pytest.raises(RuntimeError, match="offline guard blocked"):
        urllib.request.urlretrieve(BLOCKED_URL, parent_direct)
    assert not parent_direct.exists()

    parent_alias = tmp_path / "parent-alias.tar.gz"
    with pytest.raises(RuntimeError, match="offline guard blocked"):
        imported_urlretrieve(BLOCKED_URL, parent_alias)
    assert not parent_alias.exists()

    child_urllib = tmp_path / "child-urllib.tar.gz"
    child = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys, urllib.request; "
            "urllib.request.urlretrieve(sys.argv[1], sys.argv[2])",
            BLOCKED_URL,
            str(child_urllib),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=os.environ.copy(),
    )
    assert child.returncode != 0
    assert "offline guard blocked" in child.stderr
    assert not child_urllib.exists()

    wget_payload = tmp_path / "wget.tar.gz"
    wget = subprocess.run(
        ["wget", BLOCKED_URL, "-O", str(wget_payload)],
        check=False,
        capture_output=True,
        text=True,
        env=os.environ.copy(),
    )
    assert wget.returncode == 97
    assert "offline guard blocked executable: wget" in wget.stderr
    assert not wget_payload.exists()

    attempts = offline_guard["attempt_log"].read_text(encoding="utf-8").splitlines()
    assert attempts == [f"wget {BLOCKED_URL} -O {wget_payload}"]
