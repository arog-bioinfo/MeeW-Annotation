"""Exact decorated tool-version matching."""

import re


def version_contains(actual, expected):
    return bool(
        re.search(
            rf"(?<![0-9A-Za-z._+-])v?{re.escape(expected)}" r"(?![0-9A-Za-z._+-])",
            str(actual),
        )
    )
