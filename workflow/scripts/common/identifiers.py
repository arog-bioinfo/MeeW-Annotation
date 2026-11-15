"""Identifier contracts accepted by MeeW-Annotation."""

import re

USER_IDENTIFIER_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}")
GENERATED_MAG_PATTERN = re.compile(
    r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}_(?:prok|euk)_bin_[0-9]{3,}"
)
GENERATED_MAG_MAX_LENGTH = 192


def validate_annotation_identifier(value: str) -> str:
    """Accept an ordinary user ID or Binning-generated MAG ID unchanged."""
    if not isinstance(value, str) or not value.isascii():
        raise ValueError(f"Invalid Annotation identifier: {value!r}")
    if USER_IDENTIFIER_PATTERN.fullmatch(value) is not None:
        return value
    if (
        len(value) <= GENERATED_MAG_MAX_LENGTH
        and GENERATED_MAG_PATTERN.fullmatch(value) is not None
    ):
        return value
    raise ValueError(f"Invalid Annotation identifier: {value!r}")
