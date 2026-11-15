"""Helpers for validating workflow-managed argparse options."""

import shlex


def split_extra(extra, config_key):
    """Split a configured CLI fragment or raise a configuration-focused error."""
    extra = "" if extra is None else str(extra)
    try:
        return extra, shlex.split(extra, posix=True, comments=False)
    except ValueError as error:
        raise ValueError(f"{config_key} is not valid shell syntax: {error}") from error


def managed_argparse_option(
    token, *, long_options, managed_long_options, short_options
):
    """Resolve an argparse-accepted managed option represented by one token."""
    if token.startswith("--") and token != "--":
        option = token.split("=", 1)[0]
        if option in long_options:
            resolved = option
        else:
            matches = [
                candidate for candidate in long_options if candidate.startswith(option)
            ]
            resolved = matches[0] if len(matches) == 1 else None
        if resolved in managed_long_options:
            return resolved
        return None

    if not token.startswith("-") or token == "-":
        return None

    # argparse accepts attached values for one-character short options. Longer
    # short options accept a separate value or an equals-attached value.
    for option, managed in short_options.items():
        if len(option) == 2 and token.startswith(option):
            return managed
    for option, managed in short_options.items():
        if token == option or token.startswith(f"{option}="):
            return managed
    return None


def reject_managed_options(
    extra,
    *,
    config_key,
    tool,
    long_options,
    managed_long_options,
    short_options,
    explanation,
):
    """Reject every accepted short, long, equals, attached, or abbreviated form."""
    extra, tokens = split_extra(extra, config_key)
    conflicts = []
    options_enabled = True
    for token in tokens:
        if options_enabled and token == "--":
            options_enabled = False
            continue
        if not options_enabled:
            continue
        managed = managed_argparse_option(
            token,
            long_options=long_options,
            managed_long_options=managed_long_options,
            short_options=short_options,
        )
        if managed is not None:
            conflicts.append(token)

    if conflicts:
        rendered = ", ".join(repr(token) for token in dict.fromkeys(conflicts))
        raise ValueError(
            f"{config_key} must not set {tool} options managed by MeeW; "
            f"conflicting token(s): {rendered}. {explanation}"
        )
    return extra
