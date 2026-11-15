import shlex

_MANAGED_LONG_OPTIONS = ("--prefix", "--meta", "--output", "--force")
# Bakta 1.12.0 uses argparse abbreviations. These are the unambiguous spellings
# it accepts for the two output-suppressing flags, in addition to their full names.
_MANAGED_LONG_ALIASES = {
    "--skip-cd": "--skip-cds",
    "--skip-cds": "--skip-cds",
    "--skip-pl": "--skip-plot",
    "--skip-plo": "--skip-plot",
    "--skip-plot": "--skip-plot",
}
_VALUE_SHORT_OPTIONS = {"d", "m", "r", "t"}


def _managed_bakta_option(token):
    if token.startswith("--"):
        option = token.split("=", 1)[0]
        if option in _MANAGED_LONG_ALIASES:
            return _MANAGED_LONG_ALIASES[option]
        matches = [
            managed for managed in _MANAGED_LONG_OPTIONS if managed.startswith(option)
        ]
        if len(matches) == 1:
            return matches[0]
        return None

    # Bakta/argparse accepts attached values and clusters of short options.
    if token.startswith("-") and token != "-":
        for option in token[1:]:
            if option == "p":
                return "-p"
            if option == "o":
                return "-o"
            if option == "f":
                return "-f"
            if option in _VALUE_SHORT_OPTIONS:
                return None
            if option not in {"h", "v"}:
                return None
    return None


def validate_bakta_extra(extra):
    """Reject Bakta options whose behavior is managed by the workflow."""
    try:
        tokens = shlex.split(extra, posix=True, comments=False)
    except ValueError as error:
        raise ValueError(
            f"config.bakta.extra is not valid shell syntax: {error}"
        ) from error

    conflicts = []
    for token in tokens:
        managed = _managed_bakta_option(token)
        if managed is not None and token not in conflicts:
            conflicts.append(token)

    if conflicts:
        rendered = ", ".join(repr(token) for token in conflicts)
        raise ValueError(
            "config.bakta.extra must not set Bakta options managed by MeeW "
            "(--prefix/-p, --meta, --output/-o, --force/-f, --skip-cds, and "
            f"--skip-plot); conflicting token(s): {rendered}. MeeW fixes the "
            "native output prefix to 'annotation', derives --meta from "
            "input.mode, controls direct output-directory and overwrite "
            "behavior, and requires the native CDS-dependent and isolate plot "
            "outputs."
        )

    return extra
