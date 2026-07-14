"""Portable path resolution shared by MeeW-Annotation workflow rules."""

import os
from pathlib import Path


def expand_path(value):
    """Expand environment variables and ``~`` while preserving relative paths."""
    if value is None or str(value).strip() == "":
        return ""
    return os.path.expanduser(os.path.expandvars(str(value)))


def _first_environment_path(names):
    if isinstance(names, str):
        names = (names,)
    for name in names:
        value = os.environ.get(name, "")
        if value.strip():
            return expand_path(value)
    return ""


def resources_root():
    """Return the configured MeeW resource root using XDG-compatible defaults."""
    meew_resources = os.environ.get("MEEW_RESOURCES", "")
    if meew_resources.strip():
        return Path(expand_path(meew_resources))
    xdg_data_home = os.environ.get("XDG_DATA_HOME", "")
    if xdg_data_home.strip():
        return Path(expand_path(xdg_data_home)) / "meew"
    return Path(expand_path("~/.local/share/meew"))


def resolve_resource_path(explicit, tool_environment, resource_subpath):
    """Resolve explicit, tool environment, then MeeW/XDG resource locations."""
    if explicit is not None and str(explicit).strip():
        return expand_path(explicit)
    tool_path = _first_environment_path(tool_environment)
    if tool_path:
        return tool_path
    return str(resources_root() / resource_subpath)
