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


def resources_root(configured_root=None):
    """Return the configured or ambient MeeW resource root."""
    if configured_root is not None and str(configured_root).strip():
        return Path(expand_path(configured_root))
    meew_resources = os.environ.get("MEEW_RESOURCES", "")
    if meew_resources.strip():
        return Path(expand_path(meew_resources))
    xdg_data_home = os.environ.get("XDG_DATA_HOME", "")
    if xdg_data_home.strip():
        return Path(expand_path(xdg_data_home)) / "meew"
    return Path(expand_path("~/resources"))


def resolve_resource_path(
    explicit, tool_environment, resource_subpath, configured_root=None
):
    """Resolve explicit, configured-root, tool, then ambient resource paths."""
    if explicit is not None and str(explicit).strip():
        return expand_path(explicit)
    if configured_root is not None and str(configured_root).strip():
        return str(resources_root(configured_root) / resource_subpath)
    tool_path = _first_environment_path(tool_environment)
    if tool_path:
        return tool_path
    return str(resources_root() / resource_subpath)
