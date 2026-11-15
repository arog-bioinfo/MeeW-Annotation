import os
import re


def rule_count(output, rule):
    counts = re.findall(rf"^{re.escape(rule)}\s+([0-9]+)$", output, flags=re.MULTILINE)
    return int(counts[0]) if counts else 0


def file_snapshot(root):
    return {
        path.relative_to(root): (path.read_bytes(), path.stat().st_mtime_ns)
        for path in root.rglob("*")
        if path.is_file()
    }


def material_snapshot(root):
    snapshot = {}
    for path in root.rglob("*"):
        relative = path.relative_to(root)
        if path.is_symlink():
            snapshot[relative] = (
                "symlink",
                os.readlink(path),
                path.lstat().st_mtime_ns,
            )
        elif path.is_file():
            snapshot[relative] = ("file", path.read_bytes(), path.stat().st_mtime_ns)
    return snapshot
