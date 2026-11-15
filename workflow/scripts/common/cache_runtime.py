"""Failure-safe primitives for persistent, adapter-owned resource caches."""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path


class CacheIntegrityError(RuntimeError):
    """A selected cache exists but is incomplete, incompatible, or mismatched."""


@dataclass(frozen=True)
class CacheResult:
    """Describe what an adapter did with one persistent cache selection."""

    status: str
    path: str
    identity: str = "verified"
    prepared: tuple[str, ...] = field(default_factory=tuple)

    def as_dict(self):
        return {
            "status": self.status,
            "resolved_path": self.path,
            "identity": self.identity,
            "prepared": list(self.prepared),
        }


def resolved_path(path):
    """Resolve a configured path before an adapter changes its working directory."""
    return Path(path).expanduser().resolve()


def readable_nonempty_file(path):
    path = Path(path)
    try:
        return path.is_file() and path.stat().st_size > 0 and os.access(path, os.R_OK)
    except OSError:
        return False


def append_log(log_path, message):
    """Append one cache lifecycle message outside the persistent payload."""
    if not log_path:
        return
    log_path = resolved_path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(f"[MeeW cache] {message}\n")


def run_logged(argv, log_path, *, env=None, cwd=None):
    """Run an installer without a shell and append its output to the job log."""
    argv = [str(argument) for argument in argv]
    append_log(log_path, f"$ {shlex.join(argv)}")
    log_path = resolved_path(log_path)
    with log_path.open("a", encoding="utf-8") as handle:
        subprocess.run(
            argv,
            check=True,
            env=env,
            cwd=str(cwd) if cwd else None,
            stdout=handle,
            stderr=subprocess.STDOUT,
        )


def cache_lock_path(payload, namespace):
    """Return a cache-keyed lock beside, never inside, the owned payload."""
    payload = resolved_path(payload)
    key = hashlib.sha256(str(payload).encode()).hexdigest()[:20]
    return payload.parent / f".meew-{namespace}-{key}.lock"


@contextmanager
def cache_lock(payload, namespace):
    """Serialize first preparation on shared POSIX filesystems."""
    lock_path = cache_lock_path(payload, namespace)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield lock_path
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


@contextmanager
def staging_directory(payload, namespace):
    """Create and clean an isolated sibling staging directory."""
    payload = resolved_path(payload)
    payload.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".meew-{namespace}-stage-", dir=payload.parent)
    )
    try:
        yield staging
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def publish_directory(staged, destination):
    """Atomically publish one validated directory without replacing anything."""
    staged = Path(staged)
    destination = resolved_path(destination)
    if destination.exists() or destination.is_symlink():
        raise CacheIntegrityError(
            f"Refusing to replace existing persistent cache payload: {destination}"
        )
    os.replace(staged, destination)


def ensure_fresh_directory_cache(
    destination, namespace, validate, prepare, log_path, resource, prepared=(), **report
):
    """Reuse or atomically prepare one missing directory cache."""
    destination = resolved_path(destination)

    def ready(path):
        try:
            return validate(path)
        except CacheIntegrityError:
            return None

    state = ready(destination)
    status = "reused"
    if state is None:
        with cache_lock(destination, namespace):
            state = ready(destination)
            if state is not None:
                status = "waiting"
            else:
                if destination.exists() or destination.is_symlink():
                    validate(destination)
                with staging_directory(destination, namespace) as staging:
                    staged = Path(prepare(staging))
                    state = validate(staged)
                    publish_directory(staged, destination)
                status = "prepared"

    reported_path = report.get("reported_path", destination)
    detail = report.get("detail")
    suffix = f" ({detail(state)})" if detail and status != "waiting" else ""
    if status == "waiting":
        reuse = (
            f"reused {reported_path}"
            if report.get("waiting_reused_only")
            else f"reused read-only from {reported_path}"
        )
        message = f"{resource} database became ready while waiting; {reuse}."
    else:
        action = "reused read-only from" if status == "reused" else "prepared at"
        message = f"{resource} database {action} {reported_path}{suffix}."
    append_log(log_path, message)
    result_status = "prepared" if status == "prepared" else "reused"
    result_prepared = tuple(prepared) if status == "prepared" else ()
    return CacheResult(result_status, str(destination), prepared=result_prepared)


def atomic_write_json(path, value):
    """Write cache metadata atomically; callers use this only during preparation."""
    path = resolved_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()
