"""Lock-aware supervision for the original native annotation commands."""

from __future__ import annotations

import ctypes
import errno
import fcntl
import os
import signal
import subprocess
import time
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass
from pathlib import Path

from .cache_runtime import cache_lock_path, resolved_path

_TERMINATION_SIGNALS = tuple(
    candidate
    for candidate in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)
    if candidate is not None
)
_LOCK_RETRY_SECONDS = 0.05
_TERMINATION_GRACE_SECONDS = 2.0
_PR_GET_CHILD_SUBREAPER = 37
_PR_SET_CHILD_SUBREAPER = 36


class NativeExecutionError(RuntimeError):
    """A guard could not safely lock or supervise the native command."""


class NativeExecutionCancelled(NativeExecutionError):
    """Execution was cancelled before a native command could complete."""

    def __init__(self, signum):
        self.signum = signum
        super().__init__(
            f"Native execution cancelled by signal {signal.Signals(signum).name}."
        )


@dataclass(frozen=True)
class LockRequest:
    """One canonical shared-resource lock and its required access mode."""

    payload: Path
    namespace: str
    exclusive: bool
    reason: str = ""

    @property
    def lock_path(self):
        return cache_lock_path(self.payload, self.namespace)


def shared_lock(payload, namespace, reason=""):
    return LockRequest(resolved_path(payload), namespace, False, reason)


def exclusive_lock(payload, namespace, reason=""):
    return LockRequest(resolved_path(payload), namespace, True, reason)


def _normalized_lock_requests(requests):
    """Deduplicate aliases and sort lock inodes to prevent lock-order deadlocks."""
    combined = {}
    for request in requests:
        lock_path = request.lock_path.resolve()
        previous = combined.get(lock_path)
        if previous is None:
            combined[lock_path] = request
        elif request.exclusive and not previous.exclusive:
            combined[lock_path] = request
    return tuple(combined[path] for path in sorted(combined, key=str))


def ensure_lock_coverage(held_requests, required_requests):
    """Fail if a locked-state reassessment needs a lock mode not already held."""
    held = {
        request.lock_path.resolve(): request
        for request in _normalized_lock_requests(held_requests)
    }
    for required in _normalized_lock_requests(required_requests):
        path = required.lock_path.resolve()
        current = held.get(path)
        if current is None or required.exclusive and not current.exclusive:
            mode = "exclusive" if required.exclusive else "shared"
            raise NativeExecutionError(
                "Shared resource state changed while locks were being acquired; "
                f"the command now requires a new {mode} lock at {path}. Retry the "
                "job so readiness can be assessed again before native execution."
            )


class _Cancellation:
    def __init__(self):
        self.signum = None
        self.process_group = None

    @property
    def requested(self):
        return self.signum is not None

    def handler(self, signum, _frame):
        if self.signum is None:
            self.signum = signum
        if self.process_group is not None:
            try:
                os.killpg(self.process_group, signum)
            except ProcessLookupError:
                pass

    def raise_if_requested(self):
        if self.requested:
            raise NativeExecutionCancelled(self.signum)


@contextmanager
def _forwarding_signals(cancellation):
    previous = {}
    previous_mask = None
    try:
        for signum in _TERMINATION_SIGNALS:
            previous[signum] = signal.getsignal(signum)
            signal.signal(signum, cancellation.handler)
        if hasattr(signal, "pthread_sigmask"):
            previous_mask = signal.pthread_sigmask(
                signal.SIG_UNBLOCK, set(_TERMINATION_SIGNALS)
            )
        yield
    finally:
        if previous_mask is not None:
            signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
        for signum, handler in previous.items():
            signal.signal(signum, handler)


def _lock_error(lock_path, operation, error):
    return NativeExecutionError(
        f"Cannot {operation} shared resource lock {lock_path}: {error}. "
        "The resource parent must be writable and its filesystem must support "
        "POSIX advisory flock locks; the native command was not launched."
    )


@contextmanager
def _acquire_one_lock(request, cancellation):
    lock_path = request.lock_path
    try:
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        handle = lock_path.open("a+", encoding="utf-8")
    except OSError as error:
        raise _lock_error(lock_path, "open", error) from error

    operation = fcntl.LOCK_EX if request.exclusive else fcntl.LOCK_SH
    try:
        while True:
            cancellation.raise_if_requested()
            try:
                fcntl.flock(handle.fileno(), operation | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                time.sleep(_LOCK_RETRY_SECONDS)
            except OSError as error:
                if error.errno in (errno.EACCES, errno.EAGAIN):
                    time.sleep(_LOCK_RETRY_SECONDS)
                    continue
                raise _lock_error(lock_path, "acquire", error) from error
        try:
            yield lock_path
        finally:
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            except OSError as error:
                raise _lock_error(lock_path, "release", error) from error
    finally:
        handle.close()


@contextmanager
def _acquire_locks(requests, cancellation):
    with ExitStack() as stack:
        for request in _normalized_lock_requests(requests):
            stack.enter_context(_acquire_one_lock(request, cancellation))
        yield


@contextmanager
def _child_subreaper():
    """Adopt orphaned descendants so cancellation can reap the owned process tree."""
    if os.name != "posix" or not Path("/proc").is_dir():
        yield
        return
    libc = ctypes.CDLL(None, use_errno=True)
    previous = ctypes.c_int()
    if libc.prctl(_PR_GET_CHILD_SUBREAPER, ctypes.byref(previous), 0, 0, 0) != 0:
        error = ctypes.get_errno()
        raise NativeExecutionError(
            f"Cannot inspect Linux child-subreaper state: {os.strerror(error)}"
        )
    changed = previous.value == 0
    if changed and libc.prctl(_PR_SET_CHILD_SUBREAPER, 1, 0, 0, 0) != 0:
        error = ctypes.get_errno()
        raise NativeExecutionError(
            f"Cannot enable Linux child-subreaper supervision: {os.strerror(error)}"
        )
    try:
        yield
    finally:
        if changed and libc.prctl(_PR_SET_CHILD_SUBREAPER, 0, 0, 0, 0) != 0:
            error = ctypes.get_errno()
            raise NativeExecutionError(
                f"Cannot restore Linux child-subreaper state: {os.strerror(error)}"
            )


def _process_group_exists(process_group):
    try:
        os.killpg(process_group, 0)
    except ProcessLookupError:
        return False
    except PermissionError as error:
        raise NativeExecutionError(
            f"Cannot inspect owned native process group {process_group}: {error}"
        ) from error
    return True


def _reap_owned_children(process_group):
    while True:
        try:
            child, _ = os.waitpid(-process_group, os.WNOHANG)
        except ChildProcessError:
            return
        if child == 0:
            return


def _signal_owned_group(process_group, signum):
    try:
        os.killpg(process_group, signum)
    except ProcessLookupError:
        pass


def _clean_owned_process_group(process_group, *, first_signal=signal.SIGTERM):
    """Terminate and prove removal of descendants before resource locks release."""
    if not _process_group_exists(process_group):
        _reap_owned_children(process_group)
        return
    _signal_owned_group(process_group, first_signal)
    deadline = time.monotonic() + _TERMINATION_GRACE_SECONDS
    while time.monotonic() < deadline:
        _reap_owned_children(process_group)
        if not _process_group_exists(process_group):
            return
        time.sleep(0.02)
    _signal_owned_group(process_group, signal.SIGKILL)
    while True:
        _reap_owned_children(process_group)
        if not _process_group_exists(process_group):
            return
        time.sleep(0.02)


def _snakemake_shell_settings(command):
    """Apply the same Bash strict prefix/suffix used by Snakemake's shell helper."""
    from snakemake.shell import shell

    executable = shell.get_executable()
    shell._check_executable(executable)
    guarded = " ".join(
        (shell._get_process_prefix(executable), command, shell._process_suffix)
    ).strip()
    return executable, guarded


def _native_environment():
    """Retain the environment shaping of the original script-level shell call."""
    environment = dict(os.environ)
    for name in (
        "OMP_NUM_THREADS",
        "GOTO_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        environment[name] = "1"
    return environment


def _run_readonly_probe(command, cancellation, *, cwd=None):
    """Run one bounded read-only helper under the caller's locks and signals."""
    if isinstance(command, (str, bytes, os.PathLike)):
        raise NativeExecutionError(
            "Read-only native probes require an argument sequence, not a shell command."
        )
    arguments = tuple(os.fspath(argument) for argument in command)
    if not arguments:
        raise NativeExecutionError("Read-only native probes require an executable.")

    cancellation.raise_if_requested()
    with _child_subreaper():
        cancellation.raise_if_requested()
        process = subprocess.Popen(
            arguments,
            cwd=str(cwd) if cwd is not None else None,
            env=_native_environment(),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            start_new_session=True,
            restore_signals=True,
        )
        process_group = process.pid
        cancellation.process_group = process_group
        returncode = None
        try:
            while returncode is None:
                returncode = process.poll()
                if returncode is None:
                    if cancellation.requested:
                        _clean_owned_process_group(
                            process_group,
                            first_signal=cancellation.signum,
                        )
                        returncode = process.wait()
                        break
                    time.sleep(0.02)
            _reap_owned_children(process_group)
            if _process_group_exists(process_group):
                _clean_owned_process_group(process_group)
        finally:
            if process.poll() is None:
                first_signal = cancellation.signum or signal.SIGTERM
                _clean_owned_process_group(process_group, first_signal=first_signal)
                process.wait()
            elif _process_group_exists(process_group):
                _clean_owned_process_group(process_group)
            cancellation.process_group = None
        if cancellation.requested:
            raise NativeExecutionCancelled(cancellation.signum)
        if returncode:
            raise subprocess.CalledProcessError(returncode, arguments)
        return subprocess.CompletedProcess(arguments, returncode)


def run_native_command(
    command,
    lock_requests,
    *,
    cwd=None,
    prelaunch_check=None,
    prelaunch_probe_check=None,
):
    """Run one original shell command while holding all requested resource locks."""
    cancellation = _Cancellation()
    requests = _normalized_lock_requests(lock_requests)
    with _forwarding_signals(cancellation):
        with _acquire_locks(requests, cancellation):
            cancellation.raise_if_requested()
            if prelaunch_check is not None:
                prelaunch_check()
            cancellation.raise_if_requested()
            if prelaunch_probe_check is not None:
                prelaunch_probe_check(
                    lambda probe_command, *, cwd=None: _run_readonly_probe(
                        probe_command,
                        cancellation,
                        cwd=cwd,
                    )
                )
            cancellation.raise_if_requested()
            launch_command = command() if callable(command) else command
            cancellation.raise_if_requested()
            executable, guarded_command = _snakemake_shell_settings(launch_command)
            with _child_subreaper():
                cancellation.raise_if_requested()
                process = subprocess.Popen(
                    guarded_command,
                    shell=True,
                    executable=executable,
                    cwd=str(cwd) if cwd is not None else None,
                    env=_native_environment(),
                    start_new_session=True,
                    restore_signals=True,
                )
                process_group = process.pid
                cancellation.process_group = process_group
                returncode = None
                try:
                    while returncode is None:
                        returncode = process.poll()
                        if returncode is None:
                            if cancellation.requested:
                                _clean_owned_process_group(
                                    process_group,
                                    first_signal=cancellation.signum,
                                )
                                returncode = process.wait()
                                break
                            time.sleep(0.02)
                    _reap_owned_children(process_group)
                    if _process_group_exists(process_group):
                        _clean_owned_process_group(process_group)
                finally:
                    if process.poll() is None:
                        first_signal = cancellation.signum or signal.SIGTERM
                        _clean_owned_process_group(
                            process_group, first_signal=first_signal
                        )
                        process.wait()
                    elif _process_group_exists(process_group):
                        _clean_owned_process_group(process_group)
                    cancellation.process_group = None
                if cancellation.requested:
                    raise NativeExecutionCancelled(cancellation.signum)
                if returncode:
                    raise subprocess.CalledProcessError(returncode, guarded_command)
                return returncode
