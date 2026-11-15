import errno
import multiprocessing
import os
import shlex
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

from common import native_execution
from common.native_execution import (
    NativeExecutionCancelled,
    NativeExecutionError,
    exclusive_lock,
    run_native_command,
    shared_lock,
)


def python_command(source, *arguments):
    return shlex.join([sys.executable, "-c", source, *map(str, arguments)])


def run_guarded(queue, command, requests):
    started = time.monotonic()
    try:
        run_native_command(command, requests)
    except BaseException as error:  # process-safe test report
        queue.put(
            ("error", type(error).__name__, str(error), time.monotonic() - started)
        )
    else:
        queue.put(("ok", time.monotonic() - started))


def run_guarded_with_probe(queue, probe_command, command, requests):
    started = time.monotonic()

    def verify(probe):
        probe(probe_command)

    try:
        run_native_command(
            command,
            requests,
            prelaunch_probe_check=verify,
        )
    except BaseException as error:  # process-safe test report
        queue.put(
            ("error", type(error).__name__, str(error), time.monotonic() - started)
        )
    else:
        queue.put(("ok", time.monotonic() - started))


def alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


@pytest.mark.parametrize("exclusive", [False, True], ids=["shared", "exclusive"])
def test_shared_readers_overlap_and_exclusive_writers_serialize(tmp_path, exclusive):
    context = multiprocessing.get_context("fork")
    queue = context.Queue()
    payload = tmp_path / "cache"
    request = (
        exclusive_lock(payload, "test") if exclusive else shared_lock(payload, "test")
    )
    command = python_command("import time; time.sleep(0.45)")
    processes = [
        context.Process(target=run_guarded, args=(queue, command, (request,)))
        for _ in range(2)
    ]
    started = time.monotonic()
    for process in processes:
        process.start()
    for process in processes:
        process.join(5)
        assert process.exitcode == 0
    elapsed = time.monotonic() - started
    assert [queue.get(timeout=1)[0] for _ in processes] == ["ok", "ok"]
    if exclusive:
        assert elapsed >= 0.8
    else:
        assert elapsed < 0.8


def test_writer_excludes_reader_and_canonical_aliases_use_same_inode(tmp_path):
    context = multiprocessing.get_context("fork")
    queue = context.Queue()
    real = tmp_path / "real"
    real.mkdir()
    alias = tmp_path / "alias"
    alias.symlink_to(real, target_is_directory=True)
    writer = exclusive_lock(real / "cache", "alias")
    reader = shared_lock(alias / "cache", "alias")
    assert writer.lock_path == reader.lock_path

    command = python_command("import time; time.sleep(0.35)")
    first = context.Process(target=run_guarded, args=(queue, command, (writer,)))
    second = context.Process(target=run_guarded, args=(queue, command, (reader,)))
    started = time.monotonic()
    first.start()
    time.sleep(0.05)
    second.start()
    first.join(5)
    second.join(5)
    assert first.exitcode == second.exitcode == 0
    assert time.monotonic() - started >= 0.6


def test_multiple_locks_are_sorted_to_avoid_deadlock(tmp_path):
    context = multiprocessing.get_context("fork")
    queue = context.Queue()
    first = exclusive_lock(tmp_path / "a", "multi")
    second = exclusive_lock(tmp_path / "b", "multi")
    command = python_command("import time; time.sleep(0.2)")
    processes = [
        context.Process(
            target=run_guarded,
            args=(queue, command, requests),
        )
        for requests in ((first, second), (second, first))
    ]
    for process in processes:
        process.start()
    for process in processes:
        process.join(5)
        assert process.exitcode == 0
    assert [queue.get(timeout=1)[0] for _ in processes] == ["ok", "ok"]


def test_unsupported_lock_fails_before_command_launch(tmp_path, monkeypatch):
    marker = tmp_path / "launched"

    def unsupported(*_args, **_kwargs):
        raise OSError(errno.ENOLCK, "locks unavailable")

    monkeypatch.setattr(native_execution.fcntl, "flock", unsupported)
    with pytest.raises(NativeExecutionError, match="filesystem must support"):
        run_native_command(
            python_command(
                "from pathlib import Path; Path(__import__('sys').argv[1]).touch()",
                marker,
            ),
            (exclusive_lock(tmp_path / "cache", "unsupported"),),
        )
    assert not marker.exists()


def test_unwritable_lock_parent_fails_before_command_launch(tmp_path):
    parent = tmp_path / "read-only"
    parent.mkdir()
    parent.chmod(0o500)
    marker = tmp_path / "launched"
    try:
        with pytest.raises(NativeExecutionError, match="must be writable"):
            run_native_command(
                python_command(
                    "from pathlib import Path; import sys; Path(sys.argv[1]).touch()",
                    marker,
                ),
                (exclusive_lock(parent / "cache", "permissions"),),
            )
    finally:
        parent.chmod(0o700)
    assert not marker.exists()


def test_cancellation_while_waiting_never_launches_command(tmp_path):
    context = multiprocessing.get_context("fork")
    queue = context.Queue()
    payload = tmp_path / "cache"
    holder_ready = tmp_path / "holder-ready"
    waiter_launched = tmp_path / "waiter-launched"
    holder = context.Process(
        target=run_guarded,
        args=(
            queue,
            python_command(
                "from pathlib import Path; import sys,time; Path(sys.argv[1]).touch(); time.sleep(2)",
                holder_ready,
            ),
            (exclusive_lock(payload, "wait"),),
        ),
    )
    holder.start()
    deadline = time.monotonic() + 3
    while not holder_ready.exists() and time.monotonic() < deadline:
        time.sleep(0.02)
    assert holder_ready.exists()

    waiter = context.Process(
        target=run_guarded,
        args=(
            queue,
            python_command(
                "from pathlib import Path; import sys; Path(sys.argv[1]).touch()",
                waiter_launched,
            ),
            (exclusive_lock(payload, "wait"),),
        ),
    )
    waiter.start()
    time.sleep(0.15)
    os.kill(waiter.pid, signal.SIGTERM)
    waiter.join(4)
    assert waiter.exitcode == 0
    results = [queue.get(timeout=3), queue.get(timeout=3)]
    assert any(
        result[:2] == ("error", "NativeExecutionCancelled") for result in results
    )
    assert not waiter_launched.exists()
    holder.join(4)


def test_cancellation_during_prelaunch_check_never_launches(tmp_path):
    marker = tmp_path / "launched"

    def cancel():
        os.kill(os.getpid(), signal.SIGTERM)

    with pytest.raises(NativeExecutionCancelled):
        run_native_command(
            python_command(
                "from pathlib import Path; import sys; Path(sys.argv[1]).touch()",
                marker,
            ),
            (exclusive_lock(tmp_path / "cache", "prelaunch"),),
            prelaunch_check=cancel,
        )
    assert not marker.exists()


def test_readonly_probe_failure_never_launches_native_command(tmp_path):
    marker = tmp_path / "launched"

    def verify(probe):
        probe([sys.executable, "-c", "raise SystemExit(31)"])

    with pytest.raises(subprocess.CalledProcessError) as failure:
        run_native_command(
            python_command(
                "from pathlib import Path; import sys; Path(sys.argv[1]).touch()",
                marker,
            ),
            (shared_lock(tmp_path / "cache", "failed-probe"),),
            prelaunch_probe_check=verify,
        )
    assert failure.value.returncode == 31
    assert not marker.exists()


def test_cancellation_during_probe_cleans_descendants_before_unlock_and_never_launches(
    tmp_path,
):
    context = multiprocessing.get_context("fork")
    queue = context.Queue()
    child_pid_file = tmp_path / "probe-child.pid"
    launched = tmp_path / "launched"
    probe = [
        sys.executable,
        "-c",
        (
            "import signal,subprocess,sys,time; "
            "signal.signal(signal.SIGTERM, signal.SIG_IGN); "
            "child=subprocess.Popen([sys.executable,'-c','import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(30)']); "
            "open(sys.argv[1],'w').write(str(child.pid)); time.sleep(30)"
        ),
        str(child_pid_file),
    ]
    worker = context.Process(
        target=run_guarded_with_probe,
        args=(
            queue,
            probe,
            python_command(
                "from pathlib import Path; import sys; Path(sys.argv[1]).touch()",
                launched,
            ),
            (exclusive_lock(tmp_path / "cache", "probe-cancellation"),),
        ),
    )
    worker.start()
    deadline = time.monotonic() + 3
    while not child_pid_file.exists() and time.monotonic() < deadline:
        time.sleep(0.02)
    assert child_pid_file.exists()

    child_pid = int(child_pid_file.read_text())
    os.kill(worker.pid, signal.SIGTERM)
    worker.join(6)
    assert worker.exitcode == 0
    result = queue.get(timeout=1)
    assert result[:2] == ("error", "NativeExecutionCancelled")
    assert not launched.exists()
    assert not alive(child_pid)

    run_native_command(
        python_command("pass"),
        (exclusive_lock(tmp_path / "cache", "probe-cancellation"),),
    )


def test_child_exec_does_not_inherit_blocked_termination_signals(tmp_path):
    if not hasattr(signal, "pthread_sigmask"):
        pytest.skip("POSIX signal masks are unavailable")
    observed = tmp_path / "mask.txt"
    previous = signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGTERM})
    try:
        run_native_command(
            python_command(
                "import signal,sys; "
                "blocked=signal.pthread_sigmask(signal.SIG_BLOCK, set()); "
                "open(sys.argv[1],'w').write(str(int(signal.SIGTERM in blocked)))",
                observed,
            ),
            (shared_lock(tmp_path / "cache", "signal-mask"),),
        )
        restored = signal.pthread_sigmask(signal.SIG_BLOCK, set())
        assert signal.SIGTERM in restored
    finally:
        signal.pthread_sigmask(signal.SIG_SETMASK, previous)
    assert observed.read_text() == "0"


def test_active_cancellation_kills_owned_process_group_only(tmp_path):
    context = multiprocessing.get_context("fork")
    queue = context.Queue()
    child_pid_file = tmp_path / "child.pid"
    command = python_command(
        "import os,signal,subprocess,sys,time; "
        "signal.signal(signal.SIGTERM, signal.SIG_IGN); "
        "child=subprocess.Popen([sys.executable,'-c','import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(30)']); "
        "open(sys.argv[1],'w').write(str(child.pid)); time.sleep(30)",
        child_pid_file,
    )
    worker = context.Process(
        target=run_guarded,
        args=(queue, command, (exclusive_lock(tmp_path / "cache", "active"),)),
    )
    worker.start()
    deadline = time.monotonic() + 3
    while not child_pid_file.exists() and time.monotonic() < deadline:
        time.sleep(0.02)
    assert child_pid_file.exists()
    child_pid = int(child_pid_file.read_text())
    os.kill(worker.pid, signal.SIGTERM)
    worker.join(6)
    assert worker.exitcode == 0
    result = queue.get(timeout=1)
    assert result[:2] == ("error", "NativeExecutionCancelled")
    assert not alive(child_pid)


def test_parent_exit_with_surviving_grandchild_is_cleaned_before_unlock(tmp_path):
    grandchild_pid_file = tmp_path / "grandchild.pid"
    command = python_command(
        "import subprocess,sys; "
        "child=subprocess.Popen([sys.executable,'-c','import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(30)']); "
        "open(sys.argv[1],'w').write(str(child.pid))",
        grandchild_pid_file,
    )

    run_native_command(
        command,
        (exclusive_lock(tmp_path / "cache", "orphan"),),
    )

    grandchild_pid = int(grandchild_pid_file.read_text())
    assert not alive(grandchild_pid)
    # If cleanup had released early, this independent exclusive acquisition
    # would either overlap the orphan or block indefinitely.
    run_native_command(
        python_command("pass"),
        (exclusive_lock(tmp_path / "cache", "orphan"),),
    )
