"""Small bounded retry policy for explicit database transfer operations."""

from __future__ import annotations

import asyncio
import concurrent.futures
import signal
import subprocess
import time
from pathlib import Path

RETRY_DELAYS = (10, 30, 60)
INTERRUPTED_SUBPROCESS_RETURN_CODES = frozenset(
    128 + int(signum) for signum in signal.valid_signals()
)


def _subprocess_was_interrupted(error):
    return isinstance(error, subprocess.CalledProcessError) and (
        error.returncode < 0 or error.returncode in INTERRUPTED_SUBPROCESS_RETURN_CODES
    )


def _current_log_output(log_path, offset):
    if log_path is None:
        return ""
    try:
        with Path(log_path).expanduser().resolve().open("rb") as handle:
            handle.seek(offset)
            return handle.read().decode("utf-8", errors="replace")
    except OSError:
        return ""


def _log_size(log_path):
    if log_path is None:
        return 0
    try:
        return Path(log_path).expanduser().resolve().stat().st_size
    except OSError:
        return 0


def retry_download(
    operation,
    *,
    label,
    is_retryable,
    append_log,
    log_path=None,
    sleep=None,
):
    """Run one transfer up to four times and preserve the last raised exception.

    When ``log_path`` is supplied, ``is_retryable`` receives only text appended by
    the current attempt. It is called as ``is_retryable(error, attempt_output)``.
    """
    delays = RETRY_DELAYS
    attempts = len(delays) + 1
    sleeper = time.sleep if sleep is None else sleep

    for attempt in range(1, attempts + 1):
        append_log(log_path, f"{label} download attempt {attempt}/{attempts}.")
        offset = _log_size(log_path)
        try:
            return operation()
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as error:
            if isinstance(
                error, (asyncio.CancelledError, concurrent.futures.CancelledError)
            ) or _subprocess_was_interrupted(error):
                raise
            attempt_output = _current_log_output(log_path, offset)
            retryable = is_retryable(error, attempt_output)
            if not retryable:
                append_log(
                    log_path,
                    f"{label} download attempt {attempt}/{attempts} failed; "
                    "not retrying this failure.",
                )
                raise
            if attempt == attempts:
                append_log(
                    log_path,
                    f"{label} download attempt {attempt}/{attempts} failed; "
                    "retry limit exhausted.",
                )
                raise
            delay = delays[attempt - 1]
            append_log(
                log_path,
                f"{label} download attempt {attempt}/{attempts} failed; retrying "
                f"in {delay} seconds.",
            )
            sleeper(delay)

    raise AssertionError("unreachable retry state")
