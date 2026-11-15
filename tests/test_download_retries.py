import concurrent.futures
import subprocess

import pytest

from common import download_retries


def recording_log(messages):
    def append(_log_path, message):
        messages.append(message)

    return append


def test_retry_download_uses_bounded_delays_and_logs_attempts(monkeypatch):
    attempts = []
    messages = []
    sleeps = []

    def operation():
        attempts.append(len(attempts) + 1)
        if len(attempts) < 3:
            raise TimeoutError(f"transient-{len(attempts)}")
        return "ready"

    monkeypatch.setattr(download_retries.time, "sleep", sleeps.append)
    result = download_retries.retry_download(
        operation,
        label="fixture database",
        is_retryable=lambda error, output: isinstance(error, TimeoutError)
        and not output,
        append_log=recording_log(messages),
    )

    assert result == "ready"
    assert attempts == [1, 2, 3]
    assert sleeps == [10, 30]
    assert any("attempt 1/4" in message for message in messages)
    assert any("retrying in 10 seconds" in message for message in messages)
    assert any("attempt 3/4" in message for message in messages)


def test_retry_download_exhaustion_raises_original_final_exception(monkeypatch):
    errors = [TimeoutError(f"failure-{index}") for index in range(4)]
    sleeps = []
    attempts = []

    def operation():
        error = errors[len(attempts)]
        attempts.append(error)
        raise error

    monkeypatch.setattr(download_retries.time, "sleep", sleeps.append)
    with pytest.raises(TimeoutError) as raised:
        download_retries.retry_download(
            operation,
            label="fixture database",
            is_retryable=lambda error, output: True,
            append_log=lambda *args: None,
        )

    assert raised.value is errors[-1]
    assert attempts == errors
    assert sleeps == [10, 30, 60]


def test_retry_download_does_not_retry_unclassified_failure(monkeypatch):
    error = ValueError("not a transfer failure")
    sleeps = []
    monkeypatch.setattr(download_retries.time, "sleep", sleeps.append)

    with pytest.raises(ValueError) as raised:
        download_retries.retry_download(
            lambda: (_ for _ in ()).throw(error),
            label="fixture database",
            is_retryable=lambda failure, output: False,
            append_log=lambda *args: None,
        )

    assert raised.value is error
    assert sleeps == []


@pytest.mark.parametrize(
    "error",
    [KeyboardInterrupt(), SystemExit(2), concurrent.futures.CancelledError()],
    ids=["keyboard-interrupt", "system-exit", "cancelled"],
)
def test_retry_download_never_swallows_cancellation(error, monkeypatch):
    sleeps = []
    classifications = []
    monkeypatch.setattr(download_retries.time, "sleep", sleeps.append)

    with pytest.raises(type(error)) as raised:
        download_retries.retry_download(
            lambda: (_ for _ in ()).throw(error),
            label="fixture database",
            is_retryable=lambda *args: classifications.append(args) or True,
            append_log=lambda *args: None,
        )

    assert raised.value is error
    assert classifications == []
    assert sleeps == []


@pytest.mark.parametrize("returncode", [-2, -15, 130, 143])
def test_retry_download_never_retries_subprocess_signal_termination(
    returncode, monkeypatch
):
    error = subprocess.CalledProcessError(returncode, ["downloader"])
    attempts = []
    classifications = []
    sleeps = []

    def operation():
        attempts.append(1)
        raise error

    monkeypatch.setattr(download_retries.time, "sleep", sleeps.append)
    with pytest.raises(subprocess.CalledProcessError) as raised:
        download_retries.retry_download(
            operation,
            label="fixture database",
            is_retryable=lambda *args: classifications.append(args) or True,
            append_log=lambda *args: None,
        )

    assert raised.value is error
    assert attempts == [1]
    assert classifications == []
    assert sleeps == []


def test_subprocess_classification_receives_only_current_attempt_log(
    tmp_path, monkeypatch
):
    log_path = tmp_path / "download.log"
    log_path.write_text("old fatal evidence\n", encoding="utf-8")
    outputs = []
    sleeps = []
    attempts = []

    def append(path, message):
        with path.open("a", encoding="utf-8") as handle:
            handle.write(message + "\n")

    def operation():
        attempts.append(len(attempts) + 1)
        if len(attempts) == 1:
            with log_path.open("a", encoding="utf-8") as handle:
                handle.write("current transient evidence\n")
            raise RuntimeError("failed")
        return "ready"

    def classify(_error, output):
        outputs.append(output)
        return (
            "current transient evidence" in output
            and "old fatal evidence" not in output
        )

    monkeypatch.setattr(download_retries.time, "sleep", sleeps.append)
    assert (
        download_retries.retry_download(
            operation,
            label="fixture database",
            is_retryable=classify,
            append_log=append,
            log_path=log_path,
        )
        == "ready"
    )

    assert attempts == [1, 2]
    assert sleeps == [10]
    assert outputs == ["current transient evidence\n"]
