"""Fail-closed network guard inherited by Annotation test child processes."""

import os

if os.environ.get("MEEW_ANNOTATION_TEST_OFFLINE") == "1":
    import socket
    import urllib.request

    class OfflineNetworkError(RuntimeError):
        pass

    def _blocked(operation):
        def reject(*args, **kwargs):
            raise OfflineNetworkError(
                f"MeeW Annotation test offline guard blocked {operation}"
            )

        reject.__name__ = "meew_offline_reject"
        reject.__module__ = __name__
        return reject

    socket.create_connection = _blocked("socket.create_connection")
    socket.socket.connect = _blocked("socket.connect")
    urllib.request.urlopen = _blocked("urllib.request.urlopen")
    urllib.request.urlretrieve = _blocked("urllib.request.urlretrieve")
