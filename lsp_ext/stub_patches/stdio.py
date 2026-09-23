"""Process-level standard output handling."""
from __future__ import annotations

import contextlib
import ctypes
import os
import sys
from collections.abc import Iterator


@contextlib.contextmanager
def discarded_stdout() -> Iterator[None]:
    """Discard everything written to standard output, compiled code's included.

    ``contextlib.redirect_stdout`` rebinds ``sys.stdout`` alone, which misses
    extension modules writing to file descriptor 1. This points the descriptor
    itself at ``os.devnull``, flushing Python's and C's buffers on the way in and
    out, so nothing written inside surfaces after.
    """
    libc = ctypes.CDLL(None)
    sys.stdout.flush()
    libc.fflush(None)
    saved = os.dup(1)
    try:
        with open(os.devnull, "w") as null:
            os.dup2(null.fileno(), 1)
        yield
    finally:
        sys.stdout.flush()
        libc.fflush(None)
        os.dup2(saved, 1)
        os.close(saved)
