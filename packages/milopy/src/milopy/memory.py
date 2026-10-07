"""System memory management utilities for memory-safe single-cell processing.
"""

from __future__ import annotations

import ctypes
import gc


def release_system_memory() -> None:
    """Forces Python garbage collection and glibc memory arena trimming back to the OS kernel."""
    gc.collect()
    try:
        libc = ctypes.CDLL("libc.so.6")
        libc.malloc_trim(0)
    except Exception:
        pass
