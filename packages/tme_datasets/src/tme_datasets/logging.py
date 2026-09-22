"""Logging utilities for tme_datasets with seamless Jupyter notebook support."""

from __future__ import annotations

import logging
import sys
from typing import TextIO

_LOGGER_NAME = "tme_datasets"
_DEFAULT_FORMAT = "%(asctime)s [%(levelname)s] [tme_datasets] %(message)s"
_DEFAULT_DATEFMT = "%H:%M:%S"


class _TmeFlushStreamHandler(logging.StreamHandler):
    """StreamHandler that flushes immediately after every write.

    Ensures progress updates appear in real-time in Jupyter notebook output cells
    without buffering delay.
    """

    def emit(self, record: logging.LogRecord) -> None:
        super().emit(record)
        self.flush()


def configure_logging(
    level: int | str = logging.INFO,
    stream: TextIO = sys.stdout,
    propagate: bool = False,
    format_str: str | None = None,
) -> logging.Logger:
    """Configure or re-configure logging for tme_datasets.

    Idempotent: Replaces any previously attached tme_datasets handlers to avoid
    duplicated lines when cells are re-executed in Jupyter notebooks.

    Args:
        level: Logging level (e.g. logging.INFO, "INFO", "DEBUG", "WARNING").
        stream: Output stream. Defaults to sys.stdout (renders cleanly in Jupyter
            without red stderr formatting).
        propagate: Whether to propagate records to the root logger. Defaults to False.
        format_str: Custom log format string. If None, uses default format.

    Returns:
        The configured top-level package logger.
    """
    logger = logging.getLogger(_LOGGER_NAME)

    if isinstance(level, str):
        level = getattr(logging, level.upper(), logging.INFO)

    logger.setLevel(level)
    logger.propagate = propagate

    # Remove any existing handlers attached to tme_datasets to maintain idempotency
    for handler in list(logger.handlers):
        if isinstance(handler, _TmeFlushStreamHandler):
            logger.removeHandler(handler)

    fmt = format_str or _DEFAULT_FORMAT
    formatter = logging.Formatter(fmt, datefmt=_DEFAULT_DATEFMT)

    handler = _TmeFlushStreamHandler(stream)
    handler.setLevel(level)
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    return logger


def set_log_level(level: int | str) -> None:
    """Set the logging verbosity for all tme_datasets loggers.

    Args:
        level: Logging level (e.g. "DEBUG", "INFO", "WARNING", "ERROR", logging.INFO).
    """
    if isinstance(level, str):
        level = getattr(logging, level.upper(), logging.INFO)

    logger = logging.getLogger(_LOGGER_NAME)
    logger.setLevel(level)
    for handler in logger.handlers:
        handler.setLevel(level)


def get_logger(name: str = _LOGGER_NAME) -> logging.Logger:
    """Retrieve a logger scoped under the tme_datasets namespace.

    If logging has not yet been configured, initializes it with defaults.
    """
    pkg_logger = logging.getLogger(_LOGGER_NAME)
    if not any(isinstance(h, _TmeFlushStreamHandler) for h in pkg_logger.handlers):
        configure_logging(level=logging.INFO)

    if name == _LOGGER_NAME:
        return pkg_logger
    if not name.startswith(f"{_LOGGER_NAME}."):
        name = f"{_LOGGER_NAME}.{name}"
    return logging.getLogger(name)
