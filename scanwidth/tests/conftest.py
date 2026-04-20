"""Pytest configuration for local src-layout imports."""

from __future__ import annotations

import sys
from pathlib import Path


def pytest_sessionstart(session: object) -> None:
    """Add ``src`` to ``sys.path`` for test discovery and imports.

    Parameters
    ----------
    session : object
        Active pytest session object.
    """
    root_dir = Path(__file__).resolve().parents[1]
    src_dir = root_dir / "src"
    if str(src_dir) not in sys.path:
        sys.path.insert(0, str(src_dir))
