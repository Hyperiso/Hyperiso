"""Public Python package for HyperIso."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("pyhyperiso")
except PackageNotFoundError:
    # Source-tree fallback used before an editable or wheel installation.
    __version__ = "1.0.0"

__all__ = ["__version__"]
