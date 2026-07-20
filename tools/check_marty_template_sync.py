#!/usr/bin/env python3
"""Validate or synchronize the two MARTY template trees.

The Python package assets are the runtime source of truth because installed
``pyhyperiso`` wheels resolve templates from ``pyhyperiso/assets``.  The
repository-level ``Assets/template/MARTY`` directory is kept as a developer and
C++ source-tree mirror.
"""

from __future__ import annotations

import argparse
import filecmp
import shutil
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
CANONICAL_DIR = (
    REPOSITORY_ROOT
    / "Hyperiso"
    / "Hyperiso"
    / "pyhyperiso"
    / "assets"
    / "template"
    / "MARTY"
)
MIRROR_DIR = REPOSITORY_ROOT / "Assets" / "template" / "MARTY"


def relative_files(directory: Path) -> set[Path]:
    return {
        path.relative_to(directory) for path in directory.rglob("*") if path.is_file()
    }


def differences() -> list[str]:
    canonical_files = relative_files(CANONICAL_DIR)
    mirror_files = relative_files(MIRROR_DIR)
    messages: list[str] = []

    for path in sorted(canonical_files - mirror_files):
        messages.append(f"missing from mirror: {path}")
    for path in sorted(mirror_files - canonical_files):
        messages.append(f"unexpected in mirror: {path}")
    for path in sorted(canonical_files & mirror_files):
        if not filecmp.cmp(CANONICAL_DIR / path, MIRROR_DIR / path, shallow=False):
            messages.append(f"content differs: {path}")

    return messages


def synchronize() -> None:
    if MIRROR_DIR.exists():
        shutil.rmtree(MIRROR_DIR)
    shutil.copytree(CANONICAL_DIR, MIRROR_DIR)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write",
        action="store_true",
        help="replace Assets/template/MARTY with the packaged template tree",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not CANONICAL_DIR.is_dir():
        print(
            f"canonical MARTY template directory not found: {CANONICAL_DIR}",
            file=sys.stderr,
        )
        return 2
    if args.write:
        synchronize()

    mismatches = differences()
    if mismatches:
        print("MARTY template trees are not synchronized:", file=sys.stderr)
        for mismatch in mismatches:
            print(f"  - {mismatch}", file=sys.stderr)
        print(
            "Run `python tools/check_marty_template_sync.py --write` and commit both trees.",
            file=sys.stderr,
        )
        return 1

    print("MARTY template trees are synchronized.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
