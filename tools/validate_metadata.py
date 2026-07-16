#!/usr/bin/env python3
"""Validate repository-owned JSON, TOML and YAML metadata files."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tomllib
from pathlib import Path

import yaml

SKIP_PARTS = {
    ".git",
    "build",
    "dist",
    "wheelhouse",
    "third_party",
    "2HDMC",
    "node_modules",
}


def tracked_files(root: Path) -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=root,
        check=True,
        capture_output=True,
    )
    return [root / name.decode() for name in result.stdout.split(b"\0") if name]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    args = parser.parse_args()
    root = args.root.resolve()
    errors: list[str] = []

    for path in tracked_files(root):
        if not path.exists() or any(part in SKIP_PARTS for part in path.parts):
            continue
        suffix = path.suffix.lower()
        try:
            if suffix == ".json":
                json.loads(path.read_text(encoding="utf-8"))
            elif suffix == ".toml":
                with path.open("rb") as handle:
                    tomllib.load(handle)
            elif suffix in {".yaml", ".yml"}:
                yaml.safe_load(path.read_text(encoding="utf-8"))
        except Exception as exc:  # report every malformed metadata file together
            errors.append(f"{path.relative_to(root)}: {exc}")

    if errors:
        print("Metadata validation failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1

    print("Repository metadata validation OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
