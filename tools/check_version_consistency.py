#!/usr/bin/env python3
"""Fail when release version declarations disagree."""

from __future__ import annotations

import argparse
import json
import re
import sys
import tomllib
from pathlib import Path


def require(pattern: str, text: str, label: str) -> str:
    match = re.search(pattern, text, re.MULTILINE)
    if not match:
        raise RuntimeError(f"could not read {label}")
    return match.group(1)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--tag", help="Optional Git tag, e.g. v1.0.0")
    args = parser.parse_args()
    root = args.root.resolve()

    with (root / "Hyperiso/Hyperiso/pyproject.toml").open("rb") as handle:
        python_version = tomllib.load(handle)["project"]["version"]

    cmake_text = (root / "Hyperiso/Hyperiso/core/CMakeLists.txt").read_text()
    cmake_version = require(
        r"project\(Hyperiso VERSION ([0-9]+\.[0-9]+\.[0-9]+)",
        cmake_text,
        "CMake version",
    )

    cff_text = (root / "CITATION.cff").read_text()
    cff_version = require(
        r'^version:\s*["\']?([^"\'\s]+)', cff_text, "CITATION.cff version"
    )

    init_text = (root / "Hyperiso/Hyperiso/pyhyperiso/__init__.py").read_text()
    python_fallback = require(
        r'__version__\s*=\s*["\']([0-9]+\.[0-9]+\.[0-9]+)["\']',
        init_text,
        "Python fallback version",
    )
    manifest_version = json.loads((root / "reproducibility/manifest.json").read_text())[
        "hyperiso_version"
    ]

    banner_text = (
        root / "Hyperiso/Hyperiso/core/src/Core/infrastructure/HyperisoBanner.h"
    ).read_text()
    binding_text = (
        root / "Hyperiso/Hyperiso/core/src/Python/src/pyhyperiso.cpp"
    ).read_text()
    if "project_version" not in banner_text:
        raise RuntimeError(
            "Hyperiso banner is not sourced from generated project_version"
        )
    if 'm.attr("__version__") = std::string(project_version)' not in binding_text:
        raise RuntimeError(
            "Python binding version is not sourced from generated project_version"
        )

    values = {
        "pyproject.toml": python_version,
        "CMakeLists.txt": cmake_version,
        "CITATION.cff": cff_version,
        "Python fallback": python_fallback,
        "reproducibility manifest": manifest_version,
    }
    if args.tag:
        if not re.fullmatch(r"v[0-9]+\.[0-9]+\.[0-9]+", args.tag):
            print(f"Invalid release tag format: {args.tag}", file=sys.stderr)
            return 1
        values["Git tag"] = args.tag.removeprefix("v")

    expected = python_version
    mismatches = {name: value for name, value in values.items() if value != expected}
    if mismatches:
        print(f"Expected version {expected}; mismatches: {mismatches}", file=sys.stderr)
        return 1

    print(f"Version consistency OK: {expected}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
