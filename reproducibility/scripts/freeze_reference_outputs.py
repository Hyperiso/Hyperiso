#!/usr/bin/env python3
"""Copy validated outputs to the frozen reference directory and record hashes."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import platform
import subprocess
import tomllib
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    repro = root / "reproducibility"
    outputs = repro / "outputs"
    expected = repro / "expected_outputs"
    expected.mkdir(parents=True, exist_ok=True)

    manifest = json.loads((repro / "manifest.json").read_text())
    files = [example["stdout"] for example in manifest["examples"]]
    files.extend(
        example["samples_csv"]["file"]
        for example in manifest["examples"]
        if "samples_csv" in example
    )
    for name in files:
        (expected / name).write_bytes((outputs / name).read_bytes())

    with (root / "Hyperiso/Hyperiso/pyproject.toml").open("rb") as handle:
        version = tomllib.load(handle)["project"]["version"]
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"

    previous_metadata_path = expected / "reference_metadata.json"
    previous_metadata = {}
    if previous_metadata_path.exists():
        try:
            previous_metadata = json.loads(previous_metadata_path.read_text())
        except json.JSONDecodeError:
            previous_metadata = {}

    def command_version(*command: str) -> str:
        try:
            return subprocess.check_output(
                list(command), text=True, stderr=subprocess.STDOUT
            ).splitlines()[0].strip()
        except (OSError, subprocess.CalledProcessError, IndexError):
            return "unavailable"

    binary = os.environ.get("HYPERISO_BIN")
    binary_path = Path(binary).resolve() if binary else None
    metadata = {
        "suite": manifest["suite_name"],
        "hyperiso_version": version,
        "source_commit": commit,
        "source_tag": os.environ.get("GITHUB_REF_NAME", "unreleased"),
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "normalization": "startup banner and machine-specific paths removed",
        "reference_origin": previous_metadata.get(
            "reference_origin",
            "outputs regenerated from the reviewed release build",
        ),
        "provenance": {
            "operating_system": platform.platform(),
            "architecture": platform.machine(),
            "python": platform.python_version(),
            "compiler": command_version("c++", "--version"),
            "cmake": command_version("cmake", "--version"),
            "gsl": command_version("pkg-config", "--modversion", "gsl"),
            "eigen": command_version(
                "pkg-config", "--modversion", "eigen3"
            ),
            "build_type": os.environ.get("CMAKE_BUILD_TYPE", "Release"),
            "thread_count": int(os.environ.get("HYPERISO_REPRO_THREADS", "1")),
            "binary": str(binary_path) if binary_path else "not recorded",
            "binary_sha256": (
                sha256(binary_path)
                if binary_path and binary_path.is_file()
                else "not recorded"
            ),
            "marty_required": False,
            "external_spectrum_generators_required": False,
        },
        "provenance_note": (
            "R6 and R7 consume archived THDM and SUSY spectra. The release CI "
            "rebuilds the final source and must match all frozen values before publication."
        ),
        "files": {
            name: {"sha256": sha256(expected / name)} for name in sorted(set(files))
        },
    }
    (expected / "reference_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n"
    )
    print(f"Frozen {len(set(files))} reference files in {expected}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
