#!/usr/bin/env python3
"""Copy validated outputs to the frozen reference directory and record hashes."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
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

    metadata = {
        "suite": manifest["suite_name"],
        "hyperiso_version": version,
        "source_commit": commit,
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "normalization": "startup banner and machine-specific paths removed",
        "files": {name: {"sha256": sha256(expected / name)} for name in sorted(set(files))},
    }
    (expected / "reference_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n"
    )
    print(f"Frozen {len(set(files))} reference files in {expected}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
