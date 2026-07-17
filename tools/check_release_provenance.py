#!/usr/bin/env python3
"""Validate that frozen numerical references were produced by a release build."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

INVALID_MARKERS = (
    "not recorded",
    "unavailable",
    "unknown",
    "to_be",
    "to be",
    "placeholder",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git(root: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=root, text=True, stderr=subprocess.STDOUT
    ).strip()


def is_placeholder(value: object) -> bool:
    text = str(value).strip().lower()
    return not text or any(marker in text for marker in INVALID_MARKERS)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--tag", required=True, help="Release tag, for example v1.0.2")
    args = parser.parse_args()
    root = args.root.resolve()
    expected = root / "reproducibility/expected_outputs"
    metadata_path = expected / "reference_metadata.json"
    errors: list[str] = []

    try:
        metadata = json.loads(metadata_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        print(f"Cannot read {metadata_path}: {exc}", file=sys.stderr)
        return 1

    if metadata.get("source_tag") != args.tag:
        errors.append(
            f"source_tag must be {args.tag!r}, got {metadata.get('source_tag')!r}"
        )

    source_commit = str(metadata.get("source_commit", ""))
    if not re.fullmatch(r"[0-9a-f]{40}", source_commit):
        errors.append("source_commit must be a full 40-character Git commit hash")
    else:
        try:
            git(root, "cat-file", "-e", f"{source_commit}^{{commit}}")
            subprocess.run(
                ["git", "merge-base", "--is-ancestor", source_commit, "HEAD"],
                cwd=root,
                check=True,
                capture_output=True,
            )
            changed = [
                line
                for line in git(
                    root, "diff", "--name-only", f"{source_commit}..HEAD"
                ).splitlines()
                if line
            ]
            unexpected = [
                path
                for path in changed
                if not path.startswith("reproducibility/expected_outputs/")
            ]
            if unexpected:
                errors.append(
                    "non-reference files changed after the recorded build commit: "
                    + ", ".join(unexpected)
                )
        except (subprocess.CalledProcessError, OSError) as exc:
            errors.append(f"cannot validate source_commit {source_commit}: {exc}")

    provenance = metadata.get("provenance", {})
    for key in (
        "operating_system",
        "architecture",
        "python",
        "compiler",
        "cmake",
        "gsl",
        "eigen",
        "build_type",
        "binary",
        "binary_sha256",
    ):
        if is_placeholder(provenance.get(key)):
            errors.append(f"provenance.{key} is missing or still a placeholder")

    binary_sha = str(provenance.get("binary_sha256", ""))
    if (
        binary_sha
        and not is_placeholder(binary_sha)
        and not re.fullmatch(r"[0-9a-f]{64}", binary_sha)
    ):
        errors.append("provenance.binary_sha256 is not a SHA-256 digest")

    timings = metadata.get("case_timings_seconds")
    expected_case_ids = {
        example["id"]
        for example in json.loads((root / "reproducibility/manifest.json").read_text())[
            "examples"
        ]
    }
    if not isinstance(timings, dict) or set(timings) != expected_case_ids:
        errors.append(
            "case_timings_seconds must contain exactly "
            + ", ".join(sorted(expected_case_ids))
        )
    elif any(
        not isinstance(value, (int, float)) or value <= 0 for value in timings.values()
    ):
        errors.append("every reference case timing must be a positive number")

    for name, info in metadata.get("files", {}).items():
        path = expected / name
        if not path.is_file():
            errors.append(f"missing frozen file: {name}")
            continue
        if sha256(path) != info.get("sha256"):
            errors.append(f"SHA-256 mismatch for frozen file: {name}")

    inputs_dir = root / "reproducibility/inputs"
    input_metadata = metadata.get("inputs")
    expected_inputs = {path.name for path in inputs_dir.iterdir() if path.is_file()}
    if not isinstance(input_metadata, dict):
        errors.append("inputs must contain SHA-256 metadata for every frozen input")
    else:
        recorded_inputs = set(input_metadata)
        if recorded_inputs != expected_inputs:
            missing = sorted(expected_inputs - recorded_inputs)
            extra = sorted(recorded_inputs - expected_inputs)
            if missing:
                errors.append("missing frozen input metadata: " + ", ".join(missing))
            if extra:
                errors.append("unexpected frozen input metadata: " + ", ".join(extra))
        for name, info in input_metadata.items():
            path = inputs_dir / name
            if not path.is_file():
                continue
            if sha256(path) != info.get("sha256"):
                errors.append(f"SHA-256 mismatch for frozen input: {name}")

    if errors:
        print("Release provenance validation failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1

    print(f"Release provenance OK for {args.tag}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
