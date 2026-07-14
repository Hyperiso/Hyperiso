#!/usr/bin/env python3
"""Validate and compare HyperIso CPC reproducibility outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Any

FLOAT_RE = re.compile(
    r"(?<![A-Za-z_])[-+]?(?:(?:\d+\.\d*)|(?:\.\d+)|(?:\d+))(?:[eE][-+]?\d+)?(?![A-Za-z_])"
)
NONFINITE_RE = re.compile(r"(?i)(?<![A-Za-z_])(?:nan|[-+]?inf(?:inity)?)(?![A-Za-z_])")
ABSOLUTE_PATH_RE = re.compile(r"(?:^|\s)(?:/home/|/Users/|[A-Za-z]:\\)")
SUMMARY_RE = re.compile(r"^== (?:Wilson|Observable|Statistic) summary ==$")



def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_reference_metadata(expected: Path, manifest: dict[str, Any]) -> list[str]:
    metadata_path = expected / "reference_metadata.json"
    if not metadata_path.exists():
        return [f"missing frozen reference metadata {metadata_path}"]
    try:
        metadata = json.loads(metadata_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        return [f"invalid frozen reference metadata: {exc}"]

    errors: list[str] = []
    expected_version = manifest.get("hyperiso_version")
    if expected_version and metadata.get("hyperiso_version") != expected_version:
        errors.append(
            f"reference version differs: expected {expected_version}, "
            f"got {metadata.get('hyperiso_version')}"
        )

    file_metadata = metadata.get("files", {})
    required_files = [example["stdout"] for example in manifest["examples"]]
    required_files.extend(
        example["samples_csv"]["file"]
        for example in manifest["examples"]
        if "samples_csv" in example
    )
    for name in sorted(set(required_files)):
        path = expected / name
        if not path.exists():
            errors.append(f"missing frozen reference {path}")
            continue
        recorded = file_metadata.get(name, {}).get("sha256")
        if not recorded:
            errors.append(f"missing SHA-256 metadata for {name}")
            continue
        actual = sha256(path)
        if actual != recorded:
            errors.append(f"SHA-256 mismatch for {name}: expected {recorded}, got {actual}")
    return errors


def validate_text(path: Path) -> list[str]:
    errors: list[str] = []
    if not path.exists():
        return [f"missing output {path}"]
    text = path.read_text(errors="replace")
    lines = text.splitlines()
    if not lines or not SUMMARY_RE.fullmatch(lines[0]):
        errors.append("output is not normalized: expected a summary header on the first line")
    if ABSOLUTE_PATH_RE.search(text):
        errors.append("output contains a machine-specific absolute path")
    if NONFINITE_RE.search(text):
        errors.append("output contains NaN or infinity")
    return errors


def compare_text(ref_path: Path, out_path: Path, abs_tol: float, rel_tol: float) -> list[str]:
    ref_lines = ref_path.read_text(errors="replace").splitlines()
    out_lines = out_path.read_text(errors="replace").splitlines()
    errors: list[str] = []
    if len(ref_lines) != len(out_lines):
        return [f"line count differs: expected {len(ref_lines)}, got {len(out_lines)}"]

    for line_no, (ref_line, out_line) in enumerate(zip(ref_lines, out_lines), start=1):
        ref_matches = list(FLOAT_RE.finditer(ref_line))
        out_matches = list(FLOAT_RE.finditer(out_line))
        ref_shape = FLOAT_RE.sub("<NUM>", ref_line)
        out_shape = FLOAT_RE.sub("<NUM>", out_line)
        if ref_shape != out_shape or len(ref_matches) != len(out_matches):
            errors.append(
                f"line {line_no} structure differs:\n"
                f"  expected: {ref_line}\n"
                f"  got:      {out_line}"
            )
            continue

        for value_no, (ref_match, out_match) in enumerate(zip(ref_matches, out_matches), start=1):
            expected = float(ref_match.group(0))
            obtained = float(out_match.group(0))
            if not (math.isfinite(expected) and math.isfinite(obtained)):
                errors.append(
                    f"line {line_no} value {value_no} is non-finite: "
                    f"expected {expected}, got {obtained}"
                )
            elif not math.isclose(obtained, expected, rel_tol=rel_tol, abs_tol=abs_tol):
                errors.append(
                    f"line {line_no} value {value_no} differs: expected {expected:.17g}, "
                    f"got {obtained:.17g}, |diff|={abs(expected-obtained):.3g}, "
                    f"abs_tol={abs_tol:.3g}, rel_tol={rel_tol:.3g}"
                )
    return errors


def read_csv(path: Path, spec: dict[str, Any]) -> tuple[list[str], list[list[float]], list[str]]:
    errors: list[str] = []
    if not path.exists():
        return [], [], [f"missing CSV {path}"]

    with path.open(newline="") as handle:
        rows = list(csv.reader(handle))
    if not rows:
        return [], [], ["CSV is empty"]

    header = rows[0]
    expected_columns = list(spec.get("columns", []))
    if expected_columns and header != expected_columns:
        errors.append(f"CSV columns differ: expected {expected_columns}, got {header}")

    data: list[list[float]] = []
    for row_no, row in enumerate(rows[1:], start=2):
        if len(row) != len(header):
            errors.append(f"CSV row {row_no} has {len(row)} columns; expected {len(header)}")
            continue
        try:
            values = [float(value) for value in row]
        except ValueError as exc:
            errors.append(f"CSV row {row_no} contains a non-numeric field: {exc}")
            continue
        if not all(math.isfinite(value) for value in values):
            errors.append(f"CSV row {row_no} contains NaN or infinity")
            continue
        data.append(values)

    expected_rows = spec.get("rows")
    if expected_rows is not None and len(data) != int(expected_rows):
        errors.append(f"CSV row count differs: expected {expected_rows}, got {len(data)}")

    lower = spec.get("minimum")
    upper = spec.get("maximum")
    for row_no, values in enumerate(data, start=2):
        for column_no, value in enumerate(values, start=1):
            if lower is not None and value < float(lower):
                errors.append(f"CSV row {row_no}, column {column_no} is below {lower}: {value}")
            if upper is not None and value > float(upper):
                errors.append(f"CSV row {row_no}, column {column_no} is above {upper}: {value}")
            if len(errors) >= 20:
                return header, data, errors
    return header, data, errors


def compare_csv(
    ref_path: Path,
    out_path: Path,
    spec: dict[str, Any],
    abs_tol: float,
    rel_tol: float,
) -> list[str]:
    ref_header, ref_data, ref_errors = read_csv(ref_path, spec)
    out_header, out_data, out_errors = read_csv(out_path, spec)
    errors = [f"reference {error}" for error in ref_errors] + out_errors
    if errors:
        return errors
    if ref_header != out_header or len(ref_data) != len(out_data):
        return ["CSV reference and output shapes differ"]
    for row_no, (ref_row, out_row) in enumerate(zip(ref_data, out_data), start=2):
        for column_no, (expected, obtained) in enumerate(zip(ref_row, out_row), start=1):
            if not math.isclose(obtained, expected, rel_tol=rel_tol, abs_tol=abs_tol):
                errors.append(
                    f"CSV row {row_no}, column {column_no} differs: "
                    f"expected {expected:.17g}, got {obtained:.17g}"
                )
                if len(errors) >= 20:
                    return errors
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--outputs", required=True, type=Path)
    parser.add_argument("--expected", type=Path)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    if not args.validate_only and args.expected is None:
        parser.error("--expected is required unless --validate-only is used")

    manifest = json.loads(args.manifest.read_text())
    defaults = manifest.get("default_tolerances", {})
    default_abs = float(defaults.get("abs_tol", 1e-12))
    default_rel = float(defaults.get("rel_tol", 1e-8))
    all_errors: list[str] = []
    if not args.validate_only and args.expected is not None:
        all_errors.extend(
            f"metadata: {error}"
            for error in validate_reference_metadata(args.expected, manifest)
        )

    for example in manifest["examples"]:
        example_id = example["id"]
        output_path = args.outputs / example["stdout"]
        errors = validate_text(output_path)
        abs_tol = float(example.get("abs_tol", default_abs))
        rel_tol = float(example.get("rel_tol", default_rel))

        if not args.validate_only:
            reference_path = args.expected / example["stdout"]
            if not reference_path.exists():
                errors.append(f"missing frozen reference {reference_path}")
            elif output_path.exists():
                errors.extend(compare_text(reference_path, output_path, abs_tol, rel_tol))

        csv_spec = example.get("samples_csv")
        if csv_spec:
            csv_output = args.outputs / csv_spec["file"]
            _, _, csv_errors = read_csv(csv_output, csv_spec)
            errors.extend(csv_errors)
            if not args.validate_only:
                csv_reference = args.expected / csv_spec["file"]
                if not csv_reference.exists():
                    errors.append(f"missing frozen CSV reference {csv_reference}")
                elif csv_output.exists():
                    errors.extend(
                        compare_csv(
                            csv_reference,
                            csv_output,
                            csv_spec,
                            float(csv_spec.get("abs_tol", abs_tol)),
                            float(csv_spec.get("rel_tol", rel_tol)),
                        )
                    )

        if errors:
            all_errors.extend(f"{example_id}: {error}" for error in errors[:20])
        else:
            mode = "validated" if args.validate_only else "matches reference"
            print(f"[OK] {example_id} {mode}")

    if all_errors:
        for error in all_errors:
            print(f"[FAIL] {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
