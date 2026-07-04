#!/usr/bin/env python3
"""Compare HyperIso reproducibility outputs against frozen references.

The checker extracts floating-point values from each text output and compares the
current run against the corresponding reference with the tolerances declared in
``manifest.json``.  It intentionally ignores non-numeric formatting changes.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Iterable, List

FLOAT_RE = re.compile(
    r"(?<![A-Za-z_])[-+]?(?:(?:\d+\.\d*)|(?:\.\d+)|(?:\d+))(?:[eE][-+]?\d+)?(?![A-Za-z_])"
)


def read_numbers(path: Path) -> List[float]:
    text = path.read_text(errors="replace")
    return [float(match.group(0)) for match in FLOAT_RE.finditer(text)]


def compare_numbers(ref: Iterable[float], got: Iterable[float], abs_tol: float, rel_tol: float) -> list[str]:
    errors: list[str] = []
    ref_list = list(ref)
    got_list = list(got)
    if len(ref_list) != len(got_list):
        errors.append(f"number count differs: expected {len(ref_list)}, got {len(got_list)}")
        return errors
    for idx, (r, g) in enumerate(zip(ref_list, got_list)):
        tol = max(abs_tol, rel_tol * max(abs(r), 1.0))
        if not (math.isfinite(r) and math.isfinite(g)):
            errors.append(f"non-finite value at index {idx}: expected {r}, got {g}")
        elif abs(r - g) > tol:
            errors.append(
                f"value {idx} differs: expected {r:.17g}, got {g:.17g}, "
                f"|diff|={abs(r-g):.3g}, tol={tol:.3g}"
            )
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--outputs", required=True, type=Path)
    parser.add_argument("--expected", required=True, type=Path)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    default_abs = float(manifest.get("default_tolerances", {}).get("abs_tol", 1e-12))
    default_rel = float(manifest.get("default_tolerances", {}).get("rel_tol", 1e-8))

    any_missing = False
    all_errors: list[str] = []

    for example in manifest["examples"]:
        name = example["stdout"]
        ref_path = args.expected / name
        out_path = args.outputs / name
        if not ref_path.exists():
            print(f"[SKIP] {example['id']} no frozen reference output: {ref_path}")
            any_missing = True
            continue
        if not out_path.exists():
            all_errors.append(f"{example['id']}: missing current output {out_path}")
            continue

        abs_tol = float(example.get("abs_tol", default_abs))
        rel_tol = float(example.get("rel_tol", default_rel))
        errors = compare_numbers(read_numbers(ref_path), read_numbers(out_path), abs_tol, rel_tol)
        if errors:
            all_errors.extend(f"{example['id']} {err}" for err in errors[:10])
            if len(errors) > 10:
                all_errors.append(f"{example['id']} ... {len(errors)-10} more differences")
        else:
            print(f"[OK] {example['id']} {name}")

    if all_errors:
        print("\n".join(f"[FAIL] {err}" for err in all_errors), file=sys.stderr)
        return 1

    if any_missing:
        print("Some reference outputs are not frozen yet. Use --update-expected after the release setup is fixed.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
