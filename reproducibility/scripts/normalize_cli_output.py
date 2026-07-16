#!/usr/bin/env python3
"""Remove machine-specific startup diagnostics from a CLI reference output."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

ANSI_RE = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")
SUMMARY_RE = re.compile(r"^== (?:Wilson|Observable|Statistic) summary ==$")
DIAGNOSTIC_RE = re.compile(r"(?i)(?:\bWARN(?:ING)?\b|\bERROR\b|\bFATAL\b)")
NONFINITE_RE = re.compile(r"(?i)(?<![A-Za-z_])(?:nan|[-+]?inf(?:inity)?)(?![A-Za-z_])")


def normalize(text: str) -> str:
    cleaned = ANSI_RE.sub("", text)
    lines = [line.rstrip() for line in cleaned.splitlines()]

    diagnostics = [line for line in lines if DIAGNOSTIC_RE.search(line)]
    if diagnostics:
        preview = "\n".join(f"  {line}" for line in diagnostics[:10])
        suffix = "\n  ..." if len(diagnostics) > 10 else ""
        raise ValueError(
            "CLI output contains a warning, error or fatal diagnostic:\n"
            f"{preview}{suffix}"
        )

    nonfinite = [line for line in lines if NONFINITE_RE.search(line)]
    if nonfinite:
        preview = "\n".join(f"  {line}" for line in nonfinite[:10])
        suffix = "\n  ..." if len(nonfinite) > 10 else ""
        raise ValueError(f"CLI output contains NaN or infinity:\n{preview}{suffix}")
    try:
        start = next(
            index for index, line in enumerate(lines) if SUMMARY_RE.fullmatch(line)
        )
    except StopIteration as exc:
        raise ValueError(
            "CLI output does not contain a supported summary header"
        ) from exc

    stable = lines[start:]
    while stable and not stable[-1]:
        stable.pop()
    return "\n".join(stable) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(normalize(args.input.read_text(errors="replace")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
