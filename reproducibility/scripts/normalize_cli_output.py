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
    if DIAGNOSTIC_RE.search(cleaned):
        raise ValueError("CLI output contains a warning, error or fatal diagnostic")
    if NONFINITE_RE.search(cleaned):
        raise ValueError("CLI output contains NaN or infinity")
    lines = [line.rstrip() for line in cleaned.splitlines()]
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
