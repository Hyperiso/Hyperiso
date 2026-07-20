#!/usr/bin/env python3
"""Validate public release metadata against CITATION.cff.

This check intentionally covers user-visible release strings that are not all
part of the core version-consistency check: the current Zenodo DOI, release date,
README badge/citation, CPC summary, manuscript and release documentation.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


def require(pattern: str, text: str, label: str) -> str:
    match = re.search(pattern, text, re.MULTILINE)
    if not match:
        raise RuntimeError(f"could not read {label}")
    return match.group(1)


def require_text(text: str, expected: str, label: str, errors: list[str]) -> None:
    if expected not in text:
        errors.append(f"{label} does not contain {expected!r}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    args = parser.parse_args()
    root = args.root.resolve()

    cff = (root / "CITATION.cff").read_text(encoding="utf-8")
    version = require(r'^version:\s*["\']?([^"\'\s]+)', cff, "release version")
    release_date = require(r'^date-released:\s*["\']?([^"\'\s]+)', cff, "release date")
    doi = require(r'^doi:\s*["\']?([^"\'\s]+)', cff, "Zenodo DOI")

    errors: list[str] = []
    if not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", version):
        errors.append(f"invalid semantic version in CITATION.cff: {version!r}")
    if not re.fullmatch(r"[0-9]{4}-[0-9]{2}-[0-9]{2}", release_date):
        errors.append(f"invalid release date in CITATION.cff: {release_date!r}")
    if not re.fullmatch(r"10\.5281/zenodo\.[0-9]+", doi):
        errors.append(f"invalid version-specific Zenodo DOI: {doi!r}")

    readme = (root / "README.md").read_text(encoding="utf-8")
    require_text(readme, f"Zenodo v{version}", "README Zenodo badge", errors)
    require_text(readme, f"https://doi.org/{doi}", "README DOI link", errors)
    require_text(
        readme,
        f"Archived software v{version}: https://doi.org/{doi}",
        "README current citation",
        errors,
    )

    cpc = (root / "docs/cpc_program_summary.md").read_text(encoding="utf-8")
    require_text(cpc, f"**Version:** {version}", "CPC version", errors)
    require_text(cpc, f"tree/v{version}", "CPC release tag", errors)
    require_text(cpc, f"https://doi.org/{doi}", "CPC DOI", errors)

    limitations = (root / "KNOWN_LIMITATIONS.md").read_text(encoding="utf-8")
    require_text(
        limitations,
        f"# Known limitations — HyperIso {version}",
        "known-limitations heading",
        errors,
    )

    manuscript = (root / "docs/main.tex").read_text(encoding="utf-8")
    require_text(
        manuscript,
        f"\\texttt{{v{version}}}",
        "manuscript release tag",
        errors,
    )
    require_text(
        manuscript,
        f"\\texttt{{pyhyperiso=={version}}}",
        "manuscript package version",
        errors,
    )
    require_text(
        manuscript,
        f"\\url{{https://doi.org/{doi}}}",
        "manuscript DOI",
        errors,
    )

    release_doc = (root / "docs/release.md").read_text(encoding="utf-8")
    require_text(release_doc, f"v{version}", "release guide tag", errors)
    require_text(release_doc, doi, "release guide DOI", errors)

    changelog = (root / "CHANGELOG.md").read_text(encoding="utf-8")
    require_text(
        changelog,
        f"## [{version}] - {release_date}",
        "changelog release heading",
        errors,
    )

    if errors:
        print("Release metadata consistency failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1

    print(
        "Release metadata consistency OK: "
        f"version={version}, date={release_date}, doi={doi}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
