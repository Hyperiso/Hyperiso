# Security Policy

## Supported versions

HyperIso is active research software. Security fixes are provided for:

| Version | Supported |
|---|---|
| `main` branch | Yes |
| Latest tagged release | Yes |
| Older releases | Best effort |

Version 1.0.0 is the first stable release. Long-term support windows are not guaranteed.

## Reporting a vulnerability

Please do not open a public GitHub issue for vulnerabilities.

Report security issues by email to the maintainers listed in the project metadata. Include:

- affected version or commit;
- operating system and compiler/Python version;
- reproduction steps;
- impact assessment;
- whether the issue is public or under embargo.

We aim to acknowledge valid reports within 7 days and to provide a remediation plan or status update within 30 days.

## Scope

Security-sensitive areas include:

- LHA/SLHA/FLHA parsing;
- YAML/JSON input loading;
- generated code and MARTY integration;
- file-system access in examples, GUI and scripts;
- Python bindings that cross the C++/Python boundary;
- Docker images and CI artifacts.

## Implemented release hardening

The repository uses CodeQL, OpenSSF Scorecard, Dependabot, C++ sanitizer
configuration, strict pre-commit checks, frozen numerical references, SBOM
generation and provenance attestations in the release workflow.

Parser fuzzing and broader long-running numerical campaigns remain roadmap work.
