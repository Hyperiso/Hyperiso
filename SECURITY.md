# Security Policy

## Supported versions

HyperIso is currently in active research and pre-release development. Security fixes are provided for:

| Version | Supported |
|---|---|
| `main` branch | Yes |
| Latest tagged release | Yes |
| Older pre-release tags | Best effort |

Long-term support windows will be defined once the first stable release is published.

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

## Hardening roadmap

Planned security and quality improvements include:

- CodeQL analysis in CI;
- OpenSSF Scorecard monitoring;
- Dependabot updates;
- pre-commit secret detection;
- parser fuzzing for LHA/SLHA/FLHA/YAML inputs;
- sanitizers for C++ Debug CI;
- reproducible release artifacts.
