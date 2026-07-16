# HyperIso documentation

This directory contains user-facing and developer-facing documentation for HyperIso.

## Documentation map

| Page | Purpose |
|---|---|
| `../README.md` | Public project overview and quick start. |
| `installation.md` | Detailed installation and build options. |
| `reproducibility.md` | Reproducibility package for the article and benchmark figures. |
| `development.md` | Development, CI, testing and hardening notes. |
| `../Hyperiso/Hyperiso/examples_cpp/README.md` | C++ examples. |
| `../Hyperiso/Hyperiso/examples_python/README.md` | Python examples. |
| `../Hyperiso/Hyperiso/core/src/UserInterface/README.md` | CLI guide. |
| `../GHyperiso/HyperisoDashGUI/README.md` | Dash GUI guide. |

## Build API documentation

```bash
doxygen Hyperiso/Docs/Doxyfile
```

The generated HTML output should be published by the documentation workflow once the public documentation site is enabled.

## Documentation standards

Documentation should be:

- reproducible from a clean clone;
- explicit about input files and model assumptions;
- clear about optional dependencies;
- linked to examples and tests;
- updated in the same pull request as public API changes.

## Release and publication

- [Release procedure](release.md)
- [CPC Program Summary](cpc_program_summary.md)
- [Known limitations](../KNOWN_LIMITATIONS.md)
