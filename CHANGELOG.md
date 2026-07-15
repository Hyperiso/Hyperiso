# Changelog

All notable changes to HyperIso are documented in this file. The project follows
Semantic Versioning.

## [Unreleased]

### Roadmap

- Additional platform wheels after Linux support is validated.
- Expanded fuzzing and long-running numerical validation.

## [1.0.0] - 2026-07-15

### Added

- Public C++20 core, command-line interface and Python bindings.
- Standard Model, THDM, SUSY and MARTY-oriented Wilson-coefficient workflows.
- Flavour-observable and statistical-analysis interfaces.
- Fixed-seed CPC reproducibility suite with seven frozen reference cases, including archived THDM and SUSY spectra.
- Linux wheels built with cibuildwheel and Trusted Publishing workflows.
- Third-party licensing and provenance notices.

### Changed

- `mu_c_bsg` is sampled from its required flat interval `[0.9, 4.0]`.
- Release builds use optimized CMake `Release` settings.
- Python support is declared for Python 3.10–3.12, matching the syntax used by
  the package.
- Reproducibility outputs exclude machine-specific startup paths.

### Fixed

- Flat nuisance bounds are now applied exactly rather than reconstructed from a
  Gaussian standard deviation.
- Version metadata is synchronized at `1.0.0`.
- Monte-Carlo CSV output is standards-compliant and deterministic.
