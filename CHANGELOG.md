# Changelog

All notable changes to HyperIso are documented in this file. The project follows
Semantic Versioning.

## [Unreleased]

### Roadmap

- Additional platform wheels after Linux support is validated.
- Expanded fuzzing and long-running numerical validation.

## [1.0.1] - 2026-07-17

### Fixed

- Resolve statistical nuisance files through the active Core runtime path
  provider, so installed wheels honour packaged assets and pre-initialization
  path overrides instead of falling back to build-time paths.
- Add C++ adapter and installed-wheel regression tests for statistical
  initialization.

## [1.0.0] - 2026-07-16

### Added

- Public C++20 core, command-line interface and Python bindings.
- Public `DatabaseWriter`/`FileWriter` API for exporting the initialized Core
  database to JSON, YAML, LHA, SLHA and FLHA, with block and parameter filters.
- Standard Model, THDM, SUSY and MARTY-oriented Wilson-coefficient workflows.
- Flavour-observable and statistical-analysis interfaces.
- Fit-parameter sensitivity validation that blocks ill-posed fits and contours when
  selected observables are numerically insensitive to a requested parameter.
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
- Python bindings now select the const `MarginalConfigFactory::create`
  overloads correctly with pybind11.
- JSON and YAML database exports preserve complex parameter components through
  the `imaginary_value` field; LHA-family exports generate `IM...` companion
  blocks when complex values are stored directly.
- The Python extension now defines and links `BlockName::canonical()`, preventing
  import-time undefined-symbol failures in installed wheels.
- GHyperiso confidence contours use filled semi-transparent regions, darker
  boundaries and a white publication-style canvas.
