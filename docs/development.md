# Development and hardening notes

## CI strategy

HyperIso uses several GitHub Actions workflows:

| Workflow | Purpose |
|---|---|
| `ci.yml` | C++ configure, build, tests and install smoke checks. |
| `python.yml` | Python package build and tests. |
| `docs.yml` | Doxygen documentation build. |
| `docker.yml` | Container build smoke tests. |
| `codeql.yml` | Static security analysis. |
| `scorecard.yml` | OpenSSF Scorecard. |
| `release.yml` | Release artifacts on version tags. |

## Local quality checks

```bash
pre-commit run --all-files
```

## Sanitizers

Recommended local sanitizer build:

```bash
cmake -S Hyperiso/Hyperiso/core -B build-asan \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DENABLE_TESTS=ON \
  -DENABLE_ADDRESS_SANITIZER=ON \
  -DENABLE_UNDEFINED_SANITIZER=ON

cmake --build build-asan -j
ctest --test-dir build-asan --output-on-failure
```

## Fuzzing and hardening

Planned fuzzing targets:

- LHA parser;
- SLHA/FLHA parser;
- YAML override loader;
- JSON database loader;
- Wilson coefficient request construction;
- observable ID and mapper conversions.

These targets are roadmap items and are not release claims.

## Generated files

Avoid committing generated outputs unless they are intentional regression fixtures. Generated fixtures must live under a documented `testdata` or `expected` directory and include a note explaining how they were produced.
