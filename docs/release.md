# Release procedure

HyperIso releases are immutable and tag-driven. The release version must match
`pyproject.toml`, the CMake project and `CITATION.cff`.

## One-time repository configuration

Create two protected GitHub environments:

- `testpypi`: configure a TestPyPI Trusted Publisher for this repository and
  `.github/workflows/release.yml`;
- `pypi`: configure the equivalent PyPI Trusted Publisher and require a manual
  reviewer before deployment.

Protect `main` and require the CI, Python, pre-commit, documentation, Docker and
reproducibility workflows. Use a real GitHub account or organization team if a
new `CODEOWNERS` file is added; the previous dotted placeholder was invalid and
has been removed.

## Release-candidate checks

```bash
pre-commit run --all-files
python tools/check_version_consistency.py

cmake -S Hyperiso/Hyperiso/core -B build \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_TESTS=ON \
  -DBUILD_WITH_CLI=ON
cmake --build build -j2
ctest --test-dir build --output-on-failure

HYPERISO_BIN="$(find build -type f -name hyperiso-ui -perm -111 -print -quit)" \
  ./reproducibility/scripts/run_cli_suite.sh
```

Regenerate references only after a reviewed numerical change:

```bash
HYPERISO_BIN="$(find build -type f -name hyperiso-ui -perm -111 -print -quit)" \
  ./reproducibility/scripts/run_cli_suite.sh --update-expected
```

Review every numerical diff and `reference_metadata.json` before committing.

## Tagging and publishing

```bash
git tag -s v1.0.0 -m "HyperIso 1.0.0"
git push origin v1.0.0
```

The release workflow then:

1. reruns quality, C++ tests and all five frozen references;
2. builds one sdist and CPython 3.10–3.12 manylinux wheels;
3. checks the distributions and wheel platform metadata;
4. publishes those exact files to TestPyPI through OIDC;
5. installs and verifies the TestPyPI release;
6. promotes the same files to PyPI;
7. creates an SPDX SBOM and signed build-provenance attestations;
8. creates the GitHub release with the source archive, SBOM and SHA-256 checksums.

Never rebuild artifacts between TestPyPI and PyPI. PyPI versions are immutable;
a failed release must be corrected with a new version.
