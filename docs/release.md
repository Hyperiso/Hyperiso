
## Publication order

HyperIso uses a tag-driven immutable release workflow. First merge the reviewed
release commit and verify all branch checks. Then create and push the signed
`v1.0.2` tag. The tag builds the sdist and wheels once, publishes those exact
artefacts to TestPyPI, installs and verifies the TestPyPI wheel, and pauses at
the protected `pypi` environment for maintainer approval before publishing the
same files to PyPI. The GitHub Release is created only after PyPI succeeds.

Do not upload a separately rebuilt package before the tag: PyPI versions cannot
be overwritten, and rebuilding between TestPyPI and PyPI would break artefact
identity. A pre-tag release candidate may be tested with local wheels or a
distinct development version such as `1.0.2rc1`, but the final `1.0.2` files are
produced from the immutable tag.

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

## Reserve the Zenodo v1.0.2 DOI

Open the published v1.0.1 record at
`https://zenodo.org/records/21414482`, choose **New version**, update the version
to `1.0.2`, and reserve the new version-specific DOI before tagging. Insert that
DOI into `CITATION.cff`, the README/software paper and any release notes that
will be archived. Do not reuse `10.5281/zenodo.21414482`, which identifies only
v1.0.1. Keep the Zenodo draft unpublished until the final tag archive has been
created and checked.

## Freeze final numerical provenance

After all code changes are committed and the Release build has passed locally, regenerate the frozen references with the exact CLI that will be released:

```bash
GITHUB_REF_NAME=v1.0.2 \
HYPERISO_BIN=/absolute/path/to/hyperiso-ui \
CMAKE_BUILD_TYPE=Release \
./reproducibility/scripts/run_cli_suite.sh --update-expected
```

Review and commit `reproducibility/expected_outputs/` and `reference_metadata.json`. The numerical values should not change at this stage; the purpose is to record the real Git checkout, binary hash, compiler, dependency versions and per-case timings. Build the binary from the recorded commit, then make one final commit containing only the frozen files under `reproducibility/expected_outputs/`. The tag workflow rejects any other code change after the recorded build commit. Do not tag while provenance fields are incomplete.

Validate the freeze before tagging:

```bash
python tools/check_release_provenance.py --tag v1.0.2
```

## Tagging and publishing

```bash
git tag -s v1.0.2 -m "HyperIso 1.0.2"
git push origin v1.0.2
```

The release workflow then:

1. reruns quality, C++ tests and all seven frozen references;
2. builds one sdist and CPython 3.10–3.12 manylinux wheels;
3. checks the distributions and wheel platform metadata;
4. publishes those exact files to TestPyPI through OIDC;
5. installs and verifies the TestPyPI release;
6. promotes the same files to PyPI;
7. creates an SPDX SBOM and SHA-256 checksums;
8. creates the GitHub release with the source archive, SBOM and checksums.

After the tag workflow succeeds, create an exact source archive with
`git archive v1.0.2`, upload it to the reserved Zenodo new-version draft, verify
the archive checksum and metadata, and publish the record. Then add the final
version-specific DOI to the GitHub Release and project website if it was not
already present.

Never rebuild artifacts between TestPyPI and PyPI. PyPI versions are immutable;
a failed release must be corrected with a new version.
