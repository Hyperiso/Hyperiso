# HyperIso release checklist

This is the reusable release playbook for patch, minor and major releases. It is
written to be copied for every release and uses placeholders instead of a fixed
version.

The release is **tag-driven and immutable**:

1. prepare and test the release on a branch;
2. reserve the version-specific Zenodo DOI;
3. merge code and public metadata into `main`;
4. rebuild from the merged commit and freeze numerical provenance in a second PR;
5. tag the final `main` commit;
6. let GitHub Actions publish the same Python artifacts to TestPyPI and PyPI;
7. publish the exact GitHub Release source archive in the reserved Zenodo draft.

Never upload a locally rebuilt artifact after the tag. PyPI versions and
version-specific Zenodo records are immutable publication units.

## 1. Define the release variables

Choose the new version, previous version, release date and DOI before editing:

```bash
export VERSION="X.Y.Z"
export TAG="v${VERSION}"
export PREVIOUS_VERSION="A.B.C"
export PREVIOUS_TAG="v${PREVIOUS_VERSION}"
export RELEASE_DATE="YYYY-MM-DD"
export ZENODO_DOI="10.5281/zenodo.NNNNNNNN"
export PREVIOUS_ZENODO_DOI="10.5281/zenodo.MMMMMMMM"
```

Use semantic versioning:

- patch: compatible bug fixes, for example `1.0.2 -> 1.0.3`;
- minor: backward-compatible features, for example `1.0.3 -> 1.1.0`;
- major: incompatible public API or data-format changes.

Create a release branch from the latest protected `main`:

```bash
git switch main
git pull --ff-only origin main
git switch -c "release/${TAG}"
```

A focused bug-fix branch such as `fix/marty-v1.0.3` is also acceptable, but the
final release candidate must be rebased on the latest `origin/main`.

## 2. Reserve the Zenodo DOI

HyperIso embeds the DOI before tagging, so use a **manual Zenodo new-version
draft** rather than relying only on automatic GitHub-to-Zenodo archiving.

1. Open the previous published version-specific Zenodo record.
2. Choose **New version**.
3. Import or retain the previous metadata as appropriate.
4. Set the new software version and planned publication date.
5. Reserve/Get a DOI in the draft.
6. Keep the draft unpublished until the GitHub Release assets exist.
7. Record the new version DOI in `ZENODO_DOI`.

Do not confuse:

- the **version DOI**, which identifies one immutable release and belongs in
  `CITATION.cff` and the release badge;
- the **concept DOI**, which resolves to the latest version and can be shown
  separately if desired.

If Zenodo's GitHub integration is enabled, avoid publishing a second automatic
record for the same tag when using this manual reserved-DOI workflow.

## 3. Update every version declaration

The following files are part of the automated version-consistency contract:

| File | Field or text |
|---|---|
| `Hyperiso/Hyperiso/pyproject.toml` | `[project].version` |
| `Hyperiso/Hyperiso/core/CMakeLists.txt` | `project(Hyperiso VERSION ...)` |
| `Hyperiso/Hyperiso/core/src/Python/CMakeLists.txt` | Python CMake fallback |
| `Hyperiso/Hyperiso/pyhyperiso/__init__.py` | fallback `__version__` |
| `Hyperiso/Docs/Doxyfile` | `PROJECT_NUMBER` |
| `Hyperiso/Docs/pages/mainpage.md` | displayed version and install command |
| `Hyperiso/Hyperiso/README.md` | CMake, Python and Git-tag versions |
| `docs/cpc_program_summary.md` | version and repository tag |
| `docs/installation.md` | installed-version assertion |
| `reproducibility/manifest.json` | `hyperiso_version` |
| `CITATION.cff` | `version` and `date-released` |

Also review these user-facing files even when a parser does not extract a
version from every sentence:

| File | What to update |
|---|---|
| `CHANGELOG.md` | release heading, date and complete change list |
| `README.md` | current release badge, installation/citation text |
| `KNOWN_LIMITATIONS.md` | current release heading and scope |
| `docs/main.tex` | source tag, PyPI version, DOI and release statements |
| `docs/release.md` | current release-specific procedure |
| `.github/workflows/*.yml` | version-specific test names only when applicable |
| `tools/test_marty_v*.py` | release-specific regression expectations |

Historical statements such as “the FLHA migration was introduced in 1.0.2”
should remain historical. Do not blindly replace every old version string.

Run the core check immediately after editing:

```bash
python tools/check_version_consistency.py
```

## 4. Update DOI and public release metadata

Set the new version DOI in:

| File | Required update |
|---|---|
| `CITATION.cff` | `doi`, `version`, `date-released` |
| `README.md` | Zenodo badge and current archived-software citation |
| `docs/cpc_program_summary.md` | software archive DOI |
| `docs/main.tex` | program summary and software-availability DOI |
| `docs/release.md` | previous record, current reserved DOI and upload steps |

Keep the previous version DOI in a clearly labelled “Previous release” entry.

Run:

```bash
python tools/check_release_metadata_consistency.py
python tools/validate_metadata.py
```

Search for stale current-release strings:

```bash
grep -RIn \
  --exclude-dir=.git \
  --exclude-dir=Third_party \
  -E "${PREVIOUS_TAG}|${PREVIOUS_VERSION}|${PREVIOUS_ZENODO_DOI}" .
```

Classify every result:

- keep genuine history, migrations and old changelog sections;
- replace badges, current installation commands, current DOI, current tag,
  current limitations and release instructions.

## 5. Run repository-level static checks

```bash
git diff --check
python tools/check_version_consistency.py
python tools/check_release_metadata_consistency.py
python tools/validate_metadata.py
python tools/check_observable_flha_consistency.py
python tools/check_marty_template_sync.py
python tools/test_marty_v103.py --static-only
python -m unittest discover -s reproducibility/tests -v
pre-commit run --all-files
```

For a future release, rename or replace a version-specific regression script
such as `test_marty_v103.py` only when its assertions no longer represent the
current release.

## 6. Build and test the C++ release candidate

Use a clean Release build:

```bash
rm -rf build-release
cmake -S Hyperiso/Hyperiso/core -B build-release \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_TESTS=ON \
  -DBUILD_WITH_CLI=ON \
  -DHYPERISO_BUILD_COMPARISON_TESTS=ON
cmake --build build-release --parallel 2
ctest --test-dir build-release --output-on-failure
```

Check the CLI version:

```bash
HYPERISO_BIN="$(find build-release -type f -name hyperiso-ui -perm -111 -print -quit)"
test -n "${HYPERISO_BIN}"
"${HYPERISO_BIN}" --version
"${HYPERISO_BIN}" --help >/dev/null
```

The reported version must equal `VERSION`.

## 7. Test Python and package artifacts

Install development and test dependencies:

```bash
python -m pip install --upgrade pip
python -m pip install -e "./Hyperiso/Hyperiso[dev,test]"
```

Run the Python suite:

```bash
python -m pytest Hyperiso/Hyperiso/pyhyperiso/test -q --import-mode=importlib
python -m unittest discover -s reproducibility/tests -v
```

Verify the three Python version sources:

```bash
python - <<'PY'
import importlib.metadata
import pyhyperiso
from pyhyperiso.phyperiso import pyhyperiso as native

expected = importlib.metadata.version("pyhyperiso")
assert pyhyperiso.__version__ == expected
assert native.__version__ == expected
print("Python versions agree:", expected)
PY
```

Build and inspect packages:

```bash
rm -rf dist
python -m build Hyperiso/Hyperiso --outdir dist
python -m twine check --strict dist/*
```

Install the wheel in a clean virtual environment and run the same import/version
smoke test from outside the repository.

## 8. Build documentation

```bash
cd Hyperiso/Docs
rm -rf html
doxygen Doxyfile
test -s html/index.html
test "$(find html -type f | wc -l)" -gt 20
cd ../..
```

If the manuscript is part of the release, compile it with the project's normal
LaTeX workflow and inspect DOI/version rendering in the PDF.

## 9. Run model-specific and concurrency checks

For MARTY-affecting releases, test at least:

- THDM C1–C10 plus scalar and primed coefficients;
- the THDM C9 regression benchmark;
- Z-prime C9 and its large-mass decoupling;
- one repeated identical run to prove cache reuse;
- the Python multithreading test in separate processes;
- a real BSM TreeLevel model, such as a leptoquark, whenever available.

Record the exact commands and representative outputs in the pull request.

## 10. First pull request: code and public metadata

Before pushing:

```bash
git fetch origin
git rebase origin/main
git status --short
git diff --check
```

Commit and push:

```bash
git add -A
git commit -S -m "chore(release): prepare ${TAG}"
git push -u origin HEAD
```

Open a PR to `main`, require review, and wait for all required checks. Merge only
when code, documentation, package, reproducibility and security workflows are
green.

## 11. Second pull request: freeze final numerical provenance

Do this **after** the first PR is on `main`, because the frozen metadata must
record a real commit already contained in `main`.

```bash
git switch main
git pull --ff-only origin main
git switch -c "release/${TAG}-provenance"
```

Rebuild the Release CLI from this exact commit, run `ctest`, then regenerate:

```bash
HYPERISO_BIN="$(find build-release -type f -name hyperiso-ui -perm -111 -print -quit)"
GITHUB_REF_NAME="${TAG}" \
HYPERISO_BIN="${HYPERISO_BIN}" \
CMAKE_BUILD_TYPE=Release \
./reproducibility/scripts/run_cli_suite.sh --update-expected
```

Review every diff:

```bash
git status --short
git diff -- reproducibility/expected_outputs
```

Only files under `reproducibility/expected_outputs/` should change in this PR.
Check:

```bash
python tools/check_release_provenance.py --tag "${TAG}"
```

Commit and merge:

```bash
git add reproducibility/expected_outputs
git commit -S -m "chore(release): freeze ${TAG} provenance"
git push -u origin HEAD
```

## 12. Final verification on `main`

```bash
git switch main
git pull --ff-only origin main
test -z "$(git status --porcelain)"
test "$(git rev-parse HEAD)" = "$(git rev-parse origin/main)"

python tools/check_version_consistency.py --tag "${TAG}"
python tools/check_release_metadata_consistency.py
python tools/validate_metadata.py
python tools/check_release_provenance.py --tag "${TAG}"
pre-commit run --all-files
```

Repeat the clean Release build, C++ tests, reproducibility suite, Python package
build, documentation build and critical MARTY smoke tests.

Verify that the tag does not already exist:

```bash
test -z "$(git tag --list "${TAG}")"
test -z "$(git ls-remote --tags origin "refs/tags/${TAG}")"
```

## 13. Create and push the signed tag

```bash
git tag -s "${TAG}" -m "HyperIso ${VERSION}"
git tag -v "${TAG}"
git show --stat "${TAG}"
git push origin "${TAG}"
```

Do not move or replace a published release tag. A failed immutable release must
be fixed with a new version.

## 14. Monitor GitHub Actions and PyPI

The tag workflow should:

1. validate version, metadata, provenance and pre-commit;
2. build documentation;
3. build and test C++ and the frozen references;
4. build one sdist and the supported manylinux wheels;
5. assemble the source archive, SBOM, manifest and checksums;
6. publish the exact Python files to TestPyPI using OIDC;
7. install and verify the TestPyPI wheel;
8. pause for the protected PyPI environment approval;
9. publish the same artifacts to PyPI;
10. create the GitHub Release.

Before approving PyPI, inspect the TestPyPI project, filenames, wheel platforms
and installed version.

Trusted Publisher settings must match the GitHub owner, repository,
`.github/workflows/release.yml` and the `testpypi`/`pypi` environment names.

## 15. Publish the reserved Zenodo draft

After the GitHub Release succeeds:

1. download `hyperiso-${TAG}.tar.gz` and `SHA256SUMS` from the GitHub Release;
2. verify the archive against `SHA256SUMS`;
3. upload that exact archive to the reserved Zenodo draft;
4. set resource type to Software and license to `GPL-3.0-or-later`;
5. verify version, release date, creators, description, related identifiers and
   the reserved DOI;
6. preview and publish the record;
7. verify that the DOI resolves and that the Zenodo file checksum matches the
   GitHub Release asset.

Do not create another source archive locally for Zenodo.

## 16. Post-release verification

Install from PyPI in a new environment:

```bash
python -m venv /tmp/hyperiso-release-test
source /tmp/hyperiso-release-test/bin/activate
python -m pip install --upgrade pip
python -m pip install "pyhyperiso==${VERSION}"
```

Verify the installed package version and run a small calculation.

Then check:

- GitHub Release assets and checksums;
- PyPI and TestPyPI metadata;
- Zenodo DOI resolution and files;
- documentation website;
- README badges;
- release notes and changelog;
- no uncommitted release edits remain.

Finally create the next `Unreleased` section and, if desired, bump the working
development version according to the project's versioning policy.

## 17. Release record template

Keep a short signed record in the release PR or issue:

```text
Version:
Tag:
Release date:
Source commit:
Zenodo version DOI:
Previous version DOI:
C++ tests:
Python tests:
Documentation build:
Reproducibility provenance:
MARTY THDM/Z-prime tests:
Multithreading test:
TestPyPI verification:
PyPI publication:
GitHub Release:
Zenodo publication:
Maintainer:
```
