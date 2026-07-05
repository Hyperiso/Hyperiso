# Contributing to HyperIso

Thank you for considering a contribution to HyperIso. The project combines scientific software, C++ numerical code, Python bindings and user-facing interfaces. Contributions should therefore be reproducible, tested and documented.

## Development setup

```bash
git clone https://github.com/HyperIso/HyperIso.git
cd HyperIso

sudo apt-get update
sudo apt-get install -y build-essential cmake ninja-build libgsl-dev python3-dev python3-pip
python -m pip install --upgrade pip pre-commit pytest
pre-commit install
```

## Build and test before opening a pull request

```bash
cmake -S Hyperiso/Hyperiso/core -B build \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DENABLE_TESTS=ON \
  -DBUILD_WITH_CLI=ON

cmake --build build -j
ctest --test-dir build --output-on-failure
```

Python checks:

```bash
python -m pip install -e ./Hyperiso/Hyperiso
pytest Hyperiso/Hyperiso/pyhyperiso/test
```

Run formatting and lightweight checks:

```bash
pre-commit run --all-files
```

## Pull-request expectations

A pull request should include:

- a clear description of the change;
- tests or a justification if tests are not possible;
- documentation updates for public API changes;
- numerical references for physics changes;
- deterministic inputs and outputs for validation changes;
- no generated files unless they are intentional reference data.

## Commit style

Use concise, imperative commit messages, for example:

```text
Add SM Wilson coefficient regression test
Fix LHA parser handling of multi-index entries
Document Dash GUI parameter scan workflow
```

## C++ guidelines

- Use C++20 features deliberately and consistently.
- Keep public headers minimal and documented.
- Avoid hidden global state except the documented HyperIso runtime singleton.
- Prefer explicit units in variable names and comments.
- Add CTest coverage for new public behavior.
- For numerical algorithms, document references and expected precision.

## Python guidelines

- Keep wrappers thin and explicit.
- Match the C++ API naming where it improves discoverability.
- Use type hints for public wrappers.
- Add Python tests for wrapper conversions and user-facing workflows.
- Avoid writing generated files into the source tree.

## Documentation guidelines

Update the relevant README when changing:

- installation commands;
- CMake options;
- Python import paths;
- CLI options;
- Dash GUI workflows;
- example inputs or expected outputs;
- reproducibility scripts.

## Adding a new observable

A new observable should come with:

1. physics reference and formula location;
2. input parameter dependencies;
3. Wilson coefficient dependencies;
4. central-value test or regression test;
5. uncertainty/statistical behavior if relevant;
6. C++ and/or Python example when public.

## Adding a new Wilson group or MARTY workflow

Include:

1. group and coefficient identifiers;
2. matching scale convention;
3. running behavior;
4. basis conventions;
5. validation against a known model or limit;
6. example LHA/model files when possible.

## Reporting scientific discrepancies

When reporting a physics mismatch, include:

- commit hash;
- compiler and platform;
- build options;
- input file;
- observable or coefficient;
- expected value and reference;
- actual value;
- numerical tolerance.
