# HyperIso core and Python package

This directory contains the buildable HyperIso package:

- `core/` - C++20 backend, CMake build system, CLI and tests;
- `pyhyperiso/` - Python package wrapping the C++ backend;
- `examples_cpp/` - C++ examples built against an installed HyperIso package;
- `examples_python/` - Python examples using `pyhyperiso`;
- `pyproject.toml` - Python packaging entry point based on scikit-build-core.

## Build the C++ core

From the repository root:

```bash
cmake -S Hyperiso/Hyperiso/core -B build \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_WITH_APP=ON

cmake --build build -j
cmake --install build --prefix "$HOME/.local"
```

Useful options:

| Option | Default | Purpose |
|---|---:|---|
| `ENABLE_TESTS` | `OFF` | Build CTest tests. |
| `BUILD_WITH_APP` | `OFF` | Build the `hyperiso-ui` command-line executable. |
| `BUILD_WITH_PYTHON` | `OFF` | Build the Python extension from the CMake tree. |
| `BUILD_WITH_MARTY` | `OFF` | Enable MARTY integration. |
| `BUILD_WITH_2HDMC` | `OFF` | Enable 2HDMC support. |
| `BUILD_WITH_SOFTSUSY` | `OFF` | Enable SoftSusy support. |
| `BUILD_WITH_MINUIT2` | `OFF` | Enable Minuit2-related backends. |
| `ENABLE_CLANG_TIDY` | `OFF` | Run clang-tidy during the build. |
| `ENABLE_ADRESS_SANITIZER` | `OFF` | Enable AddressSanitizer. |
| `ENABLE_UNDEFINED_SANITIZER` | `OFF` | Enable UndefinedBehaviorSanitizer. |

## Run tests

```bash
cmake -S Hyperiso/Hyperiso/core -B build \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DENABLE_TESTS=ON

cmake --build build -j
ctest --test-dir build --output-on-failure
```

Labels can be used to select a subset:

```bash
ctest --test-dir build -L common --output-on-failure
ctest --test-dir build -L database --output-on-failure
ctest --test-dir build -L statistic --output-on-failure
```

Heavy validation/comparison tests are opt-in:

```bash
cmake -S Hyperiso/Hyperiso/core -B build-comparison \
  -DENABLE_TESTS=ON \
  -DHYPERISO_BUILD_COMPARISON_TESTS=ON
```

## Install the Python package

```bash
python -m pip install --upgrade pip build
python -m pip install ./Hyperiso/Hyperiso
```

Editable development install:

```bash
python -m pip install -e ./Hyperiso/Hyperiso
```

Run Python tests:

```bash
python -m pip install pytest
pytest Hyperiso/Hyperiso/pyhyperiso/test
```

## Build examples against an installed package

```bash
cmake -S Hyperiso/Hyperiso/examples_cpp -B build-examples \
  -G Ninja \
  -DCMAKE_PREFIX_PATH="$HOME/.local"

cmake --build build-examples -j
```

## Public interfaces

| Interface | Entry point |
|---|---|
| C++ API | headers and CMake targets exported by the C++ install. |
| Python API | `pyhyperiso.Common`, `pyhyperiso.Core`, `pyhyperiso.Wilson`, `pyhyperiso.Observable`, `pyhyperiso.Statistic`. |
| CLI | `hyperiso-ui`, built with `-DBUILD_WITH_APP=ON`. |
| Dash GUI | `GHyperiso/HyperisoDashGUI`. |

## Versioning

The C++ project, Python package and release tag should use the same semantic version. For example:

```text
CMake project version: 0.1.0
Python package version: 0.1.0
Git tag: v0.1.0
```

Update all version locations in the same pull request.
