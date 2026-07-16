# HyperIso C++ examples

This directory contains standalone C++ examples built against an installed HyperIso package. They are intended to serve both as user tutorials and as smoke tests for the public C++ API.

## Prerequisites

Build and install HyperIso first:

```bash
cmake -S Hyperiso/Hyperiso/core -B build \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_WITH_CLI=ON
cmake --build build -j
cmake --install build --prefix "$HOME/.local"
```

Make sure CMake can find the installed package:

```bash
export CMAKE_PREFIX_PATH="$HOME/.local:$CMAKE_PREFIX_PATH"
```

## Build all examples

From the repository root:

```bash
cmake -S Hyperiso/Hyperiso/examples_cpp -B build-examples \
  -G Ninja \
  -DCMAKE_PREFIX_PATH="$HOME/.local"

cmake --build build-examples -j
```

## Run a few smoke examples

```bash
./build-examples/core_base
./build-examples/core_database_writer
./build-examples/wilson_base
./build-examples/observable_base
```

Statistic examples require the installed package to export `Hyperiso::StatisticLib`:

```bash
./build-examples/statistic_base
```

## Organization

| Directory | Purpose |
|---|---|
| `Core/` | Initialization, LHA loading, parameter access, block inspection and QCD utilities. |
| `Wilson/` | Wilson-coefficient matching, running, MARTY-oriented examples and custom Wilson groups. |
| `Observable/` | Observable calculations, binned observables, decay configurations and parameter scans. |
| `Statistic/` | Uncertainty propagation, likelihood scans, maximum-likelihood fits and confidence contours. |
| `*/dev_examples/` | Advanced examples for developers extending HyperIso internals. |

## Recommended learning path

1. `Core/base_core_example.cpp`
2. `Core/parameter_provider_example.cpp`
3. `Core/database_writer_example.cpp`
4. `Wilson/base_wilson_example.cpp`
5. `Observable/base_observable_example.cpp`
6. `Statistic/base_statistic_example.cpp`
7. Advanced `dev_examples` only after the public API workflow is clear.

## Expected behavior

Each example should:

- compile without modifying the source tree;
- use input files from `Assets/lha` or a documented local path;
- print enough information to make the result understandable;
- avoid hidden global state except the documented HyperIso runtime singleton;
- be runnable in CI as a smoke test whenever dependencies are available.

## Adding a new C++ example

When adding a new example, include:

- a short comment at the top explaining the physics workflow;
- the required model/backend, for example SM-only, MARTY, THDM or SUSY;
- a deterministic input file;
- expected numerical output or tolerance when used as a regression test;
- a target in `CMakeLists.txt` using `add_hyperiso_example`.

## Notes on dynamic identifiers

Some examples demonstrate dynamic Wilson groups or custom observables. They are useful for extension development, but regular users should start with the built-in enum and mapper examples.
