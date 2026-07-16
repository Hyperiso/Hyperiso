# Extending HyperIso {#extending}

This page maps common contributions to the relevant source areas and validation
steps. Read @ref architecture first; new functionality should enter through the
backend rather than being implemented independently in a frontend.

## Add or modify an observable

Relevant areas:

- `BusinessLogic/domain`: decay and observable implementations;
- `BusinessLogic/adapters`: the public @ref ObservableInterface;
- shared enums and identifiers in `Common`;
- Python bindings and wrappers when the observable is public;
- C++ and Python examples for user-facing behavior.

A contribution should include unit tests for the numerical implementation and
an integration test through the public observable interface. If a published
reference changes, update the reproducibility suite only after reviewing the
numerical difference.

## Add a Wilson coefficient or group

Relevant areas:

- `PhysicalModel/domain`: coefficient definitions and group containers;
- model-specific matching and running implementations;
- @ref WilsonInterface and enum mappings;
- MARTY templates if the generic-model route supports the same group.

Document the supported perturbative order and the matching/running scales.
Never silently advertise a higher order than the implemented backend provides.

## Add an input block

Add the typed block definition and default data, then update:

- JSON schema/data assets;
- LHA/YAML readers and writers where relevant;
- block and identifier mappings;
- parameter-provider tests;
- input documentation and examples.

Use @ref HyperisoMaster::pre_init_add_block for user-defined runtime blocks that
do not belong in the distributed default database.

## Change statistical behavior

Statistical changes normally affect:

- `Statistic/domain` for likelihood, covariance and fit logic;
- `Statistic/adapters` for the public API;
- pybind11 bindings and the high-level Python `StatisticConfig`;
- the Dash service/page for graphical behavior;
- deterministic tests and fixed-seed regressions.

Keep common options in `StatisticConfig` and expert controls in
`AdvancedStatisticConfig`. Validate fit-parameter and nuisance behavior before
calling the minimizer.

## Keep interfaces synchronized

For every public C++ configuration field or enum:

1. expose it in pybind11;
2. mirror it in the Python wrapper;
3. add conversion tests;
4. update GUI controls if users can configure it graphically;
5. document defaults and failure behavior.

## Build and test locally

```bash
pre-commit run --all-files

cmake -S Hyperiso/Hyperiso/core -B build \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DENABLE_TESTS=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure

python -m pytest Hyperiso/Hyperiso/pyhyperiso/test
```

Build installed examples as well; they detect missing public headers, symbols
and transitive link dependencies that an in-tree test can hide.

## Documentation requirements

Update this Doxygen guide in the same pull request as public behavior. Public
classes should have a concise `@brief`, parameter and return documentation,
failure behavior and a minimal example where practical.

The documentation workflow builds with warnings visible but does not currently
fail on all historical warnings. New code should nevertheless avoid introducing
new unresolved references, undocumented parameters or recursive `@copydoc`
links.
