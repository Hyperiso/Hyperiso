# Architecture and data flow {#architecture}

HyperIso uses a modular, ports-and-adapters style architecture. Physics logic,
input loading, numerical utilities and user interfaces are separated so that
they can be validated independently while sharing one runtime state.

## Module map

| Module | Responsibility |
|---|---|
| `Common` | Shared identifiers, enums, mappers, logging and low-level types. |
| `DataBase` | Typed blocks, JSON/YAML/LHA readers, correlations and database export. |
| `Core` | Runtime configuration, parameter memory, initialization and QCD/flavour inputs. |
| `PhysicalModel` | Model-dependent Wilson matching and renormalization-group running. |
| `BusinessLogic` | Decay calculators, observable projections and numerical integrations. |
| `Statistic` | Nuisance distributions, uncertainty propagation, likelihoods, fits and contours. |
| `Math` | Special functions, interpolation, integration and numerical helpers. |
| `ExternalIntegration` | Optional MARTY, 2HDMC and SOFTSUSY adapters. |
| `UserInterface` | The `hyperiso-ui` command-line application. |

\image html hyperiso_data_flow.png "HyperIso modules and information flow"

## One backend, four frontends

The C++, Python, CLI and Dash layers all initialize the same core through
@ref HyperisoMaster. Consequently, the following rules do not depend on the
selected frontend:

1. default database loading;
2. YAML and LHA-family override precedence;
3. model and external-backend selection;
4. Wilson matching and running;
5. observable registration and calculation;
6. statistical nuisance and correlation handling.

This shared path is important for scientific comparison: a command reproduced
through the CLI should exercise the same core calculations as the equivalent
Python or C++ program.

## Initialization sequence

A typical initialization performs the following operations:

1. create the concrete database and correlation loaders;
2. resolve distributed asset paths and optional user path overrides;
3. load the default JSON parameter blocks;
4. apply user YAML overrides;
5. parse the LHA, SLHA or FLHA input;
6. initialize or consume the requested model spectrum;
7. cache the resulting parameter state for Wilson and observable calculations.

@ref HyperisoMaster also exposes pre-initialization hooks for additional block
prototypes, writable caches, custom database paths and external installations.
These hooks must be configured before `init()` whenever they need to affect the
first load.

## Wilson and observable pipeline

The `PhysicalModel` layer produces coefficient groups at the matching scale and
runs them to the scale requested by an observable. The `BusinessLogic` layer
then combines these coefficients with the shared flavour and QCD inputs to
calculate decay amplitudes and user-facing observables.

The origin of the Wilson coefficients is deliberately hidden from the
observable layer. A coefficient can come from a native SM/THDM/SUSY
implementation or from an optional MARTY-generated implementation; the
observable receives the same internal coefficient representation.

## Optional external integrations

### MARTY

MARTY is used for generic-model Wilson matching at the perturbative order
implemented by that workflow. It is disabled by default. A runtime installation
can be registered before initialization:

```cpp
HyperisoMaster master;
master.pre_init_set_marty_path("/path/to/MARTY_INSTALL");
```

The equivalent Python method is available on `HyperisoMaster`. HyperIso validates
the include directory and `libmarty` before entering MARTY mode.

### Spectrum calculators

2HDMC and SOFTSUSY support spectrum-oriented workflows. Release-level
reproducibility examples consume archived spectrum files where possible, so the
external generator does not need to be rerun to check the frozen numerical
outputs.

## Design implications for contributors

- Physics implementations should depend on interfaces rather than user-facing
  frontends.
- New public C++ types exposed to Python require synchronized pybind11 and
  high-level wrapper changes.
- Input schema changes require database, reader, example and reproducibility
  updates in the same pull request.
- A new observable should be tested at the decay-calculator level and through
  at least one public interface.

See @ref extending for a practical contribution map.
