# C++, Python, CLI and graphical interfaces {#interfaces}

All frontends use the same C++ calculation backend. Choose the interface based
on the surrounding workflow rather than on expected physics differences.

## C++ interface

Use C++ when integrating HyperIso into a compiled program, implementing new
physics components or requiring direct access to the native API.

The common sequence is:

1. create a `HyperisoConfig`;
2. initialize @ref HyperisoMaster with an LHA-family file;
3. create @ref WilsonInterface, @ref ObservableInterface or
   @ref StatisticInterface;
4. register the requested objects and run the calculation.

```cpp
HyperisoConfig config;
config.model = Model::SM;

HyperisoMaster master;
master.init("Assets/lha/si_input.flha", config);

ObservableInterface observables;
observables.add_observable(Observables::BR_BS_MUMU, QCDOrder::NNLO);
auto value = observables.compute_observable(Observables::BR_BS_MUMU);
```

Installed C++ examples are built against the exported HyperIso CMake package,
which also checks that public headers and link interfaces are complete.

## Python interface

Use Python for scripted scans, notebooks, plotting and integration with the
scientific Python ecosystem. High-level Python wrappers validate and convert
configuration objects before forwarding them to pybind11.

```python
from pyhyperiso.Common import Model, Observables, QCDOrder
from pyhyperiso.Core import HyperisoConfig, HyperisoMaster
from pyhyperiso.Observable import ObservableInterface

master = HyperisoMaster()
master.init(
    lha_file="Assets/lha/si_input.flha",
    config=HyperisoConfig(model=Model.SM),
)

obs = ObservableInterface()
obs.add_observable(Observables.BR_BS_MUMU, QCDOrder.NNLO)
print(obs.compute_observable(Observables.BR_BS_MUMU))
```

Create all configuration values before constructing the corresponding
interface. The C++ object receives a converted copy of the Python configuration.

## Command-line interface

Use `hyperiso-ui` for release regression tests, shell scripts and quick checks.
The command groups expose Wilson, observable and statistical summaries without
requiring a dedicated program.

```bash
hyperiso-ui observable summary \
  --model SM \
  --lha Assets/lha/si_input.flha \
  --observables BR_Bs__mu_mu \
  --order NNLO
```

Run `hyperiso-ui <group> <command> --help` for the complete option list. The CLI
is the canonical interface for the frozen release suite because its text output
can be normalized and compared directly.

## Dash graphical interface

Use the Dash application for interactive exploration. It provides pages for
configuration, Wilson coefficients, observables and statistical fits while
calling the same Python/C++ backend.

The graphical statistical page validates selected fit parameters before
running a fit or contour. Parameters to which all selected observables are
numerically insensitive are rejected with a user-facing explanation rather
than being sent to the minimizer as a flat direction.

## Optional backends

Standard calculations do not require MARTY. To use generic-model matching,
register a validated MARTY installation before `init()` and explicitly enable
the MARTY mode in the configuration. The same rule applies in C++, Python and
the GUI.

## Interface consistency checklist

When comparing two interfaces, keep all of the following equal:

- software commit or release tag;
- model, contribution and QCD order;
- LHA-family input and YAML overrides;
- matching and low-energy scales;
- selected observables and bins;
- statistical mode, random seed, draw count and thread count.
