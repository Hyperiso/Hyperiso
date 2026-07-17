# HyperIso

<div class="hyperiso-hero">
  <div class="hyperiso-kicker">HyperIso 1.0.2</div>
  <h2>Flavour observables, Wilson coefficients and statistical inference through one shared C++ backend</h2>
  <p>
    HyperIso is a modular flavour-physics framework for Standard Model and
    beyond-the-Standard-Model calculations. The same initialization,
    parameter database, Wilson-coefficient pipeline and observable calculators
    are exposed through C++, Python, the <code>hyperiso-ui</code> command-line
    interface and the Dash graphical interface.
  </p>
</div>

<div class="hyperiso-grid">
  <div class="hyperiso-card"><b>@ref getting_started "Getting started"</b><br/>Build, install and run a first calculation.</div>
  <div class="hyperiso-card"><b>@ref architecture "Architecture"</b><br/>Modules, data flow and optional backends.</div>
  <div class="hyperiso-card"><b>@ref inputs "Inputs and parameters"</b><br/>JSON defaults, YAML overrides and LHA-family files.</div>
  <div class="hyperiso-card"><b>@ref interfaces "User interfaces"</b><br/>C++, Python, CLI and Dash workflows.</div>
  <div class="hyperiso-card"><b>@ref statistics "Statistics"</b><br/>Uncertainties, chi-square fits and confidence contours.</div>
  <div class="hyperiso-card"><b>@ref reproducibility "Reproducibility"</b><br/>Frozen release references and provenance.</div>
  <div class="hyperiso-card"><b>@ref extending "Extending HyperIso"</b><br/>Where to add physics, bindings, tests and documentation.</div>
  <div class="hyperiso-card"><b><a href="annotated.html">C++ API index</a></b><br/>Classes, structures, files and call graphs.</div>
</div>

## Scientific scope

HyperIso evolves the SuperIso calculation model into a layered and extensible
software architecture. It supports:

- Wilson-coefficient matching and running for the implemented coefficient groups;
- native Standard Model, two-Higgs-doublet and supersymmetric workflows;
- optional leading-order Wilson coefficients generated for generic BSM models through MARTY;
- flavour-observable predictions in the implemented B, D and K decay families;
- Monte-Carlo uncertainty propagation, correlated experimental inputs, maximum-likelihood fits, chi-square fits and two-dimensional confidence contours;
- LHA, SLHA and FLHA input together with reproducible JSON defaults and YAML overrides.

The exact observable and coefficient catalogue is available from the generated
class and enum indices. Features outside the supported release scope are listed
in the repository's `KNOWN_LIMITATIONS.md` file.

## Data flow

The diagram below summarizes how numerical inputs enter the runtime database,
how model-dependent short-distance information is converted into Wilson
coefficients, and how observables feed the statistical layer.

\image html hyperiso_data_flow.png "HyperIso high-level data flow"

The principal user-facing C++ façades are:

- @ref HyperisoMaster for initialization, input paths and optional runtime integrations;
- @ref WilsonInterface for Wilson-coefficient matching, running and inspection;
- @ref ObservableInterface for observable selection and prediction;
- @ref StatisticInterface for uncertainty propagation, fits, scans and contours.

## First Python calculation

```python
from pyhyperiso.Common import Model, Observables, QCDOrder
from pyhyperiso.Core import HyperisoConfig, HyperisoMaster
from pyhyperiso.Observable import ObservableInterface

config = HyperisoConfig(model=Model.SM)
master = HyperisoMaster()
master.init(lha_file="Assets/lha/si_input.flha", config=config)

observables = ObservableInterface()
observables.add_observable(Observables.BR_BS_MUMU, QCDOrder.NNLO)
print(observables.compute_observable(Observables.BR_BS_MUMU))
```

## First C++ calculation

```cpp
#include <iostream>

#include "HyperisoMaster.h"
#include "Include.h"
#include "ObservableInterface.h"

int main() {
    HyperisoConfig config;
    config.model = Model::SM;

    HyperisoMaster master;
    master.init("Assets/lha/si_input.flha", config);

    ObservableInterface observables;
    observables.add_observable(Observables::BR_BS_MUMU, QCDOrder::NNLO);
    std::cout << observables.compute_observable(Observables::BR_BS_MUMU) << '\n';
}
```

## Installation routes

- Python users: `python -m pip install pyhyperiso==1.0.2` once the release is published.
- C++ and CLI users: configure and install the CMake project; see @ref getting_started.
- GUI users: install the Python package and the Dash requirements, then start `GHyperiso/HyperisoDashGUI/app.py`.

<div class="hyperiso-note">
<b>Optional backends are explicit.</b> MARTY and SOFTSUSY are not required for
standard SM calculations or for the frozen release reference suite. Enable or
register them only for workflows that need them.
</div>

## Publication and citation

The source repository, release assets, reproducibility metadata and
`CITATION.cff` are intended to identify the same immutable release. When using
HyperIso in scientific work, cite both the associated article and the software
release tag.
