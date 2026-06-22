# HyperIso Python examples

This directory contains Python examples for the public `pyhyperiso` interface. They mirror the C++ examples whenever possible so that workflows can be reproduced in both languages.

## Installation

From the repository root:

```bash
python -m pip install --upgrade pip
python -m pip install -e ./Hyperiso/Hyperiso
```

Optional plotting and notebook tools:

```bash
python -m pip install matplotlib pandas numpy jupyter
```

## Organization

| Directory | Purpose |
|---|---|
| `Core/` | Initialization, parameter access, block inspection and QCD helpers. |
| `Wilson/` | Wilson matching/running and custom Wilson groups. |
| `Observable/` | Observable calculations, decay configurations and scans. |
| `Statistic/` | Uncertainty propagation, likelihood scans and fit examples. |

## Quick start

```bash
python Hyperiso/Hyperiso/examples_python/Core/base_core_example.py
python Hyperiso/Hyperiso/examples_python/Wilson/base_wilson_example.py
python Hyperiso/Hyperiso/examples_python/Observable/base_observable_example.py
```

## Minimal Python workflow

```python
from pyhyperiso.Common import Model, QCDOrder, Observables
from pyhyperiso.Core import HyperisoConfig, HyperisoMaster
from pyhyperiso.Observable import ObservableInterface

hyp = HyperisoMaster()
hyp.init(
    lha_file="Assets/lha/si_input.flha",
    config=HyperisoConfig(model=Model.SM),
)

obs = ObservableInterface()
obs.add_observable(Observables.BR_BS_MUMU, QCDOrder.NNLO)
print(obs.compute_observable(Observables.BR_BS_MUMU))
```

## Good practice for scripts and scans

- Keep the input LHA/SLHA/FLHA file path explicit.
- Set random seeds in statistical or Monte-Carlo examples.
- Save scan metadata together with numerical outputs.
- Prefer `Assets/lha/si_input.flha` for public SM smoke tests.
- Use separate output directories and avoid writing generated files into the source tree.

## Adding a new Python example

New examples should include:

- a short module docstring explaining the workflow;
- required backend flags, if any;
- exact input paths;
- stable printed output;
- optional expected values for regression testing.

## Running tests

The package-level tests live in `Hyperiso/Hyperiso/pyhyperiso/test`:

```bash
python -m pip install pytest
pytest Hyperiso/Hyperiso/pyhyperiso/test
```
