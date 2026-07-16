# Inputs, parameter blocks and precedence {#inputs}

HyperIso stores numerical inputs in typed blocks identified by a block name and
one or more integer indices. The convention follows LHA/SLHA/FLHA whenever a
standard identifier exists and uses project-defined IDs for additional data.

## Input precedence

The runtime input system is layered from least to most specific:

1. **Distributed JSON defaults** define the reference database shipped with the
   release.
2. **User YAML overrides** replace selected central values, uncertainties or
   distribution metadata without editing the distributed assets.
3. **LHA-family input** supplies the model-specific LHA, SLHA or FLHA values and
   has the highest priority for overlapping entries.

<div class="hyperiso-note">
Do not edit the distributed JSON files for a scientific scan. Store changes in
a YAML override and archive that file together with the LHA-family input.
</div>

## Parameter identifiers

A parameter is addressed through a `ParamId`, which combines a parameter type,
block name and integer indices. This identifier is used consistently by the
core, observable dependencies and statistical layer.

Python example:

```python
from pyhyperiso.Common import ParamId, ParameterType

f_bs = ParamId(ParameterType.FLAVOR, "FCONST", [531, 1])
```

C++ example:

```cpp
ParamId f_bs(ParameterType::FLAVOR, "FCONST", {531, 1});
```

The parameter provider/setter examples show how to inspect and update values
through the public API without reaching into the internal cache.

## Path overrides and additional blocks

Before initialization, @ref HyperisoMaster can:

- register additional LHA block prototypes;
- replace selected default JSON or user YAML paths;
- set writable MARTY and spectrum-cache directories;
- register external MARTY or SOFTSUSY installations.

Path values are validated before they are accepted. Default input paths must
point to JSON files, user override paths must point to YAML files, and directory
entries must exist or be creatable where documented.

## LHA-family files

Use the input format that matches the workflow:

- LHA for generic block-oriented numerical inputs;
- SLHA for supersymmetric spectra and related parameters;
- FLHA for flavour observables, uncertainties and correlations.

A spectrum already produced by an external program can be archived and reused.
This is preferable for release regression tests because it removes dependence
on the external generator version.

## Database export

The initialized database can be exported for inspection or archival. In Python:

```python
from pyhyperiso.Core import DatabaseWriter

writer = DatabaseWriter()
writer.write("database.json")
writer.write_blocks("inputs.yaml", ["SMINPUTS", "MASS"])
```

The output suffix selects JSON, YAML, LHA, SLHA or FLHA. Full and filtered
examples are available in:

- `Hyperiso/Hyperiso/examples_python/Core/database_writer_example.py`;
- `Hyperiso/Hyperiso/examples_cpp/Core/database_writer_example.cpp`.

## Reproducible input practice

For each scientific result, retain:

- the software release tag;
- the LHA/SLHA/FLHA input;
- every YAML override;
- the model and QCD-order configuration;
- external spectrum files or the exact external-tool provenance;
- random seeds and thread settings for Monte-Carlo calculations.

The frozen release suite described in @ref reproducibility follows this model.
