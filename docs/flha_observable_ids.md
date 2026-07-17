# FLHA observable type identifiers in HyperIso 1.0.2

HyperIso extends the FLHA observable-type field for observables not assigned by
the original convention. The native C++ mapper is authoritative; the
InputHelper and Dash maps are checked against it during development.

## Canonical identifiers

| Type | Meaning |
|---|---|
| `92ij` | Polarization. `i=1`: fermion/lepton; `i=2`: longitudinal vector; `i=3`: transverse vector. `j` is the one-based position of the particle of interest in the ordered daughter list. |
| `932` | Transverse fraction `F_T`. |
| `933` | Angular coefficient `alpha_K`. |

Examples:

- `B -> D(*) tau nu`: the tau is daughter 2, so its fermion polarization uses
  `9212`;
- `B -> D* tau nu`: the `D*` is daughter 1, so its longitudinal polarization
  uses `9221`;
- a future transverse polarization of the first vector daughter would use
  `9231`.

## Migration from 1.0.0--1.0.1

| Legacy type | Canonical type | Behaviour in 1.0.2 |
|---|---|---|
| `92015` | `9212` | Accepted on input and translated automatically. |
| `921423` | `9221` | Accepted on input and translated automatically. |
| `932` used for `alpha_K` | `933` | Cannot be auto-detected because `932` is also the valid `F_T` type for the same daughters. Regenerate the file or change the type explicitly. |

HyperIso 1.0.2 always writes canonical identifiers. This makes the builtin
observable map one-to-one and guarantees deterministic reverse lookup.
