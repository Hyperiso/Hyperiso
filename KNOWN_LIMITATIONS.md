# Known limitations — HyperIso 1.0.2

This file distinguishes supported release behaviour from experimental or incomplete
paths. Unsupported operations fail explicitly where possible instead of returning
placeholder numerical values.

## Explicitly unsupported APIs

- Generic running through `EWHelper::alpha_em(scale)` is not implemented. Use the
  precomputed values populated in the `EW` block. The generic method throws
  `std::logic_error`.
- Analytical methods of `LikelihoodMarginal` (`logpdf`, derivatives, CDF, inverse
  CDF, mean and standard deviation) are not defined in 1.0.2. The class is a
  sampler; unsupported analytical calls throw `std::logic_error`.
- Direct theoretical-observable input through the LHA input flag is not implemented
  and throws instead of being silently ignored.
- `Ds -> e nu` is not implemented; the `DslnuDecay` implementation supports the
  muon and tau channels only.

## Experimental or validation-sensitive areas

- R7 validates a selected SUSY B-sector coefficient benchmark from an archived
  spectrum. Broader SUSY/NMSSM paths, especially neutral-meson mixing and less-used
  operator groups, still require additional independent benchmark coverage.
- `FWCOEF` parsing remains less extensively exercised than the standard Wilson
  input paths.
- The neutral-meson-mixing implementation contains a literature-matching shift
  inherited from the validated implementation. Its provenance must be discussed in
  the accompanying scientific documentation before changing it.
- Student-t copula and some extended covariance/correlation paths are not part of
  the seven-case 1.0.2 reproducibility gate.

## Platform support

The 1.0.2 Python binary release targets Linux x86-64 with CPython 3.10–3.12.
macOS, Windows, PyPy and musllinux are not release targets for 1.0.2.
