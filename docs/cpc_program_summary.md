# CPC Program Summary — HyperIso 1.0.0

**Program title:** HyperIso
**Version:** 1.0.0
**Licensing provisions:** GNU General Public License, version 3 or later. See
`LICENSE`, `LICENSES/` and `THIRD_PARTY_NOTICES.md`.
**Programming languages:** C++20 and Python 3.10–3.12.
**Operating systems:** Linux x86-64 for the 1.0.0 binary release.
**External libraries:** CMake, GSL, Eigen, pybind11 and the Python dependencies
listed in `pyproject.toml`. Bundled components and their licensing are recorded
in `THIRD_PARTY_NOTICES.md`.
**Nature of problem:** Calculation of flavour-physics Wilson coefficients,
flavour observables and their statistical interpretation in the Standard Model
and supported beyond-the-Standard-Model scenarios.
**Solution method:** HyperIso separates parameter ingestion, Wilson matching and
running, observable calculation and statistical inference into reusable C++
modules exposed through C++, Python and a command-line interface.
**Restrictions:** The supported and experimental paths are listed in
`KNOWN_LIMITATIONS.md`. Generic running of `alpha_em`, analytical methods of the
likelihood sampler and direct theoretical-observable LHA input are not supported
in 1.0.0 and fail explicitly.
**Typical running time:** The deterministic reference cases complete on a
standard Linux workstation; the 1,000-draw Monte-Carlo case is the longest of
the five. Exact timing is hardware-dependent and is recorded by CI rather than
claimed as a universal benchmark.
**Reproducibility:** `reproducibility/manifest.json` defines five fixed commands.
The statistical case uses 1,000 draws and seed `123456`. Normalized outputs,
sample values, hashes and strict tolerances are stored in
`reproducibility/expected_outputs/`. The release workflow rebuilds the final
source and must match all frozen references before publication.
**Related article:** *HyperIso: A general BSM calculator for flavour
observables*. Final journal identifiers should replace the current preprint
record when assigned.
