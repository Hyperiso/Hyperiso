# CPC Program Summary — HyperIso 1.0.0

**Program title:** HyperIso
**Version:** 1.0.0
**CPC Library link to program files:** To be assigned by the CPC Technical Editor.
**Developer's repository:** https://github.com/HyperIso/HyperIso/tree/v1.0.0
**Licensing provisions:** GNU General Public License 3 (GPL), with the option to use any later version.
**Programming languages:** C++20 and Python 3.10–3.12.
**Operating systems:** Linux x86-64 for the 1.0.0 binary release; source builds are tested on current Ubuntu runners.
**External libraries and optional programs:** CMake, GSL, Eigen and pybind11. The source distribution includes 2HDMC and MinuitCpp. MARTY and SOFTSUSY are optional external programs and are not required by the frozen reference suite.
**Supplementary material:** User and installation guides, API documentation, example LHA/SLHA/FLHA inputs, validation tables, fixed-seed reproducibility scripts, expected outputs and numerical metadata.
**References:** The associated HyperIso article; the SuperIso/FLHA references cited in the manuscript; the 2HDMC program article for the archived THDM spectrum.

**Nature of problem:**
Precision flavour observables constrain the Standard Model and a broad range of beyond-the-Standard-Model scenarios. Practical analyses require a consistent treatment of numerical inputs, Wilson coefficients at matching and low-energy scales, decay observables, experimental correlations, nuisance-parameter distributions and statistical fits. Existing programs often expose these tasks through model-specific or monolithic interfaces, which complicates reuse in new models and makes it difficult to keep C++, Python and command-line calculations numerically aligned. HyperIso addresses this problem by providing a common computational backend for native Standard Model, two-Higgs-doublet and supersymmetric calculations, while retaining an optional route for user-defined models through externally generated Wilson coefficients. The software also supplies a layered input system and a release-level reproducibility suite so that benchmark calculations can be repeated from archived spectrum files without rerunning external spectrum generators.

**Solution method:**
HyperIso uses a modular C++20 architecture with explicit ports and adapters for parameter loading, spectrum handling, Wilson-coefficient construction, observable evaluation and statistical inference. Distributed JSON files define the reference database, user YAML files provide controlled overrides, and LHA-family files supply model-dependent numerical inputs. Native Wilson-coefficient implementations support the Standard Model, THDM and SUSY paths at the perturbative orders available for each coefficient group. The same backend is exposed through C++, Python, a command-line interface and a graphical interface. Statistical routines combine nuisance marginals and correlations, propagate uncertainties by fixed-seed Monte Carlo sampling and provide chi-square-based profiling. The release archive includes deterministic SM cases and spectrum-file-based THDM and SUSY benchmarks; these use already generated spectra and therefore do not require 2HDMC or SOFTSUSY during reproduction.

**Additional comments, restrictions and unusual features:**
The generic MARTY backend is optional and must be installed separately. It is validated in the article against native calculations, but is deliberately excluded from the mandatory release reproducibility gate because a clean reproduction would otherwise depend on an external symbolic-software installation and its version. Native THDM and SUSY reference cases instead consume archived spectrum files with the `--spectrum` option, making the release checks independent of external spectrum generators. Generic BSM coefficients produced through MARTY are limited to the perturbative order implemented by that workflow. Features listed as experimental or unsupported in `KNOWN_LIMITATIONS.md` fail explicitly or are documented separately. The five Standard Model references include a fixed-seed 200-draw statistical regression; this is a software non-regression test rather than a Monte Carlo convergence study. Exact platform, compiler, dependency and binary hashes are recorded when the final references are frozen.

**Typical running time:** Deterministic CLI references normally complete in seconds on a contemporary Linux workstation. The fixed-seed 200-draw statistical case is the longest mandatory example; per-case timings are recorded when the final references are frozen because they depend on hardware and compiler settings.
**Example input/output:** `reproducibility/inputs/`, `reproducibility/expected_outputs/` and `reproducibility/manifest.json`.
**Classification:** 11.1, 11.3, 11.6.
**Related article:** *HyperIso: A general BSM calculator for flavour observables*.
