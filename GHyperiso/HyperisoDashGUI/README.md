# HyperIso Dash GUI

The HyperIso Dash GUI is a local graphical interface for driving the `pyhyperiso` workflow from a browser. It is intended for exploratory work, demonstrations and interactive scans, while scripted production workflows should prefer the C++ API, Python API or CLI.

## Features

The application is organized around five pages:

| Page | Purpose |
|---|---|
| Core | Initialize or switch the active LHA/SLHA/FLHA input and inspect loaded blocks by parameter type. |
| Wilson | Build Wilson interfaces, request matching/running values and scan one or two parameters. |
| Observable | Select observables or decays, handle binned observables and scan parameter dependence. |
| Stat | Run chi-square-oriented uncertainty, fit and likelihood-contour workflows. |
| QCD | Inspect `alpha_s`, running masses and QCD constants. |

## Project structure

```text
HyperisoDashGUI/
├── app.py
├── requirements.txt
├── README.md
└── pyhyperiso_dash/
    ├── callbacks.py
    ├── components.py
    ├── domain.py
    ├── figures.py
    ├── latex.py
    ├── services.py
    ├── assets/
    ├── latex_data/
    └── pages/
```

## Installation

Install HyperIso and the GUI dependencies in the same Python environment:

```bash
python -m pip install -e ../../Hyperiso/Hyperiso
python -m pip install -r requirements.txt
```

If you run from the repository root:

```bash
python -m pip install -e ./Hyperiso/Hyperiso
python -m pip install -r GHyperiso/HyperisoDashGUI/requirements.txt
```

## Run

From `GHyperiso/HyperisoDashGUI`:

```bash
python app.py
```

Open:

```text
http://127.0.0.1:8050
```

## Runtime note

HyperIso uses a C++ runtime behind the Python bindings. The Dash application keeps one process-local `HyperisoMaster` instance alive in `pyhyperiso_dash.services.RUNTIME`.

For local use, run a single Dash process. If deploying through a production server, use one worker unless you intentionally want one independent HyperIso runtime per worker.

## Statistical mode

The dashboard is intentionally locked to
`StatisticLikelihoodMode.CHI2_MC_COVARIANCE`. The profiled-nuisance likelihood
is not selectable from the GUI.

Fit parameters are listed in an editable table with one lower and upper contour
bound per parameter. Any two fitted parameters can be selected as the displayed
axes. When more than two parameters are fitted, the GUI exposes the core
`ProfilingMethod` choices (`SLICE`, `FREE_PROJECTION`, and
`PRIOR_CONSTRAINED_PROJECTION`) together with `ContourAlgorithm`, optional
fallback extraction, profiler mode, confidence levels, and one common contour
resolution. The plotted geometry comes directly from `ContourEngine.paths`, not
from a separate rectangular likelihood scan. Confidence regions use a white
publication-style theme, MathJax axis labels, semi-transparent filled regions
with darker boundaries, and an automatic viewport fitted to the computed paths
while remaining inside the bounds supplied by the selected fit parameters.

The GUI maps covariance and nuisance-pruning settings to
`StatisticConfig.advanced`, while common Monte-Carlo settings remain on the
top-level config. All C++ terminal diagnostics are disabled in the Dash path.
A non-blocking animated progress banner remains visible inside the Statistics
page while uncertainty and fit callbacks are running. The former page-wide Dash
spinner was removed because it hid the complete interface during long jobs.
Exact accepted-draw progress will still require a future callback/event API from
the C++ Monte-Carlo engine. Dynamic dropdown labels expose their complete text
through hover tooltips, and responsive grid/overflow constraints prevent cards,
tables, and selectors from overlapping at narrow or intermediate window sizes.

## Wilson identifiers

Wilson group and coefficient selectors are populated through the core mapper
API and requests are passed as dynamic string identifiers. This avoids coupling
the GUI to a potentially incomplete legacy Python enum when new built-in groups
or coefficients are added in C++.

## Parameter scans

Scan callbacks must always restore modified parameters after calculation. The current service layer follows the pattern:

1. read the original central value;
2. set a temporary value;
3. compute Wilson coefficients or observables;
4. restore the original value in a `finally` block.

This is required because the C++ runtime is global within the Python process.

## Development checklist

Before opening a pull request that changes the GUI:

```bash
python -m compileall GHyperiso/HyperisoDashGUI
python GHyperiso/HyperisoDashGUI/app.py
```

Then manually check:

- application starts;
- Core page initializes an SM input;
- Wilson page computes a built-in coefficient;
- Observable page computes at least one observable;
- QCD page renders a plot;
- no callback mutates parameters without restoring them.

## Container

A Dash-oriented container can be built from the repository root:

```bash
docker build -f docker/Dockerfile.dash -t hyperiso-dash .
docker run --rm -p 8050:8050 hyperiso-dash
```
