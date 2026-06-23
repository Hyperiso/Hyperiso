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

The dashboard focuses on the chi-square / Monte-Carlo covariance workflow for interactive use. More advanced likelihood profiling and backend choices should be exposed through the Python API first and added to the GUI only once the behavior is stable.

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
