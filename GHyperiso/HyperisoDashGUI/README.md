# PyHyperiso Dash GUI

A Dash application for driving the PyHyperiso workflow from a professional dark UI inspired by your HepMCGUI project.

The application is organized around five pages:

1. **Core** — initialize or switch the active LHA/SLHA/FLHA input and inspect loaded blocks by `ParameterType`.
2. **Wilson** — build `WilsonInterface`, run matching/running requests, and scan one or two parameters with automatic restore.
3. **Observable** — select observables or whole decays, handle binned observables, compute predictions and scan parameter dependence.
4. **Stat** — χ²-only statistical workflow: uncertainty propagation, best fits, nuisance output, correlations and 2D ΔNLL contours.
5. **QCD** — compute and plot `alpha_s`, running masses and backend QCD constants.

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
    ├── services.py
    ├── assets/
    │   └── styles.css
    ├── data/
    │   └── uploaded_lha/
    └── pages/
        ├── core.py
        ├── observable.py
        ├── qcd.py
        ├── stat.py
        └── wilson.py
```

This follows the HepMCGUI spirit: a small `app.py`, a custom CSS theme, reusable domain modules, and heavy logic kept outside the layout definitions.

## Installation

Install Dash dependencies inside the same environment where `pyhyperiso` is importable:

```bash
pip install -r requirements.txt
```

## Run

From the project root:

```bash
python app.py
```

Open:

```text
http://127.0.0.1:8050
```

## Important runtime note

Hyperiso uses a C++ singleton behind the bindings. This Dash app intentionally keeps a process-local `HyperisoMaster` alive in `pyhyperiso_dash.services.RUNTIME`.

For development and local use, run a single Dash process. If deployed with gunicorn or another server, prefer one worker unless you explicitly want one independent Hyperiso singleton per worker.

## Notes on statistical mode

The Stat page forces:

```python
StatisticLikelihoodMode.CHI2_MC_COVARIANCE
```

Minuit/Laplace-only options are intentionally not shown. The 2D contour plot is built from the available likelihood scan API rather than from the opaque C++ `Contour` wrapper, so the result can be directly displayed as a Plotly heatmap/contour.

## Parameter scans

Wilson and observable scans mutate parameters through `ParameterSetter`, read the original central values through `ParameterProvider`, and restore values in a `finally` block. This is important because the C++ runtime is global for the process.


## Patch note: model-aware ParameterType filtering

The GUI filters `ParameterType.BSM` whenever the active Hyperiso model is `SM`.
This prevents calls into the C++ singleton with a namespace that is not
registered for the current model, which can otherwise terminate the process
before Python can catch an exception.
