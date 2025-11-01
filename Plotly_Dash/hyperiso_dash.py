# -*- coding: utf-8 -*-
"""
Hyperiso Dashboard — Plotly Dash (dark theme)
-------------------------------------------------

Features (single-file app):
- Scan builder (up to 5 params; up to 100k LHA files), grid/random sampling, preview, progress, plots
- Core inspector: LHA blocks overview (bar + donut), block browser, key-values preview
- QCD panel: α_s(μ) and masses vs μ, QCD constants
- Wilson panel: running/matching coefficients, orders & contributions, parameter variation
- Observables: compute, overlay, parameter scans, uncertainty bands, leading uncertainties
- Stats: χ² calculator, 1D/2D scans, correlation matrix

Notes:
- Keeps a single C++ singleton via PyHyperisoMaster. Safe re-init when switching LHA/config.
- Robust error handling with alerts and safe fallbacks.
- Optimized with caching, vectorized sampling, and lazy compute.

Run:  python hyperiso_dash_app.py
"""

import os
import io
import re
import json
import time
import math
import itertools as it
from pathlib import Path
from typing import List, Dict, Any, Tuple

import numpy as np
import pandas as pd

import dash
from dash import dcc, html, Input, Output, State, MATCH, ALL, ctx
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import plotly.express as px

# =============== Hyperiso / wrappers (user-provided modules) ===============
# Prefer local wrappers shipped by the user. Fallback to package modules if present.
try:
    from HyperisoMaster import PyHyperisoMaster
except Exception:
    from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster  # type: ignore

try:
    from Configs import PyWilsonBuildConfig, PyWilsonRequest, PyAlphasConfig, PyMassConfig
except Exception:
    from pyhyperiso.core.PhysicalModel.WilsonInterface import PyWilsonBuildConfig, PyWilsonRequest  # type: ignore
    from pyhyperiso.core.Common.Configs import PyAlphasConfig, PyMassConfig  # type: ignore

try:
    from WilsonInterface import PyWilsonInterface
except Exception:
    from pyhyperiso.core.PhysicalModel.WilsonInterface import PyWilsonInterface  # type: ignore

try:
    from ObservableInterface import PyObservableInterface
except Exception:
    PyObservableInterface = None  # will guard below

try:
    from CorrelationProvider import PyCorrelationProvider, CorrelationType
except Exception:
    PyCorrelationProvider, CorrelationType = None, None

try:
    from ParamaterProvider import PyParameterProvider
except Exception:
    PyParameterProvider = None

try:
    from pyhyperiso.core.Core.ParameterSetter import PyParameterSetter, PyParamId
except Exception:
    PyParameterSetter, PyParamId = None, None

try:
    # Optional API adapter for LHA block-level queries
    from pyhyperiso.core.Core.APIAdapter import PyAPIAdapter
except Exception:
    PyAPIAdapter = None

from pyhyperiso.core.Common.GeneralEnum import (
    Model, QCDOrder, WGroup, WCoeff, ContributionType, WilsonBasis,
    Observables, ParameterType, MassType, ScaleType
)

try:
    from pyhyperiso.core.Common.General import PyLhaID
except Exception:
    # Minimal shim
    class PyLhaID:
        def __init__(self, *parts):
            self.parts = [int(parts[0]) if parts else 0]
        def __int__(self):
            return self.parts[0]

# ----------------- Global constants -----------------
PLOT_TEMPLATE = "plotly_dark"
DARK_BG = "#0f1116"
CARD_BG = "#151924"
ACCENT = "#04d9ff"

MAX_PARAMS = 5
MAX_FILES = 100_000

# ----------------- Helpers: figures -----------------
def make_alert(msg: str, color: str = "danger"):
    return dbc.Alert(msg, color=color, is_open=True, duration=5000, className="mt-2")

# def donut(labels, values, title: str):
#     fig = go.Figure(go.Pie(labels=labels, values=values, hole=0.55))
#     fig.update_layout(template=PLOT_TEMPLATE, title=title, legend=dict(orientation="h"))
#     return fig

def stylize(fig: go.Figure, title: str | None = None) -> go.Figure:
    fig.update_layout(
        template=PLOT_TEMPLATE,
        paper_bgcolor=CARD_BG,
        plot_bgcolor=CARD_BG,
        font=dict(color="#E8EEF8"),
        title=title if title is not None else fig.layout.title.text,
        legend=dict(bgcolor="rgba(0,0,0,0)"),
    )
    return fig

def donut(labels, values, title: str):
    fig = go.Figure(go.Pie(labels=labels, values=values, hole=0.55))
    return stylize(fig, title)

# =============== Singleton context for Hyperiso ===============
class HyperisoContext:
    _instance = None
    DEFAULT_LHA = "lha/testinput_thdm.lha"

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance.hyp = PyHyperisoMaster()
            cls._instance.config = None
            cls._instance.lha_path = None
            cls._instance.api = PyAPIAdapter() if PyAPIAdapter else None
            # tente un bootstrap très défensif
            try:
                cls._instance.ensure_ready()
            except Exception:
                pass
        return cls._instance

    def ensure_ready(self, lha_path: str | None = None):
        """
        Initialise Hyperiso si pas encore fait. À appeler en tête de TOUTE callback
        qui touche au C++ (QCD/Wilson/Obs/Core/Stats).
        """
        if self.lha_path and self.config:
            return  # déjà prêt

        from pyhyperiso.core.Core.Config import PyConfig, ExternalFlag
        # config par défaut minimaliste
        fset = {
            ExternalFlag.IS_LHA_SPECTRUM: False,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            # ExternalFlag.USE_MARTY: False,
        }
        cfg = PyConfig(flags=fset, model=Model.SM)
        # chemin LHA : priorité au paramètre, sinon valeur précédente, sinon DEFAULT
        lha = lha_path or self.lha_path or self.DEFAULT_LHA
        lha = str(lha)

        # init/switch robuste
        if self.lha_path is None:
            self.hyp.init(lha_file=lha, config=cfg)
        else:
            if hasattr(self.hyp, "switch_lha"):
                self.hyp.switch_lha(lha_file=lha, config=cfg)
            else:
                self.hyp.init(lha_file=lha, config=cfg)

        self.config = cfg
        self.lha_path = lha

    def init_or_switch(self, lha_path: str, config):
        # évite réinit coûteuse
        if self.lha_path == lha_path and repr(self.config) == repr(config):
            return
        if hasattr(self.hyp, "switch_lha"):
            self.hyp.switch_lha(lha_file=lha_path, config=config)
        else:
            self.hyp.init(lha_file=lha_path, config=config)
        self.config = config
        self.lha_path = lha_path

CTX = HyperisoContext()

# =================== Scan utilities ===================
class LhaScanPlanner:
    @staticmethod
    def grid_from_ranges(ranges: List[Dict[str, Any]]) -> pd.DataFrame:
        """Builds a full-factorial grid from param ranges.
        Each dict: {name, block, code, ptype, min, max, steps}
        """
        axes = []
        cols = []
        for r in ranges:
            vals = np.linspace(float(r["min"]), float(r["max"]), int(r["steps"]))
            axes.append(vals)
            cols.append(r["name"])
        combos = np.array(list(it.product(*axes)))
        df = pd.DataFrame(combos, columns=cols)
        return df

    @staticmethod
    def random_from_ranges(ranges: List[Dict[str, Any]], n_samples: int, seed: int = 42) -> pd.DataFrame:
        rng = np.random.default_rng(seed)
        data = {}
        for r in ranges:
            low, high = float(r["min"]), float(r["max"])
            data[r["name"]] = rng.uniform(low, high, n_samples)
        return pd.DataFrame(data)

class LhaIO:
    @staticmethod
    def read_lha(path: Path) -> str:
        return path.read_text(encoding="utf-8", errors="ignore")

    @staticmethod
    def write_lha(base_text: str, out_path: Path, updates: List[Tuple[str, int, float]]):
        """Simple, robust LHA patcher: ensure BLOCK exists and set value lines.
        updates: list of (BLOCK_NAME, code, value)
        """
        blocks = LhaIO._parse_blocks(base_text)
        for blk, code, val in updates:
            bname = blk.upper()
            if bname not in blocks:
                blocks[bname] = {}
            blocks[bname][int(code)] = float(val)
        text = LhaIO._compose_blocks(blocks)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(text)

    @staticmethod
    def _parse_blocks(text: str) -> Dict[str, Dict[int, float]]:
        blocks: Dict[str, Dict[int, float]] = {}
        cur = None
        for line in text.splitlines():
            if line.strip().upper().startswith("BLOCK "):
                cur = line.split()[1].upper()
                blocks.setdefault(cur, {})
            elif cur and line.strip() and not line.strip().startswith("#"):
                toks = line.split()
                if len(toks) >= 2 and toks[0].lstrip("-").isdigit():
                    try:
                        idx = int(toks[0])
                        val = float(toks[1])
                        blocks[cur][idx] = val
                    except Exception:
                        pass
        return blocks

    @staticmethod
    def _compose_blocks(blocks: Dict[str, Dict[int, float]]) -> str:
        out = io.StringIO()
        out.write("# Auto-generated LHA by Hyperiso Dash\n")
        for name, entries in blocks.items():
            out.write(f"BLOCK {name}\n")
            for k, v in sorted(entries.items()):
                out.write(f"  {k:5d}   {v:.10e}   # auto\n")
        return out.getvalue()

# =================== Build Dash App ===================
external_stylesheets = [dbc.themes.CYBORG]
app = dash.Dash(__name__, external_stylesheets=external_stylesheets, title="Hyperiso Dashboard", suppress_callback_exceptions=True)
server = app.server

# --------------- Reusable UI pieces ---------------

def section(title: str, children):
    return dbc.Card([
        dbc.CardHeader(html.H4(title, className="mb-0")),
        dbc.CardBody(children)
    ], className="mb-3", style={"background": CARD_BG})

# ========== Controls: Global config ==========
models_options = [{"label": m.name, "value": m.name} for m in Model]
contrib_opts = [{"label": c.name, "value": c.name} for c in ContributionType]
orders_opts = [{"label": o.name, "value": o.name} for o in QCDOrder]

global_controls = section(
    "Configuration globale",
    dbc.Row([
        dbc.Col([
            dbc.Label("Modèle"),
            dcc.Dropdown(id="cfg-model", options=models_options, value=Model.SM.name, clearable=False),
        ], md=2),
        dbc.Col([
            dbc.Checklist(
                options=[{"label": "USE_MARTY", "value": "USE_MARTY"},
                         {"label": "IS_LHA_SPECTRUM", "value": "IS_LHA_SPECTRUM"}],
                value=[], id="cfg-flags", inline=True
            ),
            dbc.Input(id="cfg-marty-name", placeholder="mty_model_name (si USE_MARTY)", style={"display": "none"}),
            dbc.Input(id="cfg-marty-path", placeholder="mty_model_path (si USE_MARTY)", style={"display": "none"}),
        ], md=4),
        dbc.Col([
            dbc.Input(id="cfg-lha-path", placeholder="Chemin vers le fichier LHA", value="lha/testinput_thdm.lha"),
        ], md=4),
        dbc.Col([
            dbc.Button("Appliquer la config", id="cfg-apply", color="info", className="mt-4"),
        ], md=2),
        html.Div(id="cfg-status", className="mt-2")
    ])
)

# ========== Tab: Scan builder ==========
param_row = lambda i: dbc.Row([
    dbc.Col(dbc.Input(id={"type": "scan-name", "index": i}, placeholder=f"Param {i+1} (ex: MASS)"), md=2),
    dbc.Col(dbc.Input(id={"type": "scan-block", "index": i}, placeholder="BLOCK (ex: MASS, BSM, ...)"), md=2),
    dbc.Col(dbc.Input(id={"type": "scan-code", "index": i}, placeholder="LHA ID (ex: 37)", type="number"), md=2),
    dbc.Col(dcc.Dropdown(id={"type": "scan-ptype", "index": i}, options=[{"label": t.name, "value": t.name} for t in ParameterType], placeholder="ParameterType"), md=2),
    dbc.Col(dbc.Input(id={"type": "scan-min", "index": i}, placeholder="min", type="number"), md=1),
    dbc.Col(dbc.Input(id={"type": "scan-max", "index": i}, placeholder="max", type="number"), md=1),
    dbc.Col(dbc.Input(id={"type": "scan-steps", "index": i}, placeholder="steps", type="number"), md=1),
    dbc.Col(dbc.Button("✖", id={"type": "scan-del", "index": i}, color="secondary"), md=1),
], className="mb-2")

scan_layout = section(
    "Création de scan LHA",
    [
        dbc.Row([
            dbc.Col(dbc.RadioItems(
                id="scan-mode", inline=True,
                options=[{"label": "Grille", "value": "grid"}, {"label": "Aléatoire", "value": "random"}],
                value="grid"
            ), md=3),
            dbc.Col(dbc.Input(id="scan-n-samples", type="number", value=1000, min=1, step=100), md=2),
            dbc.Col(dbc.Input(id="scan-outdir", placeholder="Dossier de sortie", value="./out_scan"), md=4),
            dbc.Col(dbc.Button("Aperçu", id="scan-preview", color="secondary", className="mt-1"), md=1),
            dbc.Col(dbc.Button("Générer", id="scan-generate", color="success", className="mt-1"), md=2),
        ]),
        html.Hr(),
        html.Div([param_row(i) for i in range(MAX_PARAMS)], id="scan-rows"),
        dcc.Store(id="scan-params-store"),
        html.Div(id="scan-summary"),
        dbc.Progress(id="scan-progress", striped=True, animated=True, style={"height": "20px"}, className="mt-2"),
        dbc.Row([
            dbc.Col(dcc.Graph(id="scan-hist"), md=6),
            dbc.Col(dcc.Graph(id="scan-pairs"), md=6),
        ])
    ]
)

# ========== Tab: Core inspector ==========
core_layout = section(
    "Inspection du LHA (Core)",
    [
        dbc.Row([
            dbc.Col(dcc.Graph(id="core-blocks-bar"), md=6),
            dbc.Col(dcc.Graph(id="core-blocks-donut"), md=6),
        ]),
        dbc.Row([
            dbc.Col([
                dbc.Input(id="core-block-filter", placeholder="Nom du block (ex: MASS)"),
                dcc.Loading(dcc.Graph(id="core-block-preview"), type="cube")
            ], md=6),
            dbc.Col([
                html.Div(id="core-flags"),
                html.Div(id="core-paths")
            ], md=6)
        ])
    ]
)

# ========== Tab: QCD panel ==========
qcd_layout = section(
    "QCD: α_s et masses",
    [
        dbc.Row([
            dbc.Col(dbc.Input(id="qcd-scale-min", type="number", value=1.0), md=2),
            dbc.Col(dbc.Input(id="qcd-scale-max", type="number", value=200.0), md=2),
            dbc.Col(dbc.Input(id="qcd-n", type="number", value=200), md=2),
            dbc.Col(dcc.Dropdown(id="qcd-mass-pdg", options=[{"label": f"PDG {i}", "value": i} for i in [4,5,6]], value=5), md=3),
            dbc.Col(dcc.Dropdown(id="qcd-mtype", options=[{"label": t.name, "value": t.name} for t in MassType], value=MassType.MSBAR.name), md=3),
        ]),
        dbc.Row([
            dbc.Col(dcc.Graph(id="qcd-alphas-fig"), md=6),
            dbc.Col(dcc.Graph(id="qcd-mass-fig"), md=6)
        ]),
        html.Div(id="qcd-constants")
    ]
)

# ========== Tab: Wilson ==========
wilson_layout = section(
    "Wilson coefficients",
    [
        dbc.Row([
            dbc.Col(dcc.Dropdown(id="wilson-group", options=[{"label": g.name, "value": g.name} for g in WGroup], value=WGroup.B.name), md=2),
            dbc.Col(dcc.Dropdown(id="wilson-coeff", options=[{"label": w.name, "value": w.name} for w in WCoeff], value=WCoeff.C9.name), md=3),
            dbc.Col(dcc.Dropdown(id="wilson-basis", options=[{"label": b.name, "value": b.name} for b in WilsonBasis], value=WilsonBasis.STANDARD.name), md=2),
            dbc.Col(dbc.RadioItems(id="wilson-mode", options=[{"label": "Running (FR)", "value": "running"}, {"label": "Matching (FM)", "value": "matching"}], value="running", inline=True), md=3),
            dbc.Col(dbc.Input(id="wilson-scale", type="number", value=4.18, placeholder="Scale (had/match)"), md=2),
        ]),
        dbc.Row([
            dbc.Col(dcc.Checklist(id="wilson-contrib", options=contrib_opts, value=["SM", "BSM", "TOTAL"], inline=True), md=6),
            dbc.Col(dcc.Checklist(id="wilson-orders", options=orders_opts, value=["LO", "NLO", "NNLO"], inline=True), md=6),
        ]),
        dbc.Row([
            dbc.Col(dcc.Graph(id="wilson-stacked"), md=6),
            dbc.Col(dcc.Graph(id="wilson-heat"), md=6),
        ]),
        html.Hr(),
        html.H5("Variation par paramètre"),
        dbc.Row([
            dbc.Col(dcc.Dropdown(id="wilson-var-ptype", options=[{"label": t.name, "value": t.name} for t in ParameterType], value=ParameterType.BSM.name), md=2),
            dbc.Col(dbc.Input(id="wilson-var-block", placeholder="Block (ex: MASS)", value="MASS"), md=2),
            dbc.Col(dbc.Input(id="wilson-var-code", type="number", value=37), md=2),
            dbc.Col(dbc.Input(id="wilson-var-min", type="number", value=200.0), md=2),
            dbc.Col(dbc.Input(id="wilson-var-max", type="number", value=500.0), md=2),
            dbc.Col(dbc.Input(id="wilson-var-n", type="number", value=50), md=2),
        ]),
        dbc.Row([
            dbc.Col(dbc.Button("Tracer la variation", id="wilson-var-go", color="info"), md=2)
        ]),
        dbc.Row([
            dbc.Col(dcc.Graph(id="wilson-var-fig"), md=12)
        ])
    ]
)

# ========== Tab: Observables ==========
observ_layout = section(
    "Observables",
    [
        dbc.Row([
            dbc.Col(dcc.Dropdown(id="obs-select", options=[{"label": o.name, "value": o.name} for o in Observables], multi=True), md=6),
            dbc.Col(dbc.Checklist(id="obs-opts", options=[{"label": "Afficher incertitude (aire)", "value": "uncert"}, {"label": "Sources d'incertitudes (leading)", "value": "leading"}], value=["uncert"], inline=True), md=6),
        ]),
        dbc.Button("Calculer", id="obs-calc", color="primary"),
        dbc.Row([
            dbc.Col(dcc.Graph(id="obs-fig"), md=12)
        ]),
        html.Hr(),
        html.H5("Variation vs paramètre"),
        dbc.Row([
            dbc.Col(dcc.Dropdown(id="obs-var-name", options=[{"label": o.name, "value": o.name} for o in Observables], value=Observables.BR_B_XS_GAMMA.name if hasattr(Observables, "BR_B_XS_GAMMA") else None), md=3),
            dbc.Col(dcc.Dropdown(id="obs-var-ptype", options=[{"label": t.name, "value": t.name} for t in ParameterType], value=ParameterType.BSM.name), md=2),
            dbc.Col(dbc.Input(id="obs-var-block", value="MASS"), md=2),
            dbc.Col(dbc.Input(id="obs-var-code", type="number", value=37), md=1),
            dbc.Col(dbc.Input(id="obs-var-min", type="number", value=200.0), md=2),
            dbc.Col(dbc.Input(id="obs-var-max", type="number", value=500.0), md=2),
            dbc.Col(dbc.Input(id="obs-var-n", type="number", value=40), md=2),
        ]),
        dbc.Button("Tracer variation observable", id="obs-var-go", color="info"),
        dbc.Row([dbc.Col(dcc.Graph(id="obs-var-fig"), md=12)])
    ]
)

# ========== Tab: Statistics ==========
stats_layout = section(
    "Statistiques (χ², corrélations)",
    [
        dbc.Row([
            dbc.Col(dcc.Dropdown(id="stats-obs", options=[{"label": o.name, "value": o.name} for o in Observables], multi=True), md=6),
            dbc.Col(dbc.Button("Configurer", id="stats-configure", color="secondary"), md=2),
        ]),
        dbc.Row([
            dbc.Col(dcc.Graph(id="stats-corr"), md=6),
            dbc.Col(dcc.Graph(id="stats-chi2-1d"), md=6),
        ]),
        html.Hr(),
        html.H5("Scan χ² (2 paramètres)"),
        dbc.Row([
            dbc.Col(dbc.Input(id="stats-p1-name", placeholder="Param1 name", value="MASS"), md=2),
            dbc.Col(dbc.Input(id="stats-p1-code", type="number", value=37), md=1),
            dbc.Col(dbc.Input(id="stats-p1-min", type="number", value=200.0), md=2),
            dbc.Col(dbc.Input(id="stats-p1-max", type="number", value=500.0), md=2),
            dbc.Col(dbc.Input(id="stats-p1-n", type="number", value=25), md=1),
            dbc.Col(dbc.Input(id="stats-p2-name", placeholder="Param2 name", value="ALPHA"), md=2),
            dbc.Col(dbc.Input(id="stats-p2-code", type="number", value=1), md=1),
            dbc.Col(dbc.Input(id="stats-p2-min", type="number", value=0.1), md=2),
            dbc.Col(dbc.Input(id="stats-p2-max", type="number", value=1.0), md=2),
            dbc.Col(dbc.Input(id="stats-p2-n", type="number", value=25), md=1),
        ], className="g-2"),
        dbc.Button("Tracer χ² 2D", id="stats-chi2-2d-go", color="info"),
        dbc.Row([dbc.Col(dcc.Graph(id="stats-chi2-2d"), md=12)])
    ]
)

# --------- Compose page ---------
app.layout = dbc.Container([
    html.H2("Hyperiso Dashboard", className="mt-3", style={"color": ACCENT}),
    global_controls,
    dcc.Tabs(id="tabs", value="tab-scan", children=[
        dcc.Tab(label="Scan", value="tab-scan"),
        dcc.Tab(label="Core", value="tab-core"),
        dcc.Tab(label="QCD", value="tab-qcd"),
        dcc.Tab(label="Wilson", value="tab-wilson"),
        dcc.Tab(label="Observables", value="tab-observ"),
        dcc.Tab(label="Statistiques", value="tab-stats"),
    ]),
    html.Div(id="tab-content")
], fluid=True, style={"background": DARK_BG})

# ---------------- Tab routing ----------------
@app.callback(Output("tab-content", "children"), Input("tabs", "value"))
def render_tab(tab):
    if tab == "tab-scan":
        return scan_layout
    if tab == "tab-core":
        return core_layout
    if tab == "tab-qcd":
        return qcd_layout
    if tab == "tab-wilson":
        return wilson_layout
    if tab == "tab-observ":
        return observ_layout
    if tab == "tab-stats":
        return stats_layout
    return html.Div()

# ================= Callbacks: Config =================
@app.callback(
    Output("cfg-marty-name", "style"), Output("cfg-marty-path", "style"),
    Input("cfg-flags", "value")
)
def toggle_marty(flags):
    use = "USE_MARTY" in (flags or [])
    style = {"display": "block"} if use else {"display": "none"}
    return style, style

@app.callback(
    Output("cfg-status", "children"),
    Input("cfg-apply", "n_clicks"),
    State("cfg-model", "value"), State("cfg-flags", "value"),
    State("cfg-marty-name", "value"), State("cfg-marty-path", "value"),
    State("cfg-lha-path", "value"),
    prevent_initial_call=True
)
def apply_config(_, model_name, flags, mname, mpath, lha_path):
    try:
        # Build PyConfig (user's dataclass)
        from pyhyperiso.core.Core.Config import PyConfig, ExternalFlag
        fset = {
            ExternalFlag.IS_LHA_SPECTRUM: ("IS_LHA_SPECTRUM" in (flags or [])),
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            # ExternalFlag.USE_MARTY: ("USE_MARTY" in (flags or [])),
        }
        cfg_kwargs = dict(flags=fset, model=getattr(Model, model_name))
        if fset[ExternalFlag.USE_MARTY]:
            if mname: cfg_kwargs["mty_model_name"] = mname
            if mpath: cfg_kwargs["mty_model_path"] = Path(mpath)
        config = PyConfig(**cfg_kwargs)
        CTX.init_or_switch(lha_path, config)
        return dbc.Alert(f"Config appliquée. LHA: {lha_path}", color="success", duration=4000)
    except Exception as e:
        return make_alert(f"Échec d'application de la config: {e}")

# ================= Callbacks: Scan =================

def _gather_scan_params(values_dict: Dict[str, List[Any]]) -> List[Dict[str, Any]]:
    params = []
    for i in range(MAX_PARAMS):
        name = values_dict.get("name", [None]*MAX_PARAMS)[i]
        block = values_dict.get("block", [None]*MAX_PARAMS)[i]
        code = values_dict.get("code", [None]*MAX_PARAMS)[i]
        ptype = values_dict.get("ptype", [None]*MAX_PARAMS)[i]
        vmin = values_dict.get("min", [None]*MAX_PARAMS)[i]
        vmax = values_dict.get("max", [None]*MAX_PARAMS)[i]
        steps = values_dict.get("steps", [None]*MAX_PARAMS)[i]
        if name and block and code is not None and ptype and vmin is not None and vmax is not None and steps:
            params.append(dict(name=name, block=str(block), code=int(code), ptype=ptype, min=float(vmin), max=float(vmax), steps=int(steps)))
    return params

@app.callback(
    Output("scan-params-store", "data"),
    [Input({"type": "scan-name", "index": ALL}, "value"),
     Input({"type": "scan-block", "index": ALL}, "value"),
     Input({"type": "scan-code", "index": ALL}, "value"),
     Input({"type": "scan-ptype", "index": ALL}, "value"),
     Input({"type": "scan-min", "index": ALL}, "value"),
     Input({"type": "scan-max", "index": ALL}, "value"),
     Input({"type": "scan-steps", "index": ALL}, "value")]
)
def gather_params(names, blocks, codes, ptypes, mins, maxs, steps):
    values = {"name": names, "block": blocks, "code": codes, "ptype": ptypes, "min": mins, "max": maxs, "steps": steps}
    return _gather_scan_params(values)

@app.callback(
    Output("scan-summary", "children"), Output("scan-hist", "figure"), Output("scan-pairs", "figure"),
    Input("scan-preview", "n_clicks"),
    State("scan-params-store", "data"), State("scan-mode", "value"), State("scan-n-samples", "value"),
    prevent_initial_call=True
)
def preview_scan(_, params, mode, n_samples):
    try:
        if not params:
            return make_alert("Aucun paramètre valide."), go.Figure(), go.Figure()
        if mode == "grid":
            df = LhaScanPlanner.grid_from_ranges(params)
            nfiles = len(df)
            if nfiles > MAX_FILES:
                return make_alert(f"La grille produit {nfiles} fichiers (> {MAX_FILES}). Réduisez steps."), go.Figure(), go.Figure()
        else:
            df = LhaScanPlanner.random_from_ranges(params, int(n_samples))
            nfiles = len(df)
        # Summary + plots
        desc = df.describe().T
        fig_hist = stylize(px.histogram(df, barmode="overlay", template=PLOT_TEMPLATE), "Distribution des paramètres")
        fig_pairs = stylize(px.scatter_matrix(df, template=PLOT_TEMPLATE), "Matrice de dispersion")
        return [dbc.Alert(f"Aperçu: {nfiles} points.", color="info"), dbc.Table.from_dataframe(desc.round(4), striped=True, bordered=True, hover=True)], fig_hist, fig_pairs
    except Exception as e:
        return make_alert(f"Erreur dans l'aperçu: {e}"), go.Figure(), go.Figure()

@app.callback(
    Output("scan-progress", "value"), Output("scan-progress", "label"), Output("scan-progress", "color"),
    Input("scan-generate", "n_clicks"),
    State("scan-params-store", "data"), State("scan-mode", "value"), State("scan-n-samples", "value"), State("scan-outdir", "value"), State("cfg-lha-path", "value"),
    prevent_initial_call=True
)
def generate_scan(_, params, mode, n_samples, outdir, base_lha):
    try:
        base_text = Path(base_lha).read_text()
        # Build dataframe
        if mode == "grid":
            df = LhaScanPlanner.grid_from_ranges(params)
            total = len(df)
        else:
            df = LhaScanPlanner.random_from_ranges(params, int(n_samples))
            total = len(df)
        # Write files in chunks (progressive update)
        outdir = Path(outdir)
        chk = max(1, total // 100)
        for i, row in enumerate(df.itertuples(index=False)):
            updates = []
            for r, val in zip(params, row):
                updates.append((r["block"], r["code"], float(val)))
            fname = outdir / f"scan_{i:06d}.lha"
            LhaIO.write_lha(base_text, fname, updates)
            if i % chk == 0:
                pct = int(100 * (i + 1) / total)
                yield pct, f"{pct}% ({i+1}/{total})", "success"
        yield 100, f"100% ({total}/{total}) — Terminé", "success"
    except Exception as e:
        yield 0, f"Erreur: {e}", "danger"

# ================= Core inspector =================
@app.callback(
    Output("core-blocks-bar", "figure"), Output("core-blocks-donut", "figure"), Output("core-block-preview", "figure"), Output("core-flags", "children"), Output("core-paths", "children"),
    Input("cfg-apply", "n_clicks"), Input("core-block-filter", "value"),
    prevent_initial_call=True
)
def core_refresh(_, filter_block):
    CTX.ensure_ready()
    try:
        api = CTX.api
        if api is None:
            # Fallback: parse current LHA file
            text = Path(CTX.lha_path).read_text()
            blocks = LhaIO._parse_blocks(text)
            counts = {k: len(v) for k, v in blocks.items()}
            bar = stylize(go.Figure(go.Bar(x=list(counts.keys()), y=list(counts.values()))), "Nombre d'entrées par block")
            bar.update_layout(template=PLOT_TEMPLATE, title="Blocks (fallback)")
            donut_fig = donut(list(counts.keys()), list(counts.values()), "Répartition des blocks")
            prev = go.Figure()
            if filter_block and filter_block.upper() in blocks:
                bk = blocks[filter_block.upper()]
                prev = go.Figure(go.Bar(x=list(map(str, bk.keys())), y=list(bk.values())))
                prev.update_layout(template=PLOT_TEMPLATE, title=f"{filter_block} — valeurs")
            flags_div = dbc.Alert("APIAdapter indisponible, fallback LHA parser.", color="warning")
            paths_div = html.Div()
            return bar, donut_fig, prev, flags_div, paths_div
        # With API
        blocks = [str(b) for b in api.get_all_blocks() if str(b)]
        counts = {}
        for b in blocks:
            try:
                infos = api.get_block_infos(b)
                counts[b] = len(infos)
            except Exception:
                counts[b] = 0
        bar = go.Figure(go.Bar(x=list(counts.keys()), y=list(counts.values())))
        bar.update_layout(template=PLOT_TEMPLATE, title="Nombre d'entrées par block", xaxis_tickangle=45)
        donut_fig = donut(list(counts.keys()), list(counts.values()), "Diagramme fromage — blocks")
        prev = go.Figure()
        if filter_block:
            try:
                ptypes = api.get_type_of_block(str(filter_block))
                dtype = ptypes[0].name if ptypes else "?"
                infos = api.get_block_infos(str(filter_block))
                x = [str(k) for k in infos.keys()]
                y = [float(v.real()) for v in infos.values()]
                prev = stylize(go.Figure(go.Bar(x=x, y=y)), f"{filter_block} ({dtype}) — 1ères valeurs")
                prev.update_layout(template=PLOT_TEMPLATE, title=f"{filter_block} ({dtype}) — 1ères valeurs")
            except Exception:
                pass
        # flags & paths
        from pyhyperiso.core.Core.Config import ExternalFlag as EF
        flags = [f"{ef.name}={api.check_flag(ef)}" for ef in EF]
        lha_path = api.get_path(api.__class__.__dict__.get('APIPath', None) or type('X',(object,),{}) ) if False else CTX.lha_path
        return bar, donut_fig, prev, dbc.Alert(" ".join(flags), color="info"), dbc.Alert(f"LHA: {CTX.lha_path}", color="secondary")
    except Exception as e:
        empty = go.Figure()
        empty.update_layout(template=PLOT_TEMPLATE)
        return empty, empty, empty, make_alert(f"Core error: {e}"), html.Div()

# ================= QCD panel =================
@app.callback(
    Output("qcd-alphas-fig", "figure"), Output("qcd-mass-fig", "figure"), Output("qcd-constants", "children"),
    Input("qcd-scale-min", "value"), Input("qcd-scale-max", "value"), Input("qcd-n", "value"),
    Input("qcd-mass-pdg", "value"), Input("qcd-mtype", "value")
)
def qcd_compute(mu_min, mu_max, n, pdg, mtype):
    CTX.ensure_ready()
    try:
        # Try user's QCD provider
        try:
            from pyhyperiso.phyperiso.pyhyperiso.core import QCDProvider as _CppQCDProvider
            class PyQCDProvider:
                def __init__(self):
                    self._cpp_obj = _CppQCDProvider()
                def get_alphas(self, cfg: PyAlphasConfig):
                    return self._cpp_obj.compute_alphas(cfg.to_cpp())
                def get_qcd_masses(self, cfg: PyMassConfig):
                    return self._cpp_obj.compute_mass(cfg.to_cpp())
                def get_qcd_constants(self):
                    return self._cpp_obj.get_constants()
        except Exception:
            PyQCDProvider = None
        if PyQCDProvider is None:
            return go.Figure(), go.Figure(), make_alert("QCD provider indisponible.", "warning")
        provider = PyQCDProvider()
        xs = np.linspace(float(mu_min), float(mu_max), int(n))
        al = [provider.get_alphas(PyAlphasConfig(x)).value if hasattr(provider.get_alphas(PyAlphasConfig(x)), 'value') else provider.get_alphas(PyAlphasConfig(x)) for x in xs]
        # Mass curve
        mcfg = PyMassConfig(scale=1.0, pdg_id=int(pdg), m_b_type=MassType[mtype], m_t_type=MassType[mtype])
        ms = []
        for x in xs:
            mcfg.scale = float(x)
            val = provider.get_qcd_masses(mcfg)
            ms.append(val.value if hasattr(val, 'value') else val)
        fig_a = stylize(go.Figure(go.Scatter(x=xs, y=al, mode="lines")), "α_s(μ)")
        fig_a.update_layout(template=PLOT_TEMPLATE, title="α_s(μ)")
        fig_m = stylize(go.Figure(go.Scatter(x=xs, y=ms, mode="lines")), f"m(PDG {pdg}) vs μ ({mtype})")
        fig_m.update_layout(template=PLOT_TEMPLATE, title=f"m(PDG {pdg}) vs μ ({mtype})")
        consts = provider.get_qcd_constants()
        consts_div = dbc.Alert(f"Nc={getattr(consts, 'Nc', '?')} C_F={getattr(consts, 'C_F', '?')}", color="secondary")
        return fig_a, fig_m, consts_div
    except Exception as e:
        empty = go.Figure(); empty.update_layout(template=PLOT_TEMPLATE)
        return empty, empty, make_alert(f"QCD error: {e}")

# ================= Wilson =================

def _wilson_interface(build_group: WGroup, basis: WilsonBasis, match_scale: float, had_scale: float) -> PyWilsonInterface:
    wi = PyWilsonInterface()
    cfg = PyWilsonBuildConfig(groups={build_group}, matching_scale=float(match_scale), hadronic_scale=float(had_scale), order=QCDOrder.NNLO)
    wi.build(cfg)
    return wi

@app.callback(
    Output("wilson-stacked", "figure"), Output("wilson-heat", "figure"),
    Input("wilson-group", "value"), Input("wilson-coeff", "value"), Input("wilson-basis", "value"),
    Input("wilson-mode", "value"), Input("wilson-scale", "value"), Input("wilson-contrib", "value"), Input("wilson-orders", "value")
)
def wilson_plots(gname, wname, bname, mode, scale, contribs, orders):
    CTX.ensure_ready()
    try:
        group = getattr(WGroup, gname)
        coeff = getattr(WCoeff, wname)
        basis = getattr(WilsonBasis, bname)
        # Build interface once per request (fast) at fixed scales
        wi = _wilson_interface(group, basis, match_scale=scale, had_scale=scale)
        # separate orders (with αs prefactors if needed)
        orders_en = [getattr(QCDOrder, o) for o in orders]
        contribs_en = [getattr(ContributionType, c) for c in contribs]
        # αs factor for NLO/NNLO — use hadronic scale for running, matching for matching
        try:
            from pyhyperiso.phyperiso.pyhyperiso.core import QCDProvider as _CppQCDProvider
            class _QP:  # quick shim
                def __init__(self): self._cpp = _CppQCDProvider()
                def alphas(self, mu): return self._cpp.compute_alphas(PyAlphasConfig(mu).to_cpp())
            alpha_s = _QP().alphas(scale)
            alpha_val = getattr(alpha_s, 'value', alpha_s)
        except Exception:
            alpha_val = 0.22  # fallback typical value
        fac = {QCDOrder.LO: 1.0, QCDOrder.NLO: alpha_val/(4*math.pi), QCDOrder.NNLO: (alpha_val/(4*math.pi))**2}
        # Build stacked bars per contribution
        bars = []
        for c in contribs_en:
            ys = []
            for o in orders_en:
                if mode == "running":
                    sep = wi.get_sep_order_running(group, coeff, c, basis)
                else:
                    sep = wi.get_sep_order_matching(group, coeff, c)
                val = float(sep.get(o).real()) if hasattr(sep.get(o), 'real') else float(sep.get(o))
                ys.append(val * fac[o])
            bars.append(go.Bar(name=c.name, x=[o.name for o in orders_en], y=ys))
        fig_bar = stylize(go.Figure(data=bars), f"{mode.title()} — {group.name}/{coeff.name} ({basis.name})")
        fig_bar.update_layout(template=PLOT_TEMPLATE, barmode="group", title=f"{mode.title()} — {group.name}/{coeff.name} ({basis.name})")
        # Heatmap orders x contribution of full R/M values
        z = []
        for c in contribs_en:
            row = []
            for o in orders_en:
                req = PyWilsonRequest(group, coeff, o, c, ScaleType.HADRONIC if mode=="running" else ScaleType.MATCHING, sum_qcd_orders=True)
                val = wi.get_FR(req) if mode == "running" else wi.get_FM(req)
                row.append(float(val.real()) if hasattr(val, 'real') else float(val))
            z.append(row)
        fig_heat = stylize(go.Figure(data=go.Heatmap(z=z, x=[o.name for o in orders_en], y=[c.name for c in contribs_en], colorbar=dict(title="value"))))
        fig_heat.update_layout(template=PLOT_TEMPLATE, title="Full coefficients (FR/FM)")
        return fig_bar, fig_heat
    except Exception as e:
        empty = go.Figure(); empty.update_layout(template=PLOT_TEMPLATE)
        return empty, make_alert(f"Wilson error: {e}")

@app.callback(
    Output("wilson-var-fig", "figure"),
    Input("wilson-var-go", "n_clicks"),
    State("wilson-group", "value"), State("wilson-coeff", "value"), State("wilson-basis", "value"), State("wilson-mode", "value"), State("wilson-scale", "value"),
    State("wilson-var-ptype", "value"), State("wilson-var-block", "value"), State("wilson-var-code", "value"), State("wilson-var-min", "value"), State("wilson-var-max", "value"), State("wilson-var-n", "value"),
    prevent_initial_call=True
)
def wilson_var(_, gname, wname, bname, mode, scale, ptype, block, code, vmin, vmax, n):
    CTX.ensure_ready()
    try:
        if PyParameterSetter is None:
            raise RuntimeError("PyParameterSetter indisponible")
        group = getattr(WGroup, gname)
        coeff = getattr(WCoeff, wname)
        basis = getattr(WilsonBasis, bname)
        wi = _wilson_interface(group, basis, match_scale=scale, had_scale=scale)
        xs = np.linspace(float(vmin), float(vmax), int(n))
        ys = []
        setter = PyParameterSetter()
        pid = PyParamId(getattr(ParameterType, ptype), str(block), int(code))
        for x in xs:
            setter.mutate(pid, float(x))
            req = PyWilsonRequest(group, coeff, QCDOrder.NNLO, ContributionType.TOTAL, ScaleType.HADRONIC if mode=="running" else ScaleType.MATCHING, sum_qcd_orders=True)
            val = wi.get_FR(req) if mode == "running" else wi.get_FM(req)
            ys.append(float(val.real()) if hasattr(val, 'real') else float(val))
        fig = stylize(go.Figure(go.Scatter(x=xs, y=ys, mode="lines+markers")), f"Variation {coeff.name} vs {block}[{code}]")
        fig.update_layout(template=PLOT_TEMPLATE, title=f"Variation {coeff.name} vs {block}[{code}]")
        return fig
    except Exception as e:
        empty = go.Figure(); empty.update_layout(template=PLOT_TEMPLATE)
        return empty

# ================= Observables =================
@app.callback(
    Output("obs-fig", "figure"),
    Input("obs-calc", "n_clicks"), State("obs-select", "value"),
    prevent_initial_call=True
)
def obs_compute(_, names):
    CTX.ensure_ready()
    try:
        if PyObservableInterface is None:
            raise RuntimeError("ObservableInterface indisponible")
        if not names:
            return go.Figure()
        oi = PyObservableInterface()
        for n in names:
            oi.add_observable(getattr(Observables, n), QCDOrder.NNLO, True)
        vals = {}
        for n in names:
            s = oi.compute_observable(getattr(Observables, n))
            vals[n] = float(s.real()) if hasattr(s, 'real') else float(s)
        fig = stylize(go.Figure(), "Observables (valeurs actuelles)")
        for k, v in vals.items():
            fig.add_trace(go.Bar(x=[k], y=[v], name=k))
        fig.update_layout(template=PLOT_TEMPLATE, title="Observables (valeurs actuelles)")
        return fig
    except Exception as e:
        empty = go.Figure(); empty.update_layout(template=PLOT_TEMPLATE)
        return empty

@app.callback(
    Output("obs-var-fig", "figure"),
    Input("obs-var-go", "n_clicks"),
    State("obs-var-name", "value"), State("obs-var-ptype", "value"), State("obs-var-block", "value"), State("obs-var-code", "value"),
    State("obs-var-min", "value"), State("obs-var-max", "value"), State("obs-var-n", "value"), State("obs-opts", "value"),
    prevent_initial_call=True
)
def obs_var(_, obs_name, ptype, block, code, vmin, vmax, n, opts):
    CTX.ensure_ready()
    try:
        if PyObservableInterface is None:
            raise RuntimeError("ObservableInterface indisponible")
        show_unc = "uncert" in (opts or [])
        show_lead = "leading" in (opts or [])
        xs = np.linspace(float(vmin), float(vmax), int(n))
        ys = []
        uncs = []
        setter = PyParameterSetter() if PyParameterSetter else None
        pid = PyParamId(getattr(ParameterType, ptype), str(block), int(code)) if PyParamId else None
        oi = PyObservableInterface()
        oi.add_observable(getattr(Observables, obs_name), QCDOrder.NNLO, True)
        for x in xs:
            if setter and pid:
                setter.mutate(pid, float(x))
            val = oi.compute_observable(getattr(Observables, obs_name))
            v = float(val.real()) if hasattr(val, 'real') else float(val)
            ys.append(v)
            if show_unc:
                u = oi.compute_uncertainty(getattr(Observables, obs_name))
                u = float(u.real()) if hasattr(u, 'real') else float(u)
                uncs.append(u)
        fig = stylize(go.Figure(), f"{obs_name} — variation")
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", name=obs_name))
        if show_unc and uncs:
            lo = np.array(ys) - np.array(uncs)
            hi = np.array(ys) + np.array(uncs)
            fig.add_trace(go.Scatter(x=xs, y=hi, mode="lines", line=dict(width=0), showlegend=False))
            fig.add_trace(go.Scatter(x=xs, y=lo, mode="lines", fill="tonexty", line=dict(width=0), name="±σ"))
        if show_lead:
            lead = oi.compute_leading_uncertainties()
            # lead is assumed dict-like; display in title
            fig.update_layout(title=f"{obs_name} — leading sources: {str(lead)[:120]}...")
        fig.update_layout(template=PLOT_TEMPLATE)
        return fig
    except Exception as e:
        empty = go.Figure(); empty.update_layout(template=PLOT_TEMPLATE)
        return empty

# ================= Stats =================
# @app.callback(
#     Output("stats-corr", "figure"), Output("stats-chi2-1d", "figure"),
#     Input("stats-configure", "n_clicks"), State("stats-obs", "value"),
#     prevent_initial_call=True
# )
# def stats_setup(_, names):
#     try:
#         if PyObservableInterface is None:
#             raise RuntimeError("ObservableInterface indisponible")
#         oi = PyObservableInterface()
#         if not names:
#             return go.Figure(), go.Figure()
#         for n in names:
#             oi.add_observa


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=int(os.environ.get("PORT", 8050)))