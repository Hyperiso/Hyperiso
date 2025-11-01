# -*- coding: utf-8 -*-
"""
Hyperiso Plotly Dash — Multipage Dashboard (LHA • Wilson • Observables • Stats)
-----------------------------------------------------------------------------

Notes importantes:
- HyperisoMaster reste **global** et unique (singleton de fait côté C++
  et tenu en vie par ce module).
- On **switch** simplement de LHA ou de modèle via `hyp.switch_lha(...)`.
- L'application utilise une navigation multipage (URL) dans **un seul fichier**
  pour faciliter le déploiement. Aucune structure /pages n'est nécessaire.
- De nombreuses figures sont **mémorisées** avec Flask-Caching pour accélérer.

Dépendances: dash>=2.16, plotly, flask-caching
Lancez:  `python dash_hyperiso_multipage.py`
Ouvre:   http://127.0.0.1:8055
"""
from __future__ import annotations
import json
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Tuple, Set

import dash
from dash import dcc, html, Input, Output, State, MATCH, ALL, ctx
import plotly.graph_objects as go
import plotly.express as px
from flask_caching import Cache

# ------------------------------
# Imports Hyperiso (robustes)
# ------------------------------
try:
    # Chemins "anciens" vus dans votre main initial
    from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster
    from pyhyperiso.core.Core.Config import PyConfig, ExternalFlag
    from pyhyperiso.core.Core.ParameterSetter import PyParameterSetter, PyParamId, ParameterType
    from pyhyperiso.core.PhysicalModel.WilsonInterface import (
        PyWilsonInterface,
        PyWilsonBuildConfig,
        WGroup,
        WCoeff,
        QCDOrder,
        ContributionType,
        PyWilsonRequest,
    )
    from pyhyperiso.core.Core.APIAdapter import PyAPIAdapter
except Exception:
    # Fallback vers wrappers modulaires
    from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster
    from pyhyperiso.core.Core.Config import PyConfig, ExternalFlag
    from pyhyperiso.core.Core.ParameterSetter import PyParameterSetter, PyParamId, ParameterType
    from pyhyperiso.core.Common.Configs import (
        PyWilsonBuildConfig,
        PyWilsonRequest,
    )
    from pyhyperiso.core.PhysicalModel.WilsonInterface import (
        PyWilsonInterface,
    )
    from pyhyperiso.core.Common.GeneralEnum import (
        WGroup,
        WCoeff,
        QCDOrder,
        ContributionType,
    )
    from pyhyperiso.core.Core.APIAdapter import PyAPIAdapter

# Observables + Corrélations
try:
    from pyhyperiso.core.BusinessLogic.ObservableInterface import PyObservableInterface
except Exception:
    from pyhyperiso.core.BusinessLogic.ObservableInterface import PyObservableInterface

try:
    from pyhyperiso.core.Core.CorrelationProvider import PyCorrelationProvider, CorrelationType
except Exception:
    from pyhyperiso.core.Core.CorrelationProvider import PyCorrelationProvider, CorrelationType

# Enums supplémentaires
try:
    from pyhyperiso.core.Common.GeneralEnum import Observables, Decays, WilsonBasis
except Exception:
    # Valeur par défaut si l'enum WilsonBasis est absente
    from pyhyperiso.core.Common.GeneralEnum import Observables, Decays
    class WilsonBasis:  # type: ignore
        STANDARD = 0
        TRADITIONAL = 1

# Maths / Scalar
try:
    from pyhyperiso.core.Math.scalar import Scalar
except Exception:
    from pyhyperiso.core.Math.scalar import Scalar  # chemins déjà robustes

# ------------------------------
# Global: Hyperiso singleton & interfaces
# ------------------------------
last_model_name = "THDM"
last_lha_path = "lha/testinput_thdm.lha"

hyp = PyHyperisoMaster()
# config = PyConfig(
#     flags={
#         ExternalFlag.IS_LHA_SPECTRUM: False,
#         ExternalFlag.HAS_WILSON_INPUT: False,
#         ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
#         ExternalFlag.USE_MARTY: False,
#     },
#     model=getattr(type("_Model", (), {}), "THDM", None) or None,  # modèle défini dans PyConfig côté C++
# )
# Si votre PyConfig attend un enum Model, on garde exactement votre snippet initial
from pyhyperiso.core.Common.GeneralEnum import Model as _ModelEnum
config = PyConfig(
    flags={
        ExternalFlag.IS_LHA_SPECTRUM: False,
        ExternalFlag.HAS_WILSON_INPUT: False,
        ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
        # ExternalFlag.USE_MARTY: False,
    },
    model=_ModelEnum.THDM,
)

hyp.init(lha_file=last_lha_path, config=config)

# Interfaces persistantes
wilson_interface = PyWilsonInterface()
observable_interface = PyObservableInterface()
correlation_provider = PyCorrelationProvider()
api_adapter = PyAPIAdapter()

# Build initial Wilson config (B, B', B_scalar)
DEFAULT_WILSON_GROUPS: Set[WGroup] = {WGroup.B, WGroup.BPrime, WGroup.BScalar}
DEFAULT_MATCH = 81.0
DEFAULT_HAD = 4.18
DEFAULT_ORDER = QCDOrder.NNLO

wilson_build_cfg = PyWilsonBuildConfig(
    groups=DEFAULT_WILSON_GROUPS,
    matching_scale=DEFAULT_MATCH,
    hadronic_scale=DEFAULT_HAD,
    order=DEFAULT_ORDER,
)
wilson_interface.build(wilson_build_cfg)

# ------------------------------
# Coeff sets par groupe
# ------------------------------
COEFFS_BY_GROUP: Dict[WGroup, List[WCoeff]] = {
    WGroup.B: [getattr(WCoeff, f"C{i}") for i in range(1, 11)],
    WGroup.BPrime: [getattr(WCoeff, f"CP{i}") for i in range(1, 11)] + [WCoeff.CPQ1, WCoeff.CPQ2],
    WGroup.BScalar: [WCoeff.CQ1, WCoeff.CQ2],
}

# ------------------------------
# Dash app + cache
# ------------------------------
app = dash.Dash(__name__, title="Hyperiso Dashboard", suppress_callback_exceptions=True)
server = app.server
cache = Cache(server, config={"CACHE_TYPE": "SimpleCache", "CACHE_DEFAULT_TIMEOUT": 60})

DARK = dict(plot_bgcolor="#1c1c1c", paper_bgcolor="#121212", font_color="white")
CARD_STYLE = {
    "background": "#1c1c1c",
    "borderRadius": "12px",
    "padding": "14px",
    "boxShadow": "0 2px 12px rgba(0,0,0,.25)",
}

# ------------------------------
# Utils
# ------------------------------

def _safe_scalar(x: Scalar | float | complex) -> float:
    try:
        if isinstance(x, Scalar):
            return float(x.to_double()) if hasattr(x, "to_double") else float(x.real())
        return float(x)
    except Exception:
        try:
            return x.real() if hasattr(x, "real") else 0.0
        except Exception:
            return 0.0


def _switch_lha_if_needed(model_name: str, lha_path: str):
    """Change de LHA ou de modèle si nécessaire, sans tuer le serveur.
    Conserve l'instance globale `hyp`.
    """
    global last_model_name, last_lha_path
    if model_name != last_model_name or lha_path != last_lha_path:
        from pyhyperiso.core.Common.GeneralEnum import Model
        cfg_args = dict(
            flags={
                ExternalFlag.IS_LHA_SPECTRUM: False,
                ExternalFlag.HAS_WILSON_INPUT: False,
                ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
                # ExternalFlag.USE_MARTY: False,
            },
            model=getattr(Model, model_name),
        )
        new_cfg = PyConfig(**cfg_args)
        hyp.switch_lha(lha_file=lha_path, config=new_cfg)
        last_model_name = model_name
        last_lha_path = lha_path


@dataclass
class Settings:
    model: str = "THDM"
    custom_name: str | None = None
    custom_path: str | None = None
    lha_path: str = last_lha_path
    matching_scale: float = DEFAULT_MATCH
    hadronic_scale: float = DEFAULT_HAD
    qcd_order: str = QCDOrder.NNLO.name
    basis: str = getattr(WilsonBasis, "STANDARD").__class__.__name__ if not hasattr(WilsonBasis, "name") else WilsonBasis.STANDARD.name
    contribs: List[str] = ("SM", "BSM", "TOTAL")


# ------------------------------
# Cache helpers
# ------------------------------
@cache.memoize()
def cached_lha_summary(lha_path: str) -> Dict:
    """Résumé LHA (blocs, tailles, types...)"""
    blocks = [str(b) for b in api_adapter.get_all_blocks()]  # set -> list
    info = {}
    type_of = {}
    sample_vals = {}
    for b in blocks:
        try:
            # type de block (SM/BSM/WILSON/...). On prend le premier si plusieurs
            b_types = api_adapter.get_type_of_block(b)
            type_of[b] = [t.name for t in b_types]
            # quelques clés/id présents
            b_infos = api_adapter.get_block_infos(b, b_types[0] if b_types else ParameterType.SM)
            keys = list(b_infos.keys())
            info[b] = len(keys)
            # échantillon de 3 premières valeurs
            sample_vals[b] = [str(k) for k in keys[:3]]
        except Exception:
            info[b] = 0
            type_of[b] = ["UNKWN"]
            sample_vals[b] = []
    return dict(blocks=blocks, counts=info, types=type_of, samples=sample_vals)


@cache.memoize()
def cached_wilson_values(settings_json: str, contrib: str, group_name: str) -> Tuple[List[str], List[float]]:
    s = Settings(**json.loads(settings_json))
    # Applique les échelles
    wilson_interface.set_matching_scale(float(s.matching_scale))
    wilson_interface.set_hadronic_scale(float(s.hadronic_scale))

    # Ordre + contribution
    order = getattr(QCDOrder, s.qcd_order)
    cont = getattr(ContributionType, contrib)
    basis = getattr(WilsonBasis, s.basis) if hasattr(WilsonBasis, s.basis) else WilsonBasis.STANDARD

    group = getattr(WGroup, group_name)
    coeffs = COEFFS_BY_GROUP[group]
    xs = [c.name for c in coeffs]
    ys: List[float] = []
    for c in coeffs:
        val = wilson_interface.get_FR(
            PyWilsonRequest(group, c, order, cont, wilson_basis=basis)
        )
        ys.append(_safe_scalar(val))
    return xs, ys


@cache.memoize()
def cached_wilson_order_matrix(settings_json: str, group_name: str, coeff_name: str) -> Dict[str, List[float]]:
    s = Settings(**json.loads(settings_json))
    wilson_interface.set_matching_scale(float(s.matching_scale))
    wilson_interface.set_hadronic_scale(float(s.hadronic_scale))

    basis = getattr(WilsonBasis, s.basis) if hasattr(WilsonBasis, s.basis) else WilsonBasis.STANDARD
    group = getattr(WGroup, group_name)
    coeff = getattr(WCoeff, coeff_name)
    orders = [QCDOrder.LO, QCDOrder.NLO, QCDOrder.NNLO]
    out = {}
    for cont in [ContributionType.SM, ContributionType.BSM, ContributionType.TOTAL]:
        out[cont.name] = [
            _safe_scalar(
                wilson_interface.get_FR(
                    PyWilsonRequest(group, coeff, o, cont, wilson_basis=basis)
                )
            )
            for o in orders
        ]
    return out


@cache.memoize()
def cached_scale_scan(settings_json: str, group_name: str, coeff_name: str, scale_kind: str, start: float, stop: float, steps: int) -> Tuple[List[float], Dict[str, List[float]]]:
    s = Settings(**json.loads(settings_json))
    basis = getattr(WilsonBasis, s.basis) if hasattr(WilsonBasis, s.basis) else WilsonBasis.STANDARD
    group = getattr(WGroup, group_name)
    coeff = getattr(WCoeff, coeff_name)

    xs = []
    ys: Dict[str, List[float]] = {c.name: [] for c in [ContributionType.SM, ContributionType.BSM, ContributionType.TOTAL]}
    for i in range(steps + 1):
        t = start + (stop - start) * i / steps
        if scale_kind == "matching":
            wilson_interface.set_matching_scale(t)
        else:
            wilson_interface.set_hadronic_scale(t)
        xs.append(t)
        for cont in [ContributionType.SM, ContributionType.BSM, ContributionType.TOTAL]:
            val = wilson_interface.get_FR(
                PyWilsonRequest(group, coeff, QCDOrder.NNLO, cont, wilson_basis=basis)
            )
            ys[cont.name].append(_safe_scalar(val))
    return xs, ys


@cache.memoize()
def cached_observables_compute(obs_names: Tuple[str, ...], order_name: str) -> Dict[str, Dict[str, float]]:
    # Ajoute/Met à jour les observables sélectionnées
    for o in obs_names:
        try:
            observable_interface.add_observable(getattr(Observables, o), getattr(QCDOrder, order_name), True)
        except Exception:
            pass

    results: Dict[str, Dict[str, float]] = {}
    for o in obs_names:
        On = getattr(Observables, o)
        try:
            cv = observable_interface.compute_observable(On)
            sig = observable_interface.compute_uncertainty(On)
            lead = observable_interface.compute_leading_uncertainties(On, 5)
            results[o] = {
                "central": _safe_scalar(cv),
                "sigma": _safe_scalar(sig),
                "rel_err": (abs(float(sig)) / (abs(float(cv)) + 1e-16)),
                "_lead_json": json.dumps({str(k): _safe_scalar(v) for k, v in lead.items()}),
            }
        except Exception:
            results[o] = {"central": float("nan"), "sigma": float("nan"), "rel_err": float("nan"), "_lead_json": json.dumps({})}
    return results


@cache.memoize()
def cached_observable_correlation(obs_names: Tuple[str, ...], corr_type: str) -> List[List[float]]:
    Olist = [getattr(Observables, n) for n in obs_names]
    typ = getattr(CorrelationType, corr_type)
    n = len(Olist)
    M = [[0.0 for _ in range(n)] for __ in range(n)]
    for i in range(n):
        for j in range(n):
            try:
                M[i][j] = float(correlation_provider.correlation_from_observable(Olist[i], Olist[j], typ))
            except Exception:
                M[i][j] = 0.0
    return M


# ------------------------------
# UI Helpers
# ------------------------------

def nav_link(label: str, href: str) -> html.A:
    return html.A(
        label,
        href=href,
        style={
            "padding": "8px 12px",
            "marginRight": "8px",
            "borderRadius": "10px",
            "background": "#222",
            "color": "#fff",
            "textDecoration": "none",
        },
    )


def section(title: str, *children):
    return html.Div(
        [html.H3(title, style={"marginBottom": "8px"})] + list(children),
        style={"marginBottom": "18px", **CARD_STYLE},
    )


# ------------------------------
# Layout global
# ------------------------------
app.layout = html.Div(
    style={"background": "#121212", "minHeight": "100vh", "color": "#fff", "fontFamily": "Inter, Arial"},
    children=[
        dcc.Location(id="url"),
        html.Div(
            [
                html.Div(
                    [
                        html.H1("Hyperiso Dashboard", style={"fontSize": "26px", "margin": 0}),
                        html.Div(
                            [
                                nav_link("LHA", "/lha"),
                                nav_link("Wilson", "/wilson"),
                                nav_link("Observables", "/observables"),
                                nav_link("Stats", "/stats"),
                            ],
                            style={"marginTop": "8px"},
                        ),
                    ]
                ),
                html.Div(id="right-settings", style={"display": "flex", "gap": "12px", "alignItems": "center"}, children=[
                    html.Div("Modèle:"),
                    dcc.Dropdown(
                        id="model-selector",
                        options=[
                            {"label": "THDM", "value": "THDM"},
                            {"label": "SUSY", "value": "SUSY"},
                            {"label": "SM", "value": "SM"},
                            {"label": "CUSTOM", "value": "CUSTOM"},
                        ],
                        value="THDM",
                        clearable=False,
                        style={"width": 160, "color": "#111"},
                    ),
                    dcc.Input(id="custom-model-name", type="text", placeholder="Nom modèle", style={"display": "none", "width": 150}),
                    dcc.Input(id="custom-model-path", type="text", placeholder="Chemin modèle", style={"display": "none", "width": 220}),
                    dcc.Input(id="lha-path", type="text", value=last_lha_path, placeholder="LHA path", style={"width": 280}),
                    html.Div("μW"),
                    dcc.Input(id="match-scale", type="number", value=DEFAULT_MATCH, style={"width": 100}),
                    html.Div("μh"),
                    dcc.Input(id="had-scale", type="number", value=DEFAULT_HAD, style={"width": 100}),
                    dcc.Dropdown(id="order-select", options=[{"label": o.name, "value": o.name} for o in [QCDOrder.LO, QCDOrder.NLO, QCDOrder.NNLO]], value="NNLO", clearable=False, style={"width": 120, "color": "#111"}),
                    dcc.Dropdown(id="basis-select", options=[
                        {"label": "STANDARD", "value": "STANDARD"},
                        {"label": "TRADITIONAL", "value": "TRADITIONAL"},
                    ], value="STANDARD", clearable=False, style={"width": 160, "color": "#111"}),
                    dcc.Checklist(id="contrib-check", options=["SM", "BSM", "TOTAL"], value=["SM", "BSM", "TOTAL"], inline=True),
                    html.Button("Appliquer", id="apply", n_clicks=0, style={"background": "#4a7", "border": "none", "padding": "8px 14px", "borderRadius": "8px"}),
                ]),
            ],
            style={"display": "flex", "justifyContent": "space-between", "alignItems": "center", "padding": "14px 18px", "position": "sticky", "top": 0, "zIndex": 3, "background": "#111"},
        ),
        dcc.Store(id="store-settings"),
        dcc.Store(id="store-lha-summary"),
        html.Div(id="page-content", style={"padding": "16px"}),
    ],
)


# ------------------------------
# Pages
# ------------------------------

def page_lha() -> html.Div:
    return html.Div(
        children=[
            section(
                "Résumé LHA",
                html.Div(id="lha-summary-text"),
                html.Div(
                    [
                        html.Div(dcc.Graph(id="lha-block-counts"), style={"flex": 1}),
                        html.Div(dcc.Graph(id="lha-type-pie"), style={"flex": 1}),
                    ],
                    style={"display": "flex", "gap": "12px"},
                ),
                html.Div(
                    [
                        html.Div(dcc.Graph(id="lha-sunburst"), style={"flex": 1}),
                        html.Div(dcc.Graph(id="lha-mass-spectrum"), style={"flex": 1}),
                    ],
                    style={"display": "flex", "gap": "12px"},
                ),
            ),
        ]
    )


def page_wilson() -> html.Div:
    return html.Div(
        children=[
            section(
                "Coefficients de Wilson (Full Running)",
                html.Div(
                    [
                        html.Div(dcc.Graph(id=f"wilson-bars-{g.name}"), style={"flex": 1})
                        for g in [WGroup.B, WGroup.BPrime, WGroup.BScalar]
                    ],
                    style={"display": "flex", "gap": "12px", "flexWrap": "wrap"},
                ),
            ),
            section(
                "Ordre vs Contribution (coeff sélectionné)",
                html.Div(
                    [
                        dcc.Dropdown(
                            id="wilson-group-select",
                            options=[{"label": g.name, "value": g.name} for g in COEFFS_BY_GROUP.keys()],
                            value=WGroup.B.name,
                            clearable=False,
                            style={"width": 180, "color": "#111"},
                        ),
                        dcc.Dropdown(
                            id="wilson-coeff-select",
                            options=[{"label": c.name, "value": c.name} for c in COEFFS_BY_GROUP[WGroup.B]],
                            value="C9",
                            clearable=False,
                            style={"width": 180, "color": "#111", "marginLeft": "8px"},
                        ),
                    ],
                    style={"display": "flex", "alignItems": "center", "gap": "8px", "marginBottom": "8px"},
                ),
                html.Div(
                    [
                        html.Div(dcc.Graph(id="wilson-heatmap-ovc"), style={"flex": 1}),
                        html.Div(dcc.Graph(id="wilson-basis-compare"), style={"flex": 1}),
                    ],
                    style={"display": "flex", "gap": "12px"},
                ),
            ),
            section(
                "Scan d'échelle",
                html.Div(
                    [
                        dcc.RadioItems(
                            id="scale-kind",
                            options=[
                                {"label": "Matching μW", "value": "matching"},
                                {"label": "Hadronique μh", "value": "hadronic"},
                            ],
                            value="hadronic",
                            inline=True,
                        ),
                        dcc.Input(id="scale-start", type="number", value=1.0, style={"width": 110, "marginLeft": "12px"}),
                        dcc.Input(id="scale-stop", type="number", value=10.0, style={"width": 110, "marginLeft": "8px"}),
                        dcc.Input(id="scale-steps", type="number", value=30, style={"width": 110, "marginLeft": "8px"}),
                        html.Button("Scanner", id="scan-btn", n_clicks=0, style={"marginLeft": "8px"}),
                    ],
                    style={"display": "flex", "alignItems": "center"},
                ),
                dcc.Graph(id="wilson-scale-scan"),
            ),
        ]
    )


def page_observables() -> html.Div:
    return html.Div(
        children=[
            section(
                "Sélection d'observables",
                html.Div(
                    [
                        dcc.Dropdown(
                            id="obs-select",
                            options=[{"label": o.name, "value": o.name} for o in Observables],
                            value=[Observables.BR_B_XS_GAMMA.name, Observables.BR_BS_MUMU.name, Observables.R_D.name],
                            multi=True,
                            style={"color": "#111"},
                        ),
                        dcc.Dropdown(
                            id="obs-order",
                            options=[{"label": o.name, "value": o.value} for o in [QCDOrder.LO, QCDOrder.NLO, QCDOrder.NNLO]],
                            value="NNLO",
                            clearable=False,
                            style={"width": 160, "color": "#111", "marginLeft": "8px"},
                        ),
                        html.Button("Calculer", id="obs-run", n_clicks=0, style={"marginLeft": "8px"}),
                    ],
                    style={"display": "flex", "alignItems": "center", "gap": "8px"},
                ),
                html.Div(
                    [
                        html.Div(dcc.Graph(id="obs-bar"), style={"flex": 1}),
                        html.Div(dcc.Graph(id="obs-rel-err"), style={"flex": 1}),
                    ],
                    style={"display": "flex", "gap": "12px"},
                ),
            ),
            section(
                "Leading uncertainties & Corrélations",
                html.Div(
                    [
                        dcc.Slider(id="lead-top-n", min=1, max=10, step=1, value=5, marks={i: str(i) for i in range(1, 11)}, tooltip={"placement": "bottom"}, updatemode="drag"),
                        dcc.Dropdown(id="corr-type", options=[{"label": t.name, "value": t.name} for t in CorrelationType], value="COMBINED", clearable=False, style={"width": 200, "color": "#111", "marginLeft": "8px"}),
                    ],
                    style={"display": "flex", "alignItems": "center", "gap": "12px"},
                ),
                html.Div(
                    [
                        html.Div(dcc.Graph(id="obs-leading"), style={"flex": 1}),
                        html.Div(dcc.Graph(id="obs-corr"), style={"flex": 1}),
                    ],
                    style={"display": "flex", "gap": "12px"},
                ),
            ),
        ]
    )


def page_stats() -> html.Div:
    return html.Div(
        children=[
            section(
                "χ² global et résumé incertitudes",
                html.Div(
                    [
                        html.Button("Recalculer χ²", id="chi2-run", n_clicks=0),
                        html.Div(id="chi2-val", style={"marginLeft": "12px"}),
                    ],
                    style={"display": "flex", "alignItems": "center", "gap": "8px"},
                ),
                html.Div(
                    [
                        html.Div(dcc.Graph(id="chi2-indicator"), style={"flex": 1}),
                        html.Div(dcc.Graph(id="uncertainty-summary"), style={"flex": 2}),
                    ],
                    style={"display": "flex", "gap": "12px"},
                ),
            ),
        ]
    )


# ------------------------------
# Router
# ------------------------------
@app.callback(Output("page-content", "children"), Input("url", "pathname"))
def route(pathname: str):
    if pathname in ("/", "/lha"):  
        return page_lha()
    if pathname == "/wilson":
        return page_wilson()
    if pathname == "/observables":
        return page_observables()
    if pathname == "/stats":
        return page_stats()
    return html.Div("Page inconnue", style={"padding": 24})


# ------------------------------
# Toggle champs custom model
# ------------------------------
@app.callback(
    Output("custom-model-name", "style"),
    Output("custom-model-path", "style"),
    Input("model-selector", "value"),
)
def toggle_custom_fields(model_value):
    if model_value == "CUSTOM":
        return {"display": "block"}, {"display": "block"}
    return {"display": "none"}, {"display": "none"}


# ------------------------------
# Appliquer / stocker Settings + (re)build wilson
# ------------------------------
@app.callback(
    Output("store-settings", "data"),
    Output("store-lha-summary", "data"),
    Input("apply", "n_clicks"),
    State("model-selector", "value"),
    State("custom-model-name", "value"),
    State("custom-model-path", "value"),
    State("lha-path", "value"),
    State("match-scale", "value"),
    State("had-scale", "value"),
    State("order-select", "value"),
    State("basis-select", "value"),
    State("contrib-check", "value"),
    prevent_initial_call=True,
)
def apply_settings(nc, model_value, custom_name, custom_path, lha_path, match, had, order, basis, contribs):
    # switch lha / model si besoin
    _switch_lha_if_needed(str(model_value), str(lha_path))

    # rebuild Wilson si ordre/échelles changent
    try:
        cfg = PyWilsonBuildConfig(
            groups=DEFAULT_WILSON_GROUPS,
            matching_scale=float(match),
            hadronic_scale=float(had),
            order=getattr(QCDOrder, str(order)),
        )
        wilson_interface.build(cfg)
    except Exception:
        # fallback set scales only
        wilson_interface.set_matching_scale(float(match))
        wilson_interface.set_hadronic_scale(float(had))

    s = Settings(
        model=str(model_value),
        custom_name=custom_name,
        custom_path=custom_path,
        lha_path=str(lha_path),
        matching_scale=float(match),
        hadronic_scale=float(had),
        qcd_order=str(order),
        basis=str(basis),
        contribs=list(contribs or ["SM", "BSM", "TOTAL"]),
    )

    # résumé LHA cache
    lha_summary = cached_lha_summary(s.lha_path)
    return asdict(s), lha_summary


# ------------------------------
# LHA page callbacks
# ------------------------------
@app.callback(Output("lha-summary-text", "children"), Input("store-lha-summary", "data"))
def lha_summary_text(d):
    if not d:
        return html.Div("Cliquez sur Appliquer pour charger le LHA.")
    nblocks = len(d.get("blocks", []))
    return html.Div([
        html.P(f"Fichier: {last_lha_path}"),
        html.P(f"Blocs détectés: {nblocks}"),
    ])


@app.callback(Output("lha-block-counts", "figure"), Input("store-lha-summary", "data"))
def lha_block_counts(d):
    fig = go.Figure()
    if d:
        x = list(d["counts"].keys())
        y = list(d["counts"].values())
        fig.add_bar(x=x, y=y)
    fig.update_layout(title="Nombre d'entrées par bloc", **DARK)
    return fig


@app.callback(Output("lha-type-pie", "figure"), Input("store-lha-summary", "data"))
def lha_type_pie(d):
    labels, vals = [], []
    if d:
        # Compte le premier type reporté pour chaque bloc
        from collections import Counter
        c = Counter([ (v[0] if v else "UNKWN") for v in d.get("types", {}).values() ])
        labels, vals = list(c.keys()), list(c.values())
    fig = go.Figure(go.Pie(labels=labels, values=vals, hole=0.45))
    fig.update_layout(title="Répartition par type de bloc (premier type)", **DARK)
    return fig


@app.callback(Output("lha-sunburst", "figure"), Input("store-lha-summary", "data"))
def lha_sunburst(d):
    if not d:
        return go.Figure().update_layout(**DARK)
    parents, labels, values = [], [], []
    for b, n in d["counts"].items():
        t = (d["types"].get(b, ["UNKWN"]) or ["UNKWN"])[0]
        parents += [t]
        labels += [b]
        values += [n]
    fig = go.Figure(px.sunburst(names=labels, parents=parents, values=values))
    fig.update_layout(title="Hiérarchie Type → Bloc", **DARK)
    return fig


@app.callback(Output("lha-mass-spectrum", "figure"), Input("store-lha-summary", "data"))
def lha_mass_spectrum(d):
    # essaie de récupérer un bloc MASS s'il existe
    ys, xs = [], []
    try:
        b_types = api_adapter.get_type_of_block("MASS")
        if b_types:
            infos = api_adapter.get_block_infos("MASS", b_types[0])
            for lid, val in infos.items():
                ys.append(_safe_scalar(val))
                xs.append(str(lid))
    except Exception:
        pass
    fig = go.Figure()
    if xs:
        fig.add_bar(x=xs, y=ys)
    fig.update_layout(title="Spectre de masses (bloc MASS)", xaxis_title="LHA id", yaxis_title="M", **DARK)
    return fig


# ------------------------------
# Wilson page callbacks
# ------------------------------
@app.callback(
    [Output(f"wilson-bars-{g.name}", "figure") for g in [WGroup.B, WGroup.BPrime, WGroup.BScalar]],
    Input("store-settings", "data"),
)
def wilson_bars(settings):
    figs = []
    if not settings:
        return [go.Figure().update_layout(**DARK) for _ in range(3)]
    s_json = json.dumps(settings)
    contribs = settings.get("contribs", ["SM", "BSM", "TOTAL"]) or ["SM", "BSM", "TOTAL"]
    for g in [WGroup.B, WGroup.BPrime, WGroup.BScalar]:
        bars = []
        for c in contribs:
            x, y = cached_wilson_values(s_json, c, g.name)
            bars.append(go.Bar(name=c, x=x, y=y))
        fig = go.Figure(bars)
        fig.update_layout(title=f"{g.name} — FR ({', '.join(contribs)})", barmode="group", **DARK)
        figs.append(fig)
    return figs


@app.callback(
    Output("wilson-coeff-select", "options"),
    Input("wilson-group-select", "value"),
)
def refresh_coeff_options(gn):
    g = getattr(WGroup, gn)
    return [{"label": c.name, "value": c.name} for c in COEFFS_BY_GROUP[g]]


@app.callback(
    Output("wilson-heatmap-ovc", "figure"),
    Output("wilson-basis-compare", "figure"),
    Input("store-settings", "data"),
    Input("wilson-group-select", "value"),
    Input("wilson-coeff-select", "value"),
)
def wilson_heat_and_basis(settings, gname, cname):
    if not settings:
        return go.Figure().update_layout(**DARK), go.Figure().update_layout(**DARK)
    s_json = json.dumps(settings)
    mat = cached_wilson_order_matrix(s_json, gname, cname)
    z = [mat[k] for k in ["SM", "BSM", "TOTAL"]]
    fig_hm = go.Figure(
        data=go.Heatmap(z=z, x=[o.name for o in [QCDOrder.LO, QCDOrder.NLO, QCDOrder.NNLO]], y=["SM", "BSM", "TOTAL"])  # type: ignore
    )
    fig_hm.update_layout(title=f"{gname}:{cname} — Ordre vs Contribution", **DARK)

    # Basis compare STANDARD vs TRADITIONAL à NNLO TOTAL
    xs = ["STANDARD", "TRADITIONAL"]
    ys = []
    for b in xs:
        v = wilson_interface.get_FR(
            PyWilsonRequest(getattr(WGroup, gname), getattr(WCoeff, cname), QCDOrder.NNLO, ContributionType.TOTAL, wilson_basis=getattr(WilsonBasis, b) if hasattr(WilsonBasis, b) else WilsonBasis.STANDARD)
        )
        ys.append(_safe_scalar(v))
    fig_basis = go.Figure(go.Bar(x=xs, y=ys))
    fig_basis.update_layout(title=f"{gname}:{cname} — Comparaison de base", **DARK)
    return fig_hm, fig_basis


@app.callback(
    Output("wilson-scale-scan", "figure"),
    Input("scan-btn", "n_clicks"),
    State("store-settings", "data"),
    State("wilson-group-select", "value"),
    State("wilson-coeff-select", "value"),
    State("scale-kind", "value"),
    State("scale-start", "value"),
    State("scale-stop", "value"),
    State("scale-steps", "value"),
    prevent_initial_call=True,
)
def wilson_scan(_n, settings, gname, cname, kind, start, stop, steps):
    if not settings:
        return go.Figure().update_layout(**DARK)
    xs, ys = cached_scale_scan(json.dumps(settings), gname, cname, kind, float(start), float(stop), int(steps))
    fig = go.Figure()
    for k, arr in ys.items():
        fig.add_scatter(x=xs, y=arr, mode="lines+markers", name=k)
    fig.update_layout(title=f"Scan {kind} pour {gname}:{cname}", xaxis_title=f"μ ({kind})", **DARK)
    return fig


# ------------------------------
# Observables callbacks
# ------------------------------
@app.callback(
    Output("obs-bar", "figure"),
    Output("obs-rel-err", "figure"),
    Output("obs-leading", "figure"),
    Output("obs-corr", "figure"),
    Input("obs-run", "n_clicks"),
    State("obs-select", "value"),
    State("obs-order", "value"),
    State("lead-top-n", "value"),
    State("corr-type", "value"),
    prevent_initial_call=True,
)
def compute_observables(_n, obs_list, order_name, topn, corr_type):
    obs_list = tuple(obs_list or [])
    res = cached_observables_compute(obs_list, order_name)

    # Bar avec barres d'erreur
    xs, ys, err = [], [], []
    for k in obs_list:
        xs.append(k)
        ys.append(res[k]["central"])  # type: ignore
        err.append(res[k]["sigma"])   # type: ignore
    fig_bar = go.Figure(go.Bar(x=xs, y=ys, error_y=dict(type="data", array=err, visible=True)))
    fig_bar.update_layout(title="Observables — Valeurs centrales ± σ", **DARK)

    # Erreurs relatives
    fig_rel = go.Figure(go.Bar(x=xs, y=[res[k]["rel_err"] for k in obs_list]))  # type: ignore
    fig_rel.update_layout(title="Erreur relative (σ/|central|)", **DARK)

    # Leading uncertainties (affiche top-N pids pour chaque observable
    # concaténés en sous-traces)
    fig_lead = go.Figure()
    for k in obs_list:
        lead = json.loads(res[k]["_lead_json"])  # pid->val
        # top-n
        items = sorted(lead.items(), key=lambda kv: abs(kv[1]), reverse=True)[: int(topn)]
        if items:
            fig_lead.add_bar(x=[f"{k}: {pid}" for pid, _ in items], y=[v for _, v in items], name=k)
    fig_lead.update_layout(title=f"Leading uncertainties (top-{topn})", barmode="group", **DARK)

    # Corrélations
    M = cached_observable_correlation(obs_list, corr_type)
    fig_corr = go.Figure(data=go.Heatmap(z=M, x=list(obs_list), y=list(obs_list), zmin=-1, zmax=1, colorscale="RdBu"))
    fig_corr.update_layout(title=f"Corrélations ({corr_type})", **DARK)

    return fig_bar, fig_rel, fig_lead, fig_corr


# ------------------------------
# Stats callbacks
# ------------------------------
@app.callback(
    Output("chi2-val", "children"),
    Output("chi2-indicator", "figure"),
    Output("uncertainty-summary", "figure"),
    Input("chi2-run", "n_clicks"),
    State("obs-select", "value"),
    prevent_initial_call=True,
)
def stats_compute(_n, current_obs):
    # χ² global (interne à la DB expérimentale de votre C++)
    try:
        chi2 = float(observable_interface.compute_chi2())
    except Exception:
        chi2 = float("nan")

    text = html.Div([html.B("χ² global:"), html.Span(f" {chi2:.3f}")])

    # Indicateur
    ind = go.Figure(go.Indicator(mode="gauge+number", value=chi2, title={"text": "χ²"}))
    ind.update_layout(**DARK)

    # Résumé incertitudes (si obs sélectionnées)
    fig_unc = go.Figure()
    if current_obs:
        res = cached_observables_compute(tuple(current_obs), "NNLO")
        fig_unc.add_bar(x=list(current_obs), y=[res[k]["sigma"] for k in current_obs])  # type: ignore
    fig_unc.update_layout(title="Résumé des σ (sélection d'observables)", **DARK)

    return text, ind, fig_unc


# ------------------------------
# Run
# ------------------------------
if __name__ == "__main__":
    app.run(debug=True, port=8055)
