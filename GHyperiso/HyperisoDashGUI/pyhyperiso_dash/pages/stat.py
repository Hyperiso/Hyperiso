from __future__ import annotations

from dash import dcc, html

from pyhyperiso_dash import services as svc
from pyhyperiso_dash.components import card, data_table, dropdown, enum_options, field, graph, num_input, page_title, small_note, status_box, text_input
from pyhyperiso.core.Common.GeneralEnum import QCDOrder

SELECT_MODE = [
    {"label": "Observable by observable", "value": "observable"},
    {"label": "Whole decay", "value": "decay"},
]

DEFAULT_DECAY = ["B__l_l"]
DEFAULT_OBS = ["BR_BS_MUMU"]

OBS_COLUMNS = [
    {"name": "Observable", "id": "observable_label"},
    {"name": "Raw enum", "id": "observable_raw"},
    {"name": "bin low", "id": "bin_low"},
    {"name": "bin high", "id": "bin_high"},
    {"name": "order", "id": "order"},
    {"name": "registered", "id": "registered"},
]
UNC_COLUMNS = [
    {"name": "Observable", "id": "observable_label"},
    {"name": "Raw id", "id": "raw_observable"},
    {"name": "bin low", "id": "bin_low"},
    {"name": "bin high", "id": "bin_high"},
    {"name": "central", "id": "central"},
    {"name": "sigma", "id": "sigma"},
    {"name": "sigma_minus", "id": "sigma_minus"},
    {"name": "sigma_plus", "id": "sigma_plus"},
    {"name": "mu", "id": "mu"},
    {"name": "mode", "id": "mode"},
    {"name": "skew", "id": "skew"},
    {"name": "symmetric", "id": "symmetric"},
]
PSPEC_COLUMNS = [
    {"name": "Parameter", "id": "parameter"},
    {"name": "type", "id": "type"},
    {"name": "block", "id": "block"},
    {"name": "code", "id": "code"},
]
FIT_COLUMNS = [
    {"name": "Parameter", "id": "parameter"},
    {"name": "Raw", "id": "raw_parameter"},
    {"name": "best_fit", "id": "best_fit"},
    {"name": "std", "id": "std"},
]
ETA_COLUMNS = [
    {"name": "Nuisance", "id": "nuisance"},
    {"name": "Raw", "id": "raw_nuisance"},
    {"name": "value", "id": "value"},
]


def observable_selection(prefix: str):
    return html.Div([
        field("Selection mode", dropdown(f"{prefix}-mode", SELECT_MODE, value="observable")),
        field("Decay", dropdown(f"{prefix}-decays", svc.decay_options(), value=DEFAULT_DECAY, multi=True)),
        field("Observables in selected decay", dropdown(f"{prefix}-obs", svc.observable_options_for_decays(DEFAULT_DECAY), value=DEFAULT_OBS, multi=True)),
        html.Div(className="form-grid-2", children=[
            field("QCD order", dropdown(f"{prefix}-order", enum_options(QCDOrder), value="NNLO")),
            field("Add dependencies", dcc.Checklist(id=f"{prefix}-deps", options=[{"label": "include observable dependencies", "value": "deps"}], value=["deps"], className="checklist")),
        ]),
    ])


def bin_controls(prefix: str):
    return html.Div([
        field("Bins", dcc.Checklist(id=f"{prefix}-use-bin", options=[{"label": "apply bin to selected observables", "value": "bin"}], value=[], className="checklist")),
        html.Div(id=f"{prefix}-bin-mode-wrap", style={"display": "none"}, children=[
            field("Bin grid", dcc.Checklist(id=f"{prefix}-smooth", options=[{"label": "build bins from min/max/step", "value": "smooth"}], value=[], className="checklist")),
        ]),
        html.Div(id=f"{prefix}-explicit-bin-wrap", style={"display": "none"}, className="form-grid-2", children=[
            field("Bin low", num_input(f"{prefix}-bin-low", 1.0)),
            field("Bin high", num_input(f"{prefix}-bin-high", 6.0)),
        ]),
        html.Div(id=f"{prefix}-smooth-bin-wrap", style={"display": "none"}, className="form-grid-3", children=[
            field("Smooth min", num_input(f"{prefix}-smooth-min", 1.0)),
            field("Smooth max", num_input(f"{prefix}-smooth-max", 6.0)),
            field("Step", num_input(f"{prefix}-smooth-step", 1.0)),
        ]),
    ])


def pspec_controls():
    return html.Div([
        html.Div(className="form-grid-3", children=[
            field("ParameterType", dropdown("stat-pspec-type", svc.parameter_type_options(), value=svc.default_parameter_type_name("FLAVOR"))),
            field("Block", dropdown("stat-pspec-block", [], value=None, placeholder="Choose a block...")),
            field("Code", dropdown("stat-pspec-code", [], value=None, placeholder="Choose a code...")),
        ]),
        html.Button("Add p_spec", id="stat-add-pspec-btn", n_clicks=0),
    ])


def expert_controls():
    return html.Div([
        field("Expert controls", dcc.Checklist(id="stat-advanced-mode", options=[{"label": "show decay and advanced statistic controls", "value": "advanced"}], value=[], className="checklist")),
        html.Div(id="stat-expert-panel", style={"display": "none"}, children=[
            small_note("Decay configuration is applied to the shared ObservableInterface used by both Observable and Stat pages.", tone="info"),
            field("Decay configuration target", dropdown("stat-decay-config-decay", svc.decay_options(), value="B__l_l", multi=False)),
            html.Div(id="stat-decay-config-fields"),
            html.Div(className="inline-actions", children=[
                html.Button("Apply decay configuration", id="stat-apply-decay-config-btn", n_clicks=0),
            ]),
            status_box("stat-decay-config-status", "No decay configuration applied yet."),
        ]),
    ])


def layout():
    return html.Div(children=[
        page_title("Statistics", "χ² Monte-Carlo covariance workflow: configure observables, propagate uncertainty, and fit selected parameters.", "Stat"),
        html.Div(className="page-grid", children=[
            html.Div(className="grid", children=[
                card("Statistic observables", "shared ObservableInterface", html.Div([
                    observable_selection("stat-obs"),
                    bin_controls("stat-obs"),
                    expert_controls(),
                    html.Div(className="inline-actions", children=[
                        html.Button("Configure Statistic ObservableInterface", id="stat-configure-obs-btn", n_clicks=0),
                        html.Button("Remove selected rows", id="stat-remove-obs-btn", n_clicks=0),
                    ]),
                    status_box("stat-observable-status", svc.stat_observable_status_text()),
                ])),
                card("Experiments and MC", "χ² backend only", html.Div([
                    field("Experiments", text_input("stat-experiments", "", placeholder="comma-separated; empty = all")),
                    html.Div(className="form-grid-3", children=[
                        field("MC draws", num_input("stat-mc-draws", 100)),
                        field("MC threads", num_input("stat-mc-threads", 1)),
                        field("MC seed", num_input("stat-mc-seed", 123456)),
                    ]),
                    html.Div(id="stat-advanced-options-wrap", style={"display": "none"}, children=[
                        html.Div(className="form-grid-2", children=[
                            field("Skew threshold", num_input("stat-skew-threshold", 0.2)),
                            field("Nuisance contexts", num_input("stat-nuisance-contexts", 2)),
                        ]),
                        html.Div(className="form-grid-3", children=[
                            field("Cov ridge rel", num_input("stat-ridge-rel", 1e-8)),
                            field("Cov ridge abs", num_input("stat-ridge-abs", 1e-12)),
                            field("Nuisance seed", num_input("stat-nuisance-seed", 12345)),
                        ]),
                        field("Nuisance pruning", dcc.Checklist(id="stat-nuisance-pruning", options=[{"label": "sensitivity pruning", "value": "prune"}], value=["prune"], className="checklist")),
                    ]),
                ])),
                card("Uncertainty", "GaussianSummary", html.Div([
                    field("Uncertainty display", dropdown("stat-uncertainty-mode", [{"label": "symmetric", "value": "sym"}, {"label": "asymmetric", "value": "asym"}], value="sym")),
                    html.Button("Compute uncertainties", id="stat-uncertainty-btn", n_clicks=0),
                    html.Div(
                        id="stat-uncertainty-progress-wrap",
                        className="task-progress-wrap",
                        style={"display": "none"},
                        children=[
                            html.Progress(className="task-progress", max=100),
                            html.Span("Monte-Carlo / uncertainty computation in progress…"),
                        ],
                    ),
                    status_box("stat-uncertainty-status", "No uncertainty computation yet."),
                ])),
                card("Fit and contour", "max 10 p_specs; 2D enables contour scan", html.Div([
                    pspec_controls(),
                    data_table("stat-p-specs-table", PSPEC_COLUMNS, data=[], page_size=10, row_selectable="multi"),
                    html.Button("Remove selected p_spec rows", id="stat-remove-pspec-btn", n_clicks=0),
                    html.Div(className="form-grid-4", children=[
                        field("x half width", num_input("stat-x-half-width", 1.0)),
                        field("y half width", num_input("stat-y-half-width", 1.0)),
                        field("nx", num_input("stat-nx", 25)),
                        field("ny", num_input("stat-ny", 25)),
                    ]),
                    field("Contour", dcc.Checklist(id="stat-do-contour", options=[{"label": "run 2D likelihood scan when exactly two p_specs are provided", "value": "contour"}], value=["contour"], className="checklist")),
                    html.Button("Run fit / contour", id="stat-fit-btn", n_clicks=0),
                    html.Div(
                        id="stat-fit-progress-wrap",
                        className="task-progress-wrap",
                        style={"display": "none"},
                        children=[
                            html.Progress(className="task-progress", max=100),
                            html.Span("Fit / contour computation in progress…"),
                        ],
                    ),
                    status_box("stat-fit-status", "No fit yet."),
                ])),
            ]),
            html.Div(className="grid", children=[
                card("Configured stat observables", "selection expanded from decay and bins", data_table("stat-observable-table", OBS_COLUMNS, data=svc.current_stat_observable_rows(), page_size=12, row_selectable="multi")),
                card("Uncertainty table", "GaussianSummary", data_table("stat-uncertainty-table", UNC_COLUMNS, page_size=16)),
                card("Uncertainty plot", "central value plus uncertainty band", graph("stat-uncertainty-fig", height=480), className="card graph-card"),
                html.Div(className="graph-row-2", children=[
                    card("Best-fit parameters", "p_hat", data_table("stat-fit-table", FIT_COLUMNS, page_size=12)),
                    card("Profiled nuisances", "eta_hat", data_table("stat-eta-table", ETA_COLUMNS, page_size=12)),
                ]),
                card("Fit correlations", "p_correlations", graph("stat-corr-fig", height=460), className="card graph-card"),
                card("2D likelihood contour", "ΔNLL scan around best fit", graph("stat-contour-fig", height=520), className="card graph-card"),
                small_note("LaTeX labels use the supplied observable/decay/nuisance maps. Raw enum/block/code values remain visible in companion columns."),
            ]),
        ]),
    ])
