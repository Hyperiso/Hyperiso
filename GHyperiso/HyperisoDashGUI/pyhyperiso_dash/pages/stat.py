from __future__ import annotations

from dash import dcc, html

from pyhyperiso_dash.components import card, data_table, dropdown, enum_options, field, num_input, page_title, small_note, status_box, text_input
from pyhyperiso.core.Common.GeneralEnum import Decays, Observables, ParameterType, QCDOrder

SELECT_MODE = [{"label": "Observable by observable", "value": "observable"}, {"label": "Whole decay", "value": "decay"}]
BIN_STRATEGY = [
    {"label": "No explicit bin", "value": "none"},
    {"label": "Single bin applied to selected observables", "value": "single"},
    {"label": "Smooth bin grid min/max/step", "value": "smooth"},
]


def observable_selection(prefix: str):
    return html.Div([
        field("Selection mode", dropdown(f"{prefix}-mode", SELECT_MODE, value="observable")),
        field("Observables", dropdown(f"{prefix}-obs", enum_options(Observables), value=["BR_BS_MUMU"], multi=True)),
        field("Decays", dropdown(f"{prefix}-decays", enum_options(Decays), value=[], multi=True)),
        html.Div(className="form-grid-2", children=[
            field("QCD order", dropdown(f"{prefix}-order", enum_options(QCDOrder), value="NNLO")),
            field("Add dependencies", dcc.Checklist(id=f"{prefix}-deps", options=[{"label": "include observable dependencies", "value": "deps"}], value=["deps"], className="checklist")),
        ]),
    ])


def bin_controls(prefix: str):
    return html.Div([
        field("Bin strategy", dropdown(f"{prefix}-bin-strategy", BIN_STRATEGY, value="none")),
        html.Div(className="form-grid-2", children=[field("Bin low", num_input(f"{prefix}-bin-low", 1.0)), field("Bin high", num_input(f"{prefix}-bin-high", 6.0))]),
        html.Div(className="form-grid-3", children=[field("Smooth min", num_input(f"{prefix}-smooth-min", 1.0)), field("Smooth max", num_input(f"{prefix}-smooth-max", 6.0)), field("Step", num_input(f"{prefix}-smooth-step", 1.0))]),
    ])


def stat_config_controls():
    return html.Div([
        html.Div(className="warn-box", children="χ²_MC_COVARIANCE is forced. Laplace/Minuit-only options are intentionally hidden for now."),
        html.Div(className="form-grid-3", children=[
            field("MC draws", num_input("stat-mc-draws", 100)),
            field("Skew threshold", num_input("stat-skew-threshold", 0.2)),
            field("Nuisance contexts", num_input("stat-nuisance-contexts", 2)),
        ]),
        html.Div(className="form-grid-3", children=[
            field("χ² ridge rel", num_input("stat-ridge-rel", 1e-8)),
            field("χ² ridge abs", num_input("stat-ridge-abs", 1e-12)),
            field("Sensitivity seed", num_input("stat-nuisance-seed", 12345)),
        ]),
        field("Nuisance pruning", dcc.Checklist(id="stat-nuisance-pruning", options=[{"label": "enable finite-difference sensitivity pruning", "value": "prune"}], value=["prune"], className="checklist")),
        field("Selected experiments", dcc.Textarea(id="stat-experiments", value="", placeholder="comma-separated experiment names, leave empty for all", style={"width": "100%", "height": "68px"})),
    ])


def layout():
    p_spec_columns = ["type", "block", "code"]
    default_p_specs = [{"type": "DECAY", "block": "B_ll", "code": "1"}]
    return html.Div(children=[
        page_title("Statistic interface", "χ² uncertainty propagation, best fits and two-parameter likelihood contours.", "Stat"),
        html.Div(className="page-grid", children=[
            html.Div(className="grid", children=[
                card("Observable selection", "used by StatisticInterface", html.Div([observable_selection("stat-obs"), bin_controls("stat-obs")])) ,
                card("StatisticConfig", "χ²-only controls", stat_config_controls()),
                card("Uncertainty plot", "MC propagation", html.Div([
                    field("Uncertainty display", dcc.RadioItems(id="stat-uncertainty-mode", options=[{"label": "symmetric σ", "value": "sym"}, {"label": "asymmetric σ−/σ+", "value": "asym"}], value="sym", className="checklist")),
                    html.Button("Compute uncertainties", id="stat-uncertainty-btn", n_clicks=0),
                    status_box("stat-uncertainty-status", "No uncertainty computation yet."),
                ])),
                card("Fit and contour", "max 10 p_specs; 2D enables contour scan", html.Div([
                    data_table("stat-p-specs-table", p_spec_columns, data=default_p_specs, page_size=10, editable=True),
                    html.Button("Add empty p_spec row", id="stat-add-pspec-btn", n_clicks=0),
                    html.Div(className="form-grid-4", children=[
                        field("x half width", num_input("stat-x-half-width", 1.0)),
                        field("y half width", num_input("stat-y-half-width", 1.0)),
                        field("nx", num_input("stat-nx", 25)),
                        field("ny", num_input("stat-ny", 25)),
                    ]),
                    field("Contour", dcc.Checklist(id="stat-do-contour", options=[{"label": "run 2D likelihood scan when exactly two p_specs are provided", "value": "contour"}], value=["contour"], className="checklist")),
                    html.Button("Run fit / contour", id="stat-fit-btn", n_clicks=0),
                    status_box("stat-fit-status", "No fit yet."),
                ])),
            ]),
            html.Div(className="grid", children=[
                card("Uncertainty table", "GaussianSummary", data_table("stat-uncertainty-table", ["observable", "bin_low", "bin_high", "central", "mu", "mode", "sigma", "sigma_minus", "sigma_plus", "skew", "symmetric"], page_size=16)),
                card("Uncertainty plot", "central value plus uncertainty band", dcc.Graph(id="stat-uncertainty-fig"), className="card graph-card"),
                html.Div(className="graph-row-2", children=[
                    card("Best-fit parameters", "p_hat", data_table("stat-fit-table", ["parameter", "best_fit", "std"], page_size=12)),
                    card("Profiled nuisances", "eta_hat", data_table("stat-eta-table", ["nuisance", "value"], page_size=12)),
                ]),
                card("Fit correlations", "p_correlations", dcc.Graph(id="stat-corr-fig"), className="card graph-card"),
                card("2D likelihood contour", "ΔNLL scan around best fit", dcc.Graph(id="stat-contour-fig"), className="card graph-card"),
                small_note("MC_draws above 200 can make fit+scan slow. The GUI does not expose Laplace/Minuit-specific configuration, per your current constraint."),
            ]),
        ]),
    ])
