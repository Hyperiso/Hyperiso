from __future__ import annotations

from dash import dcc, html

from pyhyperiso_dash import services as svc
from pyhyperiso_dash.components import card, data_table, dropdown, enum_options, field, graph, num_input, page_title, small_note, status_box, text_input
from pyhyperiso.core.Common.GeneralEnum import Decays, Observables, ParameterType, QCDOrder

SELECT_MODE = [{"label": "Observable by observable", "value": "observable"}, {"label": "Whole decay", "value": "decay"}]
DIM_OPTIONS = [{"label": "1D scatter", "value": "1d"}, {"label": "2D heatmap", "value": "2d"}]


def parameter_controls(prefix: str):
    return html.Div(
        className="form-grid-3",
        children=[
            field("ParameterType", dropdown(f"{prefix}-ptype", svc.parameter_type_options(), value=svc.default_parameter_type_name("SM"))),
            field("Block", text_input(f"{prefix}-block", "MASS")),
            field("Code", text_input(f"{prefix}-code", "25")),
        ],
    )


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
        field("Use one explicit bin", dcc.Checklist(id=f"{prefix}-use-bin", options=[{"label": "apply bin to selected observables", "value": "bin"}], value=[], className="checklist")),
        html.Div(className="form-grid-2", children=[
            field("Bin low", num_input(f"{prefix}-bin-low", 1.0)),
            field("Bin high", num_input(f"{prefix}-bin-high", 6.0)),
        ]),
        field("Smooth bins", dcc.Checklist(id=f"{prefix}-smooth", options=[{"label": "build bins from min/max/step", "value": "smooth"}], value=[], className="checklist")),
        html.Div(className="form-grid-3", children=[
            field("Smooth min", num_input(f"{prefix}-smooth-min", 1.0)),
            field("Smooth max", num_input(f"{prefix}-smooth-max", 6.0)),
            field("Step", num_input(f"{prefix}-smooth-step", 1.0)),
        ]),
    ])


def layout():
    return html.Div(
        children=[
            page_title("Observable interface", "Select observables or whole decays, compute predictions, and scan parameter dependence.", "Observable"),
            html.Div(className="page-grid", children=[
                html.Div(className="grid", children=[
                    card("Observable selection", "ObservableInterface", html.Div([
                        observable_selection("obs-build"),
                        bin_controls("obs-build"),
                        html.Button("Configure ObservableInterface", id="obs-build-btn", n_clicks=0),
                        status_box("obs-build-status", svc.observable_status_text()),
                    ])),
                    card("Compute", "current registered observables", html.Div([
                        html.Button("Compute current selection", id="obs-compute-btn", n_clicks=0),
                        status_box("obs-compute-status", "No computation yet."),
                    ])),
                    card("Parameter scan", "temporary parameter mutation + restore", html.Div([
                        field("Plot dimension", dropdown("obs-scan-dim", DIM_OPTIONS, value="1d")),
                        field("Observable", dropdown("obs-scan-observable", enum_options(Observables), value="BR_BS_MUMU")),
                        html.Div(className="form-grid-3", children=[
                            field("Order", dropdown("obs-scan-order", enum_options(QCDOrder), value="NNLO")),
                            field("Bin low", num_input("obs-scan-bin-low", None)),
                            field("Bin high", num_input("obs-scan-bin-high", None)),
                        ]),
                        field("Add dependencies", dcc.Checklist(id="obs-scan-deps", options=[{"label": "include dependencies", "value": "deps"}], value=["deps"], className="checklist")),
                        html.Div(className="form-grid-3", children=[
                            field("x min", num_input("obs-x-min", 0.0)),
                            field("x max", num_input("obs-x-max", 200.0)),
                            field("x points", num_input("obs-x-n", 35)),
                        ]),
                        html.Div(className="section-title", children="X parameter"),
                        parameter_controls("obs-x-param"),
                        html.Div(className="form-grid-3", children=[
                            field("y min", num_input("obs-y-min", 0.0)),
                            field("y max", num_input("obs-y-max", 200.0)),
                            field("y points", num_input("obs-y-n", 20)),
                        ]),
                        html.Div(className="section-title", children="Y parameter for 2D only"),
                        parameter_controls("obs-y-param"),
                        html.Button("Run observable scan", id="obs-scan-btn", n_clicks=0),
                        status_box("obs-scan-status", "No scan yet."),
                    ])),
                ]),
                html.Div(className="grid", children=[
                    card("Configured observables", "selection expanded from decay and bins", data_table("obs-selection-table", ["observable", "bin_low", "bin_high", "order", "registered"], page_size=12)),
                    card("Predictions", "compute_all output", data_table("obs-result-table", ["observable_id", "bin_low", "bin_high", "value"], page_size=16)),
                    card("Observable scan plot", "1D scatter or 2D heatmap", graph("obs-scan-fig", height=500), className="card graph-card"),
                    small_note("For smooth-binned plots, the backend constructs consecutive BinnedObservableId ranges from min to max with the requested step."),
                ]),
            ]),
        ]
    )
