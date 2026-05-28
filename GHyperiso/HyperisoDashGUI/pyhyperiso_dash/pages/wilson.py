from __future__ import annotations

from dash import dcc, html

from pyhyperiso_dash import services as svc
from pyhyperiso_dash.components import card, data_table, dropdown, enum_options, field, graph, num_input, page_title, small_note, status_box, text_input
from pyhyperiso.core.Common.GeneralEnum import ContributionType, ParameterType, QCDOrder, WCoeff, WGroup, WilsonBasis

METHOD_OPTIONS = [
    {"label": "M — matching at fixed QCD order", "value": "M"},
    {"label": "FM — full matching summed up to order", "value": "FM"},
    {"label": "R — running at fixed QCD order", "value": "R"},
    {"label": "FR — full running summed up to order", "value": "FR"},
]
DIM_OPTIONS = [{"label": "1D scatter", "value": "1d"}, {"label": "2D heatmap", "value": "2d"}]

SCALAR_COMPONENT_OPTIONS = [
    {"label": "Real part", "value": "real"},
    {"label": "Imaginary part", "value": "imag"},
    {"label": "Magnitude", "value": "abs"},
]


def parameter_controls(prefix: str):
    return html.Div(
        className="form-grid-3",
        children=[
            field("ParameterType", dropdown(f"{prefix}-ptype", svc.parameter_type_options(), value=svc.default_parameter_type_name("WILSON"))),
            field("Block", text_input(f"{prefix}-block", "EW_SCALE")),
            field("Code", text_input(f"{prefix}-code", "1")),
        ],
    )


def layout():
    return html.Div(
        children=[
            page_title("Wilson interface", "Build Wilson groups, query matching/running coefficients, and scan parameters.", "Wilson"),
            html.Div(
                className="page-grid",
                children=[
                    html.Div(
                        className="grid",
                        children=[
                            card(
                                "Build Wilson pipeline",
                                "WilsonBuildConfig",
                                children=html.Div([
                                    field("Groups", dropdown("wilson-groups", enum_options(WGroup), value=["B"], multi=True)),
                                    html.Div(className="form-grid-3", children=[
                                        field("Matching scale μW", num_input("wilson-matching-scale", 81.0)),
                                        field("Hadronic scale μh", num_input("wilson-hadronic-scale", 2.0)),
                                        field("QCD order", dropdown("wilson-build-order", enum_options(QCDOrder), value="LO")),
                                    ]),
                                    html.Div(className="inline-actions", children=[
                                        html.Button("Build / rebuild", id="wilson-build-btn", n_clicks=0),
                                        html.Button("Add selected groups", id="wilson-add-btn", n_clicks=0),
                                    ]),
                                    status_box("wilson-build-status", "WilsonInterface is not built yet."),
                                ]),
                            ),
                            card(
                                "Coefficient request",
                                "single value",
                                children=html.Div([
                                    html.Div(className="form-grid-2", children=[
                                        field("Method", dropdown("wilson-method", METHOD_OPTIONS, value="FM")),
                                        field("Scalar component", dropdown("wilson-component", SCALAR_COMPONENT_OPTIONS, value="real")),
                                    ]),
                                    html.Div(className="form-grid-2", children=[
                                        field("Group", dropdown("wilson-request-group", enum_options(WGroup), value="B")),
                                        field("Coefficient", dropdown("wilson-request-coeff", enum_options(WCoeff), value="C7")),
                                    ]),
                                    html.Div(className="form-grid-3", children=[
                                        field("Order", dropdown("wilson-request-order", enum_options(QCDOrder), value="NNLO")),
                                        field("Contribution", dropdown("wilson-contribution", enum_options(ContributionType), value="TOTAL")),
                                        field("Basis", dropdown("wilson-basis", enum_options(WilsonBasis), value="STANDARD")),
                                    ]),
                                    html.Button("Run request", id="wilson-query-btn", n_clicks=0),
                                    status_box("wilson-query-status", "No Wilson request yet."),
                                ]),
                            ),
                            card(
                                "Parameter scan",
                                "temporary ParameterSetter mutation + restore",
                                children=html.Div([
                                    field("Plot dimension", dropdown("wilson-scan-dim", DIM_OPTIONS, value="1d")),
                                    html.Div(className="form-grid-3", children=[
                                        field("x min", num_input("wilson-x-min", 1.0)),
                                        field("x max", num_input("wilson-x-max", 80.0)),
                                        field("x points", num_input("wilson-x-n", 40)),
                                    ]),
                                    html.Div(className="section-title", children="X parameter"),
                                    parameter_controls("wilson-x-param"),
                                    html.Div(className="form-grid-3", children=[
                                        field("y min", num_input("wilson-y-min", 1.0)),
                                        field("y max", num_input("wilson-y-max", 10.0)),
                                        field("y points", num_input("wilson-y-n", 25)),
                                    ]),
                                    html.Div(className="section-title", children="Y parameter for 2D only"),
                                    parameter_controls("wilson-y-param"),
                                    html.Button("Run Wilson scan", id="wilson-scan-btn", n_clicks=0),
                                    status_box("wilson-scan-status", "No scan yet."),
                                ]),
                            ),
                        ],
                    ),
                    html.Div(className="grid", children=[
                        card("Request result", "latest coefficient", data_table("wilson-result-table", ["method", "group", "coefficient", "order", "contribution", "basis", "component", "value", "scalar"], page_size=6)),
                        card("Wilson scan plot", "1D scatter or 2D heatmap", graph("wilson-scan-fig", height=500), className="card graph-card"),
                        small_note("The scan code stores the original central value with ParameterProvider, mutates with ParameterSetter, evaluates the coefficient, then restores parameters in a finally block."),
                    ]),
                ],
            ),
        ]
    )
