from __future__ import annotations

from dash import dcc, html

from pyhyperiso_dash import services as svc
from pyhyperiso_dash.components import card, data_table, dropdown, field, graph, num_input, page_title, small_note, status_box

QCD_RESULT_COLUMNS = [
    {"name": "Quantity", "id": "quantity"},
    {"name": "scale", "id": "scale"},
    {"name": "PDG id", "id": "pdg_id"},
    {"name": "scheme", "id": "scheme"},
    {"name": "value", "id": "value"},
]

QCD_CONSTANT_COLUMNS = [
    {"name": "Quantity", "id": "quantity"},
    {"name": "index", "id": "index"},
    {"name": "value", "id": "value"},
]


def scheme_controls(prefix: str):
    return html.Div(
        className="form-grid-2",
        children=[
            field("m_b scheme", dropdown(f"{prefix}-mb-type", svc.mass_type_options(), value="POLE")),
            field("m_t scheme", dropdown(f"{prefix}-mt-type", svc.mass_type_options(), value="POLE")),
        ],
    )


def layout():
    return html.Div(
        children=[
            page_title("QCD / QED running", "Compute and plot αs(μ), αem(μ), running masses, and constants from the initialized Hyperiso backend.", "QCD"),
            html.Div(
                className="page-grid",
                children=[
                    html.Div(
                        className="grid",
                        children=[
                            card(
                                "Single-scale values",
                                "αs and running mass",
                                html.Div([
                                    html.Div(className="form-grid-3", children=[
                                        field("Scale μ [GeV]", num_input("qcd-single-scale", 91.1876)),
                                        field("Particle / PDG id", dropdown("qcd-mass-pdg", svc.mass_particle_options(), value=5)),
                                        field("Include αem", dcc.Checklist(id="qcd-include-qed", options=[{"label": "also compute αem", "value": "qed"}], value=[], className="checklist")),
                                    ]),
                                    scheme_controls("qcd-single"),
                                    html.Button("Compute values", id="qcd-compute-btn", n_clicks=0),
                                    status_box("qcd-compute-status", "No QCD computation yet."),
                                ]),
                            ),
                            card(
                                "αs scan",
                                "plot αs as a function of μ",
                                html.Div([
                                    html.Div(className="form-grid-3", children=[
                                        field("μ min [GeV]", num_input("qcd-alpha-min", 1.0)),
                                        field("μ max [GeV]", num_input("qcd-alpha-max", 200.0)),
                                        field("points", num_input("qcd-alpha-n", 80)),
                                    ]),
                                    scheme_controls("qcd-alpha"),
                                    html.Button("Plot αs", id="qcd-alpha-scan-btn", n_clicks=0),
                                    status_box("qcd-alpha-status", "No αs scan yet."),
                                ]),
                            ),
                            card(
                                "αem scan",
                                "plot electromagnetic coupling as a function of μ",
                                html.Div([
                                    html.Div(className="form-grid-3", children=[
                                        field("μ min [GeV]", num_input("qcd-alphaem-min", 1.0)),
                                        field("μ max [GeV]", num_input("qcd-alphaem-max", 200.0)),
                                        field("points", num_input("qcd-alphaem-n", 80)),
                                    ]),
                                    scheme_controls("qcd-alphaem"),
                                    html.Button("Plot αem", id="qcd-alphaem-scan-btn", n_clicks=0),
                                    status_box("qcd-alphaem-status", "No αem scan yet."),
                                ]),
                            ),
                            card(
                                "Mass scan",
                                "plot a running mass as a function of μ",
                                html.Div([
                                    html.Div(className="form-grid-4", children=[
                                        field("Particle / PDG id", dropdown("qcd-mass-scan-pdg", svc.mass_particle_options(), value=5)),
                                        field("μ min [GeV]", num_input("qcd-mass-min", 1.0)),
                                        field("μ max [GeV]", num_input("qcd-mass-max", 200.0)),
                                        field("points", num_input("qcd-mass-n", 80)),
                                    ]),
                                    scheme_controls("qcd-mass"),
                                    html.Button("Plot mass", id="qcd-mass-scan-btn", n_clicks=0),
                                    status_box("qcd-mass-status", "No mass scan yet."),
                                ]),
                            ),
                        ],
                    ),
                    html.Div(
                        className="grid",
                        children=[
                            card("Results", "latest single-scale values", data_table("qcd-result-table", QCD_RESULT_COLUMNS, page_size=8)),
                            card("αs plot", "running coupling", graph("qcd-alpha-fig", height=480), className="card graph-card"),
                            card("αem plot", "electromagnetic coupling", graph("qcd-alphaem-fig", height=480), className="card graph-card"),
                            card("Mass plot", "running mass", graph("qcd-mass-fig", height=480), className="card graph-card"),
                            card("QCD constants", "Nc, color factors, beta/gamma coefficients", html.Div([
                                field("Constant group", dropdown("qcd-constants-kind", [
                                    {"label": "All", "value": "all"},
                                    {"label": "Color factors", "value": "color"},
                                    {"label": "β coefficients", "value": "beta"},
                                    {"label": "γ coefficients", "value": "gamma"},
                                ], value="all")),
                                html.Button("Load constants", id="qcd-constants-btn", n_clicks=0),
                                data_table("qcd-constants-table", QCD_CONSTANT_COLUMNS, page_size=12),
                            ])),
                            small_note("The QCDProvider and QEDProvider expect Hyperiso to be initialized on the Core page before use."),
                        ],
                    ),
                ],
            ),
        ]
    )
