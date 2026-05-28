from __future__ import annotations

from dash import dcc, html

from pyhyperiso_dash import services as svc
from pyhyperiso_dash.components import card, data_table, dropdown, enum_options, field, metric, num_input, page_title, small_note, status_box, text_input
from pyhyperiso.core.Common.GeneralEnum import Model, ParameterType
from pyhyperiso.core.Core.HyperisoConfig import ExternalFlag

FLAG_OPTIONS = enum_options(ExternalFlag, exclude={"HAS_TH_OBSERVABLE_INPUT"})


def layout():
    return html.Div(
        children=[
            page_title("Core configuration", "Initialize Hyperiso, keep the C++ singleton alive, and inspect parameter blocks.", "Core"),
            html.Div(
                className="page-grid",
                children=[
                    html.Div(
                        className="grid",
                        children=[
                            card(
                                "LHA and HyperisoConfig",
                                "init or switch_lha",
                                children=html.Div(
                                    children=[
                                        dcc.Upload(id="core-upload", className="upload-box", children="Drop an LHA/SLHA/FLHA file here or click to upload"),
                                        field("LHA path", text_input("core-lha-path", "lha/testinput_thdm.lha", "relative Assets/lha/... or absolute path")),
                                        field("External flags", dcc.Checklist(id="core-flags", options=FLAG_OPTIONS, value=["HYP_AS_SM_MARTY"], className="checklist")),
                                        html.Div(
                                            className="form-grid-2",
                                            children=[
                                                field("Model", dropdown("core-model", enum_options(Model), value="SM")),
                                                html.Div(id="core-marty-name-wrap", children=field("MARTY model name", text_input("core-marty-name", "MSSM_UFO"))),
                                            ],
                                        ),
                                        html.Div(id="core-marty-path-wrap", children=field("MARTY model path", text_input("core-marty-path", "/my/custom/marty/path"))),
                                        html.Button("Initialize / switch active LHA", id="core-init-btn", n_clicks=0),
                                        status_box("core-status", "Hyperiso is not initialized yet."),
                                    ]
                                ),
                            ),
                            card(
                                "Block browser",
                                "ParameterType-aware",
                                children=html.Div(
                                    children=[
                                        html.Div(
                                            className="form-grid-2",
                                            children=[
                                                field("ParameterType", dropdown("core-block-ptype", svc.parameter_type_options(), value=svc.default_parameter_type_name("SM"))),
                                                field("Block", dropdown("core-block-name", [], value=None, placeholder="Refresh after init")),
                                            ],
                                        ),
                                        html.Button("Refresh block inventory", id="core-refresh-blocks", n_clicks=0),
                                        status_box("core-block-status", "No block inventory loaded."),
                                    ]
                                ),
                            ),
                        ],
                    ),
                    html.Div(
                        className="grid",
                        children=[
                            html.Div(id="core-metrics", className="metrics-row", children=[metric("Runtime", "OFF", "Hyperiso singleton", "bad"), metric("Model", "—", "active config"), metric("LHA", "—", "input"), metric("Wilson", "not built", "current session"), metric("Observables", "not built", "current session")]),
                            card("Block inventory", "all ParameterType namespaces", dcc.Graph(id="core-block-inventory-fig"), className="card graph-card"),
                            card("Block content", "values and uncertainties when available", data_table("core-block-table", ["code", "value", "stat_std", "syst_std", "combined_std", "scale", "bin"], page_size=16)),
                            card("Block value distribution", "with combined uncertainty error bars if readable", dcc.Graph(id="core-block-values-fig"), className="card graph-card"),
                            small_note("The table reads values through BlockLogger and tries to enrich each row with ParameterProvider uncertainties. Scale/bin metadata is shown when the underlying bound parameter exposes it."),
                        ],
                    ),
                ],
            ),
        ]
    )
