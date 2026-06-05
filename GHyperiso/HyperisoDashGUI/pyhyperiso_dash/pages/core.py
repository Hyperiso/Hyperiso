from __future__ import annotations

from dash import dcc, html

from pyhyperiso_dash import services as svc
from pyhyperiso_dash.components import card, data_table, dropdown, enum_options, field, graph, metric, num_input, page_title, small_note, status_box, text_input
from pyhyperiso.core.Common.GeneralEnum import Model, ParameterType
from pyhyperiso.core.Core.HyperisoConfig import ExternalFlag

FLAG_OPTIONS = enum_options(ExternalFlag, exclude={"HAS_TH_OBSERVABLE_INPUT"})


def core_status_text():
    s = svc.runtime_summary()
    if not s["initialized"]:
        return "Hyperiso is not initialized yet."
    return f"Hyperiso is initialized. Model: {s['model']}. LHA: {s['lha_path']}"


def runtime_metrics_children():
    s = svc.runtime_summary()
    return [
        metric("Runtime", "ON" if s["initialized"] else "OFF", "Hyperiso singleton", "good" if s["initialized"] else "bad"),
        metric("Model", s["model"], "active config"),
        metric("LHA", "set" if s["lha_path"] != "—" else "—", s["lha_path"]),
        metric("Wilson", "built" if s["wilson_built"] else "not built", "current session", "good" if s["wilson_built"] else ""),
        metric("Observables", str(s.get("observable_count", 0)), "configured entries", "good" if s.get("observable_count", 0) else ""),
    ]


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
                                        status_box("core-status", core_status_text()),
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
                            card(
                                "Block dependencies",
                                "upstream/downstream and pruning",
                                children=html.Div(
                                    children=[
                                        html.Div(
                                            className="form-grid-3",
                                            children=[
                                                field("Action", dropdown("core-dep-action", [{"label": "Detach", "value": "detach"}, {"label": "Reattach", "value": "reattach"}], value="detach")),
                                                field("Scope", dropdown("core-dep-scope", [{"label": "Whole block", "value": "block"}, {"label": "Single parameter", "value": "parameter"}], value="block")),
                                                field("Parameter code", dropdown("core-dep-code", [], value=None, placeholder="Only for parameter scope")),
                                            ],
                                        ),
                                        html.Div(className="inline-actions", children=[
                                            html.Button("Refresh dependency graph", id="core-dep-refresh-btn", n_clicks=0),
                                            html.Button("Apply dependency action", id="core-dep-apply-btn", n_clicks=0),
                                        ]),
                                        status_box("core-dep-status", "Select a block and refresh to inspect dependencies."),
                                    ]
                                ),
                            ),
                        ],
                    ),
                    html.Div(
                        className="grid",
                        children=[
                            html.Div(id="core-metrics", className="metrics-row", children=runtime_metrics_children()),
                            card("Block inventory", "all ParameterType namespaces", graph("core-block-inventory-fig", height=430), className="card graph-card"),
                            card("Block content", "values and uncertainties when available", data_table("core-block-table", ["code", "name", "dependent_block", "value", "stat_std", "syst_std", "combined_std", "scale", "bin"], page_size=16)),
                            card("Dependency tree", "connected upstream/downstream component for the selected block", graph("core-dep-fig", height=470), className="card graph-card"),
                            card("Dependency links", "direct and transitive relations", data_table("core-dep-table", ["relation", "block"], page_size=12)),
                            card("Block value distribution", "with combined uncertainty error bars if readable", graph("core-block-values-fig", height=430), className="card graph-card"),
                            small_note("The table reads values through BlockLogger and tries to enrich each row with ParameterProvider uncertainties. Scale/bin metadata is shown when the underlying bound parameter exposes it."),
                        ],
                    ),
                ],
            ),
        ]
    )
