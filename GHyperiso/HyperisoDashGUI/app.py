from __future__ import annotations

from pathlib import Path

import plotly.express as px
from dash import Dash, Input, Output, dcc, html

from pyhyperiso_dash.callbacks import register_callbacks
from pyhyperiso_dash.pages import core, observable, stat, wilson

px.defaults.template = "plotly_dark"
px.defaults.color_discrete_sequence = px.colors.qualitative.Set2

PAGES = {
    "/": core.layout,
    "/core": core.layout,
    "/wilson": wilson.layout,
    "/observable": observable.layout,
    "/stat": stat.layout,
}

NAV_ITEMS = [
    ("/core", "Core", "Init, LHA, blocks"),
    ("/wilson", "Wilson", "Build, requests, scans"),
    ("/observable", "Observable", "Predictions and parameter scans"),
    ("/stat", "Stat", "χ² uncertainty, fit and contours"),
]


def make_app() -> Dash:
    app = Dash(
        __name__,
        suppress_callback_exceptions=True,
        title="PyHyperiso Dash GUI",
        assets_folder=str(Path(__file__).resolve().parent / "pyhyperiso_dash" / "assets"),
    )

    app.layout = html.Div(
        className="app-shell",
        children=[
            dcc.Location(id="url"),
            dcc.Store(id="runtime-ping", data=0),
            html.Aside(
                className="sidebar",
                children=[
                    html.Div(
                        className="brand",
                        children=[
                            html.H1("PyHyperiso GUI"),
                            html.Span("Dash", className="badge"),
                        ],
                    ),
                    html.Div(
                        className="status",
                        children=(
                            "Single-process Hyperiso runtime. The C++ singleton is kept "
                            "alive in the Python Dash process. Use Core → Switch LHA to "
                            "change the active input safely."
                        ),
                    ),
                    html.Hr(className="hr"),
                    html.Nav(
                        className="nav",
                        children=[
                            dcc.Link(
                                href=href,
                                className="nav-item",
                                children=[html.B(label), html.Span(subtitle)],
                            )
                            for href, label, subtitle in NAV_ITEMS
                        ],
                    ),
                    html.Hr(className="hr"),
                    html.Div(
                        className="small-note",
                        children=(
                            "Recommended for long fits: run with one Dash worker. "
                            "Multiple gunicorn workers each own their own C++ singleton."
                        ),
                    ),
                ],
            ),
            html.Main(
                className="content",
                children=[dcc.Loading(id="page-loading", type="circle", children=html.Div(id="page-content"))],
            ),
        ],
    )

    @app.callback(Output("page-content", "children"), Input("url", "pathname"))
    def render_page(pathname: str):
        return PAGES.get(pathname or "/core", core.layout)()

    register_callbacks(app)
    return app


app = make_app()
server = app.server


if __name__ == "__main__":
    app.run(debug=True)
