from __future__ import annotations

from typing import Any, Iterable, Sequence

from dash import dash_table, dcc, html


def page_title(title: str, subtitle: str, badge: str | None = None):
    return html.Div(
        className="page-title",
        children=[
            html.Div([html.H1(title), html.P(subtitle)]),
            html.Span(badge, className="badge") if badge else html.Span(),
        ],
    )


def card(title: str, subtitle: str | None = None, children: Any = None, className: str = "card"):
    return html.Div(
        className=className,
        children=[
            html.Div(
                className="card-title",
                children=[html.H2(title), html.Span(subtitle or "")],
            ),
            children if children is not None else html.Div(),
        ],
    )


def field(label: str, component: Any, className: str | None = None):
    return html.Div(className=className, children=[html.Label(label), component])


def text_input(id: str, value: str = "", placeholder: str = "", type: str = "text", debounce: bool = True):
    return dcc.Input(id=id, value=value, placeholder=placeholder, type=type, debounce=debounce, style={"width": "100%"})


def num_input(id: str, value: float | int | None = None, placeholder: str = "", debounce: bool = True):
    return dcc.Input(id=id, value=value, placeholder=placeholder, type="number", debounce=debounce, style={"width": "100%"})


def dropdown(id: str, options: Sequence[dict], value: Any = None, multi: bool = False, placeholder: str = "Select..."):
    return dcc.Dropdown(id=id, options=list(options), value=value, multi=multi, placeholder=placeholder, clearable=True)


def enum_options(enum_cls: type, include: Iterable[str] | None = None, exclude: Iterable[str] | None = None) -> list[dict[str, str]]:
    include_set = set(include or [])
    exclude_set = set(exclude or [])
    out = []
    for item in enum_cls:
        if include_set and item.name not in include_set:
            continue
        if item.name in exclude_set:
            continue
        out.append({"label": item.name, "value": item.name})
    return out


def status_box(id: str, text: str = "Ready."):
    return html.Pre(id=id, className="status", children=text)


def data_table(id: str, columns: Sequence[str], data: Sequence[dict] | None = None, page_size: int = 12, editable: bool = False):
    return dash_table.DataTable(
        id=id,
        columns=[{"name": c, "id": c, "editable": editable} for c in columns],
        data=list(data or []),
        page_size=page_size,
        editable=editable,
        row_deletable=editable,
        filter_action="native",
        sort_action="native",
        style_as_list_view=True,
        style_table={"overflowX": "auto"},
    )


def metric(label: str, value: Any, suffix: str = "", tone: str = ""):
    return html.Div(
        className=f"metric {tone}".strip(),
        children=[html.Div(label, className="k"), html.Div(str(value), className="v"), html.Div(suffix, className="s")],
    )


def small_note(text: str, tone: str = "warn"):
    return html.Div(text, className=f"{tone}-box")
