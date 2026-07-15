from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd
import plotly.graph_objects as go

PLOT_FONT = dict(color="rgba(229,231,235,0.92)")
DEFAULT_FIG_HEIGHT = 430


def style_fig(fig: go.Figure, title: str | None = None) -> go.Figure:
    fig.update_layout(
        title=title,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=12, r=12, t=46 if title else 20, b=12),
        height=DEFAULT_FIG_HEIGHT,
        autosize=False,
        font=PLOT_FONT,
        legend=dict(
            bgcolor="rgba(0,0,0,0)", bordercolor="rgba(148,163,184,0.16)", borderwidth=1
        ),
    )
    fig.update_xaxes(
        gridcolor="rgba(148,163,184,0.14)", zerolinecolor="rgba(148,163,184,0.22)"
    )
    fig.update_yaxes(
        gridcolor="rgba(148,163,184,0.14)", zerolinecolor="rgba(148,163,184,0.22)"
    )
    return fig


def empty_fig(title: str = "No data") -> go.Figure:
    fig = go.Figure()
    fig.add_annotation(
        text=title, x=0.5, y=0.5, xref="paper", yref="paper", showarrow=False
    )
    return style_fig(fig, title)


def block_size_bar(rows: Sequence[dict]) -> go.Figure:
    if not rows:
        return empty_fig("Block inventory")
    df = pd.DataFrame(rows)
    fig = go.Figure(
        go.Bar(
            x=df["parameter_type"],
            y=df["n_blocks"],
            text=df["n_blocks"],
            textposition="auto",
        )
    )
    fig.update_layout(xaxis_title="Parameter type", yaxis_title="Number of blocks")
    return style_fig(fig, "Blocks by ParameterType")


def block_values_scatter(rows: Sequence[dict], title: str) -> go.Figure:
    if not rows:
        return empty_fig(title)
    df = pd.DataFrame(rows)
    y = pd.to_numeric(df.get("value"), errors="coerce")
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=list(range(len(df))),
            y=y,
            mode="markers",
            text=df.get("code"),
            name="value",
        )
    )
    if "combined_std" in df:
        err = pd.to_numeric(df["combined_std"], errors="coerce")
        fig.update_traces(
            error_y=dict(type="data", array=err.fillna(0).to_numpy(), visible=True)
        )
    fig.update_layout(xaxis_title="Entry index", yaxis_title="Value")
    return style_fig(fig, title)


def series_1d(
    x: Sequence[float], y: Sequence[float], title: str, x_title: str, y_title: str
) -> go.Figure:
    if not x:
        return empty_fig(title)
    fig = go.Figure(go.Scatter(x=x, y=y, mode="markers+lines", name=y_title))
    fig.update_layout(xaxis_title=x_title, yaxis_title=y_title)
    return style_fig(fig, title)


def heatmap_2d(
    x: Sequence[float],
    y: Sequence[float],
    z: Sequence[Sequence[float]],
    title: str,
    x_title: str,
    y_title: str,
    z_title: str,
) -> go.Figure:
    if not x or not y:
        return empty_fig(title)
    fig = go.Figure(go.Heatmap(x=x, y=y, z=z, colorbar=dict(title=z_title)))
    fig.update_layout(xaxis_title=x_title, yaxis_title=y_title)
    return style_fig(fig, title)


def uncertainty_fig(rows: Sequence[dict], asymmetric: bool = False) -> go.Figure:
    if not rows:
        return empty_fig("Uncertainty")
    df = pd.DataFrame(rows)
    if (
        "bin_center" in df
        and pd.to_numeric(df["bin_center"], errors="coerce").notna().any()
    ):
        x = df["bin_center"]
    else:
        x = (
            df.get("observable_label")
            if "observable_label" in df
            else df.get("observable")
        )
    if x is None:
        x = list(range(len(df)))
    y = pd.to_numeric(df["central"], errors="coerce").to_numpy()
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x, y=y, mode="markers+lines", name="central"))
    if asymmetric:
        up = pd.to_numeric(df.get("sigma_plus"), errors="coerce").fillna(0).to_numpy()
        down = (
            pd.to_numeric(df.get("sigma_minus"), errors="coerce").fillna(0).to_numpy()
        )
    else:
        sig = pd.to_numeric(df.get("sigma"), errors="coerce").fillna(0).to_numpy()
        up, down = sig, sig
    x_vals = list(x)
    upper = y + up
    lower = y - down
    fig.add_trace(
        go.Scatter(
            x=x_vals + x_vals[::-1],
            y=list(upper) + list(lower[::-1]),
            fill="toself",
            line=dict(width=0),
            opacity=0.28,
            name="uncertainty band",
            hoverinfo="skip",
        )
    )
    fig.update_layout(xaxis_title="Observable / bin", yaxis_title="Value")
    return style_fig(fig, "Theory uncertainty")


def correlation_heatmap(rows: Sequence[dict], title: str) -> go.Figure:
    if not rows:
        return empty_fig(title)
    df = pd.DataFrame(rows)
    x = sorted(df["x"].unique())
    y = sorted(df["y"].unique())
    z = np.full((len(y), len(x)), np.nan)
    for _, r in df.iterrows():
        z[y.index(r["y"]), x.index(r["x"])] = r["corr"]
    fig = go.Figure(
        go.Heatmap(x=x, y=y, z=z, zmin=-1, zmax=1, colorbar=dict(title="ρ"))
    )
    return style_fig(fig, title)


def likelihood_contour(
    points: Sequence[dict], levels: Sequence[float], title: str = "Likelihood contour"
) -> go.Figure:
    if not points:
        return empty_fig(title)
    df = pd.DataFrame(points)
    xs = sorted(df["x"].unique())
    ys = sorted(df["y"].unique())
    z = np.full((len(ys), len(xs)), np.nan)
    for _, r in df.iterrows():
        z[ys.index(r["y"]), xs.index(r["x"])] = r["delta_nll"]
    fig = go.Figure()
    fig.add_trace(
        go.Heatmap(x=xs, y=ys, z=z, colorbar=dict(title="ΔNLL"), opacity=0.78)
    )
    if levels:
        fig.add_trace(
            go.Contour(
                x=xs,
                y=ys,
                z=z,
                contours=dict(
                    start=float(min(levels)),
                    end=float(max(levels)),
                    size=1.0,
                    coloring="none",
                ),
                line=dict(width=2),
                showscale=False,
                name="levels",
            )
        )
    fig.update_layout(xaxis_title="p₁", yaxis_title="p₂")
    return style_fig(fig, title)


def dependency_graph(data: dict, title: str = "Block dependency graph") -> go.Figure:
    """Draw the local dependency component for one selected block.

    Edges are oriented upstream → downstream.  The selected block is placed in
    the middle layer, transitive sources to the left, transitive dependents to
    the right, and unrelated nodes in the component are spread around it.
    """
    if not data or not data.get("nodes"):
        return empty_fig(title)
    block = str(data.get("block", ""))
    nodes = [str(x) for x in data.get("nodes", [])]
    sources = set(map(str, data.get("all_sources", [])))
    direct_sources = set(map(str, data.get("direct_sources", [])))
    dependents = set(map(str, data.get("all_dependents", [])))
    direct_dependents = set(map(str, data.get("direct_dependents", [])))
    edges = [(str(a), str(b)) for a, b in data.get("edges", [])]

    layers: dict[str, list[str]] = {
        "source": [],
        "selected": [],
        "dependent": [],
        "other": [],
    }
    for n in nodes:
        if n == block:
            layers["selected"].append(n)
        elif n in sources:
            layers["source"].append(n)
        elif n in dependents:
            layers["dependent"].append(n)
        else:
            layers["other"].append(n)

    xpos = {"source": -1.0, "selected": 0.0, "other": 0.0, "dependent": 1.0}
    pos: dict[str, tuple[float, float]] = {}
    for layer, ns in layers.items():
        ns = sorted(ns)
        if not ns:
            continue
        if len(ns) == 1:
            ys = [0.0]
        else:
            ys = list(np.linspace(0.82, -0.82, len(ns)))
        for n, y in zip(ns, ys):
            x = xpos[layer]
            if layer == "other":
                x = 0.0
            pos[n] = (x, float(y))

    fig = go.Figure()
    for src, dst in edges:
        if src not in pos or dst not in pos:
            continue
        x0, y0 = pos[src]
        x1, y1 = pos[dst]
        fig.add_trace(
            go.Scatter(
                x=[x0, x1],
                y=[y0, y1],
                mode="lines",
                line=dict(width=1.8),
                hoverinfo="skip",
                showlegend=False,
            )
        )
        # Small arrow marker near downstream side.
        xm = x0 + 0.78 * (x1 - x0)
        ym = y0 + 0.78 * (y1 - y0)
        fig.add_trace(
            go.Scatter(
                x=[xm],
                y=[ym],
                mode="markers",
                marker=dict(symbol="triangle-right", size=9),
                hoverinfo="skip",
                showlegend=False,
            )
        )

    for layer, ns in layers.items():
        if not ns:
            continue
        xs, ys, texts, hover = [], [], [], []
        for n in sorted(ns):
            if n not in pos:
                continue
            x, y = pos[n]
            xs.append(x)
            ys.append(y)
            texts.append(n)
            tags = []
            if n == block:
                tags.append("selected")
            if n in direct_sources:
                tags.append("direct upstream")
            elif n in sources:
                tags.append("upstream")
            if n in direct_dependents:
                tags.append("direct downstream")
            elif n in dependents:
                tags.append("downstream")
            hover.append(f"{n}<br>{', '.join(tags) if tags else layer}")
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="markers+text",
                text=texts,
                textposition="top center",
                hovertext=hover,
                hoverinfo="text",
                name=layer,
                marker=dict(
                    size=18 if layer == "selected" else 13, line=dict(width=1.2)
                ),
            )
        )

    fig.update_xaxes(visible=False, range=[-1.35, 1.35])
    fig.update_yaxes(visible=False, range=[-1.12, 1.12])
    return style_fig(fig, title)


def _rgba(hex_color: str, alpha: float) -> str:
    """Convert a ``#RRGGBB`` color to a Plotly RGBA string."""
    value = hex_color.lstrip("#")
    if len(value) != 6:
        raise ValueError("hex_color must use the #RRGGBB form")
    red, green, blue = (int(value[index : index + 2], 16) for index in (0, 2, 4))
    return f"rgba({red},{green},{blue},{alpha:.3g})"


def _confidence_contour_style(fig: go.Figure, title: str) -> go.Figure:
    """Apply the publication-style white theme used for confidence contours."""
    fig.update_layout(
        title=title,
        paper_bgcolor="white",
        plot_bgcolor="white",
        margin=dict(l=70, r=24, t=54, b=62),
        height=DEFAULT_FIG_HEIGHT,
        autosize=False,
        font=dict(color="#111827"),
        legend=dict(
            bgcolor="rgba(255,255,255,0.88)",
            bordercolor="rgba(55,65,81,0.25)",
            borderwidth=1,
        ),
    )
    axis_style = dict(
        showgrid=True,
        gridcolor="rgba(107,114,128,0.22)",
        zeroline=True,
        zerolinecolor="rgba(75,85,99,0.38)",
        linecolor="#374151",
        linewidth=1,
        mirror=True,
        ticks="outside",
        tickcolor="#374151",
    )
    fig.update_xaxes(**axis_style)
    fig.update_yaxes(**axis_style)
    return fig


def confidence_contour_paths(
    contours: Sequence[dict],
    title: str = "Confidence contours",
    x_title: str = "p₁",
    y_title: str = "p₂",
    best_fit: dict | None = None,
    bounds: Sequence[float] | None = None,
) -> go.Figure:
    """Plot filled confidence regions returned by the core ``ContourEngine``.

    Each confidence level is rendered as a semi-transparent filled region with
    a darker boundary. Larger sigma regions are drawn first, so nested smaller
    regions remain visible. Disconnected paths at the same level share their
    color and legend entry.
    """
    if not contours:
        return empty_fig(title)

    palette = ["#2563EB", "#DC2626", "#059669", "#7C3AED", "#D97706", "#0891B2"]
    fig = go.Figure()
    shown_levels: set[float] = set()
    ordered = sorted(
        contours,
        key=lambda row: (-float(row.get("sigma", 0.0)), int(row.get("path_id", 0))),
    )
    sigma_values = sorted({float(row.get("sigma", 0.0)) for row in contours})
    sigma_colors = {
        sigma: palette[index % len(palette)] for index, sigma in enumerate(sigma_values)
    }

    for item in ordered:
        points = list(item.get("points") or [])
        if len(points) < 2:
            continue
        sigma = float(item.get("sigma", 0.0))
        level = float(item.get("level", float("nan")))
        xs = [float(point["x"]) for point in points]
        ys = [float(point["y"]) for point in points]
        if xs[0] != xs[-1] or ys[0] != ys[-1]:
            xs.append(xs[0])
            ys.append(ys[0])

        color = sigma_colors[sigma]
        showlegend = sigma not in shown_levels
        shown_levels.add(sigma)
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines",
                fill="toself",
                fillcolor=_rgba(color, 0.24),
                name=f"{sigma:g}σ",
                legendgroup=f"sigma-{sigma:g}",
                showlegend=showlegend,
                line=dict(color=color, width=2.4),
                customdata=[[level]] * len(xs),
                hovertemplate=(
                    "x=%{x:.6g}<br>y=%{y:.6g}<br>"
                    "ΔNLL level=%{customdata[0]:.6g}<extra></extra>"
                ),
            )
        )

    if best_fit is not None:
        fig.add_trace(
            go.Scatter(
                x=[float(best_fit["x"])],
                y=[float(best_fit["y"])],
                mode="markers",
                marker=dict(size=12, symbol="x", color="#111827", line=dict(width=2)),
                name="best fit",
                hovertemplate="best fit<br>x=%{x:.6g}<br>y=%{y:.6g}<extra></extra>",
            )
        )

    fig.update_layout(xaxis_title=x_title, yaxis_title=y_title)
    if bounds and len(bounds) == 4:
        fig.update_xaxes(range=[float(bounds[0]), float(bounds[1])])
        fig.update_yaxes(range=[float(bounds[2]), float(bounds[3])])
    return _confidence_contour_style(fig, title)
