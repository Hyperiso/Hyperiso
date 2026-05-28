from __future__ import annotations

from typing import Iterable, Sequence

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
        legend=dict(bgcolor="rgba(0,0,0,0)", bordercolor="rgba(148,163,184,0.16)", borderwidth=1),
    )
    fig.update_xaxes(gridcolor="rgba(148,163,184,0.14)", zerolinecolor="rgba(148,163,184,0.22)")
    fig.update_yaxes(gridcolor="rgba(148,163,184,0.14)", zerolinecolor="rgba(148,163,184,0.22)")
    return fig


def empty_fig(title: str = "No data") -> go.Figure:
    fig = go.Figure()
    fig.add_annotation(text=title, x=0.5, y=0.5, xref="paper", yref="paper", showarrow=False)
    return style_fig(fig, title)


def block_size_bar(rows: Sequence[dict]) -> go.Figure:
    if not rows:
        return empty_fig("Block inventory")
    df = pd.DataFrame(rows)
    fig = go.Figure(go.Bar(x=df["parameter_type"], y=df["n_blocks"], text=df["n_blocks"], textposition="auto"))
    fig.update_layout(xaxis_title="Parameter type", yaxis_title="Number of blocks")
    return style_fig(fig, "Blocks by ParameterType")


def block_values_scatter(rows: Sequence[dict], title: str) -> go.Figure:
    if not rows:
        return empty_fig(title)
    df = pd.DataFrame(rows)
    y = pd.to_numeric(df.get("value"), errors="coerce")
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=list(range(len(df))), y=y, mode="markers", text=df.get("code"), name="value"))
    if "combined_std" in df:
        err = pd.to_numeric(df["combined_std"], errors="coerce")
        fig.update_traces(error_y=dict(type="data", array=err.fillna(0).to_numpy(), visible=True))
    fig.update_layout(xaxis_title="Entry index", yaxis_title="Value")
    return style_fig(fig, title)


def series_1d(x: Sequence[float], y: Sequence[float], title: str, x_title: str, y_title: str) -> go.Figure:
    if not x:
        return empty_fig(title)
    fig = go.Figure(go.Scatter(x=x, y=y, mode="markers+lines", name=y_title))
    fig.update_layout(xaxis_title=x_title, yaxis_title=y_title)
    return style_fig(fig, title)


def heatmap_2d(x: Sequence[float], y: Sequence[float], z: Sequence[Sequence[float]], title: str, x_title: str, y_title: str, z_title: str) -> go.Figure:
    if not x or not y:
        return empty_fig(title)
    fig = go.Figure(go.Heatmap(x=x, y=y, z=z, colorbar=dict(title=z_title)))
    fig.update_layout(xaxis_title=x_title, yaxis_title=y_title)
    return style_fig(fig, title)


def uncertainty_fig(rows: Sequence[dict], asymmetric: bool = False) -> go.Figure:
    if not rows:
        return empty_fig("Uncertainty")
    df = pd.DataFrame(rows)
    x = df.get("bin_center") if "bin_center" in df else df.get("observable")
    if x is None:
        x = list(range(len(df)))
    y = pd.to_numeric(df["central"], errors="coerce").to_numpy()
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x, y=y, mode="markers+lines", name="central"))
    if asymmetric:
        up = pd.to_numeric(df.get("sigma_plus"), errors="coerce").fillna(0).to_numpy()
        down = pd.to_numeric(df.get("sigma_minus"), errors="coerce").fillna(0).to_numpy()
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
    fig = go.Figure(go.Heatmap(x=x, y=y, z=z, zmin=-1, zmax=1, colorbar=dict(title="ρ")))
    return style_fig(fig, title)


def likelihood_contour(points: Sequence[dict], levels: Sequence[float], title: str = "Likelihood contour") -> go.Figure:
    if not points:
        return empty_fig(title)
    df = pd.DataFrame(points)
    xs = sorted(df["x"].unique())
    ys = sorted(df["y"].unique())
    z = np.full((len(ys), len(xs)), np.nan)
    for _, r in df.iterrows():
        z[ys.index(r["y"]), xs.index(r["x"])] = r["delta_nll"]
    fig = go.Figure()
    fig.add_trace(go.Heatmap(x=xs, y=ys, z=z, colorbar=dict(title="ΔNLL"), opacity=0.78))
    if levels:
        fig.add_trace(
            go.Contour(
                x=xs,
                y=ys,
                z=z,
                contours=dict(start=float(min(levels)), end=float(max(levels)), size=1.0, coloring="none"),
                line=dict(width=2),
                showscale=False,
                name="levels",
            )
        )
    fig.update_layout(xaxis_title="p₁", yaxis_title="p₂")
    return style_fig(fig, title)
