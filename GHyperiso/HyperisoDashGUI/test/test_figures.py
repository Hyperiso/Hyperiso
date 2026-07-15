"""Tests for Dash figure construction."""

from __future__ import annotations

import pytest

from pyhyperiso_dash.figures import _rgba, confidence_contour_paths


def _square(scale: float) -> list[dict[str, float]]:
    return [
        {"x": -scale, "y": -scale},
        {"x": scale, "y": -scale},
        {"x": scale, "y": scale},
        {"x": -scale, "y": scale},
    ]


def test_rgba_conversion_and_validation():
    assert _rgba("#2563EB", 0.24) == "rgba(37,99,235,0.24)"
    with pytest.raises(ValueError, match="#RRGGBB"):
        _rgba("red", 0.2)


def test_confidence_contours_are_filled_and_nested_outer_first():
    fig = confidence_contour_paths(
        [
            {"sigma": 1.0, "level": 2.30, "path_id": 0, "points": _square(1.0)},
            {"sigma": 2.0, "level": 6.18, "path_id": 0, "points": _square(2.0)},
        ],
        best_fit={"x": 0.1, "y": -0.2},
    )

    assert [trace.name for trace in fig.data] == ["2σ", "1σ", "best fit"]
    outer, inner, best_fit = fig.data
    assert outer.fill == "toself"
    assert inner.fill == "toself"
    assert outer.fillcolor != inner.fillcolor
    assert outer.line.color != inner.line.color
    assert outer.x[0] == outer.x[-1]
    assert outer.y[0] == outer.y[-1]
    assert best_fit.marker.symbol == "x"


def test_confidence_contour_uses_white_publication_style_and_bounds():
    fig = confidence_contour_paths(
        [{"sigma": 1.0, "level": 2.30, "path_id": 0, "points": _square(1.0)}],
        title="Fit",
        x_title="x parameter",
        y_title="y parameter",
        bounds=(-3.0, 3.0, -2.0, 2.0),
    )

    assert fig.layout.paper_bgcolor == "white"
    assert fig.layout.plot_bgcolor == "white"
    assert fig.layout.font.color == "#111827"
    assert tuple(fig.layout.xaxis.range) == (-3.0, 3.0)
    assert tuple(fig.layout.yaxis.range) == (-2.0, 2.0)
    assert fig.layout.xaxis.title.text == "x parameter"
    assert fig.layout.yaxis.title.text == "y parameter"
    assert fig.layout.xaxis.mirror is True
    assert fig.layout.yaxis.mirror is True


def test_disconnected_paths_share_color_and_single_legend_entry():
    fig = confidence_contour_paths(
        [
            {"sigma": 1.0, "level": 2.30, "path_id": 0, "points": _square(1.0)},
            {"sigma": 1.0, "level": 2.30, "path_id": 1, "points": _square(0.5)},
        ]
    )

    assert len(fig.data) == 2
    assert fig.data[0].line.color == fig.data[1].line.color
    assert fig.data[0].fillcolor == fig.data[1].fillcolor
    assert fig.data[0].showlegend is True
    assert fig.data[1].showlegend is False
