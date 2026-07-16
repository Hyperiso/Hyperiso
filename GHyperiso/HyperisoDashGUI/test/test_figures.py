"""Tests for Dash figure construction."""

from __future__ import annotations

import pytest

from pyhyperiso_dash.figures import (
    _confidence_axis_title,
    _padded_contour_range,
    _rgba,
    confidence_contour_paths,
)


def _square(scale: float, x0: float = 0.0, y0: float = 0.0) -> list[dict[str, float]]:
    return [
        {"x": x0 - scale, "y": y0 - scale},
        {"x": x0 + scale, "y": y0 - scale},
        {"x": x0 + scale, "y": y0 + scale},
        {"x": x0 - scale, "y": y0 + scale},
    ]


def test_rgba_conversion_and_validation():
    assert _rgba("#2563EB", 0.24) == "rgba(37,99,235,0.24)"
    with pytest.raises(ValueError, match="#RRGGBB"):
        _rgba("red", 0.2)


def test_confidence_axis_titles_preserve_latex_and_upgrade_defaults():
    assert _confidence_axis_title(r"$\delta C_9/C_9^{\rm SM}$") == (
        r"$\delta C_9/C_9^{\rm SM}$"
    )
    assert _confidence_axis_title("p₁") == r"$p_1$"
    assert _confidence_axis_title("plain label") == "plain label"


def test_confidence_contours_are_filled_and_nested_outer_first():
    fig = confidence_contour_paths(
        [
            {"sigma": 1.0, "level": 2.30, "path_id": 0, "points": _square(1.0)},
            {"sigma": 2.0, "level": 6.18, "path_id": 0, "points": _square(2.0)},
        ],
        best_fit={"x": 0.1, "y": -0.2},
    )

    assert [trace.name for trace in fig.data] == [
        "2σ (95.4% CL)",
        "1σ (68.3% CL)",
        "Best fit",
    ]
    outer, inner, best_fit = fig.data
    assert outer.fill == "toself"
    assert inner.fill == "toself"
    assert outer.fillcolor != inner.fillcolor
    assert outer.line.color != inner.line.color
    assert outer.x[0] == outer.x[-1]
    assert outer.y[0] == outer.y[-1]
    assert best_fit.marker.symbol == "x"


def test_confidence_contour_uses_tight_viewport_with_parameter_bounds():
    fig = confidence_contour_paths(
        [{"sigma": 1.0, "level": 2.30, "path_id": 0, "points": _square(1.0)}],
        title="Fit",
        x_title=r"$\delta C_9/C_9^{\rm SM}$",
        y_title=r"$\delta C_{10}/C_{10}^{\rm SM}$",
        bounds=(-10.0, 10.0, -8.0, 8.0),
    )

    assert fig.layout.paper_bgcolor == "white"
    assert fig.layout.plot_bgcolor == "white"
    assert fig.layout.font.color == "#111827"
    assert tuple(fig.layout.xaxis.range) == pytest.approx((-1.28, 1.28))
    assert tuple(fig.layout.yaxis.range) == pytest.approx((-1.28, 1.28))
    assert fig.layout.xaxis.title.text == r"$\delta C_9/C_9^{\rm SM}$"
    assert fig.layout.yaxis.title.text == r"$\delta C_{10}/C_{10}^{\rm SM}$"
    assert fig.layout.xaxis.mirror is True
    assert fig.layout.yaxis.mirror is True
    assert fig.layout.xaxis.tickformat == ".4g"
    assert fig.layout.height == 590


def test_viewport_is_clipped_to_pspec_bounds_near_an_edge():
    assert _padded_contour_range([0.75, 0.9], (0.0, 1.0)) == pytest.approx(
        [0.729, 0.921]
    )
    assert _padded_contour_range([0.95, 1.0], (0.0, 1.0)) == pytest.approx([0.943, 1.0])


def test_disconnected_paths_share_color_and_single_legend_entry():
    fig = confidence_contour_paths(
        [
            {"sigma": 1.0, "level": 2.30, "path_id": 0, "points": _square(1.0)},
            {
                "sigma": 1.0,
                "level": 2.30,
                "path_id": 1,
                "points": _square(0.25, x0=2.0),
            },
        ]
    )

    assert len(fig.data) == 2
    assert fig.data[0].line.color == fig.data[1].line.color
    assert fig.data[0].fillcolor == fig.data[1].fillcolor
    assert fig.data[0].showlegend is True
    assert fig.data[1].showlegend is False
    assert fig.layout.xaxis.range[1] > 2.25
