from __future__ import annotations

from dash import dcc, html

from pyhyperiso_dash import services as svc
from pyhyperiso_dash.components import (
    card,
    data_table,
    dropdown,
    enum_options,
    field,
    graph,
    num_input,
    page_title,
    small_note,
    status_box,
    text_input,
)
from pyhyperiso.core.Common.GeneralEnum import QCDOrder
from pyhyperiso.core.Statistic.StatisticInterface import (
    ProfilingMethod,
    ContourAlgorithm,
    ProfilerMode,
)

SELECT_MODE = [
    {"label": "Observable by observable", "value": "observable"},
    {"label": "Whole decay", "value": "decay"},
]

DEFAULT_DECAY = ["B__l_l"]
DEFAULT_OBS = ["BR_BS_MUMU"]

OBS_COLUMNS = [
    {"name": "Observable", "id": "observable_label"},
    {"name": "Raw enum", "id": "observable_raw"},
    {"name": "bin low", "id": "bin_low"},
    {"name": "bin high", "id": "bin_high"},
    {"name": "order", "id": "order"},
    {"name": "registered", "id": "registered"},
]
UNC_COLUMNS = [
    {"name": "Observable", "id": "observable_label"},
    {"name": "Raw id", "id": "raw_observable"},
    {"name": "bin low", "id": "bin_low"},
    {"name": "bin high", "id": "bin_high"},
    {"name": "central", "id": "central"},
    {"name": "sigma", "id": "sigma"},
    {"name": "sigma_minus", "id": "sigma_minus"},
    {"name": "sigma_plus", "id": "sigma_plus"},
    {"name": "mu", "id": "mu"},
    {"name": "mode", "id": "mode"},
    {"name": "skew", "id": "skew"},
    {"name": "symmetric", "id": "symmetric"},
]
PSPEC_COLUMNS = [
    {"name": "Parameter", "id": "parameter", "editable": False},
    {"name": "Source", "id": "source", "editable": False},
    {"name": "type", "id": "type", "editable": False},
    {"name": "block", "id": "block", "editable": False},
    {"name": "code", "id": "code", "editable": False},
    {"name": "fit offset", "id": "fit_offset", "type": "numeric", "editable": False},
    {"name": "Wilson coefficient", "id": "wilson_coefficient", "editable": False},
    {"name": "Wilson group", "id": "wilson_group", "editable": False},
    {"name": "Wilson scan mode", "id": "wilson_scan_mode", "editable": False},
    {
        "name": "Wilson matching scale",
        "id": "wilson_matching_scale",
        "type": "numeric",
        "editable": False,
    },
    {
        "name": "Wilson hadronic scale",
        "id": "wilson_hadronic_scale",
        "type": "numeric",
        "editable": False,
    },
    {"name": "Wilson order", "id": "wilson_order", "editable": False},
    {"name": "initial", "id": "initial", "type": "numeric", "editable": False},
    {"name": "lower bound", "id": "lower_bound", "type": "numeric", "editable": True},
    {"name": "upper bound", "id": "upper_bound", "type": "numeric", "editable": True},
]
FIT_COLUMNS = [
    {"name": "Parameter", "id": "parameter"},
    {"name": "Raw", "id": "raw_parameter"},
    {"name": "best_fit", "id": "best_fit"},
    {"name": "std", "id": "std"},
]


def _described_enum_options(
    enum_cls, descriptions: dict[str, str]
) -> list[dict[str, str]]:
    return [
        {"label": descriptions.get(item.name, item.name), "value": item.name}
        for item in enum_cls
    ]


PROFILING_OPTIONS = _described_enum_options(
    ProfilingMethod,
    {
        "SLICE": "Slice — keep other fit parameters at their best fit",
        "FREE_PROJECTION": "Free projection — profile all hidden fit parameters",
        "PRIOR_CONSTRAINED_PROJECTION": "Prior-constrained projection — Gaussian constraints from the fit",
    },
)
CONTOUR_ALGORITHM_OPTIONS = _described_enum_options(
    ContourAlgorithm,
    {
        "MINUIT": "Minuit contour",
        "AMS": "Adaptive marching squares (AMS)",
    },
)
FALLBACK_ALGORITHM_OPTIONS = [
    {"label": "No fallback", "value": "NONE"},
    *CONTOUR_ALGORITHM_OPTIONS,
]
PROFILER_MODE_OPTIONS = _described_enum_options(
    ProfilerMode,
    {
        "MINUIT": "Minuit",
        "LAPLACE_NUISANCE": "Laplace nuisance approximation (Minuit fallback)",
    },
)
CONTOUR_LEVEL_OPTIONS = [
    {"label": "1σ (68.27%)", "value": 1.0},
    {"label": "2σ (95.45%)", "value": 2.0},
    {"label": "3σ (99.73%)", "value": 3.0},
]


def task_progress(component_id: str, message: str):
    """Determinate progress banner fed by the C++ monitor through polling."""
    base = component_id.removesuffix("-progress-wrap")
    return html.Div(
        [
            dcc.Store(id=f"{base}-job", storage_type="memory"),
            dcc.Interval(
                id=f"{base}-progress-poll", interval=350, disabled=True, n_intervals=0
            ),
            html.Div(
                id=component_id,
                className="task-progress-wrap",
                style={"display": "none"},
                role="status",
                **{"aria-live": "polite", "aria-busy": "true"},
                children=[
                    html.Progress(
                        id=f"{base}-progress-bar",
                        value="0",
                        max="100",
                        className="task-progress-native",
                    ),
                    html.Div(
                        className="task-progress-copy",
                        children=[
                            html.Span(
                                message,
                                id=f"{base}-progress-message",
                                className="task-progress-message",
                            ),
                            html.Span(
                                "",
                                id=f"{base}-progress-meta",
                                className="task-progress-meta",
                            ),
                        ],
                    ),
                ],
            ),
        ]
    )


def observable_selection(prefix: str):
    return html.Div(
        [
            field(
                "Selection mode",
                dropdown(f"{prefix}-mode", SELECT_MODE, value="observable"),
            ),
            field(
                "Decay",
                dropdown(
                    f"{prefix}-decays",
                    svc.decay_options(),
                    value=DEFAULT_DECAY,
                    multi=True,
                ),
            ),
            field(
                "Observables in selected decay",
                dropdown(
                    f"{prefix}-obs",
                    svc.observable_options_for_decays(DEFAULT_DECAY),
                    value=DEFAULT_OBS,
                    multi=True,
                ),
            ),
            html.Div(
                className="form-grid-2",
                children=[
                    field(
                        "QCD order",
                        dropdown(
                            f"{prefix}-order", enum_options(QCDOrder), value="NNLO"
                        ),
                    ),
                    field(
                        "Add dependencies",
                        dcc.Checklist(
                            id=f"{prefix}-deps",
                            options=[
                                {
                                    "label": "include observable dependencies",
                                    "value": "deps",
                                }
                            ],
                            value=["deps"],
                            className="checklist",
                        ),
                    ),
                ],
            ),
        ]
    )


def bin_controls(prefix: str):
    return html.Div(
        [
            field(
                "Bins",
                dcc.Checklist(
                    id=f"{prefix}-use-bin",
                    options=[
                        {"label": "apply bin to selected observables", "value": "bin"}
                    ],
                    value=[],
                    className="checklist",
                ),
            ),
            html.Div(
                id=f"{prefix}-bin-mode-wrap",
                style={"display": "none"},
                children=[
                    field(
                        "Bin grid",
                        dcc.Checklist(
                            id=f"{prefix}-smooth",
                            options=[
                                {
                                    "label": "build bins from min/max/step",
                                    "value": "smooth",
                                }
                            ],
                            value=[],
                            className="checklist",
                        ),
                    ),
                ],
            ),
            html.Div(
                id=f"{prefix}-explicit-bin-wrap",
                style={"display": "none"},
                className="form-grid-2",
                children=[
                    field("Bin low", num_input(f"{prefix}-bin-low", 1.0)),
                    field("Bin high", num_input(f"{prefix}-bin-high", 6.0)),
                ],
            ),
            html.Div(
                id=f"{prefix}-smooth-bin-wrap",
                style={"display": "none"},
                className="form-grid-3",
                children=[
                    field("Smooth min", num_input(f"{prefix}-smooth-min", 1.0)),
                    field("Smooth max", num_input(f"{prefix}-smooth-max", 6.0)),
                    field("Step", num_input(f"{prefix}-smooth-step", 1.0)),
                ],
            ),
        ]
    )


def pspec_controls():
    return html.Div(
        [
            field(
                "Fit parameter mode",
                dcc.RadioItems(
                    id="stat-fit-parameter-mode",
                    options=[
                        {"label": "Standard parameters", "value": "STANDARD"},
                        {"label": "Wilson Scan", "value": "WILSON"},
                    ],
                    value="STANDARD",
                    className="checklist",
                ),
            ),
            html.Div(
                id="stat-standard-pspec-panel",
                children=[
                    html.Div(
                        className="form-grid-3",
                        children=[
                            field(
                                "ParameterType",
                                dropdown(
                                    "stat-pspec-type",
                                    svc.stat_parameter_type_options(),
                                    value=svc.default_parameter_type_name("FLAVOR"),
                                ),
                            ),
                            field(
                                "Block",
                                dropdown(
                                    "stat-pspec-block",
                                    [],
                                    value=None,
                                    placeholder="Choose a block...",
                                ),
                            ),
                            field(
                                "Code",
                                dropdown(
                                    "stat-pspec-code",
                                    [],
                                    value=None,
                                    placeholder="Choose a code...",
                                ),
                            ),
                        ],
                    ),
                    html.Button(
                        "Add fit parameter", id="stat-add-pspec-btn", n_clicks=0
                    ),
                ],
            ),
            html.Div(
                id="stat-wilson-scan-panel",
                style={"display": "none"},
                children=[
                    small_note(
                        "Select physical coefficients directly. The GUI builds every required Wilson group first, then converts the selection internally to the ParamId expected by the statistic core.",
                        tone="info",
                    ),
                    field(
                        "Wilson coefficients",
                        dropdown(
                            "stat-wilson-scan-coefficients",
                            svc.wilson_scan_coefficient_options(),
                            value=["C9", "C10"],
                            multi=True,
                            placeholder="Choose up to 10 coefficients...",
                            clearable=False,
                        ),
                    ),
                    field(
                        "Scan convention",
                        dcc.RadioItems(
                            id="stat-wilson-scan-kind",
                            options=[
                                {
                                    "label": "ΔC — BSM contribution only",
                                    "value": "DELTA",
                                },
                                {
                                    "label": "Full C — total coefficient",
                                    "value": "FULL",
                                },
                            ],
                            value="DELTA",
                            className="checklist",
                        ),
                    ),
                    html.Div(
                        className="form-grid-3",
                        children=[
                            field(
                                "Matching scale [GeV]",
                                num_input("stat-wilson-scan-matching-scale", 160.0),
                            ),
                            field(
                                "Hadronic scale [GeV]",
                                num_input("stat-wilson-scan-hadronic-scale", 4.8),
                            ),
                            field(
                                "Wilson evolution order",
                                dropdown(
                                    "stat-wilson-scan-order",
                                    enum_options(QCDOrder),
                                    value="NNLO",
                                    clearable=False,
                                ),
                            ),
                        ],
                    ),
                    html.Button(
                        "Use selected Wilson coefficients",
                        id="stat-apply-wilson-pspec-btn",
                        n_clicks=0,
                    ),
                    status_box(
                        "stat-wilson-scan-status", "Wilson Scan is not configured yet."
                    ),
                ],
            ),
        ]
    )


def expert_controls():
    return html.Div(
        [
            field(
                "Expert controls",
                dcc.Checklist(
                    id="stat-advanced-mode",
                    options=[
                        {
                            "label": "show decay and advanced χ² controls",
                            "value": "advanced",
                        }
                    ],
                    value=[],
                    className="checklist",
                ),
            ),
            html.Div(
                id="stat-expert-panel",
                style={"display": "none"},
                children=[
                    small_note(
                        "Decay configuration is applied to the shared ObservableInterface used by both Observable and Stat pages.",
                        tone="info",
                    ),
                    field(
                        "Decay configuration target",
                        dropdown(
                            "stat-decay-config-decay",
                            svc.decay_options(),
                            value="B__l_l",
                            multi=False,
                        ),
                    ),
                    html.Div(id="stat-decay-config-fields"),
                    html.Div(
                        className="inline-actions",
                        children=[
                            html.Button(
                                "Apply decay configuration",
                                id="stat-apply-decay-config-btn",
                                n_clicks=0,
                            ),
                        ],
                    ),
                    status_box(
                        "stat-decay-config-status",
                        "No decay configuration applied yet.",
                    ),
                ],
            ),
        ]
    )


def contour_controls():
    return html.Div(
        [
            field(
                "Contour",
                dcc.Checklist(
                    id="stat-do-contour",
                    options=[
                        {
                            "label": "compute confidence contours after the χ² fit",
                            "value": "contour",
                        }
                    ],
                    value=["contour"],
                    className="checklist",
                ),
            ),
            html.Div(
                id="stat-contour-controls",
                children=[
                    small_note(
                        "Edit the lower/upper bounds directly in the fit-parameter table. "
                        "Choose any two fitted parameters as contour axes; the remaining parameters are reduced according to ProfilingMethod.",
                        tone="info",
                    ),
                    html.Div(
                        className="form-grid-2",
                        children=[
                            field(
                                "X-axis parameter",
                                dropdown(
                                    "stat-contour-x",
                                    [],
                                    value=None,
                                    placeholder="Add at least two fit parameters",
                                ),
                            ),
                            field(
                                "Y-axis parameter",
                                dropdown(
                                    "stat-contour-y",
                                    [],
                                    value=None,
                                    placeholder="Add at least two fit parameters",
                                ),
                            ),
                        ],
                    ),
                    html.Div(
                        className="form-grid-3",
                        children=[
                            field(
                                "Confidence levels",
                                dropdown(
                                    "stat-contour-levels",
                                    CONTOUR_LEVEL_OPTIONS,
                                    value=[1.0, 2.0],
                                    multi=True,
                                    clearable=False,
                                ),
                            ),
                            field(
                                "Resolution", num_input("stat-contour-resolution", 60)
                            ),
                            field(
                                "Profiler mode",
                                dropdown(
                                    "stat-profiler-mode",
                                    PROFILER_MODE_OPTIONS,
                                    value="LAPLACE_NUISANCE",
                                    clearable=False,
                                ),
                            ),
                        ],
                    ),
                    html.Div(
                        className="form-grid-3",
                        children=[
                            field(
                                "ProfilingMethod",
                                dropdown(
                                    "stat-profiling-method",
                                    PROFILING_OPTIONS,
                                    value="SLICE",
                                    clearable=False,
                                ),
                            ),
                            field(
                                "ContourAlgorithm",
                                dropdown(
                                    "stat-contour-algorithm",
                                    CONTOUR_ALGORITHM_OPTIONS,
                                    value="MINUIT",
                                    clearable=False,
                                ),
                            ),
                            field(
                                "Fallback algorithm",
                                dropdown(
                                    "stat-contour-fallback",
                                    FALLBACK_ALGORITHM_OPTIONS,
                                    value="AMS",
                                    clearable=False,
                                ),
                            ),
                        ],
                    ),
                    small_note(
                        "ProfilingMethod is only relevant when more than two fit parameters are present; with exactly two parameters the GUI forces SLICE.",
                        tone="warn",
                    ),
                    html.Div(id="stat-contour-method-note", className="info-box"),
                ],
            ),
        ]
    )


def layout():
    return html.Div(
        children=[
            page_title(
                "Statistics",
                "χ² Monte-Carlo covariance workflow: configure observables, propagate uncertainty, fit parameters, and extract 2D confidence contours.",
                "Stat / χ²",
            ),
            html.Div(
                className="page-grid",
                children=[
                    html.Div(
                        className="grid",
                        children=[
                            card(
                                "Statistic observables",
                                "shared ObservableInterface",
                                html.Div(
                                    [
                                        observable_selection("stat-obs"),
                                        bin_controls("stat-obs"),
                                        expert_controls(),
                                        html.Div(
                                            className="inline-actions",
                                            children=[
                                                html.Button(
                                                    "Configure Statistic ObservableInterface",
                                                    id="stat-configure-obs-btn",
                                                    n_clicks=0,
                                                ),
                                                html.Button(
                                                    "Remove selected rows",
                                                    id="stat-remove-obs-btn",
                                                    n_clicks=0,
                                                ),
                                            ],
                                        ),
                                        status_box(
                                            "stat-observable-status",
                                            svc.stat_observable_status_text(),
                                        ),
                                    ]
                                ),
                            ),
                            card(
                                "Experiments and Monte Carlo",
                                "CHI2_MC_COVARIANCE only",
                                html.Div(
                                    [
                                        small_note(
                                            "The GUI is locked to the χ² MC-covariance likelihood. The profiled-nuisance likelihood mode is not exposed here.",
                                            tone="info",
                                        ),
                                        field(
                                            "Experiments",
                                            text_input(
                                                "stat-experiments",
                                                "",
                                                placeholder="comma-separated; empty = all",
                                            ),
                                        ),
                                        html.Div(
                                            className="form-grid-3",
                                            children=[
                                                field(
                                                    "MC draws",
                                                    num_input("stat-mc-draws", 100),
                                                ),
                                                field(
                                                    "MC threads",
                                                    num_input("stat-mc-threads", 1),
                                                ),
                                                field(
                                                    "MC seed",
                                                    num_input("stat-mc-seed", 123456),
                                                ),
                                            ],
                                        ),
                                        html.Div(
                                            id="stat-advanced-options-wrap",
                                            style={"display": "none"},
                                            children=[
                                                html.Div(
                                                    className="form-grid-2",
                                                    children=[
                                                        field(
                                                            "Skew threshold",
                                                            num_input(
                                                                "stat-skew-threshold",
                                                                0.2,
                                                            ),
                                                        ),
                                                        field(
                                                            "Nuisance contexts",
                                                            num_input(
                                                                "stat-nuisance-contexts",
                                                                2,
                                                            ),
                                                        ),
                                                    ],
                                                ),
                                                html.Div(
                                                    className="form-grid-3",
                                                    children=[
                                                        field(
                                                            "Cov ridge rel",
                                                            num_input(
                                                                "stat-ridge-rel", 1e-8
                                                            ),
                                                        ),
                                                        field(
                                                            "Cov ridge abs",
                                                            num_input(
                                                                "stat-ridge-abs", 1e-12
                                                            ),
                                                        ),
                                                        field(
                                                            "Nuisance seed",
                                                            num_input(
                                                                "stat-nuisance-seed",
                                                                12345,
                                                            ),
                                                        ),
                                                    ],
                                                ),
                                                field(
                                                    "MC nuisance pruning",
                                                    dcc.Checklist(
                                                        id="stat-nuisance-pruning",
                                                        options=[
                                                            {
                                                                "label": "sensitivity pruning before MC covariance propagation",
                                                                "value": "prune",
                                                            }
                                                        ],
                                                        value=["prune"],
                                                        className="checklist",
                                                    ),
                                                ),
                                            ],
                                        ),
                                    ]
                                ),
                            ),
                            card(
                                "Uncertainty",
                                "GaussianSummary",
                                html.Div(
                                    [
                                        field(
                                            "Uncertainty display",
                                            dropdown(
                                                "stat-uncertainty-mode",
                                                [
                                                    {
                                                        "label": "symmetric",
                                                        "value": "sym",
                                                    },
                                                    {
                                                        "label": "asymmetric",
                                                        "value": "asym",
                                                    },
                                                ],
                                                value="sym",
                                            ),
                                        ),
                                        html.Button(
                                            "Compute uncertainties",
                                            id="stat-uncertainty-btn",
                                            n_clicks=0,
                                        ),
                                        task_progress(
                                            "stat-uncertainty-progress-wrap",
                                            "Monte-Carlo / uncertainty computation in progress…",
                                        ),
                                        status_box(
                                            "stat-uncertainty-status",
                                            "No uncertainty computation yet.",
                                        ),
                                    ]
                                ),
                            ),
                            card(
                                "Fit and confidence contour",
                                "up to 10 fit parameters; select any pair for the 2D contour",
                                html.Div(
                                    [
                                        pspec_controls(),
                                        data_table(
                                            "stat-p-specs-table",
                                            PSPEC_COLUMNS,
                                            data=[],
                                            page_size=10,
                                            editable=True,
                                            row_selectable="multi",
                                            row_deletable=False,
                                            hidden_columns=[
                                                "type",
                                                "block",
                                                "code",
                                                "fit_offset",
                                                "wilson_coefficient",
                                                "wilson_group",
                                                "wilson_scan_mode",
                                                "wilson_matching_scale",
                                                "wilson_hadronic_scale",
                                                "wilson_order",
                                            ],
                                        ),
                                        html.Div(
                                            className="inline-actions",
                                            children=[
                                                html.Button(
                                                    "Remove selected fit parameters",
                                                    id="stat-remove-pspec-btn",
                                                    n_clicks=0,
                                                ),
                                            ],
                                        ),
                                        contour_controls(),
                                        html.Button(
                                            "Run χ² fit / contour",
                                            id="stat-fit-btn",
                                            n_clicks=0,
                                        ),
                                        task_progress(
                                            "stat-fit-progress-wrap",
                                            "χ² fit / confidence-contour computation in progress…",
                                        ),
                                        status_box("stat-fit-status", "No fit yet."),
                                    ]
                                ),
                            ),
                        ],
                    ),
                    html.Div(
                        className="grid",
                        children=[
                            card(
                                "Configured stat observables",
                                "selection expanded from decay and bins",
                                data_table(
                                    "stat-observable-table",
                                    OBS_COLUMNS,
                                    data=svc.current_stat_observable_rows(),
                                    page_size=12,
                                    row_selectable="multi",
                                ),
                            ),
                            card(
                                "Uncertainty table",
                                "GaussianSummary",
                                data_table(
                                    "stat-uncertainty-table", UNC_COLUMNS, page_size=16
                                ),
                            ),
                            card(
                                "Uncertainty plot",
                                "central value plus uncertainty band",
                                graph("stat-uncertainty-fig", height=480),
                                className="card graph-card",
                            ),
                            card(
                                "Best-fit parameters",
                                "χ² maximum-likelihood estimate",
                                data_table("stat-fit-table", FIT_COLUMNS, page_size=12),
                            ),
                            card(
                                "Fit correlations",
                                "p_correlations",
                                graph("stat-corr-fig", height=460),
                                className="card graph-card",
                            ),
                            card(
                                "2D confidence contours",
                                "core ContourEngine paths",
                                graph("stat-contour-fig", height=520),
                                className="card graph-card",
                            ),
                            small_note(
                                "LaTeX labels use the supplied observable/decay/nuisance maps. Raw enum/block/code values remain visible in companion columns."
                            ),
                        ],
                    ),
                ],
            ),
        ]
    )
