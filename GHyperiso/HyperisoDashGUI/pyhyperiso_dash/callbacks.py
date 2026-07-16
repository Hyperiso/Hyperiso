from __future__ import annotations

from dash import ALL, Input, Output, State, callback_context, html, no_update
from dash.exceptions import PreventUpdate

from pyhyperiso_dash import services as svc
from pyhyperiso_dash.components import dropdown, field, metric, num_input, small_note
from pyhyperiso_dash.figures import (
    block_size_bar,
    block_values_scatter,
    dependency_graph,
    correlation_heatmap,
    empty_fig,
    heatmap_2d,
    confidence_contour_paths,
    series_1d,
    uncertainty_fig,
)


def _err(prefix: str, exc: Exception) -> str:
    return f"{prefix}\nERROR: {type(exc).__name__}: {exc}"


def _ok(prefix: str, data) -> str:
    return f"{prefix}\n{data}"


def _progress_value(value) -> str:
    """Clamp a native HTML progress value and serialize it as a string.

    Dash's generated ``html.Progress`` prop metadata declares ``value`` as a
    string in several supported releases, even though browsers accept numeric
    attributes. Returning a string avoids client-side prop validation errors.
    """
    try:
        number = max(0.0, min(100.0, float(value)))
    except (TypeError, ValueError):
        number = 0.0
    return f"{number:.3f}".rstrip("0").rstrip(".") or "0"


def _runtime_metrics():
    s = svc.runtime_summary()
    return [
        metric(
            "Runtime",
            "ON" if s["initialized"] else "OFF",
            "Hyperiso singleton",
            "good" if s["initialized"] else "bad",
        ),
        metric("Model", s["model"], "active config"),
        metric("LHA", "set" if s["lha_path"] != "—" else "—", s["lha_path"]),
        metric(
            "Wilson",
            "built" if s["wilson_built"] else "not built",
            "current session",
            "good" if s["wilson_built"] else "",
        ),
        metric(
            "Observables",
            "built" if s["observable_built"] else "not built",
            "current session",
            "good" if s["observable_built"] else "",
        ),
    ]


def _stat_kwargs(
    mode,
    obs,
    decays,
    order,
    deps,
    bin_strategy,
    bin_low,
    bin_high,
    smooth_min,
    smooth_max,
    smooth_step,
    experiments,
    mc_draws,
    mc_threads,
    mc_seed,
    skew,
    ridge_rel,
    ridge_abs,
    prune,
    contexts,
    seed,
):
    return dict(
        mode=mode,
        obs_names=obs,
        decay_names=decays,
        order_name=order,
        add_deps="deps" in (deps or []),
        bin_strategy=bin_strategy,
        bin_low=bin_low,
        bin_high=bin_high,
        smooth_min=smooth_min,
        smooth_max=smooth_max,
        smooth_step=smooth_step,
        experiments=experiments,
        mc_draws=mc_draws,
        mc_threads=mc_threads,
        mc_seed=mc_seed,
        skew_threshold=skew,
        ridge_rel=ridge_rel,
        ridge_abs=ridge_abs,
        nuisance_pruning="prune" in (prune or []),
        nuisance_contexts=contexts,
        nuisance_seed=seed,
    )


def register_callbacks(app):
    # ---------- Cascading parameter and Wilson pickers ----------
    parameter_prefixes = [
        "wilson-x-param",
        "wilson-y-param",
        "obs-x-param",
        "obs-y-param",
    ]
    for prefix in parameter_prefixes:

        @app.callback(
            Output(f"{prefix}-block", "options"),
            Output(f"{prefix}-block", "value"),
            Input(f"{prefix}-ptype", "value"),
            Input("runtime-ping", "data"),
            State(f"{prefix}-block", "value"),
            prevent_initial_call=False,
        )
        def _update_param_blocks(ptype, _ping, current_block, _prefix=prefix):
            try:
                opts = svc.block_name_options(ptype) if ptype else []
                valid = {o["value"] for o in opts}
                value = (
                    current_block
                    if current_block in valid
                    else (opts[0]["value"] if opts else None)
                )
                return opts, value
            except Exception:
                return [], None

        @app.callback(
            Output(f"{prefix}-code", "options"),
            Output(f"{prefix}-code", "value"),
            Input(f"{prefix}-ptype", "value"),
            Input(f"{prefix}-block", "value"),
            Input("runtime-ping", "data"),
            State(f"{prefix}-code", "value"),
            prevent_initial_call=False,
        )
        def _update_param_codes(ptype, block, _ping, current_code, _prefix=prefix):
            try:
                opts = svc.block_code_options(ptype, block) if ptype and block else []
                valid = {o["value"] for o in opts}
                value = (
                    current_code
                    if current_code in valid
                    else (opts[0]["value"] if opts else None)
                )
                return opts, value
            except Exception:
                return [], None

    @app.callback(
        Output("stat-pspec-block", "options"),
        Output("stat-pspec-block", "value"),
        Input("stat-pspec-type", "value"),
        Input("runtime-ping", "data"),
        State("stat-pspec-block", "value"),
        prevent_initial_call=False,
    )
    def update_stat_pspec_blocks(ptype, _ping, current_block):
        try:
            opts = svc.block_name_options(ptype) if ptype else []
            valid = {o["value"] for o in opts}
            value = (
                current_block
                if current_block in valid
                else (opts[0]["value"] if opts else None)
            )
            return opts, value
        except Exception:
            return [], None

    @app.callback(
        Output("stat-pspec-code", "options"),
        Output("stat-pspec-code", "value"),
        Input("stat-pspec-type", "value"),
        Input("stat-pspec-block", "value"),
        Input("runtime-ping", "data"),
        State("stat-pspec-code", "value"),
        prevent_initial_call=False,
    )
    def update_stat_pspec_codes(ptype, block, _ping, current_code):
        try:
            opts = svc.block_code_options(ptype, block) if ptype and block else []
            valid = {o["value"] for o in opts}
            value = (
                current_code
                if current_code in valid
                else (opts[0]["value"] if opts else None)
            )
            return opts, value
        except Exception:
            return [], None

    @app.callback(
        Output("wilson-request-coeff", "options"),
        Output("wilson-request-coeff", "value"),
        Input("wilson-request-group", "value"),
        State("wilson-request-coeff", "value"),
        prevent_initial_call=False,
    )
    def update_wilson_coeffs(group, current_coeff):
        opts = svc.wilson_coeff_options_for_group(group)
        valid = {o["value"] for o in opts}
        value = (
            current_coeff
            if current_coeff in valid
            else (opts[0]["value"] if opts else None)
        )
        return opts, value

    # ---------- Decay -> observable cascading selectors ----------
    for prefix in ["obs-build", "stat-obs"]:

        @app.callback(
            Output(f"{prefix}-obs", "options"),
            Output(f"{prefix}-obs", "value"),
            Output(f"{prefix}-obs", "disabled"),
            Input(f"{prefix}-decays", "value"),
            Input(f"{prefix}-mode", "value"),
            State(f"{prefix}-obs", "value"),
            prevent_initial_call=False,
        )
        def _update_observables_for_decay(decays, mode, current, _prefix=prefix):
            opts = svc.observable_options_for_decays(decays or [])
            valid = {o["value"] for o in opts}
            current_values = (
                current if isinstance(current, list) else ([current] if current else [])
            )
            value = [x for x in current_values if x in valid]
            if not value and opts and mode != "decay":
                value = [opts[0]["value"]]
            disabled = mode == "decay"
            return opts, ([] if disabled else value), disabled

    @app.callback(
        Output("obs-decay-config-decay", "options"),
        Output("obs-decay-config-decay", "value"),
        Input("obs-build-decays", "value"),
        State("obs-decay-config-decay", "value"),
        prevent_initial_call=False,
    )
    def update_obs_decay_config_target(decays, current):
        selected = decays or []
        opts = [
            o for o in svc.decay_options() if not selected or o["value"] in selected
        ]
        valid = {o["value"] for o in opts}
        value = current if current in valid else (opts[0]["value"] if opts else None)
        return opts, value

    @app.callback(
        Output("stat-decay-config-decay", "options"),
        Output("stat-decay-config-decay", "value"),
        Input("stat-obs-decays", "value"),
        State("stat-decay-config-decay", "value"),
        prevent_initial_call=False,
    )
    def update_stat_decay_config_target(decays, current):
        selected = decays or []
        opts = [
            o for o in svc.decay_options() if not selected or o["value"] in selected
        ]
        valid = {o["value"] for o in opts}
        value = current if current in valid else (opts[0]["value"] if opts else None)
        return opts, value

    @app.callback(
        Output("obs-scan-observable", "options"),
        Output("obs-scan-observable", "value"),
        Input("obs-scan-decay", "value"),
        State("obs-scan-observable", "value"),
        prevent_initial_call=False,
    )
    def update_scan_observables(decay, current):
        opts = svc.observable_options_for_decays([decay] if decay else [])
        valid = {o["value"] for o in opts}
        value = current if current in valid else (opts[0]["value"] if opts else None)
        return opts, value

    # ---------- Progressive disclosure / expert controls ----------
    def _visible(flag: bool) -> dict:
        return {} if flag else {"display": "none"}

    for prefix in ["obs", "stat"]:

        @app.callback(
            Output(f"{prefix}-decay-config-fields", "children"),
            Input(f"{prefix}-decay-config-decay", "value"),
            prevent_initial_call=False,
        )
        def _render_decay_config_fields(decay_name, _prefix=prefix):
            specs = svc.decay_config_field_specs(decay_name)
            if not specs:
                return small_note(
                    "No typed decay configuration is exposed for this decay in the current binding.",
                    tone="warn",
                )
            children = []
            for spec in specs:
                cid = {
                    "type": "decay-config-field",
                    "prefix": _prefix,
                    "field": spec["name"],
                }
                if spec["kind"] == "enum":
                    component = dropdown(
                        cid,
                        spec.get("options", []),
                        value=spec.get("value"),
                        clearable=False,
                    )
                else:
                    component = num_input(cid, spec.get("value"))
                children.append(field(spec["label"], component))
            return html.Div(className="form-grid-2", children=children)

    @app.callback(
        Output("obs-decay-config-status", "children"),
        Input("obs-apply-decay-config-btn", "n_clicks"),
        State("obs-decay-config-decay", "value"),
        State({"type": "decay-config-field", "prefix": "obs", "field": ALL}, "id"),
        State({"type": "decay-config-field", "prefix": "obs", "field": ALL}, "value"),
        prevent_initial_call=True,
    )
    def apply_obs_decay_config(_clicks, decay_name, ids, values):
        try:
            field_values = {
                str(i.get("field")): v for i, v in zip(ids or [], values or [])
            }
            info = svc.apply_decay_config(decay_name, field_values)
            return f"Applied {info['config']} to {info['decay']}: {info['values']}"
        except Exception as exc:
            return _err("Decay configuration failed.", exc)

    @app.callback(
        Output("stat-decay-config-status", "children"),
        Input("stat-apply-decay-config-btn", "n_clicks"),
        State("stat-decay-config-decay", "value"),
        State({"type": "decay-config-field", "prefix": "stat", "field": ALL}, "id"),
        State({"type": "decay-config-field", "prefix": "stat", "field": ALL}, "value"),
        prevent_initial_call=True,
    )
    def apply_stat_decay_config(_clicks, decay_name, ids, values):
        try:
            field_values = {
                str(i.get("field")): v for i, v in zip(ids or [], values or [])
            }
            info = svc.apply_decay_config(decay_name, field_values)
            return f"Applied {info['config']} to {info['decay']}: {info['values']}"
        except Exception as exc:
            return _err("Decay configuration failed.", exc)

    @app.callback(
        Output("obs-expert-panel", "style"),
        Input("obs-advanced-mode", "value"),
        prevent_initial_call=False,
    )
    def toggle_obs_expert(values):
        return _visible("advanced" in (values or []))

    @app.callback(
        Output("stat-expert-panel", "style"),
        Output("stat-advanced-options-wrap", "style"),
        Input("stat-advanced-mode", "value"),
        prevent_initial_call=False,
    )
    def toggle_stat_advanced(values):
        show = "advanced" in (values or [])
        return _visible(show), _visible(show)

    @app.callback(
        Output("obs-build-bin-mode-wrap", "style"),
        Output("obs-build-explicit-bin-wrap", "style"),
        Output("obs-build-smooth-bin-wrap", "style"),
        Input("obs-build-use-bin", "value"),
        Input("obs-build-smooth", "value"),
        prevent_initial_call=False,
    )
    def toggle_obs_build_bins(use_bin, smooth):
        enabled = "bin" in (use_bin or [])
        is_smooth = "smooth" in (smooth or [])
        return (
            _visible(enabled),
            _visible(enabled and not is_smooth),
            _visible(enabled and is_smooth),
        )

    @app.callback(
        Output("stat-obs-bin-mode-wrap", "style"),
        Output("stat-obs-explicit-bin-wrap", "style"),
        Output("stat-obs-smooth-bin-wrap", "style"),
        Input("stat-obs-use-bin", "value"),
        Input("stat-obs-smooth", "value"),
        prevent_initial_call=False,
    )
    def toggle_stat_bins(use_bin, smooth):
        enabled = "bin" in (use_bin or [])
        is_smooth = "smooth" in (smooth or [])
        return (
            _visible(enabled),
            _visible(enabled and not is_smooth),
            _visible(enabled and is_smooth),
        )

    @app.callback(
        Output("obs-scan-bin-wrap", "style"),
        Input("obs-scan-use-bin", "value"),
        prevent_initial_call=False,
    )
    def toggle_obs_scan_bin(use_bin):
        return _visible("bin" in (use_bin or []))

    for prefix in ["wilson", "obs"]:

        @app.callback(
            Output(f"{prefix}-y-range-wrap", "style"),
            Output(f"{prefix}-y-param-wrap", "style"),
            Input(f"{prefix}-scan-dim", "value"),
            prevent_initial_call=False,
        )
        def _toggle_y_scan(dim, _prefix=prefix):
            show = dim == "2d"
            return _visible(show), _visible(show)

    # ---------- Core ----------
    @app.callback(
        Output("core-lha-path", "value"),
        Input("core-upload", "contents"),
        State("core-upload", "filename"),
        prevent_initial_call=True,
    )
    def upload_lha(contents, filename):
        if not contents:
            raise PreventUpdate
        return svc.save_uploaded_lha(contents, filename)

    @app.callback(
        Output("core-marty-name-wrap", "style"),
        Output("core-marty-path-wrap", "style"),
        Input("core-model", "value"),
    )
    def show_marty_fields(model):
        visible = {} if model == "MARTY" else {"display": "none"}
        return visible, visible

    @app.callback(
        Output("core-block-ptype", "options"),
        Output("core-block-ptype", "value"),
        Input("runtime-ping", "data"),
        State("core-block-ptype", "value"),
        prevent_initial_call=True,
    )
    def refresh_core_parameter_type_options(_ping, current_value):
        options = svc.parameter_type_options()
        valid_values = {opt["value"] for opt in options}
        value = (
            current_value
            if current_value in valid_values
            else svc.default_parameter_type_name("SM")
        )
        return options, value

    @app.callback(
        Output("core-status", "children"),
        Output("core-metrics", "children"),
        Output("runtime-ping", "data"),
        Input("core-init-btn", "n_clicks_timestamp"),
        State("core-lha-path", "value"),
        State("core-flags", "value"),
        State("core-model", "value"),
        State("core-marty-name", "value"),
        State("core-marty-path", "value"),
        State("runtime-ping", "data"),
        prevent_initial_call=True,
    )
    def core_init(click_ts, lha_path, flags, model, marty_name, marty_path, ping):
        if click_ts in (None, -1) or int(click_ts) <= int(
            getattr(svc.RUNTIME, "last_core_init_ts", -1)
        ):
            raise PreventUpdate
        try:
            svc.RUNTIME.last_core_init_ts = int(click_ts)
            info = svc.init_or_switch_hyperiso(
                lha_path, flags or [], model, marty_name, marty_path
            )
            return _ok("Hyperiso ready.", info), _runtime_metrics(), int(ping or 0) + 1
        except Exception as exc:
            return (
                _err("Hyperiso initialization failed.", exc),
                _runtime_metrics(),
                ping,
            )

    @app.callback(
        Output("core-block-inventory-fig", "figure"),
        Output("core-block-status", "children"),
        Output("core-block-name", "options"),
        Output("core-block-name", "value"),
        Input("core-refresh-blocks", "n_clicks"),
        Input("core-block-ptype", "value"),
        State("core-block-name", "value"),
        prevent_initial_call=True,
    )
    def refresh_blocks(_, ptype, current_block):
        try:
            inventory = svc.block_inventory_rows()
            opts = svc.block_name_options(ptype)
            valid_values = {x["value"] for x in opts}
            value = (
                current_block
                if current_block in valid_values
                else (opts[0]["value"] if opts else None)
            )
            return (
                block_size_bar(inventory),
                f"Loaded {len(opts)} blocks for {ptype}.",
                opts,
                value,
            )
        except Exception as exc:
            return (
                empty_fig("Block inventory unavailable"),
                _err("Block refresh failed.", exc),
                [],
                None,
            )

    @app.callback(
        Output("core-block-table", "data"),
        Output("core-block-values-fig", "figure"),
        Input("core-block-name", "value"),
        State("core-block-ptype", "value"),
        prevent_initial_call=True,
    )
    def load_block(block_name, ptype):
        try:
            rows = svc.block_table_rows(ptype, block_name)
            return rows, block_values_scatter(rows, f"{ptype}:{block_name}")
        except Exception as exc:
            return [], empty_fig(f"Cannot load {ptype}:{block_name}: {exc}")

    @app.callback(
        Output("core-dep-code", "options"),
        Output("core-dep-code", "value"),
        Input("core-block-ptype", "value"),
        Input("core-block-name", "value"),
        Input("runtime-ping", "data"),
        State("core-dep-code", "value"),
        prevent_initial_call=False,
    )
    def update_core_dependency_codes(ptype, block, _ping, current_code):
        try:
            opts = svc.block_code_options(ptype, block) if ptype and block else []
            valid = {o["value"] for o in opts}
            value = (
                current_code
                if current_code in valid
                else (opts[0]["value"] if opts else None)
            )
            return opts, value
        except Exception:
            return [], None

    @app.callback(
        Output("core-dep-status", "children"),
        Output("core-dep-table", "data"),
        Output("core-dep-fig", "figure"),
        Input("core-dep-refresh-btn", "n_clicks"),
        Input("core-dep-apply-btn", "n_clicks"),
        Input("core-block-name", "value"),
        State("core-block-ptype", "value"),
        State("core-dep-action", "value"),
        State("core-dep-scope", "value"),
        State("core-dep-code", "value"),
        prevent_initial_call=True,
    )
    def refresh_or_prune_dependencies(
        _refresh, _apply, block_name, ptype, action, scope, code
    ):
        try:
            trig = (
                callback_context.triggered[0]["prop_id"].split(".")[0]
                if callback_context.triggered
                else ""
            )
            if trig == "core-dep-apply-btn":
                data = svc.prune_dependency(action, scope, ptype, block_name, code)
                last = data.get("last_action", {})
                msg = f"Applied {last.get('action')} on {last.get('scope')} {ptype}:{block_name}"
                if last.get("code"):
                    msg += f":{last.get('code')}"
            else:
                data = svc.block_dependency_data(ptype, block_name)
                msg = f"Loaded dependency component for {ptype}:{block_name}."
            msg += f" Dependent block: {data.get('is_dependent')}. Nodes: {len(data.get('nodes', []))}. Edges: {len(data.get('edges', []))}."
            return (
                msg,
                data.get("rows", []),
                dependency_graph(data, f"Dependencies for {ptype}:{block_name}"),
            )
        except Exception as exc:
            return (
                _err("Dependency inspection failed.", exc),
                [],
                empty_fig("Dependency graph unavailable"),
            )

    # ---------- Wilson ----------
    @app.callback(
        Output("wilson-build-status", "children"),
        Input("wilson-build-btn", "n_clicks"),
        Input("wilson-add-btn", "n_clicks"),
        State("wilson-groups", "value"),
        State("wilson-matching-scale", "value"),
        State("wilson-hadronic-scale", "value"),
        State("wilson-build-order", "value"),
        prevent_initial_call=True,
    )
    def build_wilson(n_build, n_add, groups, matching, hadronic, order):
        try:
            triggered = callback_context.triggered[0]["prop_id"].split(".")[0]
            add = triggered == "wilson-add-btn"
            info = svc.build_wilson(groups or [], matching, hadronic, order, add=add)
            return _ok("Wilson build succeeded.", info)
        except Exception as exc:
            return _err("Wilson build failed.", exc)

    @app.callback(
        Output("wilson-result-table", "data"),
        Output("wilson-query-status", "children"),
        Input("wilson-query-btn", "n_clicks"),
        State("wilson-method", "value"),
        State("wilson-request-group", "value"),
        State("wilson-request-coeff", "value"),
        State("wilson-request-order", "value"),
        State("wilson-contribution", "value"),
        State("wilson-basis", "value"),
        State("wilson-component", "value"),
        prevent_initial_call=True,
    )
    def query_wilson(_, method, group, coeff, order, contribution, basis, component):
        try:
            row = svc.query_wilson(
                method, group, coeff, order, contribution, basis, component or "real"
            )
            return [row], _ok("Wilson request succeeded.", row)
        except Exception as exc:
            return [], _err("Wilson request failed.", exc)

    @app.callback(
        Output("wilson-scan-fig", "figure"),
        Output("wilson-scan-status", "children"),
        Input("wilson-scan-btn", "n_clicks"),
        State("wilson-scan-dim", "value"),
        State("wilson-method", "value"),
        State("wilson-request-group", "value"),
        State("wilson-request-coeff", "value"),
        State("wilson-request-order", "value"),
        State("wilson-contribution", "value"),
        State("wilson-basis", "value"),
        State("wilson-component", "value"),
        State("wilson-x-param-ptype", "value"),
        State("wilson-x-param-block", "value"),
        State("wilson-x-param-code", "value"),
        State("wilson-x-min", "value"),
        State("wilson-x-max", "value"),
        State("wilson-x-n", "value"),
        State("wilson-y-param-ptype", "value"),
        State("wilson-y-param-block", "value"),
        State("wilson-y-param-code", "value"),
        State("wilson-y-min", "value"),
        State("wilson-y-max", "value"),
        State("wilson-y-n", "value"),
        prevent_initial_call=True,
    )
    def wilson_scan(
        _,
        dim,
        method,
        group,
        coeff,
        order,
        contribution,
        basis,
        component,
        x_ptype,
        x_block,
        x_code,
        x_min,
        x_max,
        x_n,
        y_ptype,
        y_block,
        y_code,
        y_min,
        y_max,
        y_n,
    ):
        try:
            component = component or "real"
            if dim == "2d":
                xs, ys, z = svc.wilson_scan_2d(
                    method,
                    group,
                    coeff,
                    order,
                    contribution,
                    basis,
                    component,
                    (x_ptype, x_block, x_code),
                    (y_ptype, y_block, y_code),
                    x_min,
                    x_max,
                    x_n,
                    y_min,
                    y_max,
                    y_n,
                )
                title = (
                    svc.wilson_request_display_label(method, coeff, order, contribution)
                    + f" [{component}]"
                )
                return heatmap_2d(
                    xs,
                    ys,
                    z,
                    title,
                    svc.parameter_display_label(
                        x_block, x_code, f"{x_ptype}:{x_block}:{x_code}", x_ptype
                    ),
                    svc.parameter_display_label(
                        y_block, y_code, f"{y_ptype}:{y_block}:{y_code}", y_ptype
                    ),
                    f"Wilson {component}",
                ), f"2D scan complete: {len(xs)}×{len(ys)} points."
            xs, ys = svc.wilson_scan_1d(
                method,
                group,
                coeff,
                order,
                contribution,
                basis,
                component,
                x_ptype,
                x_block,
                x_code,
                x_min,
                x_max,
                x_n,
            )
            title = (
                svc.wilson_request_display_label(method, coeff, order, contribution)
                + f" [{component}]"
            )
            return series_1d(
                xs,
                ys,
                title,
                svc.parameter_display_label(
                    x_block, x_code, f"{x_ptype}:{x_block}:{x_code}", x_ptype
                ),
                f"Wilson {component}",
            ), f"1D scan complete: {len(xs)} points."
        except Exception as exc:
            return empty_fig("Wilson scan failed"), _err("Wilson scan failed.", exc)

    # ---------- Observable ----------
    @app.callback(
        Output("obs-selection-table", "data"),
        Output("obs-build-status", "children"),
        Input("obs-build-btn", "n_clicks"),
        State("obs-build-mode", "value"),
        State("obs-build-obs", "value"),
        State("obs-build-decays", "value"),
        State("obs-build-order", "value"),
        State("obs-build-deps", "value"),
        State("obs-build-use-bin", "value"),
        State("obs-build-bin-low", "value"),
        State("obs-build-bin-high", "value"),
        State("obs-build-smooth", "value"),
        State("obs-build-smooth-min", "value"),
        State("obs-build-smooth-max", "value"),
        State("obs-build-smooth-step", "value"),
        prevent_initial_call=True,
    )
    def build_observables(
        _,
        mode,
        obs,
        decays,
        order,
        deps,
        use_bin,
        low,
        high,
        smooth,
        sm_min,
        sm_max,
        sm_step,
    ):
        try:
            rows = svc.build_observables(
                mode,
                obs,
                decays,
                order,
                "deps" in (deps or []),
                "bin" in (use_bin or []),
                low,
                high,
                "smooth" in (smooth or []),
                sm_min,
                sm_max,
                sm_step,
            )
            return (
                rows,
                f"Configured {len(rows)} requested entries. {svc.observable_status_text()}",
            )
        except Exception as exc:
            return [], _err("Observable configuration failed.", exc)

    @app.callback(
        Output("obs-selection-table", "data", allow_duplicate=True),
        Output("obs-selection-table", "selected_rows"),
        Output("obs-build-status", "children", allow_duplicate=True),
        Input("obs-remove-btn", "n_clicks"),
        State("obs-selection-table", "selected_rows"),
        prevent_initial_call=True,
    )
    def remove_observable_rows(_, selected_rows):
        try:
            rows = svc.remove_observable_rows(selected_rows or [])
            return rows, [], f"Removed selected rows. {svc.observable_status_text()}"
        except Exception as exc:
            return no_update, [], _err("Observable removal failed.", exc)

    @app.callback(
        Output("obs-result-table", "data"),
        Output("obs-compute-status", "children"),
        Input("obs-compute-btn", "n_clicks"),
        prevent_initial_call=True,
    )
    def compute_observables(_):
        try:
            rows = svc.compute_current_observables()
            return rows, f"Computed {len(rows)} observable values."
        except Exception as exc:
            return [], _err("Observable computation failed.", exc)

    @app.callback(
        Output("obs-scan-fig", "figure"),
        Output("obs-scan-status", "children"),
        Input("obs-scan-btn", "n_clicks"),
        State("obs-scan-dim", "value"),
        State("obs-scan-observable", "value"),
        State("obs-scan-order", "value"),
        State("obs-scan-deps", "value"),
        State("obs-scan-use-bin", "value"),
        State("obs-scan-bin-low", "value"),
        State("obs-scan-bin-high", "value"),
        State("obs-x-param-ptype", "value"),
        State("obs-x-param-block", "value"),
        State("obs-x-param-code", "value"),
        State("obs-x-min", "value"),
        State("obs-x-max", "value"),
        State("obs-x-n", "value"),
        State("obs-y-param-ptype", "value"),
        State("obs-y-param-block", "value"),
        State("obs-y-param-code", "value"),
        State("obs-y-min", "value"),
        State("obs-y-max", "value"),
        State("obs-y-n", "value"),
        prevent_initial_call=True,
    )
    def observable_scan(
        _,
        dim,
        obs,
        order,
        deps,
        use_bin,
        bin_low,
        bin_high,
        x_ptype,
        x_block,
        x_code,
        x_min,
        x_max,
        x_n,
        y_ptype,
        y_block,
        y_code,
        y_min,
        y_max,
        y_n,
    ):
        try:
            add_deps = "deps" in (deps or [])
            if "bin" not in (use_bin or []):
                bin_low = bin_high = None
            if dim == "2d":
                xs, ys, z = svc.observable_scan_2d(
                    obs,
                    order,
                    add_deps,
                    bin_low,
                    bin_high,
                    (x_ptype, x_block, x_code),
                    (y_ptype, y_block, y_code),
                    x_min,
                    x_max,
                    x_n,
                    y_min,
                    y_max,
                    y_n,
                )
                return heatmap_2d(
                    xs,
                    ys,
                    z,
                    svc.observable_display_label(obs),
                    svc.parameter_display_label(
                        x_block, x_code, f"{x_ptype}:{x_block}:{x_code}", x_ptype
                    ),
                    svc.parameter_display_label(
                        y_block, y_code, f"{y_ptype}:{y_block}:{y_code}", y_ptype
                    ),
                    "observable",
                ), f"2D observable scan complete: {len(xs)}×{len(ys)} points."
            xs, ys = svc.observable_scan_1d(
                obs,
                order,
                add_deps,
                bin_low,
                bin_high,
                x_ptype,
                x_block,
                x_code,
                x_min,
                x_max,
                x_n,
            )
            return series_1d(
                xs,
                ys,
                svc.observable_display_label(obs),
                svc.parameter_display_label(
                    x_block, x_code, f"{x_ptype}:{x_block}:{x_code}", x_ptype
                ),
                "observable",
            ), f"1D observable scan complete: {len(xs)} points."
        except Exception as exc:
            return empty_fig("Observable scan failed"), _err(
                "Observable scan failed.", exc
            )

    # ---------- Stat ----------
    @app.callback(
        Output("stat-observable-table", "data"),
        Output("stat-observable-status", "children"),
        Input("stat-configure-obs-btn", "n_clicks"),
        State("stat-obs-mode", "value"),
        State("stat-obs-obs", "value"),
        State("stat-obs-decays", "value"),
        State("stat-obs-order", "value"),
        State("stat-obs-deps", "value"),
        State("stat-obs-use-bin", "value"),
        State("stat-obs-bin-low", "value"),
        State("stat-obs-bin-high", "value"),
        State("stat-obs-smooth", "value"),
        State("stat-obs-smooth-min", "value"),
        State("stat-obs-smooth-max", "value"),
        State("stat-obs-smooth-step", "value"),
        prevent_initial_call=True,
    )
    def configure_stat_observables(
        _,
        mode,
        obs,
        decays,
        order,
        deps,
        use_bin,
        bin_low,
        bin_high,
        smooth,
        sm_min,
        sm_max,
        sm_step,
    ):
        try:
            bin_strategy = (
                "smooth"
                if ("bin" in (use_bin or []) and "smooth" in (smooth or []))
                else ("single" if "bin" in (use_bin or []) else "none")
            )
            rows = svc.configure_stat_observables(
                mode,
                obs,
                decays,
                order,
                "deps" in (deps or []),
                bin_strategy,
                bin_low,
                bin_high,
                sm_min,
                sm_max,
                sm_step,
            )
            for row in rows:
                row["registered"] = True
            return (
                rows,
                f"Configured {len(rows)} requested entries. {svc.stat_observable_status_text()}",
            )
        except Exception as exc:
            return [], _err("Statistic observable configuration failed.", exc)

    @app.callback(
        Output("stat-observable-table", "data", allow_duplicate=True),
        Output("stat-observable-table", "selected_rows"),
        Output("stat-observable-status", "children", allow_duplicate=True),
        Input("stat-remove-obs-btn", "n_clicks"),
        State("stat-observable-table", "selected_rows"),
        prevent_initial_call=True,
    )
    def remove_stat_observable_rows(_, selected_rows):
        try:
            rows = svc.remove_stat_observable_rows(selected_rows or [])
            for row in rows:
                row["registered"] = True
            return (
                rows,
                [],
                f"Removed selected rows. {svc.stat_observable_status_text()}",
            )
        except Exception as exc:
            return no_update, [], _err("Statistic observable removal failed.", exc)

    @app.callback(
        Output("stat-standard-pspec-panel", "style"),
        Output("stat-wilson-scan-panel", "style"),
        Input("stat-fit-parameter-mode", "value"),
        prevent_initial_call=False,
    )
    def toggle_stat_fit_parameter_mode(mode):
        if str(mode or "STANDARD").upper() == "WILSON":
            return {"display": "none"}, {}
        return {}, {"display": "none"}

    @app.callback(
        Output("stat-p-specs-table", "data"),
        Output("stat-p-specs-table", "selected_rows"),
        Output("stat-wilson-scan-status", "children"),
        Input("stat-add-pspec-btn", "n_clicks"),
        Input("stat-remove-pspec-btn", "n_clicks"),
        Input("stat-apply-wilson-pspec-btn", "n_clicks"),
        State("stat-fit-parameter-mode", "value"),
        State("stat-pspec-type", "value"),
        State("stat-pspec-block", "value"),
        State("stat-pspec-code", "value"),
        State("stat-wilson-scan-coefficients", "value"),
        State("stat-wilson-scan-kind", "value"),
        State("stat-wilson-scan-matching-scale", "value"),
        State("stat-wilson-scan-hadronic-scale", "value"),
        State("stat-wilson-scan-order", "value"),
        State("stat-p-specs-table", "data"),
        State("stat-p-specs-table", "selected_rows"),
        prevent_initial_call=True,
    )
    def edit_pspec_rows(
        _add,
        _remove,
        _apply_wilson,
        fit_mode,
        ptype,
        block,
        code,
        coefficients,
        scan_kind,
        matching_scale,
        hadronic_scale,
        wilson_order,
        rows,
        selected_rows,
    ):
        rows = list(rows or [])
        triggered = callback_context.triggered[0]["prop_id"].split(".")[0]
        if triggered == "stat-remove-pspec-btn":
            remove = set(selected_rows or [])
            rows = [row for i, row in enumerate(rows) if i not in remove]
            return rows[:10], [], no_update

        if triggered == "stat-apply-wilson-pspec-btn":
            try:
                setup = svc.wilson_scan_setup(
                    coefficients,
                    scan_kind,
                    matching_scale,
                    hadronic_scale,
                    wilson_order,
                )
                rows = svc.wilson_scan_parameter_rows(setup)
                groups = sorted({row.get("wilson_group") for row in rows})
                convention = (
                    "ΔC (BSM only)"
                    if setup["scan_mode"] == "DELTA"
                    else "full C (total)"
                )
                status = (
                    f"Wilson Scan configured: {len(rows)} coefficient(s), {convention}. "
                    f"Required groups ready: {', '.join(groups)}."
                )
                return rows, [], status
            except Exception as exc:
                return rows[:10], [], _err("Wilson Scan configuration failed.", exc)

        if str(fit_mode or "STANDARD").upper() != "STANDARD":
            return (
                rows[:10],
                [],
                "Use ‘Use selected Wilson coefficients’ to configure a Wilson Scan.",
            )
        if not (ptype and block and code not in (None, "")):
            return rows[:10], [], no_update
        try:
            initial, lower, upper = svc.suggested_parameter_bounds(ptype, block, code)
        except Exception:
            initial, lower, upper = 0.0, -1.0, 1.0
        # Standard mode and Wilson mode are intentionally mutually exclusive.
        rows = [row for row in rows if str(row.get("type")) != "WILSON"]
        row = {
            "parameter": svc.parameter_display_label(
                block, code, f"{ptype}:{block}:{code}", ptype
            ),
            "source": "Standard parameter",
            "type": ptype,
            "block": block,
            "code": str(code),
            "initial": initial,
            "lower_bound": lower,
            "upper_bound": upper,
        }
        key = (row["type"], row["block"], row["code"])
        if key not in {
            (r.get("type"), r.get("block"), str(r.get("code"))) for r in rows
        }:
            rows.append(row)
        return rows[:10], [], no_update

    @app.callback(
        Output("stat-contour-x", "options"),
        Output("stat-contour-x", "value"),
        Output("stat-contour-y", "options"),
        Output("stat-contour-y", "value"),
        Output("stat-profiling-method", "disabled"),
        Output("stat-profiling-method", "value"),
        Output("stat-contour-controls", "style"),
        Output("stat-contour-method-note", "children"),
        Input("stat-p-specs-table", "data"),
        Input("stat-do-contour", "value"),
        State("stat-contour-x", "value"),
        State("stat-contour-y", "value"),
        State("stat-profiling-method", "value"),
        prevent_initial_call=False,
    )
    def update_stat_contour_controls(
        rows, do_contour, current_x, current_y, current_method
    ):
        options = svc.p_spec_axis_options(rows)
        valid = [item["value"] for item in options]
        x_value = current_x if current_x in valid else (valid[0] if valid else None)
        y_candidates = [value for value in valid if value != x_value]
        y_value = (
            current_y
            if current_y in y_candidates
            else (y_candidates[0] if y_candidates else None)
        )
        n_params = len(options)
        method_disabled = n_params <= 2
        method_value = (
            "SLICE" if method_disabled else (current_method or "FREE_PROJECTION")
        )
        style = {} if "contour" in (do_contour or []) else {"display": "none"}
        if n_params < 2:
            note = "Add at least two fit parameters to enable a contour."
        elif n_params == 2:
            note = "Exactly two fit parameters: the contour is already 2D, so ProfilingMethod is fixed to SLICE."
        else:
            note = f"{n_params} fit parameters: choose how the {n_params - 2} hidden parameter(s) are reduced to the selected 2D plane."
        return (
            options,
            x_value,
            options,
            y_value,
            method_disabled,
            method_value,
            style,
            note,
        )

    @app.callback(
        Output("stat-contour-fallback", "options"),
        Output("stat-contour-fallback", "value"),
        Input("stat-contour-algorithm", "value"),
        State("stat-contour-fallback", "value"),
        prevent_initial_call=False,
    )
    def update_contour_fallback_options(primary, current):
        labels = {
            "MINUIT": "Minuit contour",
            "AMS": "Adaptive marching squares (AMS)",
        }
        options = [{"label": "No fallback", "value": "NONE"}]
        options.extend(
            {"label": labels.get(item.name, item.name), "value": item.name}
            for item in svc.ContourAlgorithm
            if item.name != primary
        )
        valid = {item["value"] for item in options}
        if current in valid:
            value = current
        else:
            alternatives = [
                item["value"] for item in options if item["value"] != "NONE"
            ]
            value = alternatives[0] if alternatives else "NONE"
        return options, value

    @app.callback(
        Output("stat-uncertainty-table", "data"),
        Output("stat-uncertainty-fig", "figure"),
        Output("stat-uncertainty-status", "children"),
        Output("stat-uncertainty-job", "data"),
        Output("stat-uncertainty-progress-poll", "disabled"),
        Output("stat-uncertainty-btn", "disabled"),
        Output("stat-uncertainty-progress-wrap", "style"),
        Output("stat-uncertainty-progress-bar", "value"),
        Output("stat-uncertainty-progress-message", "children"),
        Output("stat-uncertainty-progress-meta", "children"),
        Input("stat-uncertainty-btn", "n_clicks"),
        Input("stat-uncertainty-progress-poll", "n_intervals"),
        State("stat-uncertainty-job", "data"),
        State("stat-obs-mode", "value"),
        State("stat-obs-obs", "value"),
        State("stat-obs-decays", "value"),
        State("stat-obs-order", "value"),
        State("stat-obs-deps", "value"),
        State("stat-obs-use-bin", "value"),
        State("stat-obs-bin-low", "value"),
        State("stat-obs-bin-high", "value"),
        State("stat-obs-smooth", "value"),
        State("stat-obs-smooth-min", "value"),
        State("stat-obs-smooth-max", "value"),
        State("stat-obs-smooth-step", "value"),
        State("stat-experiments", "value"),
        State("stat-mc-draws", "value"),
        State("stat-mc-threads", "value"),
        State("stat-mc-seed", "value"),
        State("stat-skew-threshold", "value"),
        State("stat-ridge-rel", "value"),
        State("stat-ridge-abs", "value"),
        State("stat-nuisance-pruning", "value"),
        State("stat-nuisance-contexts", "value"),
        State("stat-nuisance-seed", "value"),
        State("stat-uncertainty-mode", "value"),
        State("stat-observable-table", "data"),
        State("stat-fit-parameter-mode", "value"),
        State("stat-wilson-scan-coefficients", "value"),
        State("stat-wilson-scan-kind", "value"),
        State("stat-wilson-scan-matching-scale", "value"),
        State("stat-wilson-scan-hadronic-scale", "value"),
        State("stat-wilson-scan-order", "value"),
        prevent_initial_call=True,
    )
    def stat_uncertainty(
        _clicks,
        _ticks,
        job_data,
        mode,
        obs,
        decays,
        order,
        deps,
        use_bin,
        bin_low,
        bin_high,
        smooth,
        sm_min,
        sm_max,
        sm_step,
        experiments,
        mc_draws,
        mc_threads,
        mc_seed,
        skew,
        ridge_rel,
        ridge_abs,
        prune,
        contexts,
        seed,
        uncertainty_mode,
        configured_rows,
        fit_mode,
        coefficients,
        scan_kind,
        matching_scale,
        hadronic_scale,
        wilson_order,
    ):
        triggered = callback_context.triggered[0]["prop_id"].split(".")[0]
        if triggered == "stat-uncertainty-btn":
            try:
                bin_strategy = (
                    "smooth"
                    if ("bin" in (use_bin or []) and "smooth" in (smooth or []))
                    else ("single" if "bin" in (use_bin or []) else "none")
                )
                kwargs = _stat_kwargs(
                    mode,
                    obs,
                    decays,
                    order,
                    deps,
                    bin_strategy,
                    bin_low,
                    bin_high,
                    sm_min,
                    sm_max,
                    sm_step,
                    experiments,
                    mc_draws,
                    mc_threads,
                    mc_seed,
                    skew,
                    ridge_rel,
                    ridge_abs,
                    prune,
                    contexts,
                    seed,
                )
                if configured_rows:
                    kwargs["configured_rows"] = configured_rows
                wilson_setup = None
                if str(fit_mode or "STANDARD").upper() == "WILSON":
                    wilson_setup = svc.wilson_scan_setup(
                        coefficients,
                        scan_kind,
                        matching_scale,
                        hadronic_scale,
                        wilson_order,
                    )
                job_id = svc.start_uncertainty_job(kwargs, wilson_setup)
                return (
                    no_update,
                    no_update,
                    "Uncertainty computation started.",
                    {"job_id": job_id},
                    False,
                    True,
                    {"display": "grid"},
                    _progress_value(0),
                    "Preparing the statistic workflow…",
                    "",
                )
            except Exception as exc:
                return (
                    [],
                    empty_fig("Uncertainty failed"),
                    _err("Uncertainty computation failed.", exc),
                    no_update,
                    True,
                    False,
                    {"display": "grid"},
                    _progress_value(100),
                    "Uncertainty computation failed",
                    str(exc),
                )

        job_id = (job_data or {}).get("job_id")
        if not job_id:
            raise PreventUpdate
        try:
            state = svc.poll_statistic_job(job_id)
        except Exception as exc:
            return (
                no_update,
                no_update,
                _err("Progress polling failed.", exc),
                job_data,
                True,
                False,
                {"display": "grid"},
                _progress_value(100),
                "Progress polling failed",
                str(exc),
            )
        if not state["done"]:
            return (
                no_update,
                no_update,
                no_update,
                job_data,
                False,
                True,
                {"display": "grid"},
                _progress_value(state["percent"]),
                state["message"],
                state["meta"],
            )
        if state["error"]:
            return (
                [],
                empty_fig("Uncertainty failed"),
                f"Uncertainty computation failed.\nERROR: {state['error']}",
                job_data,
                True,
                False,
                {"display": "grid"},
                _progress_value(100),
                "Uncertainty computation failed",
                state["error"],
            )
        rows = state["result"] or []
        return (
            rows,
            uncertainty_fig(rows, asymmetric=(uncertainty_mode == "asym")),
            f"Computed uncertainties for {len(rows)} observable/bin entries.",
            job_data,
            True,
            False,
            {"display": "grid"},
            _progress_value(100),
            "Uncertainty propagation complete",
            state["meta"],
        )

    @app.callback(
        Output("stat-fit-table", "data"),
        Output("stat-corr-fig", "figure"),
        Output("stat-contour-fig", "figure"),
        Output("stat-fit-status", "children"),
        Output("stat-fit-job", "data"),
        Output("stat-fit-progress-poll", "disabled"),
        Output("stat-fit-btn", "disabled"),
        Output("stat-fit-progress-wrap", "style"),
        Output("stat-fit-progress-bar", "value"),
        Output("stat-fit-progress-message", "children"),
        Output("stat-fit-progress-meta", "children"),
        Input("stat-fit-btn", "n_clicks"),
        Input("stat-fit-progress-poll", "n_intervals"),
        State("stat-fit-job", "data"),
        State("stat-obs-mode", "value"),
        State("stat-obs-obs", "value"),
        State("stat-obs-decays", "value"),
        State("stat-obs-order", "value"),
        State("stat-obs-deps", "value"),
        State("stat-obs-use-bin", "value"),
        State("stat-obs-bin-low", "value"),
        State("stat-obs-bin-high", "value"),
        State("stat-obs-smooth", "value"),
        State("stat-obs-smooth-min", "value"),
        State("stat-obs-smooth-max", "value"),
        State("stat-obs-smooth-step", "value"),
        State("stat-experiments", "value"),
        State("stat-mc-draws", "value"),
        State("stat-mc-threads", "value"),
        State("stat-mc-seed", "value"),
        State("stat-skew-threshold", "value"),
        State("stat-ridge-rel", "value"),
        State("stat-ridge-abs", "value"),
        State("stat-nuisance-pruning", "value"),
        State("stat-nuisance-contexts", "value"),
        State("stat-nuisance-seed", "value"),
        State("stat-p-specs-table", "data"),
        State("stat-observable-table", "data"),
        State("stat-do-contour", "value"),
        State("stat-contour-x", "value"),
        State("stat-contour-y", "value"),
        State("stat-contour-levels", "value"),
        State("stat-profiling-method", "value"),
        State("stat-contour-algorithm", "value"),
        State("stat-contour-fallback", "value"),
        State("stat-profiler-mode", "value"),
        State("stat-contour-resolution", "value"),
        State("stat-fit-parameter-mode", "value"),
        State("stat-wilson-scan-coefficients", "value"),
        State("stat-wilson-scan-kind", "value"),
        State("stat-wilson-scan-matching-scale", "value"),
        State("stat-wilson-scan-hadronic-scale", "value"),
        State("stat-wilson-scan-order", "value"),
        prevent_initial_call=True,
    )
    def stat_fit(
        _clicks,
        _ticks,
        job_data,
        mode,
        obs,
        decays,
        order,
        deps,
        use_bin,
        bin_low,
        bin_high,
        smooth,
        sm_min,
        sm_max,
        sm_step,
        experiments,
        mc_draws,
        mc_threads,
        mc_seed,
        skew,
        ridge_rel,
        ridge_abs,
        prune,
        contexts,
        seed,
        p_rows_in,
        configured_rows,
        do_contour,
        x_key,
        y_key,
        levels,
        profiling_method,
        contour_algorithm,
        fallback_algorithm,
        profiler_mode,
        resolution,
        fit_mode,
        coefficients,
        scan_kind,
        matching_scale,
        hadronic_scale,
        wilson_order,
    ):
        triggered = callback_context.triggered[0]["prop_id"].split(".")[0]
        if triggered == "stat-fit-btn":
            try:
                bin_strategy = (
                    "smooth"
                    if ("bin" in (use_bin or []) and "smooth" in (smooth or []))
                    else ("single" if "bin" in (use_bin or []) else "none")
                )
                kwargs = _stat_kwargs(
                    mode,
                    obs,
                    decays,
                    order,
                    deps,
                    bin_strategy,
                    bin_low,
                    bin_high,
                    sm_min,
                    sm_max,
                    sm_step,
                    experiments,
                    mc_draws,
                    mc_threads,
                    mc_seed,
                    skew,
                    ridge_rel,
                    ridge_abs,
                    prune,
                    contexts,
                    seed,
                )
                if configured_rows:
                    kwargs["configured_rows"] = configured_rows
                wilson_setup = None
                if str(fit_mode or "STANDARD").upper() == "WILSON":
                    wilson_setup = svc.wilson_scan_setup(
                        coefficients,
                        scan_kind,
                        matching_scale,
                        hadronic_scale,
                        wilson_order,
                    )
                    selected = list(wilson_setup["coefficients"])
                    table_selected = [
                        str(row.get("wilson_coefficient"))
                        for row in (p_rows_in or [])
                        if row.get("wilson_coefficient")
                    ]
                    wilson_rows = [
                        row
                        for row in (p_rows_in or [])
                        if row.get("wilson_coefficient")
                    ]
                    table_mode = {
                        str(row.get("wilson_scan_mode")) for row in wilson_rows
                    }
                    table_matching = {
                        float(row.get("wilson_matching_scale")) for row in wilson_rows
                    }
                    table_hadronic = {
                        float(row.get("wilson_hadronic_scale")) for row in wilson_rows
                    }
                    table_order = {
                        str(row.get("wilson_order", "")).upper() for row in wilson_rows
                    }
                    setup_changed = (
                        table_selected != selected
                        or table_mode != {wilson_setup["scan_mode"]}
                        or table_matching != {float(wilson_setup["matching_scale"])}
                        or table_hadronic != {float(wilson_setup["hadronic_scale"])}
                        or table_order != {str(wilson_setup["order_name"]).upper()}
                    )
                    if setup_changed:
                        raise ValueError(
                            "Wilson Scan controls changed. Click ‘Use selected Wilson coefficients’ "
                            "before running the fit."
                        )
                contour_enabled = "contour" in (do_contour or [])
                job_id = svc.start_fit_job(
                    kwargs,
                    p_rows_in or [],
                    contour_enabled,
                    x_key,
                    y_key,
                    levels,
                    profiling_method,
                    contour_algorithm,
                    fallback_algorithm,
                    profiler_mode,
                    resolution,
                    wilson_setup,
                )
                return (
                    no_update,
                    no_update,
                    no_update,
                    "χ² fit started.",
                    {"job_id": job_id},
                    False,
                    True,
                    {"display": "grid"},
                    _progress_value(0),
                    "Preparing Wilson groups and statistic workflow…"
                    if wilson_setup
                    else "Preparing the χ² fit…",
                    "",
                )
            except Exception as exc:
                return (
                    [],
                    empty_fig("Fit correlations unavailable"),
                    empty_fig("Contour unavailable"),
                    _err("χ² fit/contour failed.", exc),
                    no_update,
                    True,
                    False,
                    {"display": "grid"},
                    _progress_value(100),
                    "χ² fit failed",
                    str(exc),
                )

        job_id = (job_data or {}).get("job_id")
        if not job_id:
            raise PreventUpdate
        try:
            state = svc.poll_statistic_job(job_id)
        except Exception as exc:
            return (
                no_update,
                no_update,
                no_update,
                _err("Progress polling failed.", exc),
                job_data,
                True,
                False,
                {"display": "grid"},
                _progress_value(100),
                "Progress polling failed",
                str(exc),
            )
        if not state["done"]:
            return (
                no_update,
                no_update,
                no_update,
                no_update,
                job_data,
                False,
                True,
                {"display": "grid"},
                _progress_value(state["percent"]),
                state["message"],
                state["meta"],
            )
        if state["error"]:
            return (
                [],
                empty_fig("Fit correlations unavailable"),
                empty_fig("Contour unavailable"),
                f"χ² fit/contour failed.\nERROR: {state['error']}",
                job_data,
                True,
                False,
                {"display": "grid"},
                _progress_value(100),
                "χ² fit failed",
                state["error"],
            )

        result = state["result"]
        corr_fig = correlation_heatmap(
            result["corr_rows"], "Fit-parameter correlations"
        )
        contour_enabled = "contour" in (do_contour or [])
        if result["contour_paths"]:
            x_label, y_label = result["axis_labels"]
            contour_fig = confidence_contour_paths(
                result["contour_paths"],
                "2D confidence contours",
                x_label,
                y_label,
                result["best_fit_point"],
                result["bounds"],
            )
        elif contour_enabled:
            contour_fig = empty_fig("Contour computation returned no usable path")
        else:
            contour_fig = empty_fig("Contour disabled")
        contour_status = ""
        if result["contour_paths"]:
            levels_done = sorted(
                {float(row["sigma"]) for row in result["contour_paths"]}
            )
            contour_status = f" Contours: {', '.join(f'{z:g}σ' for z in levels_done)}."
        if result["contour_errors"]:
            contour_status += " Contour warnings: " + " | ".join(
                result["contour_errors"]
            )
        warning = (
            "Warning: MC_draws > 200 may be slow. "
            if mc_draws and int(mc_draws) > 200
            else ""
        )
        status = (
            f"{warning}χ² fit done. fit_ok={result['fit_ok']}, "
            f"ell_hat={result['ell_hat']:.6g}, parameters={len(result['p_rows'])}."
            f"{contour_status}"
        )
        return (
            result["p_rows"],
            corr_fig,
            contour_fig,
            status,
            job_data,
            True,
            False,
            {"display": "grid"},
            _progress_value(100),
            "χ² fit / contours complete",
            state["meta"],
        )

    # ---------- QCD ----------
    @app.callback(
        Output("qcd-result-table", "data"),
        Output("qcd-compute-status", "children"),
        Input("qcd-compute-btn", "n_clicks"),
        State("qcd-single-scale", "value"),
        State("qcd-mass-pdg", "value"),
        State("qcd-single-mb-type", "value"),
        State("qcd-single-mt-type", "value"),
        State("qcd-include-qed", "value"),
        prevent_initial_call=True,
    )
    def qcd_compute(_, scale, pdg_id, mb_type, mt_type, include_qed):
        try:
            rows = svc.qcd_single_result_rows(
                scale,
                pdg_id,
                mb_type,
                mt_type,
                include_qed="qed" in (include_qed or []),
            )
            return rows, f"Computed QCD values at μ={float(scale):.6g} GeV."
        except Exception as exc:
            return [], _err("QCD computation failed.", exc)

    @app.callback(
        Output("qcd-alpha-fig", "figure"),
        Output("qcd-alpha-status", "children"),
        Input("qcd-alpha-scan-btn", "n_clicks"),
        State("qcd-alpha-min", "value"),
        State("qcd-alpha-max", "value"),
        State("qcd-alpha-n", "value"),
        State("qcd-alpha-mb-type", "value"),
        State("qcd-alpha-mt-type", "value"),
        prevent_initial_call=True,
    )
    def qcd_alpha_scan(_, scale_min, scale_max, n_points, mb_type, mt_type):
        try:
            xs, ys = svc.qcd_scan_alphas(
                scale_min, scale_max, n_points, mb_type, mt_type
            )
            return series_1d(
                xs, ys, r"$\alpha_s(\mu)$", r"$\mu\;[\mathrm{GeV}]$", r"$\alpha_s$"
            ), f"Computed αs scan with {len(xs)} points."
        except Exception as exc:
            return empty_fig("αs scan failed"), _err("αs scan failed.", exc)

    @app.callback(
        Output("qcd-alphaem-fig", "figure"),
        Output("qcd-alphaem-status", "children"),
        Input("qcd-alphaem-scan-btn", "n_clicks"),
        State("qcd-alphaem-min", "value"),
        State("qcd-alphaem-max", "value"),
        State("qcd-alphaem-n", "value"),
        State("qcd-alphaem-mb-type", "value"),
        State("qcd-alphaem-mt-type", "value"),
        prevent_initial_call=True,
    )
    def qcd_alphaem_scan(_, scale_min, scale_max, n_points, mb_type, mt_type):
        try:
            xs, ys = svc.qcd_scan_alphaem(
                scale_min, scale_max, n_points, mb_type, mt_type
            )
            return series_1d(
                xs,
                ys,
                r"$\alpha_{\mathrm{em}}(\mu)$",
                r"$\mu\;[\mathrm{GeV}]$",
                r"$\alpha_{\mathrm{em}}$",
            ), f"Computed αem scan with {len(xs)} points."
        except Exception as exc:
            return empty_fig("αem scan failed"), _err("αem scan failed.", exc)

    @app.callback(
        Output("qcd-mass-fig", "figure"),
        Output("qcd-mass-status", "children"),
        Input("qcd-mass-scan-btn", "n_clicks"),
        State("qcd-mass-min", "value"),
        State("qcd-mass-max", "value"),
        State("qcd-mass-n", "value"),
        State("qcd-mass-scan-pdg", "value"),
        State("qcd-mass-mb-type", "value"),
        State("qcd-mass-mt-type", "value"),
        prevent_initial_call=True,
    )
    def qcd_mass_scan(_, scale_min, scale_max, n_points, pdg_id, mb_type, mt_type):
        try:
            xs, ys = svc.qcd_scan_mass(
                scale_min, scale_max, n_points, pdg_id, mb_type, mt_type
            )
            return series_1d(
                xs, ys, rf"$m_{{{pdg_id}}}(\mu)$", r"$\mu\;[\mathrm{GeV}]$", r"$m(\mu)$"
            ), f"Computed mass scan with {len(xs)} points."
        except Exception as exc:
            return empty_fig("Mass scan failed"), _err("Mass scan failed.", exc)

    @app.callback(
        Output("qcd-constants-table", "data"),
        Input("qcd-constants-btn", "n_clicks"),
        State("qcd-constants-kind", "value"),
        prevent_initial_call=True,
    )
    def qcd_constants(_, kind):
        try:
            return svc.qcd_constants_rows(kind or "all")
        except Exception:
            return []
