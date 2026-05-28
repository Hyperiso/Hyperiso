from __future__ import annotations

from dash import Input, Output, State, callback_context, no_update
from dash.exceptions import PreventUpdate

from pyhyperiso_dash import services as svc
from pyhyperiso_dash.components import metric
from pyhyperiso_dash.figures import (
    block_size_bar,
    block_values_scatter,
    correlation_heatmap,
    empty_fig,
    heatmap_2d,
    likelihood_contour,
    series_1d,
    uncertainty_fig,
)


def _err(prefix: str, exc: Exception) -> str:
    return f"{prefix}\nERROR: {type(exc).__name__}: {exc}"


def _ok(prefix: str, data) -> str:
    return f"{prefix}\n{data}"


def _runtime_metrics():
    s = svc.runtime_summary()
    return [
        metric("Runtime", "ON" if s["initialized"] else "OFF", "Hyperiso singleton", "good" if s["initialized"] else "bad"),
        metric("Model", s["model"], "active config"),
        metric("LHA", "set" if s["lha_path"] != "—" else "—", s["lha_path"]),
        metric("Wilson", "built" if s["wilson_built"] else "not built", "current session", "good" if s["wilson_built"] else ""),
        metric("Observables", "built" if s["observable_built"] else "not built", "current session", "good" if s["observable_built"] else ""),
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
        skew_threshold=skew,
        ridge_rel=ridge_rel,
        ridge_abs=ridge_abs,
        nuisance_pruning="prune" in (prune or []),
        nuisance_contexts=contexts,
        nuisance_seed=seed,
    )


def register_callbacks(app):
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
        value = current_value if current_value in valid_values else svc.default_parameter_type_name("SM")
        return options, value

    @app.callback(
        Output("core-status", "children"),
        Output("core-metrics", "children"),
        Output("runtime-ping", "data"),
        Input("core-init-btn", "n_clicks"),
        State("core-lha-path", "value"),
        State("core-flags", "value"),
        State("core-model", "value"),
        State("core-marty-name", "value"),
        State("core-marty-path", "value"),
        State("runtime-ping", "data"),
        prevent_initial_call=True,
    )
    def core_init(n, lha_path, flags, model, marty_name, marty_path, ping):
        try:
            info = svc.init_or_switch_hyperiso(lha_path, flags or [], model, marty_name, marty_path)
            return _ok("Hyperiso ready.", info), _runtime_metrics(), int(ping or 0) + 1
        except Exception as exc:
            return _err("Hyperiso initialization failed.", exc), _runtime_metrics(), ping

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
            value = current_block if current_block in valid_values else (opts[0]["value"] if opts else None)
            return block_size_bar(inventory), f"Loaded {len(opts)} blocks for {ptype}.", opts, value
        except Exception as exc:
            return empty_fig("Block inventory unavailable"), _err("Block refresh failed.", exc), [], None

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
            row = svc.query_wilson(method, group, coeff, order, contribution, basis, component or "real")
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
    def wilson_scan(_, dim, method, group, coeff, order, contribution, basis, component, x_ptype, x_block, x_code, x_min, x_max, x_n, y_ptype, y_block, y_code, y_min, y_max, y_n):
        try:
            component = component or "real"
            if dim == "2d":
                xs, ys, z = svc.wilson_scan_2d(method, group, coeff, order, contribution, basis, component, (x_ptype, x_block, x_code), (y_ptype, y_block, y_code), x_min, x_max, x_n, y_min, y_max, y_n)
                return heatmap_2d(xs, ys, z, f"{method} {group}:{coeff} [{component}]", f"{x_ptype}:{x_block}:{x_code}", f"{y_ptype}:{y_block}:{y_code}", f"Wilson {component}"), f"2D scan complete: {len(xs)}×{len(ys)} points."
            xs, ys = svc.wilson_scan_1d(method, group, coeff, order, contribution, basis, component, x_ptype, x_block, x_code, x_min, x_max, x_n)
            return series_1d(xs, ys, f"{method} {group}:{coeff} [{component}]", f"{x_ptype}:{x_block}:{x_code}", f"Wilson {component}"), f"1D scan complete: {len(xs)} points."
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
    def build_observables(_, mode, obs, decays, order, deps, use_bin, low, high, smooth, sm_min, sm_max, sm_step):
        try:
            rows = svc.build_observables(mode, obs, decays, order, "deps" in (deps or []), "bin" in (use_bin or []), low, high, "smooth" in (smooth or []), sm_min, sm_max, sm_step)
            return rows, f"Configured {len(rows)} observable/bin entries."
        except Exception as exc:
            return [], _err("Observable configuration failed.", exc)

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
    def observable_scan(_, dim, obs, order, deps, bin_low, bin_high, x_ptype, x_block, x_code, x_min, x_max, x_n, y_ptype, y_block, y_code, y_min, y_max, y_n):
        try:
            add_deps = "deps" in (deps or [])
            if dim == "2d":
                xs, ys, z = svc.observable_scan_2d(obs, order, add_deps, bin_low, bin_high, (x_ptype, x_block, x_code), (y_ptype, y_block, y_code), x_min, x_max, x_n, y_min, y_max, y_n)
                return heatmap_2d(xs, ys, z, obs, f"{x_ptype}:{x_block}:{x_code}", f"{y_ptype}:{y_block}:{y_code}", "observable"), f"2D observable scan complete: {len(xs)}×{len(ys)} points."
            xs, ys = svc.observable_scan_1d(obs, order, add_deps, bin_low, bin_high, x_ptype, x_block, x_code, x_min, x_max, x_n)
            return series_1d(xs, ys, obs, f"{x_ptype}:{x_block}:{x_code}", "observable"), f"1D observable scan complete: {len(xs)} points."
        except Exception as exc:
            return empty_fig("Observable scan failed"), _err("Observable scan failed.", exc)

    # ---------- Stat ----------
    @app.callback(
        Output("stat-p-specs-table", "data"),
        Input("stat-add-pspec-btn", "n_clicks"),
        State("stat-p-specs-table", "data"),
        prevent_initial_call=True,
    )
    def add_pspec_row(_, rows):
        rows = list(rows or [])
        if len(rows) < 10:
            rows.append({"type": "DECAY", "block": "", "code": ""})
        return rows[:10]

    @app.callback(
        Output("stat-uncertainty-table", "data"),
        Output("stat-uncertainty-fig", "figure"),
        Output("stat-uncertainty-status", "children"),
        Input("stat-uncertainty-btn", "n_clicks"),
        State("stat-obs-mode", "value"),
        State("stat-obs-obs", "value"),
        State("stat-obs-decays", "value"),
        State("stat-obs-order", "value"),
        State("stat-obs-deps", "value"),
        State("stat-obs-bin-strategy", "value"),
        State("stat-obs-bin-low", "value"),
        State("stat-obs-bin-high", "value"),
        State("stat-obs-smooth-min", "value"),
        State("stat-obs-smooth-max", "value"),
        State("stat-obs-smooth-step", "value"),
        State("stat-experiments", "value"),
        State("stat-mc-draws", "value"),
        State("stat-skew-threshold", "value"),
        State("stat-ridge-rel", "value"),
        State("stat-ridge-abs", "value"),
        State("stat-nuisance-pruning", "value"),
        State("stat-nuisance-contexts", "value"),
        State("stat-nuisance-seed", "value"),
        State("stat-uncertainty-mode", "value"),
        prevent_initial_call=True,
    )
    def stat_uncertainty(_, mode, obs, decays, order, deps, bin_strategy, bin_low, bin_high, sm_min, sm_max, sm_step, experiments, mc_draws, skew, ridge_rel, ridge_abs, prune, contexts, seed, uncertainty_mode):
        try:
            kwargs = _stat_kwargs(mode, obs, decays, order, deps, bin_strategy, bin_low, bin_high, sm_min, sm_max, sm_step, experiments, mc_draws, skew, ridge_rel, ridge_abs, prune, contexts, seed)
            rows = svc.compute_uncertainty_rows(**kwargs)
            return rows, uncertainty_fig(rows, asymmetric=(uncertainty_mode == "asym")), f"Computed uncertainties for {len(rows)} observable/bin entries."
        except Exception as exc:
            return [], empty_fig("Uncertainty failed"), _err("Uncertainty computation failed.", exc)

    @app.callback(
        Output("stat-fit-table", "data"),
        Output("stat-eta-table", "data"),
        Output("stat-corr-fig", "figure"),
        Output("stat-contour-fig", "figure"),
        Output("stat-fit-status", "children"),
        Input("stat-fit-btn", "n_clicks"),
        State("stat-obs-mode", "value"),
        State("stat-obs-obs", "value"),
        State("stat-obs-decays", "value"),
        State("stat-obs-order", "value"),
        State("stat-obs-deps", "value"),
        State("stat-obs-bin-strategy", "value"),
        State("stat-obs-bin-low", "value"),
        State("stat-obs-bin-high", "value"),
        State("stat-obs-smooth-min", "value"),
        State("stat-obs-smooth-max", "value"),
        State("stat-obs-smooth-step", "value"),
        State("stat-experiments", "value"),
        State("stat-mc-draws", "value"),
        State("stat-skew-threshold", "value"),
        State("stat-ridge-rel", "value"),
        State("stat-ridge-abs", "value"),
        State("stat-nuisance-pruning", "value"),
        State("stat-nuisance-contexts", "value"),
        State("stat-nuisance-seed", "value"),
        State("stat-p-specs-table", "data"),
        State("stat-do-contour", "value"),
        State("stat-x-half-width", "value"),
        State("stat-y-half-width", "value"),
        State("stat-nx", "value"),
        State("stat-ny", "value"),
        prevent_initial_call=True,
    )
    def stat_fit(_, mode, obs, decays, order, deps, bin_strategy, bin_low, bin_high, sm_min, sm_max, sm_step, experiments, mc_draws, skew, ridge_rel, ridge_abs, prune, contexts, seed, p_rows_in, do_contour, xhw, yhw, nx, ny):
        try:
            if mc_draws and int(mc_draws) > 200:
                warning = "Warning: MC_draws > 200 may be slow. "
            else:
                warning = ""
            kwargs = _stat_kwargs(mode, obs, decays, order, deps, bin_strategy, bin_low, bin_high, sm_min, sm_max, sm_step, experiments, mc_draws, skew, ridge_rel, ridge_abs, prune, contexts, seed)
            result = svc.run_fit_and_scan(kwargs, p_rows_in or [], "contour" in (do_contour or []), xhw, yhw, nx, ny)
            corr_fig = correlation_heatmap(result["corr_rows"], "Fit-parameter correlations")
            contour_fig = likelihood_contour(result["scan_points"], [2.30, 6.18, 11.83], "2D ΔNLL contour") if result["scan_points"] else empty_fig("No 2D contour: provide exactly two p_specs and enable contour")
            status = f"{warning}Fit done. fit_ok={result['fit_ok']}, ell_hat={result['ell_hat']:.6g}. p={len(result['p_rows'])}, nuisances={len(result['eta_rows'])}."
            return result["p_rows"], result["eta_rows"], corr_fig, contour_fig, status
        except Exception as exc:
            return [], [], empty_fig("Fit correlations unavailable"), empty_fig("Contour unavailable"), _err("Fit/contour failed.", exc)
