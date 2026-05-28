from __future__ import annotations

import base64
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from pyhyperiso_dash.domain import (
    as_float,
    as_int,
    bins_from_step,
    code_to_display,
    enum_by_name,
    key_to_code,
    linspace,
    parse_code,
    scalar_to_display,
    scalar_to_float,
    split_csv,
)

from pyhyperiso.core.BusinessLogic.ObservableInterface import ObservableInterface
from pyhyperiso.core.Common.BinnedObservableId import BinnedObservableId
from pyhyperiso.core.Common.Configs import WilsonBuildConfig, WilsonRequest
from pyhyperiso.core.Common.GeneralEnum import (
    ContributionType,
    DataType,
    Decays,
    Model,
    Observables,
    ParameterType,
    QCDOrder,
    WCoeff,
    WGroup,
    WilsonBasis,
)
from pyhyperiso.core.Common.Mapper import DecayMapper, ObservableMapper
from pyhyperiso.core.Common.LhaID import LhaID
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Core.BlockProvider import BlockLogger
from pyhyperiso.core.Core.HyperisoConfig import ExternalFlag, HyperisoConfig
from pyhyperiso.core.Core.HyperisoMaster import HyperisoMaster
from pyhyperiso.core.Core.ParamaterProvider import ParameterProvider
from pyhyperiso.core.Core.ParameterSetter import ParameterSetter
from pyhyperiso.core.Statistic.Copula import CopulaKind
from pyhyperiso.core.Statistic.StatisticConfig import StatisticConfig, StatisticLikelihoodMode
from pyhyperiso.core.Statistic.StatisticInterface import StatisticInterface
from pyhyperiso.core.PhysicalModel.WilsonInterface import WilsonInterface


@dataclass
class RuntimeState:
    """Process-local state holding live wrappers around the C++ singleton."""

    hyp: HyperisoMaster | None = None
    config: HyperisoConfig | None = None
    lha_path: str | None = None
    initialized: bool = False
    block_logger: BlockLogger | None = None
    wilson: WilsonInterface | None = None
    observable: ObservableInterface | None = None
    stat: StatisticInterface | None = None
    lock: threading.RLock = field(default_factory=threading.RLock)
    last_wilson_config: WilsonBuildConfig | None = None
    last_observable_selection: list[dict] = field(default_factory=list)


RUNTIME = RuntimeState()
DATA_DIR = Path(__file__).resolve().parent / "data"
UPLOAD_DIR = DATA_DIR / "uploaded_lha"
UPLOAD_DIR.mkdir(parents=True, exist_ok=True)


# ---------- General helpers ----------

def require_initialized() -> None:
    if not RUNTIME.initialized or RUNTIME.hyp is None:
        raise RuntimeError("Hyperiso is not initialized. Go to Core and initialize an LHA file first.")


def make_param_id(param_type_name: str, block: str, code: Any) -> ParamId:
    pt = resolve_parameter_type(param_type_name, strict=True)
    return ParamId(pt, str(block), parse_code(code))


def param_to_label(pid: ParamId) -> str:
    ptype = pid.type.name if pid.type else "None"
    return f"{ptype}:{pid.block}:{pid.code.to_string()}"


def _cpp_lhaid(code: Any):
    """Return a bound C++ LhaID from a Python/C++/wrapper LHA code."""
    if hasattr(code, "_cpp_obj"):
        return code._cpp_obj
    return LhaID(code)._cpp_obj


def _raw_provider_value_by_pid(provider: ParameterProvider, pid: ParamId, dtype: DataType = DataType.VALUE):
    """Read a parameter through the raw pybind provider without Python float().

    The public ``ParameterProvider`` wrapper currently casts the returned C++
    ``scalar_t`` with ``float(...)``. That fails when pybind does not expose
    ``scalar_t`` as a native Python float. The Dash GUI needs the raw scalar so
    it can select real/imag/abs components itself.
    """
    return provider._cpp_obj(pid._cpp_obj, dtype.value)


def _raw_provider_value_by_block(provider: ParameterProvider, block: str, code: Any, dtype: DataType = DataType.VALUE):
    """Read a block entry through the raw pybind provider without Python float()."""
    return provider._cpp_obj(str(block), _cpp_lhaid(code), dtype.value)


def safe_provider_value(pid: ParamId, dtype: DataType = DataType.VALUE) -> float:
    if pid.type is None:
        raise ValueError("ParamId has no ParameterType")
    if not parameter_type_is_allowed(pid.type):
        raise ValueError(
            f"ParameterType.{pid.type.name} is not available for the active model "
            f"{active_model_name() or 'unknown'}."
        )
    provider = ParameterProvider(pid.type)
    return scalar_to_float(_raw_provider_value_by_pid(provider, pid, dtype), "real")


def save_uploaded_lha(contents: str, filename: str) -> str:
    if not contents:
        raise ValueError("No uploaded content")
    header, payload = contents.split(",", 1)
    safe_name = Path(filename or "uploaded.lha").name
    target = UPLOAD_DIR / safe_name
    target.write_bytes(base64.b64decode(payload))
    return str(target)


# ---------- Model-aware ParameterType helpers ----------

_MODEL_DEPENDENT_PARAMETER_TYPES = {"BSM"}


def active_model_name() -> str | None:
    """Return the active Hyperiso model name, if the runtime is initialized."""
    if RUNTIME.config is None:
        return None
    try:
        return RUNTIME.config.model.name
    except Exception:
        return None


def parameter_type_is_allowed(param_type_name: str | ParameterType) -> bool:
    """Return whether a ParameterType is safe for the active C++ model.

    The C++ parameter registry does not necessarily define every
    ``ParameterType`` enum value for every model. In particular,
    ``ParameterType.BSM`` is not available when the active model is ``SM``.
    Calling ``BlockLogger``/``ParameterProvider``/``ParameterSetter`` on such a
    namespace can terminate the process before Python can catch an exception, so
    every service entry point must filter it before crossing the pybind layer.
    """
    name = param_type_name.name if isinstance(param_type_name, ParameterType) else str(param_type_name)
    model = active_model_name()

    if model == "SM" and name in _MODEL_DEPENDENT_PARAMETER_TYPES:
        return False

    return True


def allowed_parameter_types() -> list[ParameterType]:
    """Return ParameterType values safe to expose for the active model."""
    return [pt for pt in ParameterType if parameter_type_is_allowed(pt)]


def allowed_parameter_type_names() -> list[str]:
    """Return safe ParameterType names for the active model."""
    return [pt.name for pt in allowed_parameter_types()]


def parameter_type_options() -> list[dict[str, str]]:
    """Return Dash dropdown options for model-safe ParameterType values."""
    return [{"label": pt.name, "value": pt.name} for pt in allowed_parameter_types()]


def default_parameter_type_name(preferred: str = "SM") -> str:
    """Return a safe default ParameterType name for dropdowns."""
    names = allowed_parameter_type_names()
    if preferred in names:
        return preferred
    return names[0] if names else preferred


def resolve_parameter_type(param_type_name: str, *, strict: bool = True) -> ParameterType | None:
    """Convert a dropdown value to a model-safe ``ParameterType``.

    Args:
        param_type_name: ParameterType enum name coming from the UI.
        strict: If ``True``, raise a clear Python error when the type is not
            available. If ``False``, return ``None`` instead.

    Returns:
        The resolved ``ParameterType`` or ``None`` when ``strict`` is false and
        the namespace is not available.

    Raises:
        ValueError: If the type is not available for the active model.
    """
    pt = enum_by_name(ParameterType, param_type_name)
    if parameter_type_is_allowed(pt):
        return pt

    message = (
        f"ParameterType.{pt.name} is not available for the active model "
        f"{active_model_name() or 'unknown'}. Switch to a BSM/MARTY model or "
        "choose another parameter namespace."
    )
    if strict:
        raise ValueError(message)
    return None


# ---------- Core ----------

def make_hyperiso_config(flag_names: Sequence[str], model_name: str, marty_name: str | None, marty_path: str | None) -> HyperisoConfig:
    selected = set(flag_names or [])
    flags = {
        ExternalFlag.IS_LHA_SPECTRUM: "IS_LHA_SPECTRUM" in selected,
        ExternalFlag.HAS_WILSON_INPUT: "HAS_WILSON_INPUT" in selected,
        # Explicitly keep the observable-input flag out of the GUI for now.
        ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
        ExternalFlag.HYP_AS_SM_MARTY: "HYP_AS_SM_MARTY" in selected,
    }
    model = enum_by_name(Model, model_name)
    return HyperisoConfig(
        flags=flags,
        model=model,
        mty_model_name=(marty_name or None) if model is Model.MARTY else None,
        mty_model_path=Path(marty_path).expanduser() if model is Model.MARTY and marty_path else None,
    )


def init_or_switch_hyperiso(lha_path: str, flag_names: Sequence[str], model_name: str, marty_name: str | None, marty_path: str | None) -> dict:
    if not lha_path:
        raise ValueError("Missing LHA path")
    cfg = make_hyperiso_config(flag_names, model_name, marty_name, marty_path)
    with RUNTIME.lock:
        if RUNTIME.hyp is None:
            RUNTIME.hyp = HyperisoMaster()
        if RUNTIME.initialized:
            RUNTIME.hyp.switch_lha(str(lha_path), cfg)
            action = "switch_lha"
        else:
            RUNTIME.hyp.init(str(lha_path), cfg)
            action = "init"
        RUNTIME.config = cfg
        RUNTIME.lha_path = str(lha_path)
        RUNTIME.initialized = True
        RUNTIME.block_logger = BlockLogger()
        # Cached interfaces are tied to the old parameter store. Rebuild after switch.
        RUNTIME.wilson = None
        RUNTIME.observable = None
        RUNTIME.stat = None
    return {
        "action": action,
        "model": cfg.model.name,
        "lha_path": str(lha_path),
        "flags": {k.name: bool(v) for k, v in cfg.flags.items()},
    }


def runtime_summary() -> dict:
    return {
        "initialized": RUNTIME.initialized,
        "model": RUNTIME.config.model.name if RUNTIME.config else "—",
        "lha_path": RUNTIME.lha_path or "—",
        "wilson_built": RUNTIME.wilson is not None,
        "observable_built": RUNTIME.observable is not None,
        "stat_built": RUNTIME.stat is not None,
    }


def block_inventory_rows() -> list[dict]:
    require_initialized()
    logger = RUNTIME.block_logger or BlockLogger()
    rows = []
    for pt in allowed_parameter_types():
        try:
            blocks = sorted(str(x) for x in logger.get_all_blocks(pt))
        except Exception:
            blocks = []
        rows.append({"parameter_type": pt.name, "n_blocks": len(blocks), "blocks": ", ".join(blocks[:15]) + ("…" if len(blocks) > 15 else "")})
    return rows


def block_name_options(param_type_name: str) -> list[dict]:
    require_initialized()
    pt = resolve_parameter_type(param_type_name, strict=False)
    if pt is None:
        return []
    blocks = sorted(str(x) for x in (RUNTIME.block_logger or BlockLogger()).get_all_blocks(pt))
    return [{"label": b, "value": b} for b in blocks]


def block_table_rows(param_type_name: str, block_name: str) -> list[dict]:
    require_initialized()
    if not block_name:
        return []
    pt = resolve_parameter_type(param_type_name, strict=False)
    if pt is None:
        return []
    logger = RUNTIME.block_logger or BlockLogger()
    values = logger.get_block(pt, block_name)
    provider = ParameterProvider(pt)
    rows = []
    for key, value in values.items():
        code = key_to_code(key)
        row = {"code": code_to_display(code), "value": scalar_to_float(value)}
        for dtype, col in [
            (DataType.STD_STAT, "stat_std"),
            (DataType.STD_SYST, "syst_std"),
            (DataType.STD_COMBINED, "combined_std"),
        ]:
            try:
                row[col] = scalar_to_float(_raw_provider_value_by_block(provider, block_name, code, dtype), "real")
            except Exception:
                row[col] = None
        # Metadata such as scale/bin can be obtained through get_parameter when available.
        try:
            pid = ParamId(pt, block_name, code)
            p = provider.get_parameter(pid)
            cpp = getattr(p, "_cpp_obj", None)
            row["scale"] = getattr(cpp, "scale", None) if cpp is not None else None
            row["bin"] = getattr(cpp, "bin", None) if cpp is not None else None
        except Exception:
            row["scale"] = None
            row["bin"] = None
        rows.append(row)
    return rows


# ---------- Wilson ----------

def build_wilson(groups: Sequence[str], matching_scale: float, hadronic_scale: float, order_name: str, add: bool = False) -> dict:
    require_initialized()
    if not groups:
        raise ValueError("Select at least one Wilson group")
    cfg = WilsonBuildConfig(
        groups={enum_by_name(WGroup, g) for g in groups},
        matching_scale=float(matching_scale),
        hadronic_scale=float(hadronic_scale),
        order=enum_by_name(QCDOrder, order_name),
    )
    with RUNTIME.lock:
        if RUNTIME.wilson is None or not add:
            RUNTIME.wilson = WilsonInterface()
            RUNTIME.wilson.build(cfg)
            action = "build"
        else:
            RUNTIME.wilson.add_wilson_group(cfg)
            action = "add_wilson_group"
        RUNTIME.last_wilson_config = cfg
    return {"action": action, "groups": list(groups), "order": order_name, "matching_scale": matching_scale, "hadronic_scale": hadronic_scale}


def query_wilson(method: str, group: str, coeff: str, order: str, contribution: str, basis: str = "STANDARD", component: str = "real") -> dict:
    require_initialized()
    if RUNTIME.wilson is None:
        raise RuntimeError("WilsonInterface is not built yet")
    req = WilsonRequest(
        group=enum_by_name(WGroup, group),
        coefficient=enum_by_name(WCoeff, coeff),
        order=enum_by_name(QCDOrder, order),
        contribution=enum_by_name(ContributionType, contribution),
        wilson_basis=enum_by_name(WilsonBasis, basis),
    )
    wi = RUNTIME.wilson
    raw_wi = getattr(wi, "_cpp_obj", None)
    if raw_wi is not None:
        if method == "M":
            value = raw_wi.get_M(req.group.value, req.coefficient.value, req.order.value, req.contribution.value)
        elif method == "FM":
            value = raw_wi.get_FM(req.group.value, req.coefficient.value, req.order.value, req.contribution.value)
        elif method == "R":
            value = raw_wi.get_R(req.group.value, req.coefficient.value, req.order.value, req.contribution.value, req.wilson_basis.value)
        elif method == "FR":
            value = raw_wi.get_FR(req.group.value, req.coefficient.value, req.order.value, req.contribution.value, req.wilson_basis.value)
        else:
            raise ValueError(f"Unknown Wilson method {method}")
    elif method == "M":
        value = wi.get_M(req)
    elif method == "FM":
        value = wi.get_FM(req)
    elif method == "R":
        value = wi.get_R(req)
    elif method == "FR":
        value = wi.get_FR(req)
    else:
        raise ValueError(f"Unknown Wilson method {method}")
    return {
        "method": method,
        "group": group,
        "coefficient": coeff,
        "order": order,
        "contribution": contribution,
        "basis": basis,
        "component": component,
        "value": scalar_to_float(value, component),
        "scalar": scalar_to_display(value),
    }


def evaluate_wilson_float(method: str, group: str, coeff: str, order: str, contribution: str, basis: str, component: str = "real") -> float:
    return scalar_to_float(query_wilson(method, group, coeff, order, contribution, basis, component)["value"], "real")


def _mutate_param_temporarily(pid: ParamId, value: float, restore_values: dict[str, float], setter: ParameterSetter) -> None:
    label = param_to_label(pid)
    if label not in restore_values:
        restore_values[label] = safe_provider_value(pid)
    setter.mutate(pid, value)


def _restore_params(pairs: Sequence[tuple[ParamId, str]], restore_values: Mapping[str, float], setter: ParameterSetter) -> None:
    for pid, label in pairs:
        if label in restore_values:
            setter.mutate(pid, restore_values[label])


def wilson_scan_1d(method: str, group: str, coeff: str, order: str, contribution: str, basis: str, component: str, ptype: str, block: str, code: Any, x_min: float, x_max: float, n_points: int) -> tuple[list[float], list[float]]:
    require_initialized()
    pid = make_param_id(ptype, block, code)
    xs = linspace(x_min, x_max, n_points)
    ys: list[float] = []
    setter = ParameterSetter()
    restore: dict[str, float] = {}
    label = param_to_label(pid)
    with RUNTIME.lock:
        try:
            for x in xs:
                _mutate_param_temporarily(pid, float(x), restore, setter)
                ys.append(evaluate_wilson_float(method, group, coeff, order, contribution, basis, component))
        finally:
            _restore_params([(pid, label)], restore, setter)
    return xs, ys


def wilson_scan_2d(method: str, group: str, coeff: str, order: str, contribution: str, basis: str, component: str, p1: tuple, p2: tuple, x_min: float, x_max: float, nx: int, y_min: float, y_max: float, ny: int) -> tuple[list[float], list[float], list[list[float]]]:
    require_initialized()
    pid1 = make_param_id(*p1)
    pid2 = make_param_id(*p2)
    xs = linspace(x_min, x_max, nx)
    ys = linspace(y_min, y_max, ny)
    z: list[list[float]] = []
    setter = ParameterSetter()
    restore: dict[str, float] = {}
    labels = [(pid1, param_to_label(pid1)), (pid2, param_to_label(pid2))]
    with RUNTIME.lock:
        try:
            for y in ys:
                row = []
                _mutate_param_temporarily(pid2, float(y), restore, setter)
                for x in xs:
                    _mutate_param_temporarily(pid1, float(x), restore, setter)
                    row.append(evaluate_wilson_float(method, group, coeff, order, contribution, basis, component))
                z.append(row)
        finally:
            _restore_params(labels, restore, setter)
    return xs, ys, z


# ---------- Observables ----------

def observable_names_from_selection(mode: str, obs_names: Sequence[str] | None, decay_names: Sequence[str] | None) -> list[str]:
    names: list[str] = []
    if mode == "decay":
        mapper = DecayMapper()
        for dn in decay_names or []:
            decay = enum_by_name(Decays, dn)
            names.extend(o.name for o in mapper.get_observables(decay))
    else:
        names.extend(str(x) for x in (obs_names or []))
    # Preserve order, drop duplicates.
    return list(dict.fromkeys(names))


def make_binned_id(obs_name: str, low: float | None, high: float | None) -> BinnedObservableId:
    obs = enum_by_name(Observables, obs_name)
    oid = ObservableMapper.to_id(obs)
    if low is None or high is None:
        return BinnedObservableId(oid)
    return BinnedObservableId(oid, (float(low), float(high)))


def build_observable_interface(mode: str, obs_names: Sequence[str] | None, decay_names: Sequence[str] | None, order_name: str, add_dependencies: bool, bin_low: float | None = None, bin_high: float | None = None, bins: Sequence[tuple[float, float]] | None = None) -> tuple[ObservableInterface, list[dict]]:
    require_initialized()
    oi = ObservableInterface()
    qcd_order = enum_by_name(QCDOrder, order_name)
    selected = observable_names_from_selection(mode, obs_names, decay_names)
    if not selected:
        raise ValueError("Select at least one observable or decay")
    rows: list[dict] = []
    for obs_name in selected:
        if bins:
            for low, high in bins:
                bid = make_binned_id(obs_name, low, high)
                oi.add_binned_observable(bid, qcd_order, bool(add_dependencies))
                rows.append({"observable": obs_name, "bin_low": low, "bin_high": high, "order": order_name})
        elif bin_low is not None and bin_high is not None:
            bid = make_binned_id(obs_name, bin_low, bin_high)
            oi.add_binned_observable(bid, qcd_order, bool(add_dependencies))
            rows.append({"observable": obs_name, "bin_low": bin_low, "bin_high": bin_high, "order": order_name})
        else:
            oi.add_observable(enum_by_name(Observables, obs_name), qcd_order, bool(add_dependencies))
            rows.append({"observable": obs_name, "bin_low": None, "bin_high": None, "order": order_name})
    return oi, rows


def build_observables(mode: str, obs_names: Sequence[str] | None, decay_names: Sequence[str] | None, order_name: str, add_dependencies: bool, use_bin: bool, bin_low: Any, bin_high: Any, smooth: bool = False, smooth_min: Any = None, smooth_max: Any = None, smooth_step: Any = None) -> list[dict]:
    bins = None
    low = high = None
    if smooth:
        bins = bins_from_step(smooth_min, smooth_max, smooth_step)
    elif use_bin:
        low = float(bin_low)
        high = float(bin_high)
    with RUNTIME.lock:
        RUNTIME.observable, rows = build_observable_interface(mode, obs_names, decay_names, order_name, add_dependencies, low, high, bins)
        RUNTIME.last_observable_selection = rows
    return rows


def _observable_rows_from_cpp_compute_all(oi: ObservableInterface) -> list[dict]:
    """Compute all observables while converting raw C++ scalar_t values safely."""
    raw = getattr(oi, "_to_cpp", lambda: None)()
    if raw is None:
        rows: list[dict] = []
        for oid, values in oi.compute_all().items():
            for value in values:
                low = high = None
                if value.bin is not None:
                    low, high = value.bin
                rows.append({
                    "observable_id": str(value.id),
                    "bin_low": low,
                    "bin_high": high,
                    "value": scalar_to_float(value.value, "real"),
                })
        return rows

    result = raw.compute_all()
    rows: list[dict] = []
    for _oid, values in result.items():
        for value in values:
            low = high = None
            b = getattr(value, "bin", None)
            if b is not None:
                low, high = float(b[0]), float(b[1])
            rows.append({
                "observable_id": str(getattr(value, "id", _oid)),
                "bin_low": low,
                "bin_high": high,
                "value": scalar_to_float(getattr(value, "value"), "real"),
            })
    return rows


def compute_current_observables() -> list[dict]:
    require_initialized()
    if RUNTIME.observable is None:
        raise RuntimeError("ObservableInterface is not configured yet")
    return _observable_rows_from_cpp_compute_all(RUNTIME.observable)


def evaluate_observable_float(obs_name: str, order_name: str, add_deps: bool, bin_low: float | None, bin_high: float | None) -> float:
    oi, _ = build_observable_interface("observable", [obs_name], None, order_name, add_deps, bin_low, bin_high)
    rows = _observable_rows_from_cpp_compute_all(oi)
    if not rows:
        raise RuntimeError("Observable computation returned no value")
    return scalar_to_float(rows[0]["value"], "real")


def observable_scan_1d(obs_name: str, order_name: str, add_deps: bool, bin_low: Any, bin_high: Any, ptype: str, block: str, code: Any, x_min: float, x_max: float, n_points: int) -> tuple[list[float], list[float]]:
    require_initialized()
    pid = make_param_id(ptype, block, code)
    low = as_float(bin_low, None)
    high = as_float(bin_high, None)
    xs = linspace(x_min, x_max, n_points)
    ys = []
    setter = ParameterSetter()
    restore: dict[str, float] = {}
    label = param_to_label(pid)
    with RUNTIME.lock:
        try:
            for x in xs:
                _mutate_param_temporarily(pid, float(x), restore, setter)
                ys.append(evaluate_observable_float(obs_name, order_name, add_deps, low, high))
        finally:
            _restore_params([(pid, label)], restore, setter)
    return xs, ys


def observable_scan_2d(obs_name: str, order_name: str, add_deps: bool, bin_low: Any, bin_high: Any, p1: tuple, p2: tuple, x_min: float, x_max: float, nx: int, y_min: float, y_max: float, ny: int) -> tuple[list[float], list[float], list[list[float]]]:
    require_initialized()
    pid1 = make_param_id(*p1)
    pid2 = make_param_id(*p2)
    low = as_float(bin_low, None)
    high = as_float(bin_high, None)
    xs = linspace(x_min, x_max, nx)
    ys = linspace(y_min, y_max, ny)
    z = []
    setter = ParameterSetter()
    restore: dict[str, float] = {}
    labels = [(pid1, param_to_label(pid1)), (pid2, param_to_label(pid2))]
    with RUNTIME.lock:
        try:
            for y in ys:
                row = []
                _mutate_param_temporarily(pid2, float(y), restore, setter)
                for x in xs:
                    _mutate_param_temporarily(pid1, float(x), restore, setter)
                    row.append(evaluate_observable_float(obs_name, order_name, add_deps, low, high))
                z.append(row)
        finally:
            _restore_params(labels, restore, setter)
    return xs, ys, z


def smooth_observable_bins(obs_name: str, order_name: str, add_deps: bool, low: Any, high: Any, step: Any) -> list[dict]:
    rows = []
    for blo, bhi in bins_from_step(low, high, step):
        val = evaluate_observable_float(obs_name, order_name, add_deps, blo, bhi)
        rows.append({"bin_low": blo, "bin_high": bhi, "bin_center": 0.5 * (blo + bhi), "value": val})
    return rows


# ---------- Statistics ----------

def make_stat_config(mc_draws: Any, skew_threshold: Any, ridge_rel: Any, ridge_abs: Any, nuisance_pruning: bool, nuisance_contexts: Any, nuisance_seed: Any) -> StatisticConfig:
    cfg = StatisticConfig()
    cfg.likelihood_mode = StatisticLikelihoodMode.CHI2_MC_COVARIANCE
    cfg.MC_draws = int(mc_draws or 100)
    cfg.skew_abs_threshold = float(skew_threshold or 0.2)
    cfg.chi2_covariance_ridge_rel = float(ridge_rel or 1e-8)
    cfg.chi2_covariance_ridge_abs = float(ridge_abs or 1e-12)
    cfg.nuisance_sensitivity_pruning = bool(nuisance_pruning)
    cfg.nuisance_sensitivity_contexts = int(nuisance_contexts or 2)
    cfg.nuisance_sensitivity_seed = int(nuisance_seed or 12345)
    # Dynamic attributes consumed by the existing StatisticInterface wrapper.
    cfg.p_specs = []
    cfg.selected_experiments = None
    return cfg


def build_stat_interface(mode: str, obs_names: Sequence[str] | None, decay_names: Sequence[str] | None, order_name: str, add_deps: bool, bin_strategy: str, bin_low: Any, bin_high: Any, smooth_min: Any, smooth_max: Any, smooth_step: Any, experiments: str | None, mc_draws: Any, skew_threshold: Any, ridge_rel: Any, ridge_abs: Any, nuisance_pruning: bool, nuisance_contexts: Any, nuisance_seed: Any, p_specs: Sequence[ParamId] | None = None) -> StatisticInterface:
    bins = None
    low = high = None
    if bin_strategy == "smooth":
        bins = bins_from_step(smooth_min, smooth_max, smooth_step)
    elif bin_strategy == "single":
        low = as_float(bin_low, None)
        high = as_float(bin_high, None)
    obs_int, _ = build_observable_interface(mode, obs_names, decay_names, order_name, add_deps, low, high, bins)
    cfg = make_stat_config(mc_draws, skew_threshold, ridge_rel, ridge_abs, nuisance_pruning, nuisance_contexts, nuisance_seed)
    cfg.p_specs = list(p_specs or [])
    exps = split_csv(experiments)
    cfg.selected_experiments = exps or None
    stat_int = StatisticInterface(cfg, observable_interface=obs_int)
    return stat_int


def compute_uncertainty_rows(**kwargs) -> list[dict]:
    require_initialized()
    with RUNTIME.lock:
        RUNTIME.stat = build_stat_interface(**kwargs)
        summaries = RUNTIME.stat.compute_uncertainties()
    rows = []
    for bid, summary in summaries.items():
        low, high = bid.p
        row = {
            "observable": str(bid.s),
            "bin_low": None if float(low) == 0.0 and float(high) == 0.0 else float(low),
            "bin_high": None if float(low) == 0.0 and float(high) == 0.0 else float(high),
            "bin_center": None if float(low) == 0.0 and float(high) == 0.0 else 0.5 * (float(low) + float(high)),
            "central": scalar_to_float(summary.mu if summary.symmetric else summary.mode, "real"),
            "mu": scalar_to_float(summary.mu, "real"),
            "mode": scalar_to_float(summary.mode, "real"),
            "sigma": scalar_to_float(summary.sigma, "real"),
            "sigma_plus": scalar_to_float(summary.sigma_p, "real"),
            "sigma_minus": scalar_to_float(summary.sigma_m, "real"),
            "skew": scalar_to_float(summary.skew, "real"),
            "symmetric": bool(summary.symmetric),
        }
        rows.append(row)
    return rows


def p_specs_from_rows(rows: Sequence[dict]) -> list[ParamId]:
    out: list[ParamId] = []
    for row in rows or []:
        if len(out) >= 10:
            break
        if not row:
            continue
        ptype, block, code = row.get("type"), row.get("block"), row.get("code")
        if ptype and block and code not in (None, ""):
            out.append(make_param_id(ptype, block, code))
    return out


def fit_result_rows(fit) -> tuple[list[dict], list[dict], list[dict]]:
    p_rows = []
    for p, v in fit.p_hat.items():
        label = param_to_label(p)
        p_rows.append({"parameter": label, "best_fit": scalar_to_float(v, "real"), "std": scalar_to_float(fit.p_hat_std.get(p, float("nan")), "real") if fit.p_hat_std else None})
    eta_rows = [{"nuisance": param_to_label(p), "value": scalar_to_float(v, "real")} for p, v in fit.eta_hat.items()]
    corr_rows = []
    for p1, inner in fit.p_correlations.items():
        for p2, corr in inner.items():
            corr_rows.append({"x": param_to_label(p1), "y": param_to_label(p2), "corr": scalar_to_float(corr, "real")})
    return p_rows, eta_rows, corr_rows


def run_fit_and_scan(stat_kwargs: dict, p_spec_rows: Sequence[dict], do_contour: bool, x_half_width: Any, y_half_width: Any, nx: Any, ny: Any) -> dict:
    require_initialized()
    p_specs = p_specs_from_rows(p_spec_rows)
    if not p_specs:
        raise ValueError("Add at least one p_spec row")
    stat_kwargs = dict(stat_kwargs)
    stat_kwargs["p_specs"] = p_specs
    with RUNTIME.lock:
        RUNTIME.stat = build_stat_interface(**stat_kwargs)
        fit = RUNTIME.stat.compute_MLE(p_specs)
        p_rows, eta_rows, corr_rows = fit_result_rows(fit)
        points = []
        if do_contour and len(p_specs) == 2:
            RUNTIME.stat.prepare_likelihood_for_scan(p_specs)
            RUNTIME.stat.set_manual_scan_point(fit.p_hat, fit.eta_hat)
            grid = RUNTIME.stat.scan_likelihood_around_current_point(
                p_specs[0],
                p_specs[1],
                float(x_half_width or 1.0),
                float(y_half_width or 1.0),
                int(nx or 25),
                int(ny or 25),
            )
            points = [{"x": scalar_to_float(p.x, "real"), "y": scalar_to_float(p.y, "real"), "nll": scalar_to_float(p.nll, "real"), "delta_nll": scalar_to_float(p.delta_nll, "real")} for p in grid.points]
    return {
        "fit_ok": bool(fit.fit_ok),
        "ell_hat": scalar_to_float(fit.ell_hat, "real"),
        "p_rows": p_rows,
        "eta_rows": eta_rows,
        "corr_rows": corr_rows,
        "scan_points": points,
    }
