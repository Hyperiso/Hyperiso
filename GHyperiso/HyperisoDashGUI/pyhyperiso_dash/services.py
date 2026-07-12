from __future__ import annotations

import base64
import threading
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from pyhyperiso_dash import latex as lx
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

try:
    from pyhyperiso.core.BusinessLogic.DecayConfig import (
        DecayConfig,
        BDlnuConfig, BDstarlnuConfig, BKllConfig, BKstarllConfig,
        BKstarGammaConfig, BsPhiConfig, BXsllConfig, KllDecayConfig, LbLllConfig,
        BDlnuBCharge, BDstarlnuBCharge, BKllBCharge, BKllLepton,
        BKstarllPowerCorrectionsImpl, BKstarllBCharge, BKstarllLepton,
        BKstarGammaBCharge, BsPhiLepton, BXsllLepton, LbLllLepton,
        BFFType, BPFFSource, BVFFSource, LbLFFSource,
    )
except Exception:  # optional on older pyhyperiso builds
    DecayConfig = None
    BDlnuConfig = BDstarlnuConfig = BKllConfig = BKstarllConfig = None
    BKstarGammaConfig = BsPhiConfig = BXsllConfig = KllDecayConfig = LbLllConfig = None
    BDlnuBCharge = BDstarlnuBCharge = BKllBCharge = BKllLepton = None
    BKstarllPowerCorrectionsImpl = BKstarllBCharge = BKstarllLepton = None
    BKstarGammaBCharge = BsPhiLepton = BXsllLepton = LbLllLepton = None
    BFFType = BPFFSource = BVFFSource = LbLFFSource = None
from pyhyperiso.core.Common.BinnedObservableId import BinnedObservableId
from pyhyperiso.core.Common.Configs import AlphasConfig, MassConfig, WilsonBuildConfig, WilsonRequest
from pyhyperiso.core.Common.GeneralEnum import (
    ContributionType,
    DataType,
    Decays,
    MassType,
    Model,
    Observables,
    ParameterType,
    QCDOrder,
    ScaleType,
    WilsonBasis,
)
from pyhyperiso.core.Common.Mapper import DecayMapper, GroupMapper, ObservableMapper, WCoefMapper
from pyhyperiso.core.Common.LhaID import LhaID
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Core.BlockProvider import BlockLogger
try:
    from pyhyperiso.core.Core.DependantBlockInfoProvider import DependantBlockInfoProvider
except Exception:  # optional on older pyhyperiso builds
    DependantBlockInfoProvider = None
try:
    from pyhyperiso.core.Core.DependencyPruner import DependencyPruner
except Exception:  # optional on older pyhyperiso builds
    DependencyPruner = None
from pyhyperiso.core.Core.HyperisoConfig import ExternalFlag, HyperisoConfig
from pyhyperiso.core.Core.HyperisoMaster import HyperisoMaster
from pyhyperiso.core.Core.ParamaterProvider import ParameterProvider
from pyhyperiso.core.Core.ParameterSetter import ParameterSetter
from pyhyperiso.core.Core.QCDProvider import QCDProvider
from pyhyperiso.core.Core.QEDProvider import QEDProvider
from pyhyperiso.core.Statistic.Copula import CopulaKind
from pyhyperiso.core.Statistic.StatisticConfig import (
    StatisticConfig, StatisticLikelihoodMode, StatisticProgressMonitor,
)
from pyhyperiso.core.Statistic.StatisticInterface import (
    StatisticInterface,
    ContourOptions,
    ContourAlgorithm,
    ProfilingMethod,
    ProfilerMode,
)
from pyhyperiso.core.PhysicalModel.WilsonInterface import WilsonInterface


# The current Python ParamId wrapper is a dataclass with value equality but no
# __hash__.  The statistic wrapper converts C++ maps into dict[ParamId, ...],
# so fits fail with ``TypeError: unhashable type: 'ParamId'`` unless the GUI
# provides a stable value hash.  Keep the hash local and deterministic.
def _param_id_hash(pid: ParamId) -> int:
    return hash((pid.type.name if pid.type else None, str(pid.block), pid.code.to_string()))

try:
    if getattr(ParamId, "__hash__", None) is None:
        ParamId.__hash__ = _param_id_hash  # type: ignore[method-assign]
except Exception:
    pass


@dataclass
class RuntimeState:
    """Process-local state holding live wrappers around the C++ singleton."""

    hyp: HyperisoMaster | None = None
    config: HyperisoConfig | None = None
    lha_path: str | None = None
    initialized: bool = False
    block_logger: BlockLogger | None = None
    dependency_info: Any | None = None
    dependency_pruner: Any | None = None
    wilson: WilsonInterface | None = None
    observable: ObservableInterface | None = None
    stat: StatisticInterface | None = None
    lock: threading.RLock = field(default_factory=threading.RLock)
    last_wilson_config: WilsonBuildConfig | None = None
    last_observable_selection: list[dict] = field(default_factory=list)
    observable_signature: tuple | None = None
    observable_registered_keys: set[tuple] = field(default_factory=set)
    observable_registered_decay_bins: set[tuple] = field(default_factory=set)
    scan_observable: ObservableInterface | None = None
    scan_observable_signature: tuple | None = None
    scan_observable_registered_keys: set[tuple] = field(default_factory=set)
    scan_observable_registered_decay_bins: set[tuple] = field(default_factory=set)
    stat_observable: ObservableInterface | None = None
    stat_observable_registered_keys: set[tuple] = field(default_factory=set)
    stat_observable_registered_decay_bins: set[tuple] = field(default_factory=set)
    last_stat_observable_selection: list[dict] = field(default_factory=list)
    wilson_status: dict = field(default_factory=dict)
    decay_configs: dict[str, dict] = field(default_factory=dict)


RUNTIME = RuntimeState()
DATA_DIR = Path(__file__).resolve().parent / "data"
UPLOAD_DIR = DATA_DIR / "uploaded_lha"
UPLOAD_DIR.mkdir(parents=True, exist_ok=True)


@dataclass
class StatisticJob:
    """One same-process statistic task with a C++ progress monitor."""

    job_id: str
    kind: str
    monitor: StatisticProgressMonitor
    created_at: float = field(default_factory=time.monotonic)
    done: bool = False
    result: Any = None
    error: str | None = None
    thread: threading.Thread | None = None
    lock: threading.RLock = field(default_factory=threading.RLock)


_STATISTIC_JOBS: dict[str, StatisticJob] = {}
_STATISTIC_JOBS_LOCK = threading.RLock()


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


def param_to_latex_label(pid: ParamId) -> str:
    ptype = pid.type.name if pid.type else None
    return lx.parameter_table_label(str(pid.block), pid.code.to_string(), param_to_label(pid), ptype)


def parameter_display_label(block: str | None, code: Any, raw: str | None = None, param_type_name: str | None = None) -> str:
    return lx.parameter_label(block, code, raw, param_type_name)


def p_spec_row_key(row: Mapping[str, Any]) -> str:
    """Return a stable GUI key for one fit-parameter row."""
    return "|".join((str(row.get("type", "")), str(row.get("block", "")), str(row.get("code", ""))))


def suggested_parameter_bounds(param_type_name: str, block: str, code: Any) -> tuple[float, float, float]:
    """Return current value and editable default contour bounds for a parameter.

    The preferred width is four combined standard deviations.  Parameters with
    no stored uncertainty receive a scale-aware fallback so the table remains
    useful for BSM coefficients initialized at zero.
    """
    require_initialized()
    ptype = enum_by_name(ParameterType, param_type_name)
    pid = make_param_id(param_type_name, block, code)
    provider = ParameterProvider(ptype)
    center = scalar_to_float(_raw_provider_value_by_pid(provider, pid, DataType.VALUE), "real")
    try:
        sigma = abs(scalar_to_float(_raw_provider_value_by_pid(provider, pid, DataType.STD_COMBINED), "real"))
    except Exception:
        sigma = 0.0
    if not (sigma > 0.0):
        sigma = 0.0
    fallback = 0.25 * max(abs(center), 1.0)
    half_width = max(4.0 * sigma, fallback, 1e-6)
    return center, center - half_width, center + half_width


def p_spec_axis_options(rows: Sequence[dict] | None) -> list[dict[str, str]]:
    """Build contour-axis dropdown options from current fit-parameter rows."""
    out: list[dict[str, str]] = []
    for row in rows or []:
        key = p_spec_row_key(row)
        if not key.strip("|"):
            continue
        out.append({"label": str(row.get("parameter") or key), "value": key})
    return out


def observable_display_label(obs_name: str | None) -> str:
    return lx.observable_label(obs_name)


def wilson_request_display_label(method: str, coeff: str, order: str, contribution: str) -> str:
    return lx.wilson_request_latex(method, coeff, order, contribution)


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


def stat_parameter_type_options() -> list[dict[str, str]]:
    """Generic fit-parameter namespaces; Wilson coefficients use Wilson Scan."""
    return [
        {"label": pt.name, "value": pt.name}
        for pt in allowed_parameter_types()
        if pt.name != "WILSON"
    ]


def default_parameter_type_name(preferred: str = "SM") -> str:
    """Return a safe default ParameterType name for dropdowns."""
    names = allowed_parameter_type_names()
    if preferred in names:
        return preferred
    return names[0] if names else preferred


def decay_options() -> list[dict]:
    """Return decay dropdown options with LaTeX labels."""
    return [lx.decay_option(dec.name) for dec in Decays]


def observable_options_for_decays(decay_names: Sequence[str] | None = None) -> list[dict]:
    """Return observable dropdown options, optionally filtered by decay."""
    names: list[str] = []
    mapper = DecayMapper()
    if decay_names:
        for decay_name in decay_names:
            try:
                decay = enum_by_name(Decays, decay_name)
                names.extend(obs.name for obs in mapper.get_observables(decay))
            except Exception:
                continue
    else:
        names = [obs.name for obs in Observables]
    return [lx.observable_option(name) for name in dict.fromkeys(names)]


def observable_is_in_decay(obs_name: str | None, decay_names: Sequence[str] | None) -> bool:
    if not obs_name or not decay_names:
        return True
    valid = {opt["value"] for opt in observable_options_for_decays(decay_names)}
    return str(obs_name) in valid


def _enum_choice_options(enum_cls: Any) -> list[dict[str, str]]:
    """Return dropdown options for a Python Enum class."""
    if enum_cls is None:
        return []
    return [{"label": item.name, "value": item.name} for item in enum_cls]


_DECAY_CONFIG_SPECS: dict[str, dict[str, Any]] = {}

def _build_decay_config_specs() -> dict[str, dict[str, Any]]:
    """Build lazy decay-config metadata from the available wrapper classes."""
    if _DECAY_CONFIG_SPECS:
        return _DECAY_CONFIG_SPECS
    if DecayConfig is None:
        return _DECAY_CONFIG_SPECS

    def add(decay: str, cls: Any, fields: list[tuple[str, str, str, Any, str]]):
        if cls is not None:
            _DECAY_CONFIG_SPECS[decay] = {"class": cls, "fields": fields}

    add("B__D_l_nu", BDlnuConfig, [("charge", "B charge", "enum", BDlnuBCharge, "B-meson charge convention")])
    add("B__Dstar_l_nu", BDstarlnuConfig, [("charge", "B charge", "enum", BDstarlnuBCharge, "B-meson charge convention")])
    add("B__K_l_l", BKllConfig, [
        ("ff_src", "Form-factor source", "enum", BPFFSource, "B → K form-factor source"),
        ("ff_type", "Form-factor type", "enum", BFFType, "Full or soft form-factor treatment"),
        ("charge", "B charge", "enum", BKllBCharge, "B-meson charge convention"),
        ("gen", "Lepton", "enum", BKllLepton, "Lepton generation"),
        ("n_threads", "Threads", "int", None, "Requested worker threads"),
    ])
    add("B__Kstar_l_l", BKstarllConfig, [
        ("ff_src", "Form-factor source", "enum", BVFFSource, "B → K* form-factor source"),
        ("ff_type", "Form-factor type", "enum", BFFType, "Full or soft form-factor treatment"),
        ("power_corr_impl", "Power corrections", "enum", BKstarllPowerCorrectionsImpl, "Non-factorisable power-correction prescription"),
        ("charge", "B charge", "enum", BKstarllBCharge, "B-meson charge convention"),
        ("gen", "Lepton", "enum", BKstarllLepton, "Lepton generation"),
        ("n_threads", "Threads", "int", None, "Requested worker threads"),
    ])
    add("B__Kstar_gamma", BKstarGammaConfig, [
        ("ff_src", "Form-factor source", "enum", BVFFSource, "B → K* form-factor source"),
        ("charge", "B charge", "enum", BKstarGammaBCharge, "B-meson charge convention"),
    ])
    add("Bs__phi_l_l", BsPhiConfig, [
        ("ff_src", "Form-factor source", "enum", BVFFSource, "Bs → φ form-factor source"),
        ("ff_type", "Form-factor type", "enum", BFFType, "Full or soft form-factor treatment"),
        ("gen", "Lepton", "enum", BsPhiLepton, "Lepton generation"),
        ("n_threads", "Threads", "int", None, "Requested worker threads"),
    ])
    add("B__Xs_l_l", BXsllConfig, [("gen", "Lepton", "enum", BXsllLepton, "Lepton generation")])
    add("K__l_l", KllDecayConfig, [
        ("N_L_sign", "N_L sign", "int", None, "Long-distance sign convention"),
        ("gen", "Lepton generation", "int", None, "C++ integer lepton-generation convention"),
    ])
    add("Lambda_b__Lambda_l_l", LbLllConfig, [
        ("ff_src", "Form-factor source", "enum", LbLFFSource, "Λb → Λ form-factor source"),
        ("gen", "Lepton", "enum", LbLllLepton, "Lepton generation"),
    ])
    return _DECAY_CONFIG_SPECS


def decay_config_field_specs(decay_name: str | None) -> list[dict[str, Any]]:
    """Return UI metadata for the decay configuration associated with a decay."""
    specs = _build_decay_config_specs()
    entry = specs.get(str(decay_name or ""))
    if not entry:
        return []
    cfg = entry["class"]()
    out: list[dict[str, Any]] = []
    for name, label, kind, enum_cls, help_text in entry["fields"]:
        value = getattr(cfg, name)
        if kind == "enum":
            current = value.name if hasattr(value, "name") else str(value)
            out.append({
                "name": name,
                "label": label,
                "kind": kind,
                "value": current,
                "options": _enum_choice_options(enum_cls),
                "help": help_text,
            })
        else:
            out.append({
                "name": name,
                "label": label,
                "kind": kind,
                "value": int(value),
                "options": [],
                "help": help_text,
            })
    return out


def apply_decay_config(decay_name: str | None, values: Mapping[str, Any]) -> dict[str, Any]:
    """Apply a typed decay configuration to the shared ObservableInterface."""
    require_initialized()
    decay_name = str(decay_name or "")
    specs = _build_decay_config_specs()
    entry = specs.get(decay_name)
    if not entry:
        raise ValueError(f"Decay {decay_name or '—'} has no configurable options exposed by the current binding.")
    kwargs: dict[str, Any] = {}
    for name, _label, kind, enum_cls, _help in entry["fields"]:
        raw = values.get(name)
        if kind == "enum":
            if enum_cls is None:
                continue
            kwargs[name] = enum_cls[str(raw)]
        elif kind == "int":
            kwargs[name] = int(raw)
        else:
            kwargs[name] = raw
    cfg = entry["class"](**kwargs)
    with RUNTIME.lock:
        if RUNTIME.observable is None:
            RUNTIME.observable = ObservableInterface()
        RUNTIME.observable.set_decay_config(enum_by_name(Decays, decay_name), cfg)
        # The statistic manager caches observables/nuisances; invalidate it so
        # the next uncertainty/fit sees the new decay configuration.
        RUNTIME.stat = None
        RUNTIME.decay_configs[decay_name] = {"class": type(cfg).__name__, "values": dict(values)}
    return {"decay": decay_name, "config": type(cfg).__name__, "values": dict(values)}


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
        RUNTIME.dependency_info = DependantBlockInfoProvider() if DependantBlockInfoProvider is not None else None
        RUNTIME.dependency_pruner = DependencyPruner() if DependencyPruner is not None else None
        # Cached interfaces are tied to the old parameter store. Rebuild after switch.
        RUNTIME.wilson = None
        RUNTIME.last_wilson_config = None
        RUNTIME.wilson_status = {}
        # Keep a long-lived ObservableInterface as soon as Hyperiso is initialized.
        # The C++ ObsManager/DecayParent layer is stateful: repeated add_obs/add_bin
        # calls on fresh or duplicated interfaces can leave decay-level bins and
        # Wilson builders in inconsistent states.  Dash therefore adds only missing
        # observables/bins to this persistent interface.
        RUNTIME.observable = ObservableInterface()
        RUNTIME.observable_signature = tuple()
        RUNTIME.observable_registered_keys.clear()
        RUNTIME.observable_registered_decay_bins.clear()
        RUNTIME.last_observable_selection = []
        # Scan and statistic workflows deliberately reuse the same persistent
        # ObservableInterface as the Observable page. Creating a second C++
        # ObservableInterface for the same observable can trigger a second
        # ObsManager::add_obs / DecayParent::enable sequence and lead to
        # std::map::at from stale C++ observable/decay state.
        RUNTIME.scan_observable = None
        RUNTIME.scan_observable_signature = tuple()
        RUNTIME.scan_observable_registered_keys.clear()
        RUNTIME.scan_observable_registered_decay_bins.clear()
        RUNTIME.stat_observable = None
        RUNTIME.stat_observable_registered_keys.clear()
        RUNTIME.stat_observable_registered_decay_bins.clear()
        RUNTIME.last_stat_observable_selection = []
        RUNTIME.stat = None
        RUNTIME.decay_configs.clear()
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
        "wilson_status": dict(RUNTIME.wilson_status or {}),
        "observable_built": bool(RUNTIME.observable_registered_keys),
        "observable_count": len(RUNTIME.observable_registered_keys),
        "stat_observable_count": len(RUNTIME.stat_observable_registered_keys or RUNTIME.observable_registered_keys),
        "stat_built": RUNTIME.stat is not None,
    }


def wilson_status_text() -> str:
    """Return a user-facing status string for the persistent WilsonInterface."""
    if RUNTIME.wilson is None:
        return "WilsonInterface is not built yet."
    info = RUNTIME.wilson_status or {}
    groups = ", ".join(info.get("groups", [])) or "configured"
    action = info.get("action", "built")
    order = info.get("order", "—")
    return f"WilsonInterface {action}. Groups: {groups}. Order: {order}."


def observable_status_text() -> str:
    """Return a user-facing status string for the persistent ObservableInterface."""
    if not RUNTIME.initialized:
        return "Initialize Hyperiso on the Core page first."
    if not RUNTIME.observable_registered_keys:
        return "ObservableInterface is alive but no observable has been added yet."
    return f"ObservableInterface ready with {len(RUNTIME.observable_registered_keys)} observable/bin entries."


def stat_observable_status_text() -> str:
    """Return a status string for the statistic observable selection.

    Statistics intentionally reuse the main persistent ObservableInterface.
    This avoids constructing another C++ ObsManager for the same selected
    observables, which can duplicate add_obs/add_bin calls and provoke
    std::map::at failures in the observable layer.
    """
    if not RUNTIME.initialized:
        return "Initialize Hyperiso on the Core page first."
    n = len(RUNTIME.stat_observable_registered_keys or RUNTIME.observable_registered_keys)
    if not n:
        return "Statistic workflow will reuse the main ObservableInterface; no observable has been added yet."
    return f"Statistic workflow ready with {n} observable/bin entries on the shared ObservableInterface."


def current_observable_rows() -> list[dict]:
    """Rows currently known by the Observable page table.

    Dash rebuilds page layouts when navigating.  The live C++ ObservableInterface
    stays in ``RUNTIME.observable``; these rows repopulate the table from the
    process-local registry instead of showing an empty GUI after page changes.
    """
    return list(RUNTIME.last_observable_selection or [])


def current_stat_observable_rows() -> list[dict]:
    """Rows currently known by the Statistic observable table.

    The statistic workflow reuses the main ObservableInterface.  If the stat page
    has no separate selection yet, show the main observable rows so the GUI
    reflects the true shared C++ state when navigating Observable ↔ Stat.
    """
    if RUNTIME.last_stat_observable_selection:
        return list(RUNTIME.last_stat_observable_selection)
    return [dict(row, registered=True) for row in (RUNTIME.last_observable_selection or [])]


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


def block_code_options(param_type_name: str, block_name: str | None) -> list[dict]:
    """Return code dropdown options for a block using BlockLogger.

    The values are the canonical LHA strings accepted by ``LhaID`` such as
    ``"25"`` or ``"511_1"``.  This lets the UI constrain fit/scan parameters
    to existing entries instead of free-typing invalid ids.
    """
    require_initialized()
    if not block_name:
        return []
    pt = resolve_parameter_type(param_type_name, strict=False)
    if pt is None:
        return []
    try:
        values = (RUNTIME.block_logger or BlockLogger()).get_block(pt, str(block_name))
    except Exception:
        return []
    seen: set[str] = set()
    codes: list[str] = []
    for key in values.keys():
        code = code_to_display(key_to_code(key))
        if code in seen:
            continue
        seen.add(code)
        codes.append(code)

    def _code_sort_key(code: str):
        parts = str(code).replace(",", "_").split("_")
        parsed = []
        for part in parts:
            try:
                parsed.append((0, int(part)))
            except Exception:
                parsed.append((1, part))
        return parsed

    out: list[dict] = []
    for code in sorted(codes, key=_code_sort_key):
        try:
            option = lx.parameter_option(str(block_name), code, pt.name if pt else param_type_name)
        except Exception:
            # Never let a display-label problem make the dropdown look empty.
            raw = f"{block_name}:{code}"
            option = {"label": str(code), "value": str(code), "title": raw, "search": raw}
        option.setdefault("value", str(code))
        option.setdefault("title", f"{block_name}:{code}")
        option.setdefault("search", f"{block_name}:{code}")
        out.append(option)
    return out


def parameter_picker_options(param_type_name: str | None, block_name: str | None = None) -> tuple[list[dict], str | None, list[dict], str | None]:
    """Return block/code options plus safe current values for cascaded pickers."""
    ptype = param_type_name or default_parameter_type_name("SM")
    block_opts = block_name_options(ptype)
    block_values = {opt["value"] for opt in block_opts}
    block_value = block_name if block_name in block_values else (block_opts[0]["value"] if block_opts else None)
    code_opts = block_code_options(ptype, block_value) if block_value else []
    code_value = code_opts[0]["value"] if code_opts else None
    return block_opts, block_value, code_opts, code_value


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
        code_s = code_to_display(code)
        row = {
            "code": code_s,
            "name": lx.parameter_label(block_name, code_s, f"{block_name}:{code_s}", pt.name if pt else param_type_name),
            "dependent_block": is_dependent_block(pt.name if pt else param_type_name, block_name),
            "value": scalar_to_float(value),
        }
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



# ---------- Block dependencies ----------

def _dependency_info_provider():
    require_initialized()
    if DependantBlockInfoProvider is None:
        raise RuntimeError("DependantBlockInfoProvider is not available in this pyhyperiso build.")
    if RUNTIME.dependency_info is None:
        RUNTIME.dependency_info = DependantBlockInfoProvider()
    return RUNTIME.dependency_info


def _dependency_pruner():
    require_initialized()
    if DependencyPruner is None:
        raise RuntimeError("DependencyPruner is not available in this pyhyperiso build.")
    if RUNTIME.dependency_pruner is None:
        RUNTIME.dependency_pruner = DependencyPruner()
    return RUNTIME.dependency_pruner


def is_dependent_block(param_type_name: str, block_name: str | None) -> bool:
    """Safely report whether a block is a dependent block."""
    if not block_name:
        return False
    try:
        pt = resolve_parameter_type(param_type_name, strict=False)
        if pt is None:
            return False
        return bool(_dependency_info_provider().is_dependent_block(pt, str(block_name)))
    except Exception:
        return False


def block_dependency_data(param_type_name: str, block_name: str | None) -> dict:
    """Return dependency rows and graph edges for one selected block.

    The provider exposes direct and transitive upstream/downstream relations. To
    draw the local connected tree, this helper takes the transitive component of
    the selected block and then queries direct source/dependent links inside that
    component.
    """
    require_initialized()
    if not block_name:
        raise ValueError("Select a block first.")
    pt = resolve_parameter_type(param_type_name, strict=True)
    provider = _dependency_info_provider()
    block = str(block_name)

    def _safe_list(fn, name):
        try:
            return sorted({str(x) for x in fn(pt, name)})
        except Exception:
            return []

    direct_sources = _safe_list(provider.get_source_blocks, block)
    direct_dependents = _safe_list(provider.get_dependent_blocks, block)
    all_sources = _safe_list(provider.get_all_source_blocks, block)
    all_dependents = _safe_list(provider.get_all_dependent_blocks, block)
    dependent = bool(provider.is_dependent_block(pt, block))

    nodes = {block, *direct_sources, *direct_dependents, *all_sources, *all_dependents}
    edges: set[tuple[str, str]] = set()
    for node in list(nodes):
        for src in _safe_list(provider.get_source_blocks, node):
            if src in nodes:
                edges.add((src, node))
        for dst in _safe_list(provider.get_dependent_blocks, node):
            if dst in nodes:
                edges.add((node, dst))

    # Ensure the direct one-hop relations are present even if the implementation
    # does not expose the reverse direction for a detached block.
    for src in direct_sources:
        edges.add((src, block))
    for dst in direct_dependents:
        edges.add((block, dst))

    rows: list[dict] = []
    for relation, values in [
        ("direct upstream", direct_sources),
        ("direct downstream", direct_dependents),
        ("all upstream", all_sources),
        ("all downstream", all_dependents),
    ]:
        if values:
            rows.extend({"relation": relation, "block": value} for value in values)
        else:
            rows.append({"relation": relation, "block": "—"})

    return {
        "parameter_type": pt.name,
        "block": block,
        "is_dependent": dependent,
        "direct_sources": direct_sources,
        "direct_dependents": direct_dependents,
        "all_sources": all_sources,
        "all_dependents": all_dependents,
        "nodes": sorted(nodes),
        "edges": sorted(edges),
        "rows": rows,
    }


def prune_dependency(action: str, scope: str, param_type_name: str, block_name: str | None, code: Any | None = None) -> dict:
    """Detach or reattach a block/parameter dependency and return fresh graph data."""
    require_initialized()
    if not block_name:
        raise ValueError("Select a block first.")
    pt = resolve_parameter_type(param_type_name, strict=True)
    pruner = _dependency_pruner()
    action = str(action or "detach")
    scope = str(scope or "block")
    block = str(block_name)

    if scope == "parameter":
        if code in (None, ""):
            raise ValueError("Select a parameter code for parameter-level pruning.")
        raw = getattr(pruner, "_cpp_obj", None)
        if raw is not None:
            # The standalone Python LhaID wrapper may not expose to_cpp(); use
            # the bound C++ object directly so parameter-level pruning works with
            # both old and new pyhyperiso wrappers.
            if action == "reattach":
                raw.reattach_parameter(pt.value, block, _cpp_lhaid(code))
            else:
                raw.detach_parameter(pt.value, block, _cpp_lhaid(code))
        else:
            lhaid = LhaID(str(code))
            if action == "reattach":
                pruner.reattach_parameter(pt, block, lhaid)
            else:
                pruner.detach_parameter(pt, block, lhaid)
    else:
        if action == "reattach":
            pruner.reattach_block(pt, block)
        else:
            pruner.detach_block(pt, block)

    # Some providers cache C++ handles internally; recreate the info provider so
    # the next graph reflects the current dependency state.
    RUNTIME.dependency_info = DependantBlockInfoProvider() if DependantBlockInfoProvider is not None else None
    data = block_dependency_data(pt.name, block)
    data["last_action"] = {"action": action, "scope": scope, "block": block, "code": str(code) if code not in (None, "") else None}
    return data

# ---------- Wilson ----------

# The current core resolves Wilson groups/coefficients through dynamic symbol
# ids.  Keep GUI values as strings instead of coercing them back to the legacy
# Python enums, whose surface may lag behind the C++ enum after a core update.
_WILSON_GROUP_ALIASES: dict[str, tuple[str, ...]] = {
    "B": ("B", "BCoefficients"),
    "BPrime": ("BPrime", "BPrimeCoefficients"),
    "BScalar": ("BScalar", "BScalarCoefficients"),
    "CC_bc": ("CC_bc", "BcChargedCurrentCoefficients"),
    "CC_bu": ("CC_bu", "BuChargedCurrentCoefficients"),
    "CC_cs": ("CC_cs", "DsChargedCurrentCoefficients"),
    "CC_cd": ("CC_cd", "DdChargedCurrentCoefficiens"),
    "CC_su": ("CC_su", "KuChargedCurrentCoefficients"),
    "CC_du": ("CC_du", "PIuChargedCurrentCoefficients"),
    "MESON_MIXING": ("MESON_MIXING", "MesonMixing"),
    "K": ("K",),
}

_WILSON_GROUP_COEFFS: dict[str, list[str]] = {
    "B": [f"C{i}" for i in range(1, 11)],
    "BPrime": [
        *(f"CP{i}" for i in range(1, 11)),
        "CPQ1", "CPQ2",
        "CPQ1_E", "CPQ1_MU", "CPQ1_TA", "CPQ2_E", "CPQ2_MU", "CPQ2_TA",
    ],
    "BScalar": ["CQ1", "CQ2", "CQ1_E", "CQ1_MU", "CQ1_TA", "CQ2_E", "CQ2_MU", "CQ2_TA"],
    "CC_bc": ["C_V1_bc", "C_V2_bc", "C_S1_bc", "C_S2_bc", "C_T_bc"],
    "CC_bu": ["C_V1_bu", "C_V2_bu", "C_S1_bu", "C_S2_bu", "C_T_bu"],
    "CC_cs": ["C_V1_cs", "C_V2_cs", "C_S1_cs", "C_S2_cs", "C_T_cs"],
    "CC_cd": ["C_V1_cd", "C_V2_cd", "C_S1_cd", "C_S2_cd", "C_T_cd"],
    "CC_su": ["C_V1_su", "C_V2_su", "C_S1_su", "C_S2_su", "C_T_su"],
    "CC_du": ["C_V1_du", "C_V2_du", "C_S1_du", "C_S2_du", "C_T_du"],
    "MESON_MIXING": [
        "C_BD_1", "CT_BD_1", "C_BD_2", "CT_BD_2", "C_BD_3", "CT_BD_3", "C_BD_4", "C_BD_5",
        "C_BS_1", "CT_BS_1", "C_BS_2", "CT_BS_2", "C_BS_3", "CT_BS_3", "C_BS_4", "C_BS_5",
        "C_SD_1", "CT_SD_1", "C_SD_2", "CT_SD_2", "C_SD_3", "CT_SD_3", "C_SD_4", "C_SD_5",
        "C_CU_1", "CT_CU_1", "C_CU_2", "CT_CU_2", "C_CU_3", "CT_CU_3", "C_CU_4", "C_CU_5",
    ],
    "K": ["CK9", "CPK9", "CK10", "CPK10", "CKQ1", "CKQ2", "CPKQ1", "CPKQ2", "CK_L"],
}


def _mapper_strings(mapper: Any) -> list[str]:
    """Return canonical mapper strings without depending on enum completeness."""
    raw = mapper.get_str()
    values = raw.values() if isinstance(raw, Mapping) else raw
    return list(dict.fromkeys(str(value) for value in values))


def _wilson_group_key(group_name: str | None) -> str | None:
    value = str(group_name or "")
    folded = value.casefold()
    for key, aliases in _WILSON_GROUP_ALIASES.items():
        if any(folded == alias.casefold() for alias in aliases):
            return key
    return None


def wilson_group_options() -> list[dict[str, str]]:
    """Return all builtin groups advertised by the current core mapper."""
    try:
        canonical_names = _mapper_strings(GroupMapper())
    except Exception:
        canonical_names = [aliases[-1] for aliases in _WILSON_GROUP_ALIASES.values()]

    label_by_canonical = {
        aliases[-1]: key for key, aliases in _WILSON_GROUP_ALIASES.items()
    }
    options = []
    for canonical in canonical_names:
        short = label_by_canonical.get(canonical, canonical)
        label = short if short == canonical else f"{short} — {canonical}"
        options.append({"label": label, "value": canonical})
    return options


def wilson_coeff_options_for_group(group_name: str | None) -> list[dict[str, str]]:
    """Return coefficient options compatible with the selected Wilson group."""
    names = _WILSON_GROUP_COEFFS.get(_wilson_group_key(group_name) or "")
    try:
        valid = set(_mapper_strings(WCoefMapper()))
    except Exception:
        valid = {name for group_names in _WILSON_GROUP_COEFFS.values() for name in group_names}
    if not names:
        names = sorted(valid)
    return [lx.wilson_option(name) for name in names if name in valid]


def build_wilson(groups: Sequence[str], matching_scale: float, hadronic_scale: float, order_name: str, add: bool = False) -> dict:
    require_initialized()
    if not groups:
        raise ValueError("Select at least one Wilson group")
    cfg = WilsonBuildConfig(
        groups={str(g) for g in groups},
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
        existing = set(RUNTIME.wilson_status.get("groups", [])) if RUNTIME.wilson_status else set()
        if add:
            existing.update(str(g) for g in groups)
        else:
            existing = {str(g) for g in groups}
        RUNTIME.wilson_status = {
            "action": action,
            "groups": sorted(existing),
            "order": order_name,
            "matching_scale": float(matching_scale),
            "hadronic_scale": float(hadronic_scale),
        }
    return {"action": action, "groups": list(groups), "order": order_name, "matching_scale": matching_scale, "hadronic_scale": hadronic_scale}


def wilson_scan_coefficient_options() -> list[dict]:
    """All Wilson coefficients supported by the current core, with LaTeX labels."""
    try:
        valid = set(_mapper_strings(WCoefMapper()))
    except Exception:
        valid = {name for names in _WILSON_GROUP_COEFFS.values() for name in names}
    ordered: list[str] = []
    for names in _WILSON_GROUP_COEFFS.values():
        for name in names:
            if name in valid and name not in ordered:
                ordered.append(name)
    ordered.extend(sorted(valid.difference(ordered)))
    # Keep this large multi-select text-only. Dash 2.x ships an old
    # react-virtualized-select memoizer that can throw ``undefined.join`` when
    # component labels are used in a virtualized multi dropdown.
    return [
        {
            "label": f"{lx.compact_math_text(lx.wilson_latex(name), name)}  —  {name}",
            "value": name,
            "title": f"{lx.compact_math_text(lx.wilson_latex(name), name)} ({name})",
            "search": f"{name} {lx.text_latex(lx.wilson_latex(name), '')}",
        }
        for name in ordered
    ]


def _wilson_group_for_coefficient(coefficient: str) -> str:
    for key, names in _WILSON_GROUP_COEFFS.items():
        if coefficient in names:
            return _WILSON_GROUP_ALIASES[key][-1]
    raise ValueError(f"No Wilson group is known for coefficient '{coefficient}'.")


def wilson_scan_setup(
    coefficients: Sequence[str] | None,
    scan_mode: str | None,
    matching_scale: Any,
    hadronic_scale: Any,
    order_name: str | None,
) -> dict:
    """Normalize the high-level Wilson Scan controls."""
    selected = list(dict.fromkeys(str(c) for c in (coefficients or []) if c))
    if len(selected) > 10:
        raise ValueError("Wilson Scan supports at most 10 fitted coefficients.")
    normalized_mode = str(scan_mode or "DELTA").upper()
    if normalized_mode not in {"DELTA", "FULL"}:
        raise ValueError("Wilson Scan convention must be DELTA or FULL.")
    matching = float(matching_scale if matching_scale not in (None, "") else 160.0)
    hadronic = float(hadronic_scale if hadronic_scale not in (None, "") else 4.8)
    if matching <= 0.0 or hadronic <= 0.0:
        raise ValueError("Wilson matching and hadronic scales must be strictly positive.")
    normalized_order = str(order_name or "NNLO").upper()
    valid_orders = {item.name for item in QCDOrder}
    if normalized_order not in valid_orders:
        raise ValueError(f"Unknown Wilson evolution order '{normalized_order}'.")
    return {
        "coefficients": selected,
        "scan_mode": normalized_mode,
        "matching_scale": matching,
        "hadronic_scale": hadronic,
        "order_name": normalized_order,
    }


def ensure_wilson_scan_ready(setup: Mapping[str, Any], monitor: StatisticProgressMonitor | None = None) -> dict:
    """Build every Wilson group required by a scan before any stat operation."""
    coefficients = list(setup.get("coefficients") or [])
    if not coefficients:
        raise ValueError("Select at least one Wilson coefficient for Wilson Scan.")
    groups = sorted({_wilson_group_for_coefficient(name) for name in coefficients})
    matching = float(setup.get("matching_scale", 160.0))
    hadronic = float(setup.get("hadronic_scale", 4.8))
    order_name = str(setup.get("order_name", "NNLO")).upper()
    if monitor is not None:
        monitor.set_progress("wilson_build", "Building the required Wilson groups", 0.02)

    with RUNTIME.lock:
        status = dict(RUNTIME.wilson_status or {})
        existing = set(str(g) for g in status.get("groups", []))
        same_config = (
            RUNTIME.wilson is not None
            and str(status.get("order", "")).upper() == order_name
            and abs(float(status.get("matching_scale", matching)) - matching) < 1e-12
            and abs(float(status.get("hadronic_scale", hadronic)) - hadronic) < 1e-12
        )
        missing = [group for group in groups if group not in existing]
        if RUNTIME.wilson is None:
            result = build_wilson(groups, matching, hadronic, order_name, add=False)
        elif same_config and missing:
            result = build_wilson(missing, matching, hadronic, order_name, add=True)
        elif same_config:
            result = {"action": "reuse", "groups": groups, "order": order_name}
        else:
            # Preserve already selected groups when changing scales/order, while
            # guaranteeing all scan coefficients exist before statistics starts.
            result = build_wilson(sorted(existing.union(groups)), matching, hadronic, order_name, add=False)

    if monitor is not None:
        monitor.set_progress("wilson_build", f"Wilson groups ready: {', '.join(groups)}", 0.05)
    return result


def wilson_scan_parameter_rows(setup: Mapping[str, Any]) -> list[dict]:
    """Translate user-facing Wilson choices into internal statistic ParamIds."""
    ensure_wilson_scan_ready(setup)
    scan_mode = str(setup.get("scan_mode", "DELTA")).upper()
    if scan_mode not in {"DELTA", "FULL"}:
        raise ValueError("Wilson Scan mode must be DELTA or FULL.")
    rows: list[dict] = []
    internal_targets: dict[tuple[str, str], str] = {}
    for coefficient in list(setup.get("coefficients") or [])[:10]:
        coefficient = str(coefficient)
        group = _wilson_group_for_coefficient(coefficient)
        final_block = GroupMapper.block_name(group, ScaleType.HADRONIC, WilsonBasis.STANDARD)
        bsm_block = final_block + "__BSM_INTERMEDIATE"
        bsm_code = WCoefMapper.flha_full(coefficient, QCDOrder.LO, ContributionType.BSM).to_string()
        target_key = (bsm_block, bsm_code)
        previous = internal_targets.get(target_key)
        if previous is not None:
            raise ValueError(
                f"Wilson coefficients '{previous}' and '{coefficient}' resolve to the same "
                "internal FLHA parameter. Select only one of these aliases."
            )
        internal_targets[target_key] = coefficient

        # The dependency graph is driven from the BSM intermediate block.  Even
        # for a full-C scan, fitting the final TOTAL block directly would let an
        # upstream SM refresh overwrite the trial value.  The statistic core
        # therefore receives the stable BSM ParamId plus an affine offset so its
        # minimizer coordinate is the physical full coefficient C = C_SM + ΔC.
        try:
            bsm_initial, bsm_lower, bsm_upper = suggested_parameter_bounds(
                "WILSON", bsm_block, bsm_code
            )
        except Exception:
            bsm_initial, bsm_lower, bsm_upper = 0.0, -2.0, 2.0

        fit_offset = 0.0
        if scan_mode == "FULL":
            total_code = WCoefMapper.flha_full(
                coefficient, QCDOrder.LO, ContributionType.TOTAL
            ).to_string()
            try:
                total_initial, total_lower, total_upper = suggested_parameter_bounds(
                    "WILSON", final_block, total_code
                )
            except Exception:
                total_initial = float(bsm_initial)
                total_lower, total_upper = total_initial - 2.0, total_initial + 2.0
            fit_offset = float(total_initial) - float(bsm_initial)
            initial, lower, upper = total_initial, total_lower, total_upper
            contribution = ContributionType.TOTAL
        else:
            initial, lower, upper = bsm_initial, bsm_lower, bsm_upper
            half_width = max(2.0, abs(float(initial)) + 1.0)
            lower, upper = float(initial) - half_width, float(initial) + half_width
            contribution = ContributionType.BSM

        label = lx.wilson_request_table_latex(
            "R" if scan_mode == "DELTA" else "FM",
            coefficient,
            "LO",
            contribution.name,
        )
        rows.append({
            "parameter": label,
            "source": "Wilson ΔC" if scan_mode == "DELTA" else "Wilson full C",
            # Internal model target stays the BSM intermediate parameter.
            "type": "WILSON",
            "block": bsm_block,
            "code": bsm_code,
            "initial": float(initial),
            "lower_bound": float(lower),
            "upper_bound": float(upper),
            "fit_offset": float(fit_offset),
            "wilson_coefficient": coefficient,
            "wilson_group": group,
            "wilson_scan_mode": scan_mode,
            "wilson_matching_scale": float(setup.get("matching_scale", 160.0)),
            "wilson_hadronic_scale": float(setup.get("hadronic_scale", 4.8)),
            "wilson_order": str(setup.get("order_name", "NNLO")).upper(),
        })
    return rows


def query_wilson(method: str, group: str, coeff: str, order: str, contribution: str, basis: str = "STANDARD", component: str = "real") -> dict:
    require_initialized()
    if RUNTIME.wilson is None:
        raise RuntimeError("WilsonInterface is not built yet")
    req = WilsonRequest(
        group=str(group),
        coefficient=str(coeff),
        order=enum_by_name(QCDOrder, order),
        contribution=enum_by_name(ContributionType, contribution),
        wilson_basis=enum_by_name(WilsonBasis, basis),
    )
    wi = RUNTIME.wilson
    if method == "M":
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
        "coefficient_latex": lx.wilson_request_table_latex(method, coeff, order, contribution),
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


def _normalize_optional_float(value: Any) -> float | None:
    """Normalize UI float-like values for stable observable signatures."""
    if value in (None, ""):
        return None
    return float(value)


def _observable_requested_rows(
    mode: str,
    obs_names: Sequence[str] | None,
    decay_names: Sequence[str] | None,
    order_name: str,
    add_dependencies: bool,
    bin_low: float | None = None,
    bin_high: float | None = None,
    bins: Sequence[tuple[float, float]] | None = None,
) -> list[dict]:
    """Build a pure-Python description of the requested observable selection.

    This function intentionally does *not* touch the C++ ObservableInterface.
    It lets the Dash layer compare the requested selection to the currently
    registered one before calling ``add_observable`` again. Re-adding the same
    observable/bin can leave the C++ manager in an inconsistent state for some
    decay engines.
    """
    selected = observable_names_from_selection(mode, obs_names, decay_names)
    if not selected:
        raise ValueError("Select at least one observable or decay")

    rows: list[dict] = []
    for obs_name in selected:
        if bins:
            for low, high in bins:
                rows.append({
                    "observable": obs_name,
                    "bin_low": float(low),
                    "bin_high": float(high),
                    "order": str(order_name),
                    "dependencies": bool(add_dependencies),
                })
        elif bin_low is not None and bin_high is not None:
            rows.append({
                "observable": obs_name,
                "bin_low": float(bin_low),
                "bin_high": float(bin_high),
                "order": str(order_name),
                "dependencies": bool(add_dependencies),
            })
        else:
            rows.append({
                "observable": obs_name,
                "bin_low": None,
                "bin_high": None,
                "order": str(order_name),
                "dependencies": bool(add_dependencies),
            })

    # Deduplicate while preserving order. This also prevents duplicate bins from
    # being sent to C++ when a decay-level selection resolves to repeated obs ids.
    unique: list[dict] = []
    seen: set[tuple] = set()
    for row in rows:
        key = (
            row["observable"],
            row["bin_low"],
            row["bin_high"],
            row["order"],
            row["dependencies"],
        )
        if key in seen:
            continue
        seen.add(key)
        row["observable_label"] = lx.observable_table_label(row["observable"])
        row["observable_raw"] = row["observable"]
        unique.append(row)
    return unique


def _observable_signature(rows: Sequence[dict]) -> tuple:
    """Return a stable immutable key for an observable selection."""
    return tuple(
        (
            str(row.get("observable")),
            None if row.get("bin_low") is None else float(row.get("bin_low")),
            None if row.get("bin_high") is None else float(row.get("bin_high")),
            str(row.get("order")),
            bool(row.get("dependencies", False)),
        )
        for row in rows
    )


def _registry_signature(keys: set[tuple]) -> tuple:
    """Return a deterministic, sortable signature for registry keys."""
    return tuple(sorted(
        (str(k[0]), "" if k[1] is None else f"{float(k[1]):.17g}", "" if k[2] is None else f"{float(k[2]):.17g}")
        for k in keys
    ))


def _observable_add_key(row: Mapping[str, Any]) -> tuple:
    """Key used to prevent repeated C++ add_observable/add_bin calls.

    The C++ ObsManager stores bins at the decay level and ``DecayParent.add_bin``
    blindly appends to a vector.  The key deliberately ignores QCD order and
    dependency flags: once an observable/bin exists in a C++ manager, re-adding
    it only to change metadata is unsafe.  Use a fresh Hyperiso/LHA session if
    you need a different order for an already registered decay.
    """
    return (
        str(row.get("observable")),
        None if row.get("bin_low") is None else float(row.get("bin_low")),
        None if row.get("bin_high") is None else float(row.get("bin_high")),
    )


def _decay_name_for_observable(obs_name: str) -> str:
    """Return the decay enum name associated with an observable."""
    obs = enum_by_name(Observables, obs_name)
    try:
        decay = DecayMapper.get_decay(obs)
    except Exception:
        decay = DecayMapper().get_decay(obs)
    return decay.name


def _decay_bin_key(row: Mapping[str, Any]) -> tuple | None:
    low = row.get("bin_low")
    high = row.get("bin_high")
    if low is None or high is None:
        return None
    return (_decay_name_for_observable(str(row.get("observable"))), float(low), float(high))


def _add_observable_rows_once(
    oi: ObservableInterface,
    rows: Sequence[dict],
    registered_keys: set[tuple],
    registered_decay_bins: set[tuple],
) -> dict:
    """Add only missing observables/bins to a long-lived C++ interface.

    Important C++ details handled here:
    * ``ObsManager.add_obs(ObservableId, ...)`` does not check whether the
      observable was already present before touching the decay/order state.
    * ``ObsManager.add_obs(BinnedObservableId, ...)`` always calls
      ``DecayParent.add_bin``; bins are stored per decay, not per observable.
      Therefore, for several observables in the same decay and same bin, only
      the first one may go through the binned overload.  The following ones are
      registered as regular observables and inherit the decay-level bins.
    """
    added: list[tuple] = []
    skipped: list[tuple] = []

    for row in rows:
        key = _observable_add_key(row)
        if key in registered_keys:
            skipped.append(key)
            continue

        obs_name = str(row["observable"])
        obs_enum = enum_by_name(Observables, obs_name)
        qcd_order = enum_by_name(QCDOrder, str(row["order"]))
        add_dependencies = bool(row.get("dependencies", False))
        low = row.get("bin_low")
        high = row.get("bin_high")

        if low is not None and high is not None:
            dkey = _decay_bin_key(row)
            if dkey not in registered_decay_bins:
                oi.add_binned_observable(
                    make_binned_id(obs_name, float(low), float(high)),
                    qcd_order,
                    add_dependencies,
                )
                registered_decay_bins.add(dkey)
            else:
                # The decay already owns this bin.  Register only the observable
                # to avoid appending the same bin again to DecayParent::bins.
                oi.add_observable(obs_enum, qcd_order, add_dependencies)
        else:
            oi.add_observable(obs_enum, qcd_order, add_dependencies)

        registered_keys.add(key)
        added.append(key)

    return {"added": added, "skipped": skipped, "total": len(registered_keys)}


def _build_observable_interface_from_rows(rows: Sequence[dict]) -> ObservableInterface:
    """Create a fresh C++ ObservableInterface and add rows safely once."""
    oi = ObservableInterface()
    _add_observable_rows_once(oi, rows, set(), set())
    return oi


def build_observable_interface(
    mode: str,
    obs_names: Sequence[str] | None,
    decay_names: Sequence[str] | None,
    order_name: str,
    add_dependencies: bool,
    bin_low: float | None = None,
    bin_high: float | None = None,
    bins: Sequence[tuple[float, float]] | None = None,
) -> tuple[ObservableInterface, list[dict]]:
    require_initialized()
    rows = _observable_requested_rows(
        mode,
        obs_names,
        decay_names,
        order_name,
        add_dependencies,
        bin_low,
        bin_high,
        bins,
    )
    return _build_observable_interface_from_rows(rows), rows


def build_observables(
    mode: str,
    obs_names: Sequence[str] | None,
    decay_names: Sequence[str] | None,
    order_name: str,
    add_dependencies: bool,
    use_bin: bool,
    bin_low: Any,
    bin_high: Any,
    smooth: bool = False,
    smooth_min: Any = None,
    smooth_max: Any = None,
    smooth_step: Any = None,
) -> list[dict]:
    bins = None
    low = high = None
    if smooth:
        bins = bins_from_step(smooth_min, smooth_max, smooth_step)
    elif use_bin:
        low = _normalize_optional_float(bin_low)
        high = _normalize_optional_float(bin_high)
        if low is None or high is None:
            raise ValueError("Both bin low and bin high are required when binning is enabled")

    rows = _observable_requested_rows(
        mode,
        obs_names,
        decay_names,
        order_name,
        add_dependencies,
        low,
        high,
        bins,
    )

    with RUNTIME.lock:
        if RUNTIME.observable is None:
            RUNTIME.observable = ObservableInterface()
        info = _add_observable_rows_once(
            RUNTIME.observable,
            rows,
            RUNTIME.observable_registered_keys,
            RUNTIME.observable_registered_decay_bins,
        )
        RUNTIME.last_observable_selection = list(rows)
        RUNTIME.observable_signature = _registry_signature(RUNTIME.observable_registered_keys)

    for row in rows:
        row["registered"] = _observable_add_key(row) in RUNTIME.observable_registered_keys
    return rows


def _cpp_observable_values(oi: ObservableInterface, obs_name: str):
    """Compute one observable from the raw C++ interface when possible."""
    obs = enum_by_name(Observables, obs_name)
    raw = getattr(oi, "_to_cpp", lambda: None)()
    if raw is not None:
        return list(raw.compute_observable(obs.value))
    return list(oi.compute_observable(obs))


def _bin_tuple_from_cpp_value(value: Any) -> tuple[float | None, float | None]:
    b = getattr(value, "bin", None)
    if b is None:
        return None, None
    return float(b[0]), float(b[1])


def _observable_rows_for_target(
    oi: ObservableInterface,
    obs_name: str,
    bin_low: float | None = None,
    bin_high: float | None = None,
) -> list[dict]:
    """Compute rows for one target observable, optionally filtering one bin."""
    values = _cpp_observable_values(oi, obs_name)
    rows: list[dict] = []
    for value in values:
        low, high = _bin_tuple_from_cpp_value(value)
        if bin_low is not None and bin_high is not None:
            if low is None or high is None:
                continue
            if abs(low - float(bin_low)) > 1e-12 or abs(high - float(bin_high)) > 1e-12:
                continue
        raw_obs = str(getattr(value, "id", obs_name))
        rows.append({
            "observable_id": raw_obs,
            "observable_label": lx.observable_table_label(obs_name),
            "raw_observable": obs_name,
            "bin_low": low,
            "bin_high": high,
            "value": scalar_to_float(getattr(value, "value"), "real"),
        })

    if bin_low is not None and bin_high is not None and not rows:
        available = [
            _bin_tuple_from_cpp_value(v) for v in values
        ]
        raise RuntimeError(
            f"Observable {obs_name} did not return requested bin "
            f"[{bin_low}, {bin_high}]. Available bins: {available}"
        )
    return rows


def _observable_rows_for_selection(oi: ObservableInterface, rows: Sequence[dict]) -> list[dict]:
    """Compute configured rows one target at a time.

    This avoids relying on ``compute_all`` for interactive scans/computes, and
    prevents one stale observable inside the C++ manager from breaking all UI
    computations with ``map::at``.
    """
    out: list[dict] = []
    seen: set[tuple] = set()
    for row in rows:
        key = _observable_add_key(row)
        if key in seen:
            continue
        seen.add(key)
        out.extend(
            _observable_rows_for_target(
                oi,
                str(row["observable"]),
                row.get("bin_low"),
                row.get("bin_high"),
            )
        )
    return out


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
                raw_obs = str(value.id)
                rows.append({
                    "observable_id": raw_obs,
                    "observable_label": lx.observable_table_label(raw_obs),
                    "raw_observable": raw_obs,
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
            raw_obs = str(getattr(value, "id", _oid))
            rows.append({
                "observable_id": raw_obs,
                "observable_label": lx.observable_table_label(raw_obs),
                "raw_observable": raw_obs,
                "bin_low": low,
                "bin_high": high,
                "value": scalar_to_float(getattr(value, "value"), "real"),
            })
    return rows


def _rows_union(*collections: Sequence[dict]) -> list[dict]:
    """Return an ordered union of observable rows by C++ add key."""
    out: list[dict] = []
    seen: set[tuple] = set()
    for rows in collections:
        for row in rows or []:
            key = _observable_add_key(row)
            if key in seen:
                continue
            seen.add(key)
            out.append(dict(row))
    return out


def _rebuild_shared_observable_interface(rows: Sequence[dict]) -> None:
    """Replace the shared ObservableInterface and register the provided rows once."""
    RUNTIME.observable = ObservableInterface()
    RUNTIME.observable_registered_keys.clear()
    RUNTIME.observable_registered_decay_bins.clear()
    if rows:
        _add_observable_rows_once(
            RUNTIME.observable,
            rows,
            RUNTIME.observable_registered_keys,
            RUNTIME.observable_registered_decay_bins,
        )
    RUNTIME.observable_signature = _registry_signature(RUNTIME.observable_registered_keys)
    # Keep legacy aliases pointing to the shared registry.
    RUNTIME.scan_observable = RUNTIME.observable
    RUNTIME.stat_observable = RUNTIME.observable
    RUNTIME.scan_observable_registered_keys = RUNTIME.observable_registered_keys
    RUNTIME.scan_observable_registered_decay_bins = RUNTIME.observable_registered_decay_bins
    RUNTIME.stat_observable_registered_keys = RUNTIME.observable_registered_keys
    RUNTIME.stat_observable_registered_decay_bins = RUNTIME.observable_registered_decay_bins
    RUNTIME.stat = None


def remove_observable_rows(indices: Sequence[int] | None) -> list[dict]:
    """Remove rows from the main observable selection and rebuild safely.

    Since the C++ API can remove observables but not individual decay-level
    bins, rebuilding the single shared ObservableInterface from the remaining
    Python selection is the safest way to implement UI removal.
    """
    require_initialized()
    if not indices:
        return list(RUNTIME.last_observable_selection)
    remove = {int(i) for i in indices}
    with RUNTIME.lock:
        RUNTIME.last_observable_selection = [
            dict(row) for i, row in enumerate(RUNTIME.last_observable_selection) if i not in remove
        ]
        # Statistic rows may still be present; preserve them when rebuilding.
        union_rows = _rows_union(RUNTIME.last_observable_selection, RUNTIME.last_stat_observable_selection)
        _rebuild_shared_observable_interface(union_rows)
    return list(RUNTIME.last_observable_selection)


def remove_stat_observable_rows(indices: Sequence[int] | None) -> list[dict]:
    """Remove rows from the statistic observable selection and rebuild safely."""
    require_initialized()
    if not indices:
        return list(RUNTIME.last_stat_observable_selection)
    remove = {int(i) for i in indices}
    with RUNTIME.lock:
        RUNTIME.last_stat_observable_selection = [
            dict(row) for i, row in enumerate(RUNTIME.last_stat_observable_selection) if i not in remove
        ]
        union_rows = _rows_union(RUNTIME.last_observable_selection, RUNTIME.last_stat_observable_selection)
        _rebuild_shared_observable_interface(union_rows)
    return list(RUNTIME.last_stat_observable_selection)


def compute_current_observables() -> list[dict]:
    require_initialized()
    if RUNTIME.observable is None or not RUNTIME.last_observable_selection:
        raise RuntimeError("ObservableInterface has no configured observables yet")
    with RUNTIME.lock:
        return _observable_rows_for_selection(RUNTIME.observable, RUNTIME.last_observable_selection)


def _refresh_observable_interface(oi: ObservableInterface) -> None:
    """Compatibility hook kept for older code paths.

    Do not call ``ObservableInterface.reload_params()`` here: in the current C++
    implementation it loops over *all* decays and calls ``load_params()`` even
    for decays that were never enabled.  ``compute_all()`` already calls
    ``enable_obs()``, and enabled decays reload their parameters there, so scans
    should simply mutate the global parameter store and then compute.
    """
    return None


def _compute_observable_float_from_interface(oi: ObservableInterface, obs_name: str, bin_low: float | None, bin_high: float | None) -> float:
    """Compute one target value from an already configured observable interface."""
    rows = _observable_rows_for_target(oi, obs_name, bin_low, bin_high)
    if not rows:
        raise RuntimeError("Observable computation returned no value")
    return scalar_to_float(rows[0]["value"], "real")


def _get_cached_scan_observable_interface(
    obs_name: str,
    order_name: str,
    add_deps: bool,
    bin_low: float | None,
    bin_high: float | None,
) -> ObservableInterface:
    """Return the shared persistent ObservableInterface for scans.

    The scan target is added to the main interface only if it is not already
    registered.  This prevents a second C++ ObservableInterface from emitting a
    second ``Adding observable ...`` for the same target.
    """
    rows = _observable_requested_rows(
        "observable",
        [obs_name],
        None,
        order_name,
        add_deps,
        bin_low,
        bin_high,
        None,
    )
    oi = _ensure_rows_on_shared_observable_interface(rows)
    RUNTIME.scan_observable = oi
    RUNTIME.scan_observable_signature = _registry_signature(RUNTIME.observable_registered_keys)
    # Keep aliases for UI/debug only; they are no longer independent registries.
    RUNTIME.scan_observable_registered_keys = RUNTIME.observable_registered_keys
    RUNTIME.scan_observable_registered_decay_bins = RUNTIME.observable_registered_decay_bins
    return oi


def evaluate_observable_float(obs_name: str, order_name: str, add_deps: bool, bin_low: float | None, bin_high: float | None) -> float:
    require_initialized()
    with RUNTIME.lock:
        oi = _get_cached_scan_observable_interface(obs_name, order_name, add_deps, bin_low, bin_high)
        return _compute_observable_float_from_interface(oi, obs_name, bin_low, bin_high)


def observable_scan_1d(obs_name: str, order_name: str, add_deps: bool, bin_low: Any, bin_high: Any, ptype: str, block: str, code: Any, x_min: float, x_max: float, n_points: int) -> tuple[list[float], list[float]]:
    require_initialized()
    pid = make_param_id(ptype, block, code)
    low = as_float(bin_low, None)
    high = as_float(bin_high, None)
    xs = linspace(x_min, x_max, n_points)
    ys: list[float] = []
    setter = ParameterSetter()
    restore: dict[str, float] = {}
    label = param_to_label(pid)
    with RUNTIME.lock:
        oi = _get_cached_scan_observable_interface(obs_name, order_name, add_deps, low, high)
        try:
            for x in xs:
                _mutate_param_temporarily(pid, float(x), restore, setter)
                ys.append(_compute_observable_float_from_interface(oi, obs_name, low, high))
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
    z: list[list[float]] = []
    setter = ParameterSetter()
    restore: dict[str, float] = {}
    labels = [(pid1, param_to_label(pid1)), (pid2, param_to_label(pid2))]
    with RUNTIME.lock:
        oi = _get_cached_scan_observable_interface(obs_name, order_name, add_deps, low, high)
        try:
            for y in ys:
                row: list[float] = []
                _mutate_param_temporarily(pid2, float(y), restore, setter)
                for x in xs:
                    _mutate_param_temporarily(pid1, float(x), restore, setter)
                    row.append(_compute_observable_float_from_interface(oi, obs_name, low, high))
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

def make_stat_config(
    mc_draws: Any,
    mc_threads: Any,
    mc_seed: Any,
    skew_threshold: Any,
    ridge_rel: Any,
    ridge_abs: Any,
    nuisance_pruning: bool,
    nuisance_contexts: Any,
    nuisance_seed: Any,
    progress_monitor: StatisticProgressMonitor | None = None,
) -> StatisticConfig:
    def supplied(value: Any, default: Any) -> Any:
        return default if value is None or value == "" else value

    cfg = StatisticConfig()
    cfg.MC_draws = max(1, int(supplied(mc_draws, 100)))
    cfg.MC_threads = max(1, int(supplied(mc_threads, 1)))
    cfg.MC_seed = int(supplied(mc_seed, 123456))
    cfg.skew_abs_threshold = float(supplied(skew_threshold, 0.2))

    # A GUI should never emit terminal diagnostics.  Progress is represented by
    # the Dash progress element instead of the C++ stdout reporters.
    cfg.print_mc_progress = False
    cfg.print_chi2_pipeline_progress = False
    cfg.print_mc_config = False
    cfg.print_fit_summary = False
    cfg.print_scan_summary = False
    cfg.print_cache_summary = False
    cfg.print_debug = False

    advanced = cfg.advanced
    advanced.likelihood_mode = StatisticLikelihoodMode.CHI2_MC_COVARIANCE
    advanced.chi2_covariance_ridge_rel = float(supplied(ridge_rel, 1e-8))
    advanced.chi2_covariance_ridge_abs = float(supplied(ridge_abs, 1e-12))
    advanced.nuisance_sensitivity_pruning = bool(nuisance_pruning)
    advanced.nuisance_sensitivity_contexts = int(supplied(nuisance_contexts, 2))
    advanced.nuisance_sensitivity_seed = int(supplied(nuisance_seed, 12345))
    cfg.progress_monitor = progress_monitor
    # Dynamic attributes consumed by the existing StatisticInterface wrapper.
    cfg.p_specs = []
    cfg.selected_experiments = None
    return cfg


def _stat_rows_from_kwargs(
    mode: str,
    obs_names: Sequence[str] | None,
    decay_names: Sequence[str] | None,
    order_name: str,
    add_deps: bool,
    bin_strategy: str,
    bin_low: Any,
    bin_high: Any,
    smooth_min: Any,
    smooth_max: Any,
    smooth_step: Any,
) -> list[dict]:
    bins = None
    low = high = None
    if bin_strategy == "smooth":
        bins = bins_from_step(smooth_min, smooth_max, smooth_step)
    elif bin_strategy == "single":
        low = as_float(bin_low, None)
        high = as_float(bin_high, None)
    return _observable_requested_rows(mode, obs_names, decay_names, order_name, add_deps, low, high, bins)


def _ensure_rows_on_shared_observable_interface(rows: Sequence[dict]) -> ObservableInterface:
    """Add rows once to the single process-wide ObservableInterface.

    All observable consumers -- direct compute, scans and statistics -- must use
    this same interface.  The C++ ObservableInterface constructor builds a fresh
    ObsManager and fresh decay/Wilson ports; for the current Hyperiso singleton,
    duplicating those objects for the same observable can re-enter add_obs and
    leave the underlying manager/proxy state inconsistent.
    """
    if RUNTIME.observable is None:
        RUNTIME.observable = ObservableInterface()
    _add_observable_rows_once(
        RUNTIME.observable,
        rows,
        RUNTIME.observable_registered_keys,
        RUNTIME.observable_registered_decay_bins,
    )
    RUNTIME.observable_signature = _registry_signature(RUNTIME.observable_registered_keys)
    RUNTIME.scan_observable = RUNTIME.observable
    RUNTIME.stat_observable = RUNTIME.observable
    RUNTIME.scan_observable_registered_keys = RUNTIME.observable_registered_keys
    RUNTIME.scan_observable_registered_decay_bins = RUNTIME.observable_registered_decay_bins
    RUNTIME.stat_observable_registered_keys = RUNTIME.observable_registered_keys
    RUNTIME.stat_observable_registered_decay_bins = RUNTIME.observable_registered_decay_bins
    return RUNTIME.observable


def configure_stat_observables(
    mode: str,
    obs_names: Sequence[str] | None,
    decay_names: Sequence[str] | None,
    order_name: str,
    add_deps: bool,
    bin_strategy: str,
    bin_low: Any,
    bin_high: Any,
    smooth_min: Any,
    smooth_max: Any,
    smooth_step: Any,
) -> list[dict]:
    """Add statistic observables once to the persistent stat ObservableInterface."""
    rows = _stat_rows_from_kwargs(
        mode, obs_names, decay_names, order_name, add_deps,
        bin_strategy, bin_low, bin_high, smooth_min, smooth_max, smooth_step,
    )
    with RUNTIME.lock:
        oi = _ensure_rows_on_shared_observable_interface(rows)
        RUNTIME.stat_observable = oi
        RUNTIME.stat_observable_registered_keys = RUNTIME.observable_registered_keys
        RUNTIME.stat_observable_registered_decay_bins = RUNTIME.observable_registered_decay_bins
        RUNTIME.last_stat_observable_selection = list(rows)
    return rows


def build_stat_interface(mode: str, obs_names: Sequence[str] | None, decay_names: Sequence[str] | None, order_name: str, add_deps: bool, bin_strategy: str, bin_low: Any, bin_high: Any, smooth_min: Any, smooth_max: Any, smooth_step: Any, experiments: str | None, mc_draws: Any, mc_threads: Any, mc_seed: Any, skew_threshold: Any, ridge_rel: Any, ridge_abs: Any, nuisance_pruning: bool, nuisance_contexts: Any, nuisance_seed: Any, p_specs: Sequence[ParamId] | None = None, configured_rows: Sequence[dict] | None = None, progress_monitor: StatisticProgressMonitor | None = None, fit_parameter_bounds: Mapping[ParamId, tuple[float, float]] | None = None, fit_parameter_offsets: Mapping[ParamId, float] | None = None) -> StatisticInterface:
    if configured_rows:
        rows = [dict(row) for row in configured_rows]
        with RUNTIME.lock:
            oi = _ensure_rows_on_shared_observable_interface(rows)
            RUNTIME.stat_observable = oi
            RUNTIME.stat_observable_registered_keys = RUNTIME.observable_registered_keys
            RUNTIME.stat_observable_registered_decay_bins = RUNTIME.observable_registered_decay_bins
            RUNTIME.last_stat_observable_selection = list(rows)
    else:
        rows = configure_stat_observables(
            mode, obs_names, decay_names, order_name, add_deps,
            bin_strategy, bin_low, bin_high, smooth_min, smooth_max, smooth_step,
        )
    if RUNTIME.observable is None:
        raise RuntimeError("Shared ObservableInterface is not available")
    cfg = make_stat_config(
        mc_draws, mc_threads, mc_seed, skew_threshold, ridge_rel, ridge_abs,
        nuisance_pruning, nuisance_contexts, nuisance_seed, progress_monitor,
    )
    cfg.p_specs = list(p_specs or [])
    cfg.fit_parameter_bounds = dict(fit_parameter_bounds or {})
    cfg.fit_parameter_offsets = dict(fit_parameter_offsets or {})
    exps = split_csv(experiments)
    cfg.selected_experiments = exps or None
    stat_int = StatisticInterface(cfg, observable_interface=RUNTIME.observable)
    if exps:
        stat_int.select_experiments(exps)
    else:
        stat_int.select_experiments_all()
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
            "observable_label": lx.observable_table_label(str(bid.s)),
            "raw_observable": str(bid.s),
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


def p_spec_rows_and_ids(rows: Sequence[dict] | None) -> list[tuple[dict, ParamId]]:
    """Return validated fit-parameter rows paired with their ``ParamId``."""
    out: list[tuple[dict, ParamId]] = []
    for row in rows or []:
        if len(out) >= 10:
            break
        if not row:
            continue
        ptype, block, code = row.get("type"), row.get("block"), row.get("code")
        if ptype and block and code not in (None, ""):
            out.append((row, make_param_id(ptype, block, code)))
    return out


def p_specs_from_rows(rows: Sequence[dict]) -> list[ParamId]:
    return [pid for _, pid in p_spec_rows_and_ids(rows)]


def fit_result_rows(fit, display_rows: Sequence[dict] | None = None) -> tuple[list[dict], list[dict], list[dict]]:
    """Format fit output while preserving high-level Wilson Scan labels."""
    display_by_key = {p_spec_row_key(row): row for row in (display_rows or [])}

    def display(pid: ParamId) -> tuple[str, str]:
        raw = param_to_label(pid)
        key = "|".join((pid.type.name if pid.type else "", str(pid.block), pid.code.to_string()))
        row = display_by_key.get(key, {})
        label = str(row.get("parameter") or lx.parameter_table_label(
            str(pid.block), pid.code.to_string(), raw, pid.type.name if pid.type else None
        ))
        return label, raw

    p_rows = []
    for p, v in fit.p_hat.items():
        label, raw = display(p)
        p_rows.append({
            "parameter": label,
            "raw_parameter": raw,
            "best_fit": scalar_to_float(v, "real"),
            "std": scalar_to_float(fit.p_hat_std.get(p, float("nan")), "real") if fit.p_hat_std else None,
        })
    eta_rows = []
    for p, v in fit.eta_hat.items():
        raw = param_to_label(p)
        eta_rows.append({
            "nuisance": lx.parameter_table_label(str(p.block), p.code.to_string(), raw, p.type.name if p.type else None),
            "raw_nuisance": raw,
            "value": scalar_to_float(v, "real"),
        })
    corr_rows = []
    for p1, inner in fit.p_correlations.items():
        label1, raw1 = display(p1)
        for p2, corr in inner.items():
            label2, raw2 = display(p2)
            corr_rows.append({
                "x": label1, "y": label2, "raw_x": raw1, "raw_y": raw2,
                "corr": scalar_to_float(corr, "real"),
            })
    return p_rows, eta_rows, corr_rows


def _enum_member(enum_cls, name: str, fallback: str):
    key = str(name or fallback).upper()
    try:
        return enum_cls[key]
    except KeyError as exc:
        valid = ", ".join(item.name for item in enum_cls)
        raise ValueError(f"Unknown {enum_cls.__name__} '{name}'. Valid values: {valid}") from exc


def _row_bounds(row: Mapping[str, Any]) -> tuple[float, float]:
    low = float(row.get("lower_bound"))
    high = float(row.get("upper_bound"))
    if not low < high:
        raise ValueError(f"Invalid bounds for {row.get('parameter')}: lower bound must be smaller than upper bound.")
    return low, high


def run_fit_and_contours(
    stat_kwargs: dict,
    p_spec_rows: Sequence[dict],
    do_contour: bool,
    x_key: str | None,
    y_key: str | None,
    confidence_levels: Sequence[Any] | None,
    profiling_method: str | None,
    contour_algorithm: str | None,
    fallback_algorithm: str | None,
    profiler_mode: str | None,
    resolution: Any,
    progress_monitor: StatisticProgressMonitor | None = None,
    wilson_setup: Mapping[str, Any] | None = None,
) -> dict:
    """Run the χ² fit and optional core confidence contours.

    Bounds are read from the two selected rows in the fit-parameter table.  For
    fits with more than two parameters, ``ProfilingMethod`` determines how the
    hidden fit coordinates are reduced to the displayed plane.
    """
    require_initialized()
    if wilson_setup is not None:
        ensure_wilson_scan_ready(wilson_setup, progress_monitor)
    rows = list(p_spec_rows or [])
    row_pid_pairs = p_spec_rows_and_ids(rows)
    p_specs = [pid for _, pid in row_pid_pairs]
    if not p_specs:
        raise ValueError("Add at least one fit parameter.")

    row_by_key = {p_spec_row_key(row): row for row, _ in row_pid_pairs}
    pid_by_key = {p_spec_row_key(row): pid for row, pid in row_pid_pairs}

    stat_kwargs = dict(stat_kwargs)
    stat_kwargs["p_specs"] = p_specs
    stat_kwargs["progress_monitor"] = progress_monitor
    stat_kwargs["fit_parameter_bounds"] = {
        pid: _row_bounds(row) for row, pid in row_pid_pairs
    }
    stat_kwargs["fit_parameter_offsets"] = {
        pid: float(row.get("fit_offset", 0.0) or 0.0)
        for row, pid in row_pid_pairs
        if float(row.get("fit_offset", 0.0) or 0.0) != 0.0
    }
    with RUNTIME.lock:
        RUNTIME.stat = build_stat_interface(**stat_kwargs)
        fit = RUNTIME.stat.compute_MLE(p_specs)
        p_rows, eta_rows, corr_rows = fit_result_rows(fit, rows)

        contour_paths: list[dict] = []
        contour_errors: list[str] = []
        axis_labels = ("p₁", "p₂")
        best_fit_point = None
        bounds = None

        if do_contour:
            if len(p_specs) < 2:
                raise ValueError("A confidence contour requires at least two fit parameters.")
            if not x_key or not y_key or x_key not in pid_by_key or y_key not in pid_by_key:
                raise ValueError("Choose two valid contour-axis parameters.")
            if x_key == y_key:
                raise ValueError("The X and Y contour axes must be different parameters.")
            if not fit.fit_ok:
                raise RuntimeError("The χ² fit did not converge; no contour can be computed.")

            x_pid = pid_by_key[x_key]
            y_pid = pid_by_key[y_key]
            x_row = row_by_key[x_key]
            y_row = row_by_key[y_key]
            x_low, x_high = _row_bounds(x_row)
            y_low, y_high = _row_bounds(y_row)
            bounds = [x_low, x_high, y_low, y_high]
            axis_labels = (str(x_row.get("parameter") or x_key), str(y_row.get("parameter") or y_key))

            x_hat = scalar_to_float(fit.p_hat[x_pid], "real")
            y_hat = scalar_to_float(fit.p_hat[y_pid], "real")
            if not (x_low <= x_hat <= x_high):
                raise ValueError(f"The X best fit ({x_hat:.6g}) lies outside its bounds [{x_low:.6g}, {x_high:.6g}].")
            if not (y_low <= y_hat <= y_high):
                raise ValueError(f"The Y best fit ({y_hat:.6g}) lies outside its bounds [{y_low:.6g}, {y_high:.6g}].")
            best_fit_point = {"x": x_hat, "y": y_hat}

            method = ProfilingMethod.SLICE if len(p_specs) == 2 else _enum_member(ProfilingMethod, profiling_method or "SLICE", "SLICE")
            primary = _enum_member(ContourAlgorithm, contour_algorithm or "MINUIT", "MINUIT")
            fallback = None
            fallback_name = str(fallback_algorithm or "NONE").upper()
            if fallback_name != "NONE" and fallback_name != primary.name:
                fallback = _enum_member(ContourAlgorithm, fallback_name, "AMS")
            backend = _enum_member(ProfilerMode, profiler_mode or "LAPLACE_NUISANCE", "LAPLACE_NUISANCE")
            contour_options = ContourOptions(
                profiling_method=method,
                profile_backend=backend,
                primary_contour_method=primary,
                fallback_contour_method=fallback,
                resolution=max(8, min(500, int(resolution or 60))),
            )

            levels = sorted({float(z) for z in (confidence_levels or [1.0]) if float(z) > 0.0})
            if not levels:
                raise ValueError("Select at least one positive confidence level.")
            if progress_monitor is not None:
                progress_monitor.set_progress("contours", "Computing confidence contours", 0.0, total=len(levels))
            for level_index, z in enumerate(levels):
                if progress_monitor is not None:
                    progress_monitor.set_progress(
                        "contours", f"Computing the {z:g}σ contour",
                        level_index / max(1, len(levels)),
                        completed=level_index, total=len(levels),
                    )
                try:
                    contour = RUNTIME.stat.compute_confidence_contour(
                        x_pid, y_pid, z, bounds, contour_options
                    )
                    for path_id, path in enumerate(contour.paths):
                        if len(path) < 2:
                            continue
                        contour_paths.append({
                            "sigma": z,
                            "level": contour.level,
                            "path_id": path_id,
                            "success": contour.success,
                            "points": [{"x": x, "y": y} for x, y in path],
                        })
                    if not contour.success:
                        contour_errors.append(f"{z:g}σ contour returned success=false")
                    elif not contour.paths:
                        contour_errors.append(f"{z:g}σ contour returned no path")
                except Exception as exc:
                    contour_errors.append(f"{z:g}σ: {type(exc).__name__}: {exc}")

    if progress_monitor is not None:
        progress_monitor.set_progress("complete", "Statistic workflow complete", 1.0, finished=True)

    return {
        "fit_ok": bool(fit.fit_ok),
        "ell_hat": scalar_to_float(fit.ell_hat, "real"),
        "p_rows": p_rows,
        "eta_rows": eta_rows,
        "corr_rows": corr_rows,
        "contour_paths": contour_paths,
        "contour_errors": contour_errors,
        "axis_labels": axis_labels,
        "best_fit_point": best_fit_point,
        "bounds": bounds,
    }


def _prune_statistic_jobs(max_age_seconds: float = 3600.0) -> None:
    now = time.monotonic()
    with _STATISTIC_JOBS_LOCK:
        stale = [
            job_id for job_id, job in _STATISTIC_JOBS.items()
            if job.done and now - job.created_at > max_age_seconds
        ]
        for job_id in stale:
            _STATISTIC_JOBS.pop(job_id, None)


def _start_statistic_job(kind: str, target) -> str:
    _prune_statistic_jobs()
    job_id = uuid.uuid4().hex
    monitor = StatisticProgressMonitor()
    monitor.reset("queued", "Statistic task queued")
    job = StatisticJob(job_id=job_id, kind=kind, monitor=monitor)

    def runner() -> None:
        try:
            result = target(monitor)
            with job.lock:
                job.result = result
                job.done = True
        except Exception as exc:
            monitor.set_progress("failed", f"{type(exc).__name__}: {exc}", 1.0, finished=True)
            with job.lock:
                job.error = f"{type(exc).__name__}: {exc}"
                job.done = True

    thread = threading.Thread(target=runner, name=f"hyperiso-{kind}-{job_id[:8]}", daemon=True)
    job.thread = thread
    with _STATISTIC_JOBS_LOCK:
        _STATISTIC_JOBS[job_id] = job
    thread.start()
    return job_id


def start_uncertainty_job(
    stat_kwargs: Mapping[str, Any],
    wilson_setup: Mapping[str, Any] | None = None,
) -> str:
    """Start uncertainty propagation in the same process as the Hyperiso singleton."""
    kwargs = dict(stat_kwargs)

    def target(monitor: StatisticProgressMonitor):
        if wilson_setup is not None:
            ensure_wilson_scan_ready(wilson_setup, monitor)
        kwargs["progress_monitor"] = monitor
        rows = compute_uncertainty_rows(**kwargs)
        monitor.set_progress("complete", "Uncertainty propagation complete", 1.0, finished=True)
        return rows

    return _start_statistic_job("uncertainty", target)


def start_fit_job(
    stat_kwargs: Mapping[str, Any],
    p_spec_rows: Sequence[dict],
    do_contour: bool,
    x_key: str | None,
    y_key: str | None,
    confidence_levels: Sequence[Any] | None,
    profiling_method: str | None,
    contour_algorithm: str | None,
    fallback_algorithm: str | None,
    profiler_mode: str | None,
    resolution: Any,
    wilson_setup: Mapping[str, Any] | None = None,
) -> str:
    """Start a χ² fit/contour task with live C++ progress snapshots."""
    kwargs = dict(stat_kwargs)
    rows = [dict(row) for row in (p_spec_rows or [])]

    def target(monitor: StatisticProgressMonitor):
        return run_fit_and_contours(
            kwargs, rows, do_contour, x_key, y_key, confidence_levels,
            profiling_method, contour_algorithm, fallback_algorithm,
            profiler_mode, resolution, progress_monitor=monitor,
            wilson_setup=wilson_setup,
        )

    return _start_statistic_job("fit", target)


def _statistic_job_percent(kind: str, phase: str, fraction: float, done: bool) -> float:
    frac = max(0.0, min(1.0, float(fraction)))
    if done or phase == "complete":
        return 100.0
    if phase == "failed":
        return 100.0
    if phase in {"queued", "preparing", "wilson_build"}:
        return min(5.0, max(1.0, frac * 100.0))
    if phase == "monte_carlo":
        return 5.0 + (90.0 if kind == "uncertainty" else 70.0) * frac
    if phase == "chi2_pipeline":
        # Continue from the completed fit-MC segment instead of jumping back to
        # the beginning when the covariance/likelihood pipeline starts.
        return 75.0 + 20.0 * frac
    if phase == "contours":
        return 95.0 + 5.0 * frac
    return min(99.0, max(1.0, frac * 100.0))


def _duration_label(seconds: float) -> str:
    if seconds < 0.0:
        return ""
    seconds_i = max(0, int(round(seconds)))
    if seconds_i < 60:
        return f"{seconds_i}s"
    minutes, sec = divmod(seconds_i, 60)
    if minutes < 60:
        return f"{minutes}m {sec:02d}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h {minutes:02d}m"


def poll_statistic_job(job_id: str | None) -> dict:
    """Return the latest progress snapshot and final result, when available."""
    if not job_id:
        raise ValueError("Missing statistic job id")
    with _STATISTIC_JOBS_LOCK:
        job = _STATISTIC_JOBS.get(str(job_id))
    if job is None:
        raise KeyError("Statistic job is no longer available")
    snapshot = job.monitor.snapshot()
    with job.lock:
        done = bool(job.done)
        result = job.result if done else None
        error = job.error
    percent = _statistic_job_percent(job.kind, snapshot.phase, snapshot.fraction, done)
    meta_parts = []
    if snapshot.total > 0:
        meta_parts.append(f"{snapshot.completed}/{snapshot.total}")
    if snapshot.elapsed_seconds > 0:
        meta_parts.append(f"elapsed {_duration_label(snapshot.elapsed_seconds)}")
    if snapshot.eta_seconds >= 0 and not done:
        meta_parts.append(f"ETA {_duration_label(snapshot.eta_seconds)}")
    if snapshot.attempts > 0 and snapshot.attempts != snapshot.completed:
        meta_parts.append(f"attempts {snapshot.attempts}")
    if snapshot.failures > 0:
        meta_parts.append(f"rejected {snapshot.failures}")
    return {
        "job_id": job.job_id,
        "kind": job.kind,
        "done": done,
        "error": error,
        "result": result,
        "phase": snapshot.phase,
        "message": snapshot.message or "Statistic computation in progress",
        "percent": round(percent, 2),
        "meta": " · ".join(meta_parts),
    }


# Backward-compatible service name for external callers of the Dash helper.
def run_fit_and_scan(stat_kwargs: dict, p_spec_rows: Sequence[dict], do_contour: bool, x_half_width: Any, y_half_width: Any, nx: Any, ny: Any) -> dict:
    """Deprecated compatibility wrapper using table bounds and core contours."""
    rows = list(p_spec_rows or [])
    options = p_spec_axis_options(rows)
    x_key = options[0]["value"] if len(options) > 0 else None
    y_key = options[1]["value"] if len(options) > 1 else None
    for row, half_width in zip(rows[:2], (x_half_width, y_half_width)):
        center = float(row.get("initial", 0.0) or 0.0)
        width = abs(float(half_width or 1.0))
        row.setdefault("lower_bound", center - width)
        row.setdefault("upper_bound", center + width)
    return run_fit_and_contours(
        stat_kwargs, rows, do_contour, x_key, y_key, [1.0, 2.0],
        "SLICE", "AMS", "NONE", "MINUIT", max(int(nx or 25), int(ny or 25)),
    )


# ---------- QCD / QED running ----------

_PARTICLE_LABELS = {
    1: "d quark",
    2: "u quark",
    3: "s quark",
    4: "c quark",
    5: "b quark",
    6: "t quark",
    11: "electron",
    13: "muon",
    15: "tau",
}


def mass_particle_options() -> list[dict[str, str | int]]:
    """Return common PDG ids for the QCD mass calculator."""
    return [{"label": f"{name} ({pid})", "value": pid} for pid, name in _PARTICLE_LABELS.items()]


def mass_type_options() -> list[dict[str, str]]:
    """Return MassType options usable by the running providers."""
    return [{"label": mt.name, "value": mt.name} for mt in MassType]


def _mass_type(name: str | None, fallback: str = "POLE") -> MassType:
    return enum_by_name(MassType, name or fallback)


def _qcd_provider() -> QCDProvider:
    require_initialized()
    return QCDProvider()


def qcd_constants_rows(kind: str = "all") -> list[dict]:
    """Return QCD constants as table rows."""
    qcd = _qcd_provider()
    consts = qcd.get_qcd_constants()
    kind = str(kind or "all").lower()
    rows = []
    if kind in {"all", "color"}:
        rows.extend([
            {"quantity": "N_c", "index": "", "value": consts.Nc},
            {"quantity": "C_F", "index": "", "value": consts.C_F},
            {"quantity": "C_A", "index": "", "value": consts.C_A},
        ])
    if kind in {"all", "beta"}:
        try:
            for i, value in enumerate(list(consts.beta)):
                rows.append({"quantity": "beta", "index": i, "value": value})
        except Exception:
            pass
    if kind in {"all", "gamma"}:
        try:
            for i, value in enumerate(list(consts.gamma)):
                rows.append({"quantity": "gamma", "index": i, "value": value})
        except Exception:
            pass
    return rows


def qcd_alphas(scale: float, mb_type: str = "POLE", mt_type: str = "POLE") -> float:
    """Compute alpha_s at one scale."""
    qcd = _qcd_provider()
    value = qcd.get_alphas(AlphasConfig(float(scale), _mass_type(mb_type), _mass_type(mt_type)))
    return scalar_to_float(value, "real")


def qcd_alphaem(scale: float, mb_type: str = "POLE", mt_type: str = "POLE") -> float:
    """Compute alpha_em at one scale through QEDProvider."""
    require_initialized()
    qed = QEDProvider()
    value = qed.get_alpha_em(AlphasConfig(float(scale), _mass_type(mb_type), _mass_type(mt_type)))
    return scalar_to_float(value, "real")


def qcd_mass(scale: float, pdg_id: int, mb_type: str = "POLE", mt_type: str = "POLE") -> float:
    """Compute the running mass for one PDG id and scale."""
    qcd = _qcd_provider()
    cfg = MassConfig(float(scale), _mass_type(mb_type), _mass_type(mt_type), int(pdg_id))
    value = qcd.get_qcd_masses(cfg)
    return scalar_to_float(value, "real")


def qcd_scan_alphas(scale_min: float, scale_max: float, n_points: int, mb_type: str = "POLE", mt_type: str = "POLE") -> tuple[list[float], list[float]]:
    """Compute alpha_s on a scale grid."""
    xs = linspace(float(scale_min), float(scale_max), int(n_points or 1))
    ys = [qcd_alphas(x, mb_type, mt_type) for x in xs]
    return xs, ys


def qcd_scan_alphaem(scale_min: float, scale_max: float, n_points: int, mb_type: str = "POLE", mt_type: str = "POLE") -> tuple[list[float], list[float]]:
    """Compute alpha_em on a scale grid."""
    xs = linspace(float(scale_min), float(scale_max), int(n_points or 1))
    ys = [qcd_alphaem(x, mb_type, mt_type) for x in xs]
    return xs, ys


def qcd_scan_mass(scale_min: float, scale_max: float, n_points: int, pdg_id: int, mb_type: str = "POLE", mt_type: str = "POLE") -> tuple[list[float], list[float]]:
    """Compute a running mass on a scale grid."""
    xs = linspace(float(scale_min), float(scale_max), int(n_points or 1))
    ys = [qcd_mass(x, int(pdg_id), mb_type, mt_type) for x in xs]
    return xs, ys


def qcd_single_result_rows(scale: float, pdg_id: int, mb_type: str = "POLE", mt_type: str = "POLE", include_qed: bool = False) -> list[dict]:
    """Return a compact table containing alpha_s and one running mass."""
    rows = [
        {"quantity": r"$\alpha_s(\mu)$", "scale": float(scale), "pdg_id": "", "scheme": f"m_b={mb_type}, m_t={mt_type}", "value": qcd_alphas(scale, mb_type, mt_type)},
        {"quantity": r"$m(\mu)$", "scale": float(scale), "pdg_id": int(pdg_id), "scheme": f"m_b={mb_type}, m_t={mt_type}", "value": qcd_mass(scale, int(pdg_id), mb_type, mt_type)},
    ]
    if include_qed:
        try:
            rows.append({"quantity": r"$\alpha_{\rm em}(\mu)$", "scale": float(scale), "pdg_id": "", "scheme": f"m_b={mb_type}, m_t={mt_type}", "value": qcd_alphaem(scale, mb_type, mt_type)})
        except Exception as exc:
            rows.append({"quantity": r"$\alpha_{\rm em}(\mu)$", "scale": float(scale), "pdg_id": "", "scheme": "unavailable", "value": str(exc)})
    return rows
