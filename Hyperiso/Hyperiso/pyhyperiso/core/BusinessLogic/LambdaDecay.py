"""Python wrappers for lambda-backed custom decays and observables.

The C++ binding exposes the low-level ``LambdaDecay`` callback context.  This
module keeps the public Python API wrapper-only: user callbacks receive
``LambdaDecayContext`` and normal Python ids/enums instead of raw pybind objects.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Sequence, Set

from pyhyperiso.phyperiso.pyhyperiso import observable as _obs
from pyhyperiso.core.Common.GeneralEnum import QCDOrder, ContributionType, DataType
from pyhyperiso.core.Common.LhaID import LhaID
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Common.SymbolId import WGroupId, ObservableId
from pyhyperiso.core.Common.Mapper import GroupMapper
from pyhyperiso.core.PhysicalModel.WilsonInterface import CustomWilsonGroupConfig
from pyhyperiso.core.Common.Configs import _cpp_group_id, _cpp_coef_id
from pyhyperiso.core.Math.Scalar import Scalar


def _enum_value(value):
    """Return the bound C++ enum stored in the Python enum wrapper.

    The public wrappers use Python ``Enum`` objects whose ``.value`` is the
    pybind enum.  If a caller already passes the bound enum, it is returned as-is.
    """
    return getattr(value, "value", value)


def _to_observable_id(obs_id) -> ObservableId:
    """Convert a C++/Python observable id to the Python wrapper id."""
    if isinstance(obs_id, ObservableId):
        return obs_id
    return ObservableId(str(obs_id))


class LambdaDecayContext:
    """Wrapper around the C++ ``LambdaDecay`` callback context.

    User-defined observable callbacks receive this object.  It hides the raw
    pybind context and accepts the normal Python wrappers for Wilson ids,
    parameters and enums.
    """

    def __init__(self, cpp_obj):
        self._cpp_obj = cpp_obj

    def get_M(self, group, coeff, order: QCDOrder, contribution: ContributionType) -> Scalar:
        """Return a matching coefficient as ``Scalar``."""
        return Scalar.from_cpp(
            self._cpp_obj.get_M(
                _cpp_group_id(group),
                _cpp_coef_id(coeff),
                _enum_value(order),
                _enum_value(contribution),
            )
        )

    def get_FM(self, group, coeff, order: QCDOrder, contribution: ContributionType) -> Scalar:
        """Return a full matching coefficient as ``Scalar``."""
        return Scalar.from_cpp(
            self._cpp_obj.get_FM(
                _cpp_group_id(group),
                _cpp_coef_id(coeff),
                _enum_value(order),
                _enum_value(contribution),
            )
        )

    def get_R(self, group, coeff, order: QCDOrder, contribution: ContributionType) -> Scalar:
        """Return a running coefficient as ``Scalar``."""
        return Scalar.from_cpp(
            self._cpp_obj.get_R(
                _cpp_group_id(group),
                _cpp_coef_id(coeff),
                _enum_value(order),
                _enum_value(contribution),
            )
        )

    def get_FR(self, group, coeff, order: QCDOrder, contribution: ContributionType) -> Scalar:
        """Return a full running coefficient as ``Scalar``."""
        return Scalar.from_cpp(
            self._cpp_obj.get_FR(
                _cpp_group_id(group),
                _cpp_coef_id(coeff),
                _enum_value(order),
                _enum_value(contribution),
            )
        )

    def get_sm_param(self, pid: ParamId, data_type: DataType = DataType.VALUE) -> Scalar:
        """Return an SM parameter as ``Scalar``."""
        return Scalar.from_cpp(self._cpp_obj.get_sm_param(pid.to_cpp(), _enum_value(data_type)))

    def get_flavor_param(self, pid: ParamId, data_type: DataType = DataType.VALUE) -> Scalar:
        """Return a FLAVOR parameter as ``Scalar``."""
        return Scalar.from_cpp(self._cpp_obj.get_flavor_param(pid.to_cpp(), _enum_value(data_type)))

    def current_bins(self):
        """Return the currently requested bins."""
        return self._cpp_obj.current_bins()


class LambdaObservableConfig:
    """Runtime observable computed by a Python callable.

    Args:
        canonical: Canonical observable name registered in the mapper layer.
        compute: Callback receiving ``(ctx, observable_id)`` for scalar
            observables, or ``(ctx, bin, observable_id)`` for binned observables.
        aliases: Optional aliases accepted by ``ObservableMapper.id_of``.
        flha: Optional FLHA id for experimental input/output conventions.
            Required when the observable is passed to ``StatisticInterface``.
        dependencies: Parameter dependencies visible to Statistic.
        binned: Whether ``compute`` should be called with the current bin.
    """

    def __init__(
        self,
        canonical: str,
        compute: Callable,
        aliases: Sequence[str] | None = None,
        flha: LhaID | None = None,
        dependencies: Sequence[ParamId] | Set[ParamId] | None = None,
        binned: bool = False,
    ):
        if binned:

            def wrapped_compute(ctx, bin_range, obs_id):
                return compute(LambdaDecayContext(ctx), bin_range, _to_observable_id(obs_id))

            factory = _obs.LambdaObservableConfig.binned_scalar
        else:

            def wrapped_compute(ctx, obs_id):
                return compute(LambdaDecayContext(ctx), _to_observable_id(obs_id))

            factory = _obs.LambdaObservableConfig.scalar

        self._cpp_obj = factory(canonical, wrapped_compute)
        self._cpp_obj.aliases = list(aliases or [])
        self._cpp_obj.flha = None if flha is None else flha.to_cpp()
        self._cpp_obj.dependencies = {
            p.to_cpp() if isinstance(p, ParamId) else p for p in (dependencies or [])
        }

    @classmethod
    def scalar(
        cls,
        canonical: str,
        compute: Callable,
        aliases: Sequence[str] | None = None,
        flha: LhaID | None = None,
        dependencies: Sequence[ParamId] | Set[ParamId] | None = None,
    ) -> "LambdaObservableConfig":
        """Create an unbinned scalar observable.

        ``compute`` receives ``(ctx, observable_id)`` and can return ``float`` or
        any object implementing ``__float__`` such as ``Scalar``.
        """
        return cls(
            canonical, compute, aliases=aliases, flha=flha, dependencies=dependencies, binned=False
        )

    @classmethod
    def binned_scalar(
        cls,
        canonical: str,
        compute: Callable,
        aliases: Sequence[str] | None = None,
        flha: LhaID | None = None,
        dependencies: Sequence[ParamId] | Set[ParamId] | None = None,
    ) -> "LambdaObservableConfig":
        """Create a binned scalar observable.

        ``compute`` receives ``(ctx, (q2_min, q2_max), observable_id)`` and can
        return ``float`` or any object implementing ``__float__``.
        """
        return cls(
            canonical, compute, aliases=aliases, flha=flha, dependencies=dependencies, binned=True
        )

    def to_cpp(self):
        """Return the bound C++ ``LambdaObservableConfig``."""
        return self._cpp_obj


@dataclass
class LambdaDecayConfig:
    """Runtime decay backed by Python observable lambdas.

    This config is the Python counterpart of C++ ``LambdaDecayConfig``. It can
    declare builtin Wilson groups, lambda-backed custom Wilson groups and one or
    more custom observables. Declared dependencies are propagated to Statistic.
    """

    canonical: str
    observables: Sequence[LambdaObservableConfig]
    aliases: Sequence[str] = field(default_factory=list)
    matching_scale: float = 81.0
    hadronic_scale: float = 4.8
    order: QCDOrder = QCDOrder.LO
    max_order: QCDOrder = QCDOrder.NNLO
    wilson_groups: Sequence[WGroupId | str] = field(default_factory=list)
    custom_wilson_groups: Sequence[CustomWilsonGroupConfig] = field(default_factory=list)
    propagate_custom_wilson_dependencies: bool = True

    def to_cpp(self):
        """Convert to the bound C++ ``LambdaDecayConfig``."""
        cpp = _obs.LambdaDecayConfig()
        cpp.canonical = self.canonical
        cpp.aliases = list(self.aliases)
        cpp.matching_scale = float(self.matching_scale)
        cpp.hadronic_scale = float(self.hadronic_scale)
        cpp.order = _enum_value(self.order)
        cpp.max_order = _enum_value(self.max_order)
        cpp.wilson_groups = {
            g._to_cpp() if isinstance(g, WGroupId) else GroupMapper.id_of(g)._to_cpp()
            for g in self.wilson_groups
        }
        cpp.custom_wilson_groups = [
            g.to_cpp() if isinstance(g, CustomWilsonGroupConfig) else g
            for g in self.custom_wilson_groups
        ]
        cpp.observables = [
            o.to_cpp() if isinstance(o, LambdaObservableConfig) else o for o in self.observables
        ]
        cpp.propagate_custom_wilson_dependencies = bool(self.propagate_custom_wilson_dependencies)
        return cpp


__all__ = ["LambdaObservableConfig", "LambdaDecayConfig", "LambdaDecayContext"]
