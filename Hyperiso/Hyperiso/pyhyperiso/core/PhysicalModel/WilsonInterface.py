"""Python wrapper for Wilson coefficients, including dynamic lambda groups.

The C++ ``WilsonInterface`` computes builtin Wilson groups and now also accepts
runtime groups/coefficient ids. This module keeps the old enum-based API working
while exposing the new ``WGroupId``/``WCoefId`` and custom-lambda workflow.
"""

from __future__ import annotations

from typing import Dict, Mapping, Sequence, Set

from pyhyperiso.phyperiso.pyhyperiso.wilson.wilson_interface import (
    WilsonInterface as _CppWilsonInterface,
    CustomWilsonCoefficientConfig as _CppCustomWilsonCoefficientConfig,
    CustomWilsonGroupConfig as _CppCustomWilsonGroupConfig,
)
from pyhyperiso.core.Common.Configs import (
    WilsonBuildConfig,
    WilsonRequest,
    WilsonGroupLike,
    WilsonCoefLike,
    _cpp_group_id,
    _cpp_coef_id,
)
from pyhyperiso.core.Common.GeneralEnum import (
    QCDOrder,
    WCoeff,
    WGroup,
    ContributionType,
    WilsonBasis,
    ParameterType,
)
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Common.SymbolId import WCoefId
from pyhyperiso.core.Math.Scalar import Scalar


def _cpp_order(order: QCDOrder):
    return order.value if isinstance(order, QCDOrder) else order


def _cpp_contribution(contribution: ContributionType):
    return contribution.value if isinstance(contribution, ContributionType) else contribution


def _cpp_basis(basis: WilsonBasis):
    return basis.value if isinstance(basis, WilsonBasis) else basis


def _cpp_sources(sources: Sequence[ParamId] | Set[ParamId]):
    return {p.to_cpp() if isinstance(p, ParamId) else p for p in sources}


def _to_cpp_scalar(value):
    """Convert float/complex/Scalar to what the binding can pass as scalar_t."""
    if isinstance(value, Scalar):
        return value._cpp_obj
    return value


class ParamSrcView:
    """Python view over the C++ ``ParamSrc`` callback object.

    Users should receive this wrapper in custom Wilson matching callbacks instead
    of the raw pybind object.
    """

    def __init__(self, cpp_obj):
        self._cpp_obj = cpp_obj

    def has(self, pid: ParamId) -> bool:
        """Return whether a parameter is available."""
        return bool(self._cpp_obj.has(pid.to_cpp() if isinstance(pid, ParamId) else pid))

    def get_val(self, pid: ParamId) -> Scalar:
        """Return a parameter value as a ``Scalar`` wrapper."""
        return Scalar.from_cpp(
            self._cpp_obj.get_val(pid.to_cpp() if isinstance(pid, ParamId) else pid)
        )

    def size(self) -> int:
        """Return the number of parameters in this view."""
        return int(self._cpp_obj.size())


class BlockSrcView:
    """Python view over the C++ ``BlockSrc`` callback object."""

    def __init__(self, cpp_obj):
        self._cpp_obj = cpp_obj

    def has_block(self, block: str) -> bool:
        """Return whether a block exists in the view."""
        return bool(self._cpp_obj.has_block(str(block)))

    def get_val(self, block: str, code) -> Scalar:
        """Return a block value as a ``Scalar`` wrapper."""
        if hasattr(code, "to_cpp"):
            code = code.to_cpp()
        return Scalar.from_cpp(self._cpp_obj.get_val(str(block), code))

    def size(self) -> int:
        """Return the number of blocks in this view."""
        return int(self._cpp_obj.size())


class CustomWilsonCoefficientConfig:
    """Python wrapper for one lambda-backed Wilson coefficient.

    Args:
        coefficient: Dynamic coefficient id, builtin enum, or name/alias.

    A matching lambda receives a C++ ``ParamSrc`` object. Use
    ``src.get_val(pid.to_cpp())`` or ``src.get_val(ParameterType.SM.value,
    "SMINPUTS", 6)`` inside the callback.
    """

    def __init__(self, coefficient: WilsonCoefLike):
        self._cpp_obj = _CppCustomWilsonCoefficientConfig(_cpp_coef_id(coefficient))

    def set_matching(
        self,
        order: QCDOrder,
        sources: Sequence[ParamId] | Set[ParamId],
        compute,
        contribution: ContributionType = ContributionType.SM,
    ) -> "CustomWilsonCoefficientConfig":
        """Attach a matching lambda for one QCD order.

        Args:
            order: QCD order of the matching contribution.
            sources: Parameter dependencies required by ``compute``.
            compute: Callable ``compute(src) -> float | complex``.
            contribution: SM/BSM component represented by this lambda.
        """

        def wrapped_compute(src):
            return _to_cpp_scalar(compute(ParamSrcView(src)))

        self._cpp_obj.set_matching(
            _cpp_order(order),
            _cpp_sources(sources),
            wrapped_compute,
            _cpp_contribution(contribution),
        )
        return self

    def to_cpp(self):
        """Return the bound C++ config object."""
        return self._cpp_obj


class CustomWilsonGroupConfig:
    """Python wrapper for a lambda-backed Wilson group.

    Args:
        group: Dynamic group id, builtin enum, or group name/alias.
        matching_scale: Matching scale in GeV.
        hadronic_scale: Hadronic/running scale in GeV.
        order: Maximum QCD order supplied by the group.
        contribution: Contribution component for the group.
    """

    def __init__(
        self,
        group: WilsonGroupLike,
        matching_scale: float = 81.0,
        hadronic_scale: float = 4.8,
        order: QCDOrder = QCDOrder.LO,
        contribution: ContributionType = ContributionType.SM,
        display_name: str = "",
    ):
        self._cpp_obj = _CppCustomWilsonGroupConfig(_cpp_group_id(group))
        self._cpp_obj.matching_scale = float(matching_scale)
        self._cpp_obj.hadronic_scale = float(hadronic_scale)
        self._cpp_obj.order = _cpp_order(order)
        self._cpp_obj.contribution = _cpp_contribution(contribution)
        self._cpp_obj.display_name = display_name

    def add_coefficient(
        self, coefficient: CustomWilsonCoefficientConfig
    ) -> "CustomWilsonGroupConfig":
        """Append a coefficient config to this group."""
        if isinstance(coefficient, CustomWilsonCoefficientConfig):
            coefficient = coefficient.to_cpp()
        self._cpp_obj.add_coefficient(coefficient)
        return self

    def set_running(
        self,
        basis: WilsonBasis,
        order: QCDOrder,
        sources: Mapping[ParameterType, Sequence[str]],
        compute,
    ) -> "CustomWilsonGroupConfig":
        """Attach a running lambda for one basis/order.

        The callable receives ``(matching, block_src)`` and must return a mapping
        ``{WCoefId_cpp: scalar}``. For simple dev use cases, keep
        ``install_identity_running_if_empty=True`` and skip this method.
        """
        cpp_sources = {
            (k.value if isinstance(k, ParameterType) else k): list(v)
            for k, v in dict(sources).items()
        }

        def wrapped_running(matching, block_src):
            py_matching = {
                QCDOrder(o): {WCoefId(str(k)): Scalar.from_cpp(v) for k, v in vals.items()}
                for o, vals in matching.items()
            }
            result = compute(py_matching, BlockSrcView(block_src))
            return {
                (k._to_cpp() if isinstance(k, WCoefId) else _cpp_coef_id(k)): _to_cpp_scalar(v)
                for k, v in dict(result).items()
            }

        self._cpp_obj.set_running(
            _cpp_basis(basis), _cpp_order(order), cpp_sources, wrapped_running
        )
        return self

    @property
    def install_identity_running_if_empty(self) -> bool:
        """Whether C++ installs identity running if no running lambda is given."""
        return bool(self._cpp_obj.install_identity_running_if_empty)

    @install_identity_running_if_empty.setter
    def install_identity_running_if_empty(self, value: bool) -> None:
        """Enable or disable automatic identity running for an empty group."""
        self._cpp_obj.install_identity_running_if_empty = bool(value)

    def to_cpp(self):
        """Return the bound C++ config object."""
        return self._cpp_obj


class WilsonInterface:
    """User-facing wrapper for the C++ Wilson-coefficient interface."""

    def __init__(self) -> None:
        """Create an unbuilt Wilson interface."""
        self._cpp_obj = _CppWilsonInterface()

    def build(self, config: WilsonBuildConfig) -> None:
        """Build the C++ Wilson pipeline from a ``WilsonBuildConfig``."""
        self._cpp_obj.build(config.to_cpp())

    def add_wilson_group(self, config: WilsonBuildConfig) -> None:
        """Add builtin or already-registered dynamic Wilson groups."""
        self._cpp_obj.add_wilson_group(config.to_cpp())

    def add_custom_group(self, config: CustomWilsonGroupConfig) -> "WilsonInterface":
        """Add a lambda-backed custom Wilson group."""
        if isinstance(config, CustomWilsonGroupConfig):
            config = config.to_cpp()
        self._cpp_obj.add_custom_group(config)
        return self

    def set_matching_scale(self, mu_W: float) -> None:
        """Set the matching scale ``mu_W``."""
        self._cpp_obj.set_matching_scale(mu_W)

    def set_hadronic_scale(self, mu_h: float) -> None:
        """Set the hadronic running scale ``mu_h``."""
        self._cpp_obj.set_hadronic_scale(mu_h)

    def _req(
        self, req_or_group, coeff=None, order=None, contribution=None, basis=None
    ) -> WilsonRequest:
        if isinstance(req_or_group, WilsonRequest):
            return req_or_group
        return WilsonRequest(
            req_or_group,
            coeff,
            order or QCDOrder.LO,
            contribution or ContributionType.TOTAL,
            wilson_basis=basis or WilsonBasis.STANDARD,
        )

    def get_M(self, req_or_group, coeff=None, order=None, contribution=None) -> Scalar:
        """Return one matching coefficient at the requested QCD order."""
        req = self._req(req_or_group, coeff, order, contribution)
        return Scalar.from_cpp(
            self._cpp_obj.get_M(
                _cpp_group_id(req.group),
                _cpp_coef_id(req.coefficient),
                _cpp_order(req.order),
                _cpp_contribution(req.contribution),
            )
        )

    def get_FM(self, req_or_group, coeff=None, order=None, contribution=None) -> Scalar:
        """Return one full matching coefficient summed up to ``order``."""
        req = self._req(req_or_group, coeff, order, contribution)
        return Scalar.from_cpp(
            self._cpp_obj.get_FM(
                _cpp_group_id(req.group),
                _cpp_coef_id(req.coefficient),
                _cpp_order(req.order),
                _cpp_contribution(req.contribution),
            )
        )

    def get_R(
        self, req_or_group, coeff=None, order=None, contribution=None, basis=WilsonBasis.STANDARD
    ) -> Scalar:
        """Return one running coefficient at the requested QCD order."""
        req = self._req(req_or_group, coeff, order, contribution, basis)
        return Scalar.from_cpp(
            self._cpp_obj.get_R(
                _cpp_group_id(req.group),
                _cpp_coef_id(req.coefficient),
                _cpp_order(req.order),
                _cpp_contribution(req.contribution),
                _cpp_basis(req.wilson_basis),
            )
        )

    def get_FR(
        self, req_or_group, coeff=None, order=None, contribution=None, basis=WilsonBasis.STANDARD
    ) -> Scalar:
        """Return one full running coefficient summed up to ``order``."""
        req = self._req(req_or_group, coeff, order, contribution, basis)
        return Scalar.from_cpp(
            self._cpp_obj.get_FR(
                _cpp_group_id(req.group),
                _cpp_coef_id(req.coefficient),
                _cpp_order(req.order),
                _cpp_contribution(req.contribution),
                _cpp_basis(req.wilson_basis),
            )
        )

    def get_sep_order_matching(
        self, group: WilsonGroupLike, coeff: WilsonCoefLike, contribution: ContributionType
    ) -> Dict[QCDOrder, Scalar]:
        """Return matching coefficients separated by QCD order."""
        cpp_map = self._cpp_obj.get_sep_order_matching_coefficient(
            _cpp_group_id(group), _cpp_coef_id(coeff), _cpp_contribution(contribution)
        )
        return {QCDOrder(order): Scalar.from_cpp(val) for order, val in cpp_map.items()}

    def get_sep_order_running(
        self,
        group: WilsonGroupLike,
        coeff: WilsonCoefLike,
        contribution: ContributionType,
        basis: WilsonBasis = WilsonBasis.STANDARD,
    ) -> Dict[QCDOrder, Scalar]:
        """Return running coefficients separated by QCD order."""
        cpp_map = self._cpp_obj.get_sep_order_run_coefficient(
            _cpp_group_id(group),
            _cpp_coef_id(coeff),
            _cpp_contribution(contribution),
            _cpp_basis(basis),
        )
        return {QCDOrder(order): Scalar.from_cpp(val) for order, val in cpp_map.items()}

    # Builtin-group map helpers remain enum-only because the C++ return key is WCoeff.
    def get_all_matching(
        self, group: WGroup, order: QCDOrder, contribution: ContributionType
    ) -> Dict[WCoeff, Scalar]:
        """Return all builtin matching coefficients in one static group."""
        cpp_map = self._cpp_obj.get_all_matching_coefficient(
            group.value, order.value, contribution.value
        )
        return {WCoeff(k): Scalar.from_cpp(v) for k, v in cpp_map.items()}

    def get_all_running(
        self,
        group: WGroup,
        order: QCDOrder,
        contribution: ContributionType,
        basis: WilsonBasis = WilsonBasis.STANDARD,
    ) -> Dict[WCoeff, Scalar]:
        """Return all builtin running coefficients in one static group."""
        cpp_map = self._cpp_obj.get_all_run_coefficient(
            group.value, order.value, contribution.value, basis.value
        )
        return {WCoeff(k): Scalar.from_cpp(v) for k, v in cpp_map.items()}


__all__ = [
    "WilsonInterface",
    "CustomWilsonCoefficientConfig",
    "CustomWilsonGroupConfig",
    "ParamSrcView",
    "BlockSrcView",
]
