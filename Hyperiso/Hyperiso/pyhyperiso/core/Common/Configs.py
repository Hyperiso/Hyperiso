from dataclasses import dataclass, field
from typing import Iterable, Set, Union

from pyhyperiso.phyperiso.pyhyperiso.common import (
    WilsonBuildConfig as _CppWilsonBuildConfig,
    WilsonRequest as _CppWilsonRequest,
    AlphasConfig as _CppAlphasConfig,
    MassConfig as _CppMassConfig,
)
from pyhyperiso.core.Common.Mapper import GroupMapper, WCoefMapper
from pyhyperiso.core.Common.SymbolId import WGroupId, WCoefId
from pyhyperiso.core.Common.GeneralEnum import (
    QCDOrder,
    WGroup,
    WCoeff,
    ContributionType,
    ScaleType,
    MassType,
    WilsonBasis,
)

WilsonGroupLike = Union[WGroup, WGroupId, str]
WilsonCoefLike = Union[WCoeff, WCoefId, str]


def _cpp_group_id(group: WilsonGroupLike):
    """Convert a group enum/string/id to the bound C++ ``WGroupId``."""
    if isinstance(group, WGroupId):
        return group._to_cpp()
    return GroupMapper.id_of(group)._to_cpp()


def _cpp_coef_id(coef: WilsonCoefLike):
    """Convert a coefficient enum/string/id to the bound C++ ``WCoefId``."""
    if isinstance(coef, WCoefId):
        return coef._to_cpp()
    return WCoefMapper.id_of(coef)._to_cpp()


@dataclass
class WilsonBuildConfig:
    """Configuration for building Wilson coefficient groups.

    The C++ layer stores groups as dynamic ``WGroupId`` values. Python accepts
    legacy :class:`WGroup` enums, strings/aliases, or :class:`WGroupId` objects.

    Args:
        groups: Wilson groups to build. Examples: ``WGroup.B``,
            ``"BCoefficients"`` or ``GroupMapper.id_of("DEV_GROUP")``.
        matching_scale: Matching scale in GeV.
        hadronic_scale: Hadronic/running scale in GeV.
        order: Maximum QCD order requested for matching/running.
    """

    groups: Set[WilsonGroupLike] = field(default_factory=set)
    matching_scale: float = 81.0
    hadronic_scale: float = 4.8
    order: QCDOrder = QCDOrder.LO

    def to_cpp(self) -> _CppWilsonBuildConfig:
        """Convert this Python config to the bound C++ config object."""
        cpp = _CppWilsonBuildConfig()
        cpp.groups = {_cpp_group_id(g) for g in self.groups}
        cpp.matching_scale = float(self.matching_scale)
        cpp.hadronic_scale = float(self.hadronic_scale)
        cpp.order = self.order.value
        return cpp


@dataclass
class WilsonRequest:
    """Request for one Wilson coefficient.

    ``group`` and ``coefficient`` accept both legacy enums and dynamic ids. This
    means the same wrapper can query builtin Wilson coefficients and custom
    lambda-backed coefficients.
    """

    group: WilsonGroupLike
    coefficient: WilsonCoefLike
    order: QCDOrder = QCDOrder.LO
    contribution: ContributionType = ContributionType.TOTAL
    scale_type: ScaleType = ScaleType.HADRONIC
    wilson_basis: WilsonBasis = WilsonBasis.STANDARD
    sum_qcd_orders: bool = False

    def to_cpp(self) -> _CppWilsonRequest:
        """Convert this request to the bound C++ ``WilsonRequest``."""
        return _CppWilsonRequest(
            _cpp_group_id(self.group),
            _cpp_coef_id(self.coefficient),
            self.order.value,
            self.contribution.value,
            self.scale_type.value,
            bool(self.sum_qcd_orders),
        )

    def to_matching_args(self):
        """Return keyword arguments suitable for matching-scale getters."""
        return {
            "group": self.group,
            "coeff": self.coefficient,
            "order": self.order,
            "cont_type": self.contribution,
        }

    def to_running_args(self):
        """Return keyword arguments suitable for running-scale getters."""
        return {
            "group": self.group,
            "coeff": self.coefficient,
            "order": self.order,
            "cont_type": self.contribution,
            "basis": self.wilson_basis,
        }


@dataclass
class AlphasConfig:
    """Configuration for evaluating the strong coupling constant ``alpha_s``."""

    scale: float
    m_b_type: MassType = field(default=MassType.POLE)
    m_t_type: MassType = field(default=MassType.POLE)

    def to_cpp(self) -> _CppAlphasConfig:
        """Convert this Python config to the bound C++ config object."""
        return _CppAlphasConfig(self.scale, self.m_b_type.value, self.m_t_type.value)


@dataclass
class MassConfig(AlphasConfig):
    """Configuration for computing a particle mass at a given scale."""

    pdg_id: int = field(default=6)

    def to_cpp(self) -> _CppMassConfig:
        """Convert this Python config to the bound C++ config object."""
        return _CppMassConfig(self.pdg_id, self.scale, self.m_b_type.value, self.m_t_type.value)


__all__ = [
    "WilsonBuildConfig",
    "WilsonRequest",
    "AlphasConfig",
    "MassConfig",
    "WilsonGroupLike",
    "WilsonCoefLike",
    "_cpp_group_id",
    "_cpp_coef_id",
]
    
    
if __name__ == "__main__":
    py_alpha_config = AlphasConfig(scale=91.1876, m_b_type=MassType.POLE, m_t_type=MassType.MSBAR)
    py_mass_config = MassConfig(pdg_id=5, scale=91.1876, m_b_type=MassType.MSBAR, m_t_type=MassType.MSBAR)
    py_wilson = WilsonRequest(WGroup.B, WCoeff.C7, QCDOrder.LO, ContributionType.SM, ScaleType.HADRONIC)
    print(py_alpha_config)
    print(py_mass_config)
    print(py_wilson)

    cpp_alpha_config = py_alpha_config.to_cpp()
    cpp_mass_config = py_mass_config.to_cpp()
    py_wilson.to_cpp()
    