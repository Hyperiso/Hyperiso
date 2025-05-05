from dataclasses import dataclass, field
from typing import Set
from pyhyperiso.phyperiso.pyhyperiso.common import (
    WilsonBuildConfig as _CppWilsonBuildConfig,
    WilsonRequest as _CppWilsonRequest,
    AlphasConfig as _CppAlphasConfig,
    MassConfig as _CppMassConfig,
)

from pyhyperiso.core.Common.GeneralEnum import (
    QCDOrder, WGroup, WCoeff, ContributionType, ScaleType, MassType
)


@dataclass
class PyWilsonBuildConfig:
    """Python wrapper for the WilsonBuildConfig C++ struct."""
    groups: Set[WGroup] = field(default_factory=set)
    matching_scale: float = 0.0
    hadronic_scale: float = 0.0
    order: QCDOrder = QCDOrder.LO

    def to_cpp(self) -> _CppWilsonBuildConfig:
        """Converts the Python wrapper to a native C++ WilsonBuildConfig object."""
        cpp = _CppWilsonBuildConfig()
        cpp.groups = {g.value for g in self.groups}
        cpp.matching_scale = self.matching_scale
        cpp.hadronic_scale = self.hadronic_scale
        cpp.order = self.order.value
        return cpp


@dataclass
class PyWilsonRequest:
    """Python wrapper for the WilsonRequest C++ struct."""
    group: WGroup
    coefficient: WCoeff
    order: QCDOrder = QCDOrder.LO
    contribution: ContributionType = ContributionType.TOTAL
    scale_type: ScaleType = ScaleType.HADRONIC
    sum_qcd_orders: bool = False

    def to_cpp(self) -> _CppWilsonRequest:
        """Converts the Python wrapper to a native C++ WilsonRequest object."""
        cpp = _CppWilsonRequest()
        cpp.group = self.group.value
        cpp.coefficient = self.coefficient.value
        cpp.order = self.order.value
        cpp.contribution = self.contribution.value
        cpp.scale_type = self.scale_type.value
        cpp.sum_qcd_orders = self.sum_qcd_orders
        return cpp


@dataclass
class PyAlphasConfig:
    """Python wrapper for the AlphasConfig C++ struct."""
    scale: float
    m_b_type: MassType
    m_t_type: MassType

    def to_cpp(self) -> _CppAlphasConfig:
        """Converts the Python wrapper to a native C++ AlphasConfig object."""
        return _CppAlphasConfig(self.scale, self.m_b_type.value, self.m_t_type.value)


@dataclass
class PyMassConfig(PyAlphasConfig):
    """Python wrapper for the MassConfig C++ struct, extending AlphasConfig."""
    pdg_id: int

    def to_cpp(self) -> _CppMassConfig:
        """Converts the Python wrapper to a native C++ MassConfig object."""
        return _CppMassConfig(self.pdg_id, self.scale, self.m_b_type.value, self.m_t_type.value)