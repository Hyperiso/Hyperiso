from dataclasses import dataclass, field
from typing import Set
from pyhyperiso.phyperiso.pyhyperiso.common import (
    WilsonBuildConfig as _CppWilsonBuildConfig,
    WilsonRequest as _CppWilsonRequest,
    AlphasConfig as _CppAlphasConfig,
    MassConfig as _CppMassConfig,
)
from pyhyperiso.core.Common.Mapper import GroupMapper
from pyhyperiso.core.Common.GeneralEnum import (
    QCDOrder, WGroup, WCoeff, ContributionType, ScaleType, MassType, WilsonBasis
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
        cpp.groups = {GroupMapper().id_of(g) for g in self.groups}
        cpp.matching_scale = self.matching_scale
        cpp.hadronic_scale = self.hadronic_scale
        cpp.order = self.order.value
        print(cpp)
        print("cpp.order ", cpp.order)
        print("cpp.groups ", cpp.groups)
        return cpp


@dataclass
class PyWilsonRequest:
    """Python wrapper for the WilsonRequest C++ struct."""
    group: WGroup
    coefficient: WCoeff
    order: QCDOrder = QCDOrder.LO
    contribution: ContributionType = ContributionType.TOTAL
    scale_type: ScaleType = ScaleType.HADRONIC
    wilson_basis : WilsonBasis = WilsonBasis.STANDARD
    sum_qcd_orders: bool = False

    def to_cpp(self) -> _CppWilsonRequest:
        """Converts the Python wrapper to a native C++ WilsonRequest object."""
        cpp = _CppWilsonRequest(self.group.value, self.coefficient.value,
                                self.order.value, self.contribution.value,
                                self.scale_type.value, self.sum_qcd_orders)

        return cpp

    def to_matching_args(self):
        return {
            "group": self.group.value,
            "coeff": self.coefficient.value,
            "order": self.order.value,
            "cont_type": self.contribution.value
        }

    def to_running_args(self):
        return {
            "group": self.group.value,
            "coeff": self.coefficient.value,
            "order": self.order.value,
            "cont_type": self.contribution.value,
            "basis": self.scale_type.value
        }

@dataclass
class PyAlphasConfig:
    """Python wrapper for the AlphasConfig C++ struct."""
    scale: float
    m_b_type: MassType = field(default=MassType.POLE)
    m_t_type: MassType = field(default=MassType.POLE)

    def to_cpp(self) -> _CppAlphasConfig:
        """Converts the Python wrapper to a native C++ AlphasConfig object."""
        return _CppAlphasConfig(self.scale, self.m_b_type.value, self.m_t_type.value)


@dataclass
class PyMassConfig(PyAlphasConfig):
    """Python wrapper for the MassConfig C++ struct, extending AlphasConfig."""
    pdg_id: int = field(default=6)

    def to_cpp(self) -> _CppMassConfig:
        """Converts the Python wrapper to a native C++ MassConfig object."""
        return _CppMassConfig(self.pdg_id, self.scale, self.m_b_type.value, self.m_t_type.value)
    
    
if __name__ == "__main__":
    py_alpha_config = PyAlphasConfig(scale=91.1876, m_b_type=MassType.POLE, m_t_type=MassType.MSBAR)
    py_mass_config = PyMassConfig(pdg_id=5, scale=91.1876, m_b_type=MassType.MSBAR, m_t_type=MassType.MSBAR)
    py_wilson = PyWilsonRequest(WGroup.B, WCoeff.C7, QCDOrder.LO, ContributionType.SM, ScaleType.HADRONIC)
    print(py_alpha_config)
    print(py_mass_config)
    print(py_wilson)

    cpp_alpha_config = py_alpha_config.to_cpp()
    cpp_mass_config = py_mass_config.to_cpp()
    py_wilson.to_cpp()
    