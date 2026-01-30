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
class WilsonBuildConfig:
    """Configuration for building (computing) sets of Wilson coefficients.

    This is a Python wrapper around the C++ ``WilsonBuildConfig`` struct. It
    specifies:
      - which Wilson operator groups must be included,
      - the matching scale (UV -> EFT matching) in GeV,
      - the hadronic (low-energy) evaluation scale in GeV,
      - the perturbative QCD order (LO, NLO, ...).

    In C++, the group set is stored as a set of ``WGroupId``. In Python, this
    wrapper exposes a set of ``WGroup`` enums and maps them to identifiers using
    ``GroupMapper``.

    Attributes:
        groups (Set[WGroup]): Set of Wilson operator groups to include.
        matching_scale (float): Matching scale in GeV (high scale where matching
            is performed).
        hadronic_scale (float): Hadronic (low-energy) scale in GeV where Wilson
            coefficients are evaluated in the EFT.
        order (QCDOrder): Perturbative QCD order used for matching/running.

    Notes:
        - The C++ default constructor leaves scales uninitialized; this Python
          wrapper defaults them to ``0.0``.
        - ``to_cpp()`` currently prints debug information (``print(cpp)`` etc.).
          If this is not desired in production, consider removing those prints.
    """
    groups: Set[WGroup] = field(default_factory=set)
    matching_scale: float = 0.0
    hadronic_scale: float = 0.0
    order: QCDOrder = QCDOrder.LO

    def to_cpp(self) -> _CppWilsonBuildConfig:
        """Convert this Python config to the bound C++ config object.

        Returns:
            _CppWilsonBuildConfig: A native/bound C++ ``WilsonBuildConfig`` with
            mapped group identifiers and numeric enum values.
        """
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
class WilsonRequest:
    """Request for a specific Wilson coefficient.

    This is a Python wrapper around the C++ ``WilsonRequest`` struct. It models
    a single query for one Wilson coefficient:
      - which operator group (e.g. b->s sector),
      - which coefficient within that group,
      - the perturbative QCD order,
      - the contribution type (TOTAL, SM, BSM, ...),
      - the scale type at which to evaluate it (matching vs hadronic),
      - whether to sum over QCD orders,
      - and an optional Wilson basis.

    Attributes:
        group (WGroup): Operator group for the requested coefficient.
        coefficient (WCoeff): Specific coefficient within the group.
        order (QCDOrder): Perturbative QCD order for the request.
        contribution (ContributionType): Contribution to extract (TOTAL/SM/BSM...).
        scale_type (ScaleType): Scale choice (e.g. HADRONIC vs MATCHING).
        wilson_basis (WilsonBasis): Basis in which the coefficient is expressed.
        sum_qcd_orders (bool): If ``True``, contributions from multiple QCD orders
            are summed (implementation-dependent). If ``False``, only ``order`` is used.

    Notes:
        - In the C++ struct, the basis is an ``std::optional<WilsonBasis>``.
          Your Python wrapper stores ``wilson_basis`` but **does not currently
          forward it** to the C++ object in ``to_cpp()`` (the C++ constructor you
          call does not include a basis argument). If the C++ side supports basis
          elsewhere, you may want to extend the binding/constructor accordingly.
        - ``to_matching_args()`` and ``to_running_args()`` return dictionaries
          for internal plumbing. Pay attention that ``to_running_args()`` currently
          sets ``"basis"`` to ``self.scale_type.value`` (likely a bug/typo if the
          intent was to pass ``wilson_basis``).
    """
    group: WGroup
    coefficient: WCoeff
    order: QCDOrder = QCDOrder.LO
    contribution: ContributionType = ContributionType.TOTAL
    scale_type: ScaleType = ScaleType.HADRONIC
    wilson_basis : WilsonBasis = WilsonBasis.STANDARD
    sum_qcd_orders: bool = False

    def to_cpp(self) -> _CppWilsonRequest:
        """Convert this Python request to the bound C++ request object.

        Returns:
            _CppWilsonRequest: A native/bound C++ ``WilsonRequest`` created from
            the enum numeric values.

        Notes:
            The returned C++ object is constructed without an explicit Wilson
            basis (see class Notes).
        """
        cpp = _CppWilsonRequest(self.group.value, self.coefficient.value,
                                self.order.value, self.contribution.value,
                                self.scale_type.value, self.sum_qcd_orders)

        return cpp

    def to_matching_args(self):
        """Build a dict of arguments commonly used for matching computations.

        Returns:
            Dict[str, Any]: Dictionary with keys ``group``, ``coeff``, ``order``,
            and ``cont_type`` containing numeric enum values.
        """
        return {
            "group": self.group.value,
            "coeff": self.coefficient.value,
            "order": self.order.value,
            "cont_type": self.contribution.value
        }

    def to_running_args(self):
        """Build a dict of arguments commonly used for running computations.

        Returns:
            Dict[str, Any]: Dictionary with keys ``group``, ``coeff``, ``order``,
            ``cont_type``, and ``basis`` (see Notes about the current value of
            ``basis`` in this wrapper).
        """
        return {
            "group": self.group.value,
            "coeff": self.coefficient.value,
            "order": self.order.value,
            "cont_type": self.contribution.value,
            "basis": self.scale_type.value
        }

@dataclass
class AlphasConfig:
    """Configuration for evaluating the strong coupling constant α_s.

    This is a Python wrapper around the C++ ``AlphasConfig`` struct. It specifies:
      - the renormalization scale (in GeV) at which α_s is evaluated,
      - the mass schemes used for the bottom and top quarks.

    Attributes:
        scale (float): Renormalization scale (GeV) for α_s evaluation.
        m_b_type (MassType): Mass scheme for the bottom quark (e.g. POLE, MSBAR).
        m_t_type (MassType): Mass scheme for the top quark (e.g. POLE, MSBAR).
    """
    scale: float
    m_b_type: MassType = field(default=MassType.POLE)
    m_t_type: MassType = field(default=MassType.POLE)

    def to_cpp(self) -> _CppAlphasConfig:
        """Convert this Python config to the bound C++ config object.

        Returns:
            _CppAlphasConfig: A native/bound C++ ``AlphasConfig`` created from
            scalar values and numeric enum values.
        """
        return _CppAlphasConfig(self.scale, self.m_b_type.value, self.m_t_type.value)


@dataclass
class MassConfig(AlphasConfig):
    """Configuration for computing a particle mass at a given scale.

    This is a Python wrapper around the C++ ``MassConfig`` struct. It extends
    :class:`AlphasConfig` by adding the PDG identifier of the particle whose mass
    is requested.

    Attributes:
        pdg_id (int): PDG identifier of the particle (e.g. 5 for b-quark, 6 for t-quark).
        scale (float): Renormalization scale in GeV at which the mass is evaluated.
        m_b_type (MassType): Bottom-quark mass scheme used internally (passed through).
        m_t_type (MassType): Top-quark mass scheme used internally (passed through).
    """
    pdg_id: int = field(default=6)

    def to_cpp(self) -> _CppMassConfig:
        """Convert this Python config to the bound C++ config object.

        Returns:
            _CppMassConfig: A native/bound C++ ``MassConfig`` created from the PDG id,
            scale, and numeric enum values.
        """
        return _CppMassConfig(self.pdg_id, self.scale, self.m_b_type.value, self.m_t_type.value)
    
    
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
    