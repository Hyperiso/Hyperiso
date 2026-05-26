"""Python wrapper for Wilson-coefficient construction and queries.

The C++ ``WilsonInterface`` builds Wilson-coefficient groups and exposes
coefficients at the matching scale and at the hadronic, running scale. This
module provides a Python-facing API based on ``WilsonBuildConfig`` and
``WilsonRequest`` wrapper objects.

Typical workflow:
    1. Build the interface with a ``WilsonBuildConfig``.
    2. Query matching coefficients with :meth:`WilsonInterface.get_M` or
       :meth:`WilsonInterface.get_FM`.
    3. Query running coefficients with :meth:`WilsonInterface.get_R` or
       :meth:`WilsonInterface.get_FR`.

Example:
    >>> from pyhyperiso.core.Common.Configs import WilsonBuildConfig, WilsonRequest
    >>> from pyhyperiso.core.Common.GeneralEnum import WGroup, WCoeff, QCDOrder, ContributionType
    >>> wi = WilsonInterface()
    >>> wi.build(WilsonBuildConfig(groups={WGroup.B}, order=QCDOrder.NNLO))
    >>> req = WilsonRequest(WGroup.B, WCoeff.C9, QCDOrder.NNLO, ContributionType.TOTAL)
    >>> c9 = wi.get_FM(req)
"""

from __future__ import annotations

from typing import Dict

from pyhyperiso.phyperiso.pyhyperiso.wilson.wilson_interface import WilsonInterface as _CppWilsonInterface
from pyhyperiso.core.Common.GeneralEnum import QCDOrder, WCoeff, WGroup, ContributionType, WilsonBasis, Model, ParameterType
from pyhyperiso.core.Common.Configs import WilsonBuildConfig, WilsonRequest
from pyhyperiso.core.Math.scalar import Scalar


class WilsonInterface:
    """User-facing wrapper for the C++ Wilson-coefficient interface.

    The interface must be built before any coefficient query is made. The build
    step instantiates the C++ Wilson pipeline for the requested groups, QCD
    order, contribution model and scales. Additional groups can be added later
    with :meth:`add_wilson_group`.

    Coefficient naming mirrors the C++ API:
        * ``M``: matching coefficient at the matching scale, fixed QCD order;
        * ``FM``: full matching coefficient summed up to the requested order;
        * ``R``: running coefficient at the hadronic scale, fixed QCD order;
        * ``FR``: full running coefficient summed up to the requested order.

    Notes:
        When the MARTY backend is active, the underlying C++ layer may restrict
        Wilson calculations to LO QCD even if a higher order is requested.
    """

    def __init__(self) -> None:
        """Create an unbuilt Wilson interface."""
        self._cpp_obj = _CppWilsonInterface()

    def build(self, config: WilsonBuildConfig) -> None:
        """Build the C++ Wilson pipeline.

        Args:
            config: Build configuration containing Wilson groups, requested QCD
                order, matching scale, hadronic scale and backend options.

        Raises:
            TypeError: Propagated from ``config.to_cpp()`` if ``config`` is not
                a valid ``WilsonBuildConfig`` wrapper.
        """
        self._cpp_obj.build(config.to_cpp())

    def add_wilson_group(self, config: WilsonBuildConfig) -> None:
        """Add Wilson groups to an already built interface.

        Args:
            config: Build configuration whose ``groups`` field specifies the
                additional Wilson groups to add. Scale and order fields are
                forwarded to the C++ builder.

        Raises:
            RuntimeError: Propagated from C++ if the interface has not been
                built before adding groups.
        """
        self._cpp_obj.add_wilson_group(config.to_cpp())

    def set_matching_scale(self, mu_W: float) -> None:
        """Set the matching scale ``mu_W``.

        Args:
            mu_W: Matching scale used for matching coefficients.

        Notes:
            This forwards to the C++ scale setter. Dependent coefficient blocks
            are updated by the underlying C++ provider according to its caching
            rules.
        """
        self._cpp_obj.set_matching_scale(mu_W)

    def set_hadronic_scale(self, mu_h: float) -> None:
        """Set the hadronic running scale ``mu_h``.

        Args:
            mu_h: Hadronic scale used for running coefficients.
        """
        self._cpp_obj.set_hadronic_scale(mu_h)

    def get_M(self, req: WilsonRequest) -> Scalar:
        """Return one matching coefficient at the requested QCD order.

        Args:
            req: Wilson request specifying group, coefficient, QCD order and
                contribution component.

        Returns:
            Matching coefficient as a ``Scalar`` wrapper.
        """
        return Scalar.from_cpp(
            self._cpp_obj.get_M(
                req.group.value,
                req.coefficient.value,
                req.order.value,
                req.contribution.value,
            )
        )

    def get_FM(self, req: WilsonRequest) -> Scalar:
        """Return one full matching coefficient summed up to ``req.order``.

        Args:
            req: Wilson request specifying group, coefficient, QCD order and
                contribution component.

        Returns:
            Full matching coefficient as a ``Scalar`` wrapper.
        """
        return Scalar.from_cpp(
            self._cpp_obj.get_FM(
                req.group.value,
                req.coefficient.value,
                req.order.value,
                req.contribution.value,
            )
        )

    def get_R(self, req: WilsonRequest) -> Scalar:
        """Return one running coefficient at the requested QCD order.

        Args:
            req: Wilson request specifying group, coefficient, QCD order,
                contribution component and Wilson basis.

        Returns:
            Hadronic-scale running coefficient as a ``Scalar`` wrapper.
        """
        return Scalar.from_cpp(
            self._cpp_obj.get_R(
                req.group.value,
                req.coefficient.value,
                req.order.value,
                req.contribution.value,
                req.wilson_basis.value,
            )
        )

    def get_FR(self, req: WilsonRequest) -> Scalar:
        """Return one full running coefficient summed up to ``req.order``.

        Args:
            req: Wilson request specifying group, coefficient, QCD order,
                contribution component and Wilson basis.

        Returns:
            Full hadronic-scale running coefficient as a ``Scalar`` wrapper.
        """
        return Scalar.from_cpp(
            self._cpp_obj.get_FR(
                req.group.value,
                req.coefficient.value,
                req.order.value,
                req.contribution.value,
                req.wilson_basis.value,
            )
        )

    def get_sep_order_matching(
        self,
        group: WGroup,
        coeff: WCoeff,
        contribution: ContributionType,
    ) -> Dict[QCDOrder, Scalar]:
        """Return matching coefficients separated by QCD order.

        Args:
            group: Wilson group, for example ``WGroup.B``.
            coeff: Wilson coefficient identifier inside the group.
            contribution: Contribution component to query, typically ``SM``,
                ``BSM`` or ``TOTAL``.

        Returns:
            Dictionary mapping ``QCDOrder`` to matching coefficients.
        """
        cpp_map = self._cpp_obj.get_sep_order_matching_coefficient(group.value, coeff.value, contribution.value)
        return {QCDOrder(order): Scalar.from_cpp(val) for order, val in cpp_map.items()}

    def get_sep_order_running(
        self,
        group: WGroup,
        coeff: WCoeff,
        contribution: ContributionType,
        basis: WilsonBasis = WilsonBasis.STANDARD,
    ) -> Dict[QCDOrder, Scalar]:
        """Return running coefficients separated by QCD order.

        Args:
            group: Wilson group, for example ``WGroup.B``.
            coeff: Wilson coefficient identifier inside the group.
            contribution: Contribution component to query.
            basis: Operator basis used at the hadronic scale.

        Returns:
            Dictionary mapping ``QCDOrder`` to running coefficients.
        """
        cpp_map = self._cpp_obj.get_sep_order_run_coefficient(group.value, coeff.value, contribution.value, basis.value)
        return {QCDOrder(order): Scalar.from_cpp(val) for order, val in cpp_map.items()}

    def get_all_matching(
        self,
        group: WGroup,
        order: QCDOrder,
        contribution: ContributionType,
    ) -> Dict[WCoeff, Scalar]:
        """Return all matching coefficients in a Wilson group.

        Args:
            group: Wilson group to iterate over.
            order: QCD order to query.
            contribution: Contribution component to query.

        Returns:
            Dictionary mapping each coefficient in ``group`` to its matching
            coefficient.
        """
        cpp_map = self._cpp_obj.get_all_matching_coefficient(group.value, order.value, contribution.value)
        return {WCoeff(k): Scalar.from_cpp(v) for k, v in cpp_map.items()}

    def get_all_running(
        self,
        group: WGroup,
        order: QCDOrder,
        contribution: ContributionType,
        basis: WilsonBasis = WilsonBasis.STANDARD,
    ) -> Dict[WCoeff, Scalar]:
        """Return all running coefficients in a Wilson group.

        Args:
            group: Wilson group to iterate over.
            order: QCD order to query.
            contribution: Contribution component to query.
            basis: Operator basis used at the hadronic scale.

        Returns:
            Dictionary mapping each coefficient in ``group`` to its running
            coefficient.
        """
        cpp_map = self._cpp_obj.get_all_run_coefficient(group.value, order.value, contribution.value, basis.value)
        return {WCoeff(k): Scalar.from_cpp(v) for k, v in cpp_map.items()}

    def get_all_full_matching(
        self,
        group: WGroup,
        order: QCDOrder,
        contribution: ContributionType,
    ) -> Dict[WCoeff, Scalar]:
        """Return all full matching coefficients in a Wilson group.

        Args:
            group: Wilson group to iterate over.
            order: Maximum QCD order included in the perturbative sum.
            contribution: Contribution component to query.

        Returns:
            Dictionary mapping each coefficient in ``group`` to its full
            matching coefficient.
        """
        cpp_map = self._cpp_obj.get_all_full_matching_coefficient(group.value, order.value, contribution.value)
        return {WCoeff(k): Scalar.from_cpp(v) for k, v in cpp_map.items()}

    def get_all_full_running(
        self,
        group: WGroup,
        order: QCDOrder,
        contribution: ContributionType,
        basis: WilsonBasis = WilsonBasis.STANDARD,
    ) -> Dict[WCoeff, Scalar]:
        """Return all full running coefficients in a Wilson group.

        Args:
            group: Wilson group to iterate over.
            order: Maximum QCD order included in the perturbative sum.
            contribution: Contribution component to query.
            basis: Operator basis used at the hadronic scale.

        Returns:
            Dictionary mapping each coefficient in ``group`` to its full running
            coefficient.
        """
        cpp_map = self._cpp_obj.get_all_full_run_coefficient(group.value, order.value, contribution.value, basis.value)
        return {WCoeff(k): Scalar.from_cpp(v) for k, v in cpp_map.items()}


__all__ = ["WilsonInterface"]

if __name__ == "__main__":
    from pyhyperiso.core.Core.HyperisoMaster import HyperisoMaster
    from pathlib import Path
    from pyhyperiso.core.Core.HyperisoConfig import HyperisoConfig, ExternalFlag
    from pyhyperiso.core.Core.ParamaterProvider import ParameterProvider
    from pyhyperiso.core.Core.BlockProvider import BlockLogger
    
    print("🔧 Initializing PyHyperisoMaster with custom PyHyperisoConfig...")

    config = HyperisoConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: False,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            ExternalFlag.HYP_AS_SM_MARTY: True,
        },
        model=Model.SM,
        mty_model_name="MSSM_UFO",
        mty_model_path=Path("/my/custom/marty/path")
    )

    print("🔧 PyHyperisoConfig content:")
    print(config)

    hyp = HyperisoMaster()
    lha_file_path = "lha/testinput_thdm.lha" 
    # lha_file_path = "/home/theo/hyperiso/Hyperiso/Hyperiso/core/Test/InputFiles/testInput.slha"
    print("\n🚀 Calling init with config...")
    hyp.init(lha_file=lha_file_path, config=config)
    
    
    config = WilsonBuildConfig(
        groups={WGroup.B, WGroup.BScalar},
        matching_scale=81.0,
        hadronic_scale=2.0,
        order=QCDOrder.LO
    )
    # config = WilsonBuildConfig(
    #     groups={WGroup.B},
    #     matching_scale=81.0,
    #     hadronic_scale=2.0,
    #     order=QCDOrder.LO
    # )
    # config = WilsonBuildConfig(
    #     groups={},
    #     matching_scale=81.0,
    #     hadronic_scale=2.0,
    #     order=QCDOrder.LO
    # )


        
    BlockLogger().log_all_blocks(ParameterType.WILSON)
    print("trying to build wilsoninterface")
    interface = WilsonInterface()
    interface.build(config)
    print("build successful")
    # interface.set_matching_scale(81.0)

    config.groups = [WGroup.BPrime]
    interface.add_wilson_group(config)
     
    req = WilsonRequest(
        group=WGroup.B,
        coefficient=WCoeff.C9,
        order=QCDOrder.NNLO,
        contribution=ContributionType.TOTAL
    )
    coefs = {WCoeff.C1, WCoeff.C2, WCoeff.C3, WCoeff.C4, WCoeff.C5, WCoeff.C6, WCoeff.C7, WCoeff.C8, WCoeff.C9, WCoeff.C10}
    coefs_primes = {WCoeff.CP1, WCoeff.CP2, WCoeff.CP3, WCoeff.CP4, WCoeff.CP5, WCoeff.CP6, WCoeff.CP7, WCoeff.CP8, WCoeff.CP9, WCoeff.CP10, WCoeff.CPQ1, WCoeff.CPQ2}
    coefs_scalar = {WCoeff.CQ1, WCoeff.CQ2}
    for coef in coefs:
        print(coef.name, " : ", interface.get_sep_order_matching(WGroup.B, coef, ContributionType.TOTAL))
    
    print("\n\n\n")

    for coef in coefs_primes:
        print(coef.name, " : ", interface.get_sep_order_matching(WGroup.BPrime, coef, ContributionType.TOTAL))

    print("\n\n\n")

    for coef in coefs_scalar:
        print(coef.name, " : ", interface.get_sep_order_matching(WGroup.BScalar, coef, ContributionType.TOTAL))

    
    # print(interface.get_sep_order_matching(WGroup.B, WCoeff.C7, ContributionType.BSM))
    value = interface.get_M(req)
    print(value)  # Scalar(...)
    
    test_values = []
    from pyhyperiso.core.Core.ParameterSetter import ParameterSetter, ParamId, ParameterType
    py_set = ParameterSetter()
    for i in range(1, 81):
        py_set.mutate(ParamId(ParameterType.WILSON, "EW_SCALE", 1), i)
        test_values.append(interface.get_FM(WilsonRequest(WGroup.B, WCoeff.C7, QCDOrder.NNLO, ContributionType.TOTAL)))
    
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
    
    plt.scatter(range(1, 81), test_values)
    
    plt.show()