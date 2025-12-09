from pyhyperiso.phyperiso.pyhyperiso.wilson.wilson_interface import WilsonInterface as _CppWilsonInterface
from pyhyperiso.core.Common.General import PyParamId
from pyhyperiso.core.Common.GeneralEnum import Model, QCDOrder, WCoeff, WGroup, ContributionType, WilsonBasis
from pyhyperiso.core.Common.Configs import PyWilsonBuildConfig, PyWilsonRequest
from pyhyperiso.core.Math.scalar import Scalar

class PyWilsonInterface:
    """User-facing Python wrapper for the C++ WilsonInterface class."""

    def __init__(self):
        self._cpp_obj = _CppWilsonInterface()

    def build(self, config: PyWilsonBuildConfig):
        """Initializes the WilsonInterface with a build config."""
        self._cpp_obj.build(config.to_cpp())

    def add_wilson_group(self, config: PyWilsonBuildConfig):
        self._cpp_obj.add_wilson_group(config.to_cpp())
        
    def set_matching_scale(self, mu_W: float):
        """Sets the matching scale (μ_W)."""
        self._cpp_obj.set_matching_scale(mu_W)

    def set_hadronic_scale(self, mu_h: float):
        """Sets the hadronic scale (μ_h)."""
        self._cpp_obj.set_hadronic_scale(mu_h)

    def get_M(self, req: PyWilsonRequest) -> Scalar:
        """Gets the matching coefficient (alias: getM)."""
        return Scalar.from_cpp(
            self._cpp_obj.get_M(
                req.group.value,
                req.coefficient.value,
                req.order.value,
                req.contribution.value
            )
        )

    def get_FM(self, req: PyWilsonRequest) -> Scalar:
        """Gets the full matching coefficient (alias: getFM)."""
        return Scalar.from_cpp(
            self._cpp_obj.get_FM(
                req.group.value,
                req.coefficient.value,
                req.order.value,
                req.contribution.value
            )
        )

    def get_R(self, req: PyWilsonRequest) -> Scalar:
        """Gets the running coefficient (alias: getR)."""
        return Scalar.from_cpp(
            self._cpp_obj.get_R(
                req.group.value,
                req.coefficient.value,
                req.order.value,
                req.contribution.value,
                req.wilson_basis.value
            )
        )

    def get_FR(self, req: PyWilsonRequest) -> Scalar:
        """Gets the full running coefficient (alias: getFR)."""
        return Scalar.from_cpp(
            self._cpp_obj.get_FR(
                req.group.value,
                req.coefficient.value,
                req.order.value,
                req.contribution.value,
                req.wilson_basis.value
            )
        )

    def get_sep_order_matching(self, group: WGroup, coeff: WCoeff, contribution: ContributionType) -> dict:
        """Returns {QCDOrder: Scalar} for separate matching orders."""
        cpp_map = self._cpp_obj.get_sep_order_matching_coefficient(group.value, coeff.value, contribution.value)
        return {QCDOrder(order): Scalar.from_cpp(val) for order, val in cpp_map.items()}

    def get_sep_order_running(self, group: WGroup, coeff: WCoeff, contribution: ContributionType, basis: WilsonBasis = WilsonBasis.STANDARD) -> dict:
        """Returns {QCDOrder: Scalar} for separate running orders."""
        cpp_map = self._cpp_obj.get_sep_order_run_coefficient(group.value, coeff.value, contribution.value, basis.value)
        return {QCDOrder(order): Scalar.from_cpp(val) for order, val in cpp_map.items()}

    def get_all_matching(self, group: WGroup, order: QCDOrder, contribution: ContributionType) -> dict:
        """Returns {WCoeff: Scalar} for all matching coefficients."""
        cpp_map = self._cpp_obj.get_all_matching_coefficient(group.value, order.value, contribution.value)
        return {WCoeff(k): Scalar.from_cpp(v) for k, v in cpp_map.items()}

    def get_all_running(self, group: WGroup, order: QCDOrder, contribution: ContributionType, basis: WilsonBasis = WilsonBasis.STANDARD) -> dict:
        """Returns {WCoeff: Scalar} for all running coefficients."""
        cpp_map = self._cpp_obj.get_all_run_coefficient(group.value, order.value, contribution.value, basis.value)
        return {WCoeff(k): Scalar.from_cpp(v) for k, v in cpp_map.items()}

    def get_all_full_matching(self, group: WGroup, order: QCDOrder, contribution: ContributionType) -> dict:
        """Returns {WCoeff: Scalar} for full matching coefficients summed over QCD orders."""
        cpp_map = self._cpp_obj.get_all_full_matching_coefficient(group.value, order.value, contribution.value)
        return {WCoeff(k): Scalar.from_cpp(v) for k, v in cpp_map.items()}

    def get_all_full_running(self, group: WGroup, order: QCDOrder, contribution: ContributionType, basis: WilsonBasis = WilsonBasis.STANDARD) -> dict:
        """Returns {WCoeff: Scalar} for full running coefficients summed over QCD orders."""
        cpp_map = self._cpp_obj.get_all_full_run_coefficient(group.value, order.value, contribution.value, basis.value)
        return {WCoeff(k): Scalar.from_cpp(v) for k, v in cpp_map.items()}


if __name__ == "__main__":
    from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster
    from pathlib import Path
    from pyhyperiso.core.Core.Config import PyConfig, ExternalFlag
    from pyhyperiso.core.Core.ParamaterProvider import PyParameterProvider
    print("🔧 Initializing PyHyperisoMaster with custom PyConfig...")

    config = PyConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: False,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            # ExternalFlag.USE_MARTY: False
        },
        model=Model.SM,
        mty_model_name="MSSM_UFO",
        mty_model_path=Path("/my/custom/marty/path")
    )

    print("🔧 PyConfig content:")
    print(config)

    hyp = PyHyperisoMaster()
    lha_file_path = "lha/testinput_thdm.lha" 
    # lha_file_path = "/home/theo/hyperiso/Hyperiso/Hyperiso/core/Test/InputFiles/testInput.slha"
    print("\n🚀 Calling init with config...")
    hyp.init(lha_file=lha_file_path, config=config)
    
    
    config = PyWilsonBuildConfig(
        groups={WGroup.B, WGroup.BScalar},
        matching_scale=81.0,
        hadronic_scale=2.0,
        order=QCDOrder.NNLO
    )
    print("trying to build wilsoninterface")
    interface = PyWilsonInterface()
    interface.build(config)
    print("build successful")
    # interface.set_matching_scale(81.0)

    config.groups = [WGroup.BPrime]
    interface.add_wilson_group(config)
     
    req = PyWilsonRequest(
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
    
    # test_values = []
    # from pyhyperiso.core.Core.ParameterSetter import PyParameterSetter, PyParamId, ParameterType
    # py_set = PyParameterSetter()
    # for i in range(1, 81):
    #     py_set.mutate(PyParamId(ParameterType.WILSON, "B_SCALE", 1), i)
    #     test_values.append(interface.get_FR(PyWilsonRequest(WGroup.B, WCoeff.C9, QCDOrder.NNLO, ContributionType.TOTAL)))
    
    # import matplotlib
    # matplotlib.use("TkAgg")
    # import matplotlib.pyplot as plt
    
    # plt.scatter(range(1, 81), test_values)
    
    # plt.show()