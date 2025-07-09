from pyhyperiso.phyperiso.pyhyperiso.observable import ObservableInterface as _CppObservableInterface
from pyhyperiso.core.Common.GeneralEnum import Observables, Decays, QCDOrder
from pyhyperiso.core.Common.General import PyParamId
from pyhyperiso.core.Math.scalar import Scalar
from pyhyperiso.core.BusinessLogic.compounds import Estimate
from typing import Dict, Set

class PyObservableInterface:
    def __init__(self):
        self._cpp_obj = _CppObservableInterface()

    def add_observable(self, obs_name: Observables, qcd_order: QCDOrder, add_dependencies: bool = False):
        self._cpp_obj.add_observable(obs_name.value, qcd_order.value, add_dependencies)

    def add_observables(self, obs_names: Dict[Observables, QCDOrder], add_dependencies: bool = False):
        self._cpp_obj.add_observables({k.value: v.value for k, v in obs_names.items()}, add_dependencies)

    def add_observables_from_decay(self, decay: Decays, qcd_order: QCDOrder, add_dependencies: bool = False):
        self._cpp_obj.add_observables(decay.value, qcd_order.value, add_dependencies)

    def add_observable_parameter(self, obs_name: Observables, pid: PyParamId):
        self._cpp_obj.add_observable_parameter(obs_name.value, pid._cpp_obj)

    def add_observable_parameters(self, obs_name: Observables, pids: Set[PyParamId]):
        self._cpp_obj.add_observable_parameters(obs_name.value, {p._cpp_obj for p in pids})

    def compute_observable(self, obs_name: Observables):
        return Scalar.from_cpp(self._cpp_obj.compute_observable(obs_name.value))

    def compute_uncertainty(self, obs_name: Observables):
        return Scalar.from_cpp(self._cpp_obj.compute_uncertainty(obs_name.value))

    def compute_leading_uncertainties(self, obs_name: Observables, n: int):
        return {PyParamId(k): Scalar.from_cpp(v) for k, v in self._cpp_obj.compute_leading_uncertainties(obs_name.value, n).items()}

    def compute_all_uncertainties(self):
        return {Observables(k): Scalar.from_cpp(v) for k, v in self._cpp_obj.compute_all_uncertainties().items()}

    def compute_all(self):
        return {Observables(k): v for k, v in self._cpp_obj.compute_all().items()}

    def compute_chi2(self):
        return self._cpp_obj.compute_chi2()
    
    
if __name__ == "__main__":
    from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster
    from pathlib import Path
    from pyhyperiso.core.Core.Config import PyConfig, ExternalFlag
    from pyhyperiso.core.Common.GeneralEnum import Model
    from pyhyperiso.core.Core.ParamaterProvider import PyParameterProvider
    print("🔧 Initializing PyHyperisoMaster with custom PyConfig...")

    config = PyConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: True,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            ExternalFlag.USE_MARTY: False
        },
        model=Model.SM,
        mty_model_name="MSSM_UFO",
        mty_model_path=Path("/my/custom/marty/path")
    )

    print("🔧 PyConfig content:")
    print(config)

    hyp = PyHyperisoMaster()
    lha_file_path = "lha/camilia.flha" 

    print("\n🚀 Calling init with config...")
    hyp.init(lha_file=lha_file_path, config=config)
    


    interface = PyObservableInterface()
    interface.add_observable(Observables.BR_B_XS_GAMMA, QCDOrder.LO)

    print(interface.compute_observable(Observables.BR_B_XS_GAMMA))  # Scalar(...)
    
    # test_values = []
    # from pyhyperiso.core.Core.ParameterSetter import PyParameterSetter, PyParamId, ParameterType
    # py_set = PyParameterSetter()
    # for i in range(1, 81):
    #     py_set.mutate(PyParamId(ParameterType.WILSON, "B_SCALE", 1), i)
    #     test_values.append(interface.get_R(PyWilsonRequest(WGroup.B, WCoeff.C9, QCDOrder.LO, ContributionType.TOTAL)))
    
    # import matplotlib
    # matplotlib.use("TkAgg")
    # import matplotlib.pyplot as plt
    
    # plt.scatter(range(1, 81), test_values)
    
    # plt.show()