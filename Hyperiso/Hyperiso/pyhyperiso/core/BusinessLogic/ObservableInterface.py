from pyhyperiso.phyperiso.pyhyperiso.observable import ObservableInterface as _CppObservableInterface
from pyhyperiso.core.Common.GeneralEnum import Observables, Decays, QCDOrder
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Math.scalar import Scalar
from pyhyperiso.core.BusinessLogic.ObservableValue import PyObservableValue, _key_to_cpp_obs_arg, ObsLike, _cpp_id_to_python_key
# from pyhyperiso.core.BusinessLogic.compounds import Estimate
from typing import Dict, Set, List

class PyObservableInterface:
    def __init__(self):
        self._cpp_obj = _CppObservableInterface()

    def add_observable(self, obs_name: Observables, qcd_order: QCDOrder, add_dependencies: bool = False):
        self._cpp_obj.add_observable(obs_name.value, qcd_order.value, add_dependencies)

    def add_observables(self, obs_names: Dict[Observables, QCDOrder], add_dependencies: bool = False):
        self._cpp_obj.add_observables({k.value: v.value for k, v in obs_names.items()}, add_dependencies)

    def add_observables_from_decay(self, decay: Decays, qcd_order: QCDOrder, add_dependencies: bool = False):
        self._cpp_obj.add_observables(decay.value, qcd_order.value, add_dependencies)

    def add_observable_parameter(self, obs_name: Observables, pid: ParamId):
        self._cpp_obj.add_observable_parameter(obs_name.value, pid._cpp_obj)

    def add_observable_parameters(self, obs_name: Observables, pids: Set[ParamId]):
        self._cpp_obj.add_observable_parameters(obs_name.value, {p._cpp_obj for p in pids})


    def compute_observable(self, obs: ObsLike) -> List[PyObservableValue]:
        print("truc : ", obs)
        print(_key_to_cpp_obs_arg(obs))
        print(type(_key_to_cpp_obs_arg(obs)))
        # cpp_vals = self._cpp_obj.compute_observable(_key_to_cpp_obs_arg(obs))  # vector<ObservableValue>
        cpp_vals = self._cpp_obj.compute_observable(obs.value)  # vector<ObservableValue>
        return [PyObservableValue.from_cpp(v) for v in cpp_vals]
    
    
    def compute_observable_central(self, obs: ObsLike) -> float:
        """Convenience: pour les non-binnés typiques (vector de taille 1)."""
        vals = self.compute_observable(obs)
        if not vals:
            raise RuntimeError("compute_observable a renvoyé une liste vide.")
        if len(vals) != 1:
            raise RuntimeError("Observable binné: utilisez compute_observable() pour obtenir toutes les bins.")
        return vals[0].value

    def compute_all(self) -> Dict[ObsLike, List[PyObservableValue]]:
        """Nécessite que tu exposes .def("compute_all", &ObservableInterface::compute_all) côté binding."""
        cpp_map = self._cpp_obj.compute_all()  # unordered_map<ObservableId, vector<ObservableValue>>
        out: Dict[ObsLike, List[PyObservableValue]] = {}
        for k, vec in cpp_map.items():
            out[_cpp_id_to_python_key(k)] = [PyObservableValue.from_cpp(v) for v in vec]
        return out

    
    
if __name__ == "__main__":
    from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster
    from pathlib import Path
    from pyhyperiso.core.Core.HyperisoConfig import PyHyperisoConfig, ExternalFlag
    from pyhyperiso.core.Common.GeneralEnum import Model
    from pyhyperiso.core.Core.ParamaterProvider import PyParameterProvider
    print("🔧 Initializing PyHyperisoMaster with custom PyHyperisoConfig...")

    config = PyHyperisoConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: True,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            ExternalFlag.HYP_AS_SM_MARTY: True,
            # ExternalFlag.USE_MARTY: False
        },
        model=Model.SM,
        mty_model_name="MSSM_UFO",
        mty_model_path=Path("/my/custom/marty/path")
    )

    print("🔧 PyHyperisoConfig content:")
    print(config)

    hyp = PyHyperisoMaster()
    lha_file_path = "lha/si_input.flha" 

    print("\n🚀 Calling init with config...")
    hyp.init(lha_file=lha_file_path, config=config)
    


    interface = PyObservableInterface()
    interface.add_observable(Observables.BR_B_XS_GAMMA, QCDOrder.NNLO, True)
    interface.add_observable(Observables.BR_BS_MUMU, QCDOrder.NNLO, True)
    interface.add_observable(Observables.BR_BD_MUMU, QCDOrder.NNLO, True)
    interface.add_observable(Observables.BR_B0__D_TAU_NU, QCDOrder.NNLO, True)
    interface.add_observable(Observables.BR_B__KSTAR_GAMMA, QCDOrder.NNLO, True)
    
    print(interface.compute_observable(Observables.BR_B_XS_GAMMA))  # Scalar(...)
    print(interface.compute_observable(Observables.BR_BS_MUMU))  # Scalar(...)
    print(interface.compute_all())
    # test_values = []
    # from pyhyperiso.core.Core.ParameterSetter import PyParameterSetter, ParamId, ParameterType
    # py_set = PyParameterSetter()
    # for i in range(1, 81):
    #     py_set.mutate(ParamId(ParameterType.WILSON, "B_SCALE", 1), i)
    #     test_values.append(interface.get_R(WilsonRequest(WGroup.B, WCoeff.C9, QCDOrder.LO, ContributionType.TOTAL)))
    
    # import matplotlib
    # matplotlib.use("TkAgg")
    # import matplotlib.pyplot as plt
    
    # plt.scatter(range(1, 81), test_values)
    
    # plt.show()