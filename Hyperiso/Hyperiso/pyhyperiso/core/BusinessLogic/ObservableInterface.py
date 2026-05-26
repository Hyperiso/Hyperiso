from __future__ import annotations

from typing import Dict, List, Mapping, Sequence, Set

from pyhyperiso.phyperiso.pyhyperiso.observable import ObservableInterface as _CppObservableInterface
from pyhyperiso.core.BusinessLogic.ObservableValue import ObservableValue
from pyhyperiso.core.Common.BinnedObservableId import BinnedObservableId
from pyhyperiso.core.Common.GeneralEnum import Decays, Observables, ParameterType, QCDOrder, UncertaintyType
from pyhyperiso.core.Common.LhaID import LhaID
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Common.SymbolId import ObservableId


def _require(value, typ, name: str):
    if not isinstance(value, typ):
        raise TypeError(f"{name} doit être {typ.__name__}, reçu {type(value)!r}.")
    return value


def _cpp_observable_enum(obs: Observables):
    return _require(obs, Observables, "obs").value


def _cpp_observable_id(obs: ObservableId):
    return _require(obs, ObservableId, "obs")._to_cpp()


def _cpp_binned_observable_id(obs: BinnedObservableId):
    return _require(obs, BinnedObservableId, "obs").to_cpp()


def _cpp_qcd_order(order: QCDOrder):
    return _require(order, QCDOrder, "qcd_order").value


def _cpp_decay(decay: Decays):
    return _require(decay, Decays, "decay").value


def _cpp_uncertainty_type(u_type: UncertaintyType):
    return _require(u_type, UncertaintyType, "u_type").value


def _cpp_param_id(pid: ParamId):
    return _require(pid, ParamId, "pid").to_cpp()


def _param_from_cpp(cpp_obj) -> ParamId:
    return ParamId.from_cpp(cpp_obj)


def _single_lha_code(code: LhaID) -> int:
    _require(code, LhaID, "pid.code")
    parts = code.get_parts()
    if len(parts) != 1:
        raise ValueError("set_param/get_param via ObservableInterface exige un LhaID à une seule entrée.")
    return int(parts[0])


class ObservableInterface:
    """API observable côté Python.

    Les entrées publiques sont exclusivement les wrappers Python du package.
    Les conversions vers pybind11 restent confinées dans cette classe.
    """

    def __init__(self) -> None:
        self._cpp_obj = _CppObservableInterface()

    @classmethod
    def from_cpp(cls, cpp_obj) -> "ObservableInterface":
        inst = cls.__new__(cls)
        inst._cpp_obj = cpp_obj
        return inst

    def _to_cpp(self):
        return self._cpp_obj

    def add_observable(
        self,
        obs: Observables,
        qcd_order: QCDOrder,
        add_dependencies: bool = False,
    ) -> "ObservableInterface":
        self._cpp_obj.add_observable(_cpp_observable_enum(obs), _cpp_qcd_order(qcd_order), bool(add_dependencies))
        return self

    def add_observable_id(
        self,
        obs: ObservableId,
        qcd_order: QCDOrder,
        add_dependencies: bool = False,
    ) -> "ObservableInterface":
        self._cpp_obj.add_observable(_cpp_observable_id(obs), _cpp_qcd_order(qcd_order), bool(add_dependencies))
        return self

    def add_binned_observable(
        self,
        obs: BinnedObservableId,
        qcd_order: QCDOrder,
        add_dependencies: bool = False,
    ) -> "ObservableInterface":
        self._cpp_obj.add_observable(_cpp_binned_observable_id(obs), _cpp_qcd_order(qcd_order), bool(add_dependencies))
        return self

    def add_observables(
        self,
        obs_names: Mapping[Observables, QCDOrder],
        add_dependencies: bool = False,
    ) -> "ObservableInterface":
        cpp_map = {_cpp_observable_enum(k): _cpp_qcd_order(v) for k, v in obs_names.items()}
        self._cpp_obj.add_observables(cpp_map, bool(add_dependencies))
        return self

    def add_observable_ids(
        self,
        obs_names: Mapping[ObservableId, QCDOrder],
        add_dependencies: bool = False,
    ) -> "ObservableInterface":
        cpp_map = {_cpp_observable_id(k): _cpp_qcd_order(v) for k, v in obs_names.items()}
        self._cpp_obj.add_observables(cpp_map, bool(add_dependencies))
        return self

    def add_observables_from_decay(
        self,
        decay: Decays,
        qcd_order: QCDOrder,
        add_dependencies: bool = False,
    ) -> "ObservableInterface":
        self._cpp_obj.add_observables(_cpp_decay(decay), _cpp_qcd_order(qcd_order), bool(add_dependencies))
        return self

    def add_observable_parameter(self, obs: Observables, pid: ParamId) -> "ObservableInterface":
        self._cpp_obj.add_observable_parameter(_cpp_observable_enum(obs), _cpp_param_id(pid))
        return self

    def add_observable_id_parameter(self, obs: ObservableId, pid: ParamId) -> "ObservableInterface":
        self._cpp_obj.add_observable_parameter(_cpp_observable_id(obs), _cpp_param_id(pid))
        return self

    def add_observable_parameters(self, obs: Observables, pids: Set[ParamId]) -> "ObservableInterface":
        self._cpp_obj.add_observable_parameters(_cpp_observable_enum(obs), {_cpp_param_id(p) for p in pids})
        return self

    def add_observable_id_parameters(self, obs: ObservableId, pids: Set[ParamId]) -> "ObservableInterface":
        self._cpp_obj.add_observable_parameters(_cpp_observable_id(obs), {_cpp_param_id(p) for p in pids})
        return self

    def remove_observable(self, obs: Observables) -> "ObservableInterface":
        self._cpp_obj.remove_observable(_cpp_observable_enum(obs))
        return self

    def remove_observable_id(self, obs: ObservableId) -> "ObservableInterface":
        self._cpp_obj.remove_observable(_cpp_observable_id(obs))
        return self

    def remove_observables(self, ids: Set[Observables]) -> "ObservableInterface":
        self._cpp_obj.remove_observables({_cpp_observable_enum(x) for x in ids})
        return self

    def remove_observables_from_decay(self, decay: Decays) -> "ObservableInterface":
        self._cpp_obj.remove_observables(_cpp_decay(decay))
        return self

    def compute_observable(self, obs: Observables) -> List[ObservableValue]:
        return [ObservableValue.from_cpp(v) for v in self._cpp_obj.compute_observable(_cpp_observable_enum(obs))]

    def compute_observable_id(self, obs: ObservableId) -> List[ObservableValue]:
        return [ObservableValue.from_cpp(v) for v in self._cpp_obj.compute_observable(_cpp_observable_id(obs))]

    def compute_observable_central(self, obs: Observables) -> float:
        vals = self.compute_observable(obs)
        if not vals:
            raise RuntimeError("compute_observable a renvoyé une liste vide.")
        if len(vals) != 1:
            raise RuntimeError("Observable binné : utilisez compute_observable() pour récupérer toutes les bins.")
        return vals[0].value

    def compute_observable_id_central(self, obs: ObservableId) -> float:
        vals = self.compute_observable_id(obs)
        if not vals:
            raise RuntimeError("compute_observable_id a renvoyé une liste vide.")
        if len(vals) != 1:
            raise RuntimeError("Observable binné : utilisez compute_observable_id() pour récupérer toutes les bins.")
        return vals[0].value

    def compute_all(self) -> Dict[ObservableId, List[ObservableValue]]:
        cpp_map = self._cpp_obj.compute_all()
        return {
            ObservableId(str(cpp_id)): [ObservableValue.from_cpp(v) for v in values]
            for cpp_id, values in cpp_map.items()
        }

    def get_exp_value(self, obs: Observables) -> float:
        return float(self._cpp_obj.get_exp_value(_cpp_observable_enum(obs)))

    def get_exp_value_id(self, obs: ObservableId) -> float:
        return float(self._cpp_obj.get_exp_value(_cpp_observable_id(obs)))

    def get_exp_uncertainty(
        self,
        obs: Observables,
        u_type: UncertaintyType = UncertaintyType.COMBINED,
    ) -> float:
        return float(self._cpp_obj.get_exp_uncertainty(_cpp_observable_enum(obs), _cpp_uncertainty_type(u_type)))

    def get_exp_uncertainty_id(
        self,
        obs: ObservableId,
        u_type: UncertaintyType = UncertaintyType.COMBINED,
    ) -> float:
        return float(self._cpp_obj.get_exp_uncertainty(_cpp_observable_id(obs), _cpp_uncertainty_type(u_type)))

    def get_current_observables(self) -> List[BinnedObservableId]:
        return [BinnedObservableId.from_cpp(x) for x in self._cpp_obj.get_current_observables()]

    def get_all_ops_deps(self, obs: Observables) -> Set[ParamId]:
        return {_param_from_cpp(x) for x in self._cpp_obj.get_all_ops_deps(_cpp_observable_enum(obs))}

    def get_all_ops_deps_id(self, obs: ObservableId) -> Set[ParamId]:
        return {_param_from_cpp(x) for x in self._cpp_obj.get_all_ops_deps(_cpp_observable_id(obs))}

    def set_param(self, pid: ParamId, value: float) -> None:
        _require(pid, ParamId, "pid")
        if pid.type is None:
            raise ValueError("pid.type doit être défini pour set_param().")
        _require(pid.type, ParameterType, "pid.type")
        self._cpp_obj.set_param(str(pid.block), _single_lha_code(pid.code), float(value), pid.type.value)

    def get_param(self, pid: ParamId):
        _require(pid, ParamId, "pid")
        if pid.type is None:
            raise ValueError("pid.type doit être défini pour get_param().")
        _require(pid.type, ParameterType, "pid.type")
        return self._cpp_obj.get_param(str(pid.block), _single_lha_code(pid.code), pid.type.value)

    def reload_params(self) -> None:
        self._cpp_obj.reload_params()

    def enable_obs(self) -> None:
        self._cpp_obj.enable_obs()

    def set_bkstarll_threads(self, n_threads: int) -> None:
        self._cpp_obj.set_bkstarll_threads(int(n_threads))

    def set_bkll_threads(self, n_threads: int) -> None:
        self._cpp_obj.set_bkll_threads(int(n_threads))

    def set_bsphi_threads(self, n_threads: int) -> None:
        self._cpp_obj.set_bsphi_threads(int(n_threads))


__all__ = ["ObservableInterface"]

    
if __name__ == "__main__":
    from pyhyperiso.core.Core.HyperisoMaster import HyperisoMaster
    from pathlib import Path
    from pyhyperiso.core.Core.HyperisoConfig import HyperisoConfig, ExternalFlag
    from pyhyperiso.core.Common.GeneralEnum import Model
    from pyhyperiso.core.Core.ParamaterProvider import ParameterProvider
    print("🔧 Initializing PyHyperisoMaster with custom PyHyperisoConfig...")

    config = HyperisoConfig(
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

    hyp = HyperisoMaster()
    lha_file_path = "lha/si_input.flha" 

    print("\n🚀 Calling init with config...")
    hyp.init(lha_file=lha_file_path, config=config)
    


    interface = ObservableInterface()
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