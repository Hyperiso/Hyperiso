from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Sequence, Tuple

from pyhyperiso.core.Statistic.StatisticConfig import StatisticConfig
from pyhyperiso.phyperiso.pyhyperiso import statistic as st

from pyhyperiso.core.Common.BinnedObservableId import BinnedObservableId
from pyhyperiso.core.Common.ParamId import ParamId

from pyhyperiso.core.Statistic.MCResult import MCResult
from pyhyperiso.core.Statistic.GaussianSummary import GaussianSummary


_CppStatisticInterface = st.StatisticInterface


def _unwrap_cpp(obj: Any) -> Any:
    if hasattr(obj, "to_cpp") and callable(getattr(obj, "to_cpp")):
        return obj.to_cpp()
    if hasattr(obj, "_cpp"):
        return obj._cpp
    if hasattr(obj, "_cpp_obj"):
        return obj._cpp_obj
    return obj


def _extract_observable_interface(obj: Any) -> Any:
    """
    Accepte soit un ObservableInterface wrapper, soit un objet plus haut niveau
    qui possède observable_interface / get_observable_interface / _observable_interface.
    """
    if obj is None:
        raise TypeError(
            "StatisticInterface attend maintenant un ObservableInterface "
            "ou un wrapper qui permet de le récupérer."
        )

    return _unwrap_cpp(obj)


def _param_from_cpp(x: Any) -> ParamId:
    if isinstance(x, ParamId):
        return x
    return ParamId.from_cpp(x)


def _param_to_cpp(x: Any) -> Any:
    if hasattr(x, "to_cpp") and callable(getattr(x, "to_cpp")):
        return x.to_cpp()
    if hasattr(x, "_cpp_obj"):
        return x._cpp_obj
    return x


def _param_float_map_from_cpp(m: Any) -> Dict[ParamId, float]:
    return {
        _param_from_cpp(k): float(v)
        for k, v in dict(m).items()
    }


def _param_nested_float_map_from_cpp(m: Any) -> Dict[ParamId, Dict[ParamId, float]]:
    return {
        _param_from_cpp(k): _param_float_map_from_cpp(v)
        for k, v in dict(m).items()
    }


@dataclass(frozen=True)
class FitResultWithMaps:
    fit_ok: bool = False
    p_hat: Dict[ParamId, float] = field(default_factory=dict)
    eta_hat: Dict[ParamId, float] = field(default_factory=dict)
    p_hat_std: Dict[ParamId, float] = field(default_factory=dict)
    p_correlations: Dict[ParamId, Dict[ParamId, float]] = field(default_factory=dict)
    ell_hat: float = 0.0

    @classmethod
    def from_cpp(cls, cpp: Any) -> "FitResultWithMaps":
        return cls(
            fit_ok=bool(cpp.fit_ok),
            p_hat=_param_float_map_from_cpp(cpp.p_hat),
            eta_hat=_param_float_map_from_cpp(cpp.eta_hat),
            p_hat_std=_param_float_map_from_cpp(cpp.p_hat_std),
            p_correlations=_param_nested_float_map_from_cpp(cpp.p_correlations),
            ell_hat=float(cpp.ell_hat),
        )

    def to_cpp(self) -> Any:
        cpp = st.FitResultWithMaps()
        cpp.fit_ok = bool(self.fit_ok)
        cpp.p_hat = {_param_to_cpp(k): float(v) for k, v in self.p_hat.items()}
        cpp.eta_hat = {_param_to_cpp(k): float(v) for k, v in self.eta_hat.items()}
        cpp.p_hat_std = {_param_to_cpp(k): float(v) for k, v in self.p_hat_std.items()}
        cpp.p_correlations = {
            _param_to_cpp(k): {_param_to_cpp(k2): float(v2) for k2, v2 in row.items()}
            for k, row in self.p_correlations.items()
        }
        cpp.ell_hat = float(self.ell_hat)
        return cpp


class ProfilingMethod(Enum):
    SLICE = "SLICE"

    def to_cpp(self) -> Any:
        return getattr(st.ProfilingMethod, self.value)


class ContourAlgorithm(Enum):
    MINUIT = "MINUIT"

    def to_cpp(self) -> Any:
        return getattr(st.ContourAlgorithm, self.value)


def _to_cpp_enum(x: Any) -> Any:
    if x is None:
        return None
    if hasattr(x, "to_cpp") and callable(getattr(x, "to_cpp")):
        return x.to_cpp()
    return x


@dataclass(frozen=True)
class ContourOptions:
    profiling_method: Any = ProfilingMethod.SLICE
    primary_contour_method: Any = ContourAlgorithm.MINUIT
    fallback_contour_method: Optional[Any] = None
    resolution: int = 100

    def to_cpp(self) -> Any:
        cpp = st.ContourOptions()
        cpp.profiling_method = _to_cpp_enum(self.profiling_method)
        cpp.primary_contour_method = _to_cpp_enum(self.primary_contour_method)
        cpp.fallback_contour_method = _to_cpp_enum(self.fallback_contour_method)
        cpp.resolution = int(self.resolution)
        return cpp


class Contour:
    """
    Wrapper minimal.

    Pour faire mieux, il faut binder les champs publics réels de `Contour`
    dans `ContourEngine.h` / `IContourExtractor.h`.
    """
    __slots__ = ("_cpp_obj",)

    def __init__(self, cpp_obj: Any):
        self._cpp_obj = cpp_obj

    @classmethod
    def from_cpp(cls, cpp_obj: Any) -> "Contour":
        return cls(cpp_obj)

    def to_cpp(self) -> Any:
        return self._cpp_obj

    def __repr__(self) -> str:
        return "Contour(<cpp>)"


class StatisticInterface:
    """
    High-level Python interface to the underlying C++ statistical engine.

    Le wrapper cache les objets pybind autant que possible.
    """

    def __init__(self, config: StatisticConfig, observable_interface: Any):
        self.config = config
        oi_cpp = _extract_observable_interface(observable_interface)
        self._cpp = _CppStatisticInterface(config.to_cpp(), oi_cpp)

        # Optionnel : ne marche que si tu exposes ces méthodes côté binding.
        if config.selected_experiments is not None and hasattr(self._cpp, "select_experiments"):
            self._cpp.select_experiments([str(x) for x in config.selected_experiments])

    def compute_uncertainties(self) -> Dict[BinnedObservableId, GaussianSummary]:
        cpp_res = self._cpp.compute_uncertainties()
        return {
            BinnedObservableId.from_cpp(x): GaussianSummary.from_cpp(y)
            for x, y in cpp_res.items()
        }

    def compute_uncertainties_and_sampling(self) -> MCResult:
        cpp_res = self._cpp.compute_uncertainties_and_sampling()
        return MCResult.from_cpp(cpp_res)

    def compute_MLE(self, p_specs: Optional[Sequence[ParamId]] = None) -> FitResultWithMaps:
        if p_specs is None:
            p_specs = self.config.p_specs

        p_specs_cpp = [_param_to_cpp(p) for p in p_specs]
        cpp_res = self._cpp.compute_MLE(p_specs_cpp)
        return FitResultWithMaps.from_cpp(cpp_res)

    def compute_confidence_contour(
        self,
        p1: ParamId,
        p2: ParamId,
        z: float,
        bounds: Sequence[float],
        options: Optional[ContourOptions] = None,
    ) -> Contour:
        if len(bounds) != 4:
            raise ValueError("bounds doit contenir exactement 4 valeurs: [xmin, xmax, ymin, ymax].")

        cpp_options = options.to_cpp() if options is not None else st.ContourOptions()

        cpp_res = self._cpp.compute_confidence_contour(
            _param_to_cpp(p1),
            _param_to_cpp(p2),
            float(z),
            [float(x) for x in bounds],
            cpp_options,
        )
        return Contour.from_cpp(cpp_res)

    def select_experiment(self, experiment: str) -> None:
        if not hasattr(self._cpp, "select_experiment"):
            raise RuntimeError("select_experiment n'est pas exposé dans le binding C++.")
        self._cpp.select_experiment(str(experiment))

    def select_experiments(self, experiments: Sequence[str]) -> None:
        if not hasattr(self._cpp, "select_experiments"):
            raise RuntimeError("select_experiments n'est pas exposé dans le binding C++.")
        self._cpp.select_experiments([str(x) for x in experiments])

    def select_experiments_all(self) -> None:
        if not hasattr(self._cpp, "select_experiments_all"):
            raise RuntimeError("select_experiments_all n'est pas exposé dans le binding C++.")
        self._cpp.select_experiments_all()

    def selected_experiments(self) -> List[str]:
        if not hasattr(self._cpp, "selected_experiments"):
            raise RuntimeError("selected_experiments n'est pas exposé dans le binding C++.")
        return [str(x) for x in self._cpp.selected_experiments()]


__all__ = [
    "StatisticInterface",
    "FitResultWithMaps",
    "Contour",
    "ContourOptions",
    "ProfilingMethod",
    "ContourAlgorithm",
]
    

if __name__ == "__main__":
    
    from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster
    from pyhyperiso.core.Core.HyperisoConfig import PyHyperisoConfig, ExternalFlag
    from pyhyperiso.core.Common.GeneralEnum import Model, Observables, QCDOrder, ParameterType
    from pyhyperiso.core.Common.ParamId import ParamId
    from pathlib import Path
    config = PyHyperisoConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: False,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            ExternalFlag.HYP_AS_SM_MARTY: True
        },
        model=Model.SM,
        mty_model_name="MSSM_UFO",
        mty_model_path=Path("/my/custom/marty/path")
    )
    
    hyp = PyHyperisoMaster()
    lha_file_path = "/home/theo/hyperiso/Assets/lha/si_input.flha" 
    
    hyp.init(lha_file=lha_file_path, config=config)
    
    obs = {Observables.BR_BS_MUMU : QCDOrder.LO, Observables.BR_BS_MUMU_UNTAG : QCDOrder.LO, Observables.BR_BD_MUMU : QCDOrder.LO}

    from pyhyperiso.core.BusinessLogic.ObservableInterface import PyObservableInterface
    
    obs_int = PyObservableInterface()
    
    obs_int.add_observables(obs, True)
    
    config_stat = StatisticConfig()
    config_stat.p_specs = [ParamId(ParameterType.DECAY, "B_ll", 1)]

    si = StatisticInterface(config_stat, observable_interface=obs_int)
    fit = si.compute_MLE()

    print(fit.fit_ok)
    print(fit.p_hat)
    print(fit.p_hat_std)