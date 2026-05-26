from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

from pyhyperiso.phyperiso.pyhyperiso import statistic as st
from pyhyperiso.core.BusinessLogic.ObservableInterface import ObservableInterface
from pyhyperiso.core.Common.BinnedObservableId import BinnedObservableId
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Statistic.ExperimentObs import ExperimentObs
from pyhyperiso.core.Statistic.GaussianSummary import GaussianSummary
from pyhyperiso.core.Statistic.MCResult import MCResult
from pyhyperiso.core.Statistic.StatisticConfig import StatisticConfig


class ProfilingMethod(Enum):
    SLICE = st.ProfilingMethod.SLICE

    def to_cpp(self):
        return self.value


class ContourAlgorithm(Enum):
    MINUIT = st.ContourAlgorithm.MINUIT

    def to_cpp(self):
        return self.value


class ProfileBackend(Enum):
    MINUIT = st.ProfileBackend.MINUIT
    LAPLACE_NUISANCE = st.ProfileBackend.LAPLACE_NUISANCE

    def to_cpp(self):
        return self.value


def _require(value, typ, name: str):
    if not isinstance(value, typ):
        raise TypeError(f"{name} doit être {typ.__name__}, reçu {type(value)!r}.")
    return value


def _cpp_param_id(pid: ParamId):
    return _require(pid, ParamId, "ParamId").to_cpp()


def _param_from_cpp(cpp_obj) -> ParamId:
    return ParamId.from_cpp(cpp_obj)


def _param_float_map_from_cpp(cpp_map) -> Dict[ParamId, float]:
    return {_param_from_cpp(k): float(v) for k, v in dict(cpp_map).items()}


def _param_float_map_to_cpp(values: Mapping[ParamId, float]):
    return {_cpp_param_id(k): float(v) for k, v in values.items()}


def _param_nested_float_map_from_cpp(cpp_map) -> Dict[ParamId, Dict[ParamId, float]]:
    return {_param_from_cpp(k): _param_float_map_from_cpp(v) for k, v in dict(cpp_map).items()}


def _experiment_float_map_from_cpp(cpp_map) -> Dict[ExperimentObs, float]:
    return {ExperimentObs.from_cpp(k): float(v) for k, v in dict(cpp_map).items()}


def _experiment_nested_float_map_from_cpp(cpp_map) -> Dict[ExperimentObs, Dict[ExperimentObs, float]]:
    return {ExperimentObs.from_cpp(k): _experiment_float_map_from_cpp(v) for k, v in dict(cpp_map).items()}


def _cpp_profiling_method(value: ProfilingMethod):
    return _require(value, ProfilingMethod, "profiling_method").to_cpp()


def _cpp_profile_backend(value: ProfileBackend):
    return _require(value, ProfileBackend, "profile_backend").to_cpp()


def _cpp_contour_algorithm(value: ContourAlgorithm):
    return _require(value, ContourAlgorithm, "contour_algorithm").to_cpp()


@dataclass(frozen=True)
class FitResultWithMaps:
    fit_ok: bool = False
    p_hat: Dict[ParamId, float] = field(default_factory=dict)
    eta_hat: Dict[ParamId, float] = field(default_factory=dict)
    p_hat_std: Dict[ParamId, float] = field(default_factory=dict)
    p_correlations: Dict[ParamId, Dict[ParamId, float]] = field(default_factory=dict)
    ell_hat: float = 0.0

    @classmethod
    def from_cpp(cls, cpp_obj) -> "FitResultWithMaps":
        return cls(
            fit_ok=bool(cpp_obj.fit_ok),
            p_hat=_param_float_map_from_cpp(cpp_obj.p_hat),
            eta_hat=_param_float_map_from_cpp(cpp_obj.eta_hat),
            p_hat_std=_param_float_map_from_cpp(cpp_obj.p_hat_std),
            p_correlations=_param_nested_float_map_from_cpp(cpp_obj.p_correlations),
            ell_hat=float(cpp_obj.ell_hat),
        )

    def to_cpp(self):
        cpp = st.FitResultWithMaps()
        cpp.fit_ok = bool(self.fit_ok)
        cpp.p_hat = _param_float_map_to_cpp(self.p_hat)
        cpp.eta_hat = _param_float_map_to_cpp(self.eta_hat)
        cpp.p_hat_std = _param_float_map_to_cpp(self.p_hat_std)
        cpp.p_correlations = {_cpp_param_id(k): _param_float_map_to_cpp(v) for k, v in self.p_correlations.items()}
        cpp.ell_hat = float(self.ell_hat)
        return cpp


@dataclass(frozen=True)
class MLFitOptions:
    run_hesse: bool = True
    request_minos: bool = False
    verbose: bool = False
    strategy: int = 0
    max_fcn: int = 0
    tolerance: float = 0.0
    allow_profile_hessian_fallback: bool = True
    profile_hessian_step_scale: float = 1.0
    profile_hessian_eig_floor_rel: float = 1e-8
    trace_first_evals: bool = False
    trace_max_evals: int = 25

    @classmethod
    def from_cpp(cls, cpp_obj) -> "MLFitOptions":
        return cls(
            run_hesse=bool(cpp_obj.run_hesse),
            request_minos=bool(cpp_obj.request_minos),
            verbose=bool(cpp_obj.verbose),
            strategy=int(cpp_obj.strategy),
            max_fcn=int(cpp_obj.max_fcn),
            tolerance=float(cpp_obj.tolerance),
            allow_profile_hessian_fallback=bool(cpp_obj.allow_profile_hessian_fallback),
            profile_hessian_step_scale=float(cpp_obj.profile_hessian_step_scale),
            profile_hessian_eig_floor_rel=float(cpp_obj.profile_hessian_eig_floor_rel),
            trace_first_evals=bool(cpp_obj.trace_first_evals),
            trace_max_evals=int(cpp_obj.trace_max_evals),
        )

    def to_cpp(self):
        cpp = st.MLFitOptions()
        cpp.run_hesse = bool(self.run_hesse)
        cpp.request_minos = bool(self.request_minos)
        cpp.verbose = bool(self.verbose)
        cpp.strategy = int(self.strategy)
        cpp.max_fcn = int(self.max_fcn)
        cpp.tolerance = float(self.tolerance)
        cpp.allow_profile_hessian_fallback = bool(self.allow_profile_hessian_fallback)
        cpp.profile_hessian_step_scale = float(self.profile_hessian_step_scale)
        cpp.profile_hessian_eig_floor_rel = float(self.profile_hessian_eig_floor_rel)
        cpp.trace_first_evals = bool(self.trace_first_evals)
        cpp.trace_max_evals = int(self.trace_max_evals)
        return cpp


@dataclass(frozen=True)
class ContourOptions:
    profiling_method: ProfilingMethod = ProfilingMethod.SLICE
    profile_backend: ProfileBackend = ProfileBackend.LAPLACE_NUISANCE
    primary_contour_method: ContourAlgorithm = ContourAlgorithm.MINUIT
    fallback_contour_method: Optional[ContourAlgorithm] = None
    resolution: int = 40

    @classmethod
    def from_cpp(cls, cpp_obj) -> "ContourOptions":
        return cls(
            profiling_method=ProfilingMethod(cpp_obj.profiling_method),
            profile_backend=ProfileBackend(cpp_obj.profile_backend),
            primary_contour_method=ContourAlgorithm(cpp_obj.primary_contour_method),
            fallback_contour_method=None if cpp_obj.fallback_contour_method is None else ContourAlgorithm(cpp_obj.fallback_contour_method),
            resolution=int(cpp_obj.resolution),
        )

    def to_cpp(self):
        cpp = st.ContourOptions()
        cpp.profiling_method = _cpp_profiling_method(self.profiling_method)
        cpp.profile_backend = _cpp_profile_backend(self.profile_backend)
        cpp.primary_contour_method = _cpp_contour_algorithm(self.primary_contour_method)
        cpp.fallback_contour_method = None if self.fallback_contour_method is None else _cpp_contour_algorithm(self.fallback_contour_method)
        cpp.resolution = int(self.resolution)
        return cpp


class Contour:
    __slots__ = ("_cpp_obj",)

    def __init__(self, cpp_obj) -> None:
        self._cpp_obj = cpp_obj

    @classmethod
    def from_cpp(cls, cpp_obj) -> "Contour":
        return cls(cpp_obj)

    def _to_cpp(self):
        return self._cpp_obj

    def __repr__(self) -> str:
        return "Contour(<opaque>)"


@dataclass(frozen=True)
class LikelihoodScanPoint:
    x: float = 0.0
    y: float = 0.0
    nll: float = 0.0
    delta_nll: float = 0.0

    @classmethod
    def from_cpp(cls, cpp_obj) -> "LikelihoodScanPoint":
        return cls(x=float(cpp_obj.x), y=float(cpp_obj.y), nll=float(cpp_obj.nll), delta_nll=float(cpp_obj.delta_nll))

    def to_cpp(self):
        cpp = st.LikelihoodScanPoint()
        cpp.x = float(self.x)
        cpp.y = float(self.y)
        cpp.nll = float(self.nll)
        cpp.delta_nll = float(self.delta_nll)
        return cpp


@dataclass(frozen=True)
class LikelihoodScanGrid:
    x_param: ParamId
    y_param: ParamId
    x_center: float = 0.0
    y_center: float = 0.0
    nx: int = 0
    ny: int = 0
    points: List[LikelihoodScanPoint] = field(default_factory=list)

    def __post_init__(self) -> None:
        _require(self.x_param, ParamId, "x_param")
        _require(self.y_param, ParamId, "y_param")

    @classmethod
    def from_cpp(cls, cpp_obj) -> "LikelihoodScanGrid":
        return cls(
            x_param=_param_from_cpp(cpp_obj.x_param),
            y_param=_param_from_cpp(cpp_obj.y_param),
            x_center=float(cpp_obj.x_center),
            y_center=float(cpp_obj.y_center),
            nx=int(cpp_obj.nx),
            ny=int(cpp_obj.ny),
            points=[LikelihoodScanPoint.from_cpp(p) for p in cpp_obj.points],
        )

    def to_cpp(self):
        cpp = st.LikelihoodScanGrid()
        cpp.x_param = _cpp_param_id(self.x_param)
        cpp.y_param = _cpp_param_id(self.y_param)
        cpp.x_center = float(self.x_center)
        cpp.y_center = float(self.y_center)
        cpp.nx = int(self.nx)
        cpp.ny = int(self.ny)
        cpp.points = [p.to_cpp() for p in self.points]
        return cpp


class StatisticInterface:
    def __init__(self, config: StatisticConfig, observable_interface: ObservableInterface) -> None:
        _require(config, StatisticConfig, "config")
        _require(observable_interface, ObservableInterface, "observable_interface")
        self.config = config
        self._cpp = st.StatisticInterface(config.to_cpp(), observable_interface._to_cpp())
        if config.selected_experiments is not None:
            self.select_experiments(config.selected_experiments)

    def _to_cpp(self):
        return self._cpp

    def select_experiment(self, experiment: str) -> None:
        self._cpp.select_experiment(str(experiment))

    def select_experiments(self, experiments: Sequence[str]) -> None:
        self._cpp.select_experiments([str(x) for x in experiments])

    def select_experiments_all(self) -> None:
        self._cpp.select_experiments_all()

    def has_experiment_selection(self) -> bool:
        return bool(self._cpp.has_experiment_selection())

    def selected_experiments(self) -> List[str]:
        return [str(x) for x in self._cpp.selected_experiments()]

    def reload_nuisance_specs(self) -> None:
        self._cpp.reload_nuisance_specs()

    def set_nuisance_user_file(self, user_yaml_path: str | Path) -> None:
        self._cpp.set_nuisance_user_file(str(user_yaml_path))

    def clear_nuisance_user_file(self) -> None:
        self._cpp.clear_nuisance_user_file()

    def update_cache(self, p_specs: Optional[Sequence[ParamId]] = None) -> None:
        p_specs = self.config.p_specs if p_specs is None else p_specs
        self._cpp.update_cache([_cpp_param_id(p) for p in p_specs])

    def compute_uncertainties(self) -> Dict[BinnedObservableId, GaussianSummary]:
        return {BinnedObservableId.from_cpp(k): GaussianSummary.from_cpp(v) for k, v in self._cpp.compute_uncertainties().items()}

    def compute_uncertainties_and_sampling(self) -> MCResult:
        return MCResult.from_cpp(self._cpp.compute_uncertainties_and_sampling())

    def compute_MLE(self, p_specs: Optional[Sequence[ParamId]] = None) -> FitResultWithMaps:
        p_specs = self.config.p_specs if p_specs is None else p_specs
        return FitResultWithMaps.from_cpp(self._cpp.compute_MLE([_cpp_param_id(p) for p in p_specs]))

    def compute_confidence_contour(
        self,
        p1: ParamId,
        p2: ParamId,
        z: float,
        bounds: Sequence[float],
        options: Optional[ContourOptions] = None,
    ) -> Contour:
        if len(bounds) != 4:
            raise ValueError("bounds doit contenir exactement 4 valeurs : [xmin, xmax, ymin, ymax].")
        cpp_options = ContourOptions().to_cpp() if options is None else _require(options, ContourOptions, "options").to_cpp()
        return Contour.from_cpp(
            self._cpp.compute_confidence_contour(
                _cpp_param_id(p1),
                _cpp_param_id(p2),
                float(z),
                [float(x) for x in bounds],
                cpp_options,
            )
        )

    def prepare_likelihood_for_scan(self, p_specs: Optional[Sequence[ParamId]] = None) -> None:
        p_specs = self.config.p_specs if p_specs is None else p_specs
        self._cpp.prepare_likelihood_for_scan([_cpp_param_id(p) for p in p_specs])

    def set_manual_scan_point(self, p_hat: Mapping[ParamId, float], eta_hat: Mapping[ParamId, float]) -> None:
        self._cpp.set_manual_scan_point(_param_float_map_to_cpp(p_hat), _param_float_map_to_cpp(eta_hat))

    def scan_likelihood_around_current_point(
        self,
        p1: ParamId,
        p2: ParamId,
        x_half_width: float,
        y_half_width: float,
        nx: int,
        ny: int,
    ) -> LikelihoodScanGrid:
        return LikelihoodScanGrid.from_cpp(
            self._cpp.scan_likelihood_around_current_point(
                _cpp_param_id(p1),
                _cpp_param_id(p2),
                float(x_half_width),
                float(y_half_width),
                int(nx),
                int(ny),
            )
        )

    def save_likelihood_scan_csv(self, path: str | Path, grid: LikelihoodScanGrid) -> None:
        self._cpp.save_likelihood_scan_csv(str(path), _require(grid, LikelihoodScanGrid, "grid").to_cpp())

    def get_all_obss_deps(self) -> Dict[ParamId, float]:
        return _param_float_map_from_cpp(self._cpp.get_all_obss_deps())

    def get_p_specs(self, p_specs: Optional[Sequence[ParamId]] = None) -> Dict[ParamId, float]:
        p_specs = self.config.p_specs if p_specs is None else p_specs
        return _param_float_map_from_cpp(self._cpp.get_p_specs([_cpp_param_id(p) for p in p_specs]))

    def get_all_correlations(self) -> Dict[ParamId, Dict[ParamId, float]]:
        return _param_nested_float_map_from_cpp(self._cpp.get_all_correlations())

    def get_all_obs_correlations(self) -> Dict[ExperimentObs, Dict[ExperimentObs, float]]:
        return _experiment_nested_float_map_from_cpp(self._cpp.get_all_obs_correlations())

    def get_obs_exp(self) -> Dict[ExperimentObs, float]:
        return _experiment_float_map_from_cpp(self._cpp.get_obs_exp())

    def print_cache(self) -> None:
        self._cpp.print_cache()


__all__ = [
    "StatisticInterface",
    "FitResultWithMaps",
    "Contour",
    "ContourOptions",
    "ProfilingMethod",
    "ContourAlgorithm",
    "ProfileBackend",
    "MLFitOptions",
    "LikelihoodScanPoint",
    "LikelihoodScanGrid",
]


if __name__ == "__main__":
    
    from pyhyperiso.core.Core.HyperisoMaster import HyperisoMaster
    from pyhyperiso.core.Core.HyperisoConfig import HyperisoConfig, ExternalFlag
    from pyhyperiso.core.Common.GeneralEnum import Model, Observables, QCDOrder, ParameterType
    from pyhyperiso.core.Common.ParamId import ParamId
    from pathlib import Path
    config = HyperisoConfig(
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
    
    hyp = HyperisoMaster()
    lha_file_path = "/home/theo/hyperiso/Assets/lha/si_input.flha" 
    
    hyp.init(lha_file=lha_file_path, config=config)
    
    obs = {Observables.BR_BS_MUMU : QCDOrder.LO, Observables.BR_BS_MUMU_UNTAG : QCDOrder.LO, Observables.BR_BD_MUMU : QCDOrder.LO}

    from pyhyperiso.core.BusinessLogic.ObservableInterface import ObservableInterface
    
    obs_int = ObservableInterface()
    
    obs_int.add_observables(obs, True)
    
    config_stat = StatisticConfig()
    config_stat.p_specs = [ParamId(ParameterType.DECAY, "B_ll", 1)]

    si = StatisticInterface(config_stat, observable_interface=obs_int)
    fit = si.compute_MLE()

    print(fit.fit_ok)
    print(fit.p_hat)
    print(fit.p_hat_std)