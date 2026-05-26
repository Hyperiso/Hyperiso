from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Sequence

from pyhyperiso.phyperiso.pyhyperiso import statistic as _cpp_stat
from pyhyperiso.core.Common.BinnedObservableId import BinnedObservableId
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Math.real_matrix import Matrix
from pyhyperiso.core.Statistic.GaussianSummary import GaussianSummary


ObsSample = Dict[BinnedObservableId, float]
ObsSamples = List[ObsSample]
NuisanceSample = Dict[ParamId, float]
NuisanceSamples = List[NuisanceSample]


def _require(value, typ, name: str):
    if not isinstance(value, typ):
        raise TypeError(f"{name} doit être {typ.__name__}, reçu {type(value)!r}.")
    return value


def _cpp_binned_observable_id(obs: BinnedObservableId):
    return _require(obs, BinnedObservableId, "BinnedObservableId").to_cpp()


def _cpp_param_id(pid: ParamId):
    return _require(pid, ParamId, "ParamId").to_cpp()


def _obs_sample_from_cpp(row) -> ObsSample:
    return {BinnedObservableId.from_cpp(k): float(v) for k, v in dict(row).items()}


def _obs_sample_to_cpp(row: Mapping[BinnedObservableId, float]):
    return {_cpp_binned_observable_id(k): float(v) for k, v in row.items()}


def _nuisance_sample_from_cpp(row) -> NuisanceSample:
    return {ParamId.from_cpp(k): float(v) for k, v in dict(row).items()}


def _nuisance_sample_to_cpp(row: Mapping[ParamId, float]):
    return {_cpp_param_id(k): float(v) for k, v in row.items()}


@dataclass
class MCConfig:
    draws: int = 10000
    skew_abs_threshold: float = 0.2
    covariance_ridge_rel: float = 1e-8
    covariance_ridge_abs: float = 1e-12
    retry_failed_predictions: bool = True
    max_prediction_failures: int = 20000

    @classmethod
    def from_cpp(cls, cpp_obj) -> "MCConfig":
        return cls(
            draws=int(cpp_obj.draws),
            skew_abs_threshold=float(cpp_obj.skew_abs_threshold),
            covariance_ridge_rel=float(cpp_obj.covariance_ridge_rel),
            covariance_ridge_abs=float(cpp_obj.covariance_ridge_abs),
            retry_failed_predictions=bool(cpp_obj.retry_failed_predictions),
            max_prediction_failures=int(cpp_obj.max_prediction_failures),
        )

    def to_cpp(self):
        cpp = _cpp_stat.MCConfig()
        cpp.draws = int(self.draws)
        cpp.skew_abs_threshold = float(self.skew_abs_threshold)
        cpp.covariance_ridge_rel = float(self.covariance_ridge_rel)
        cpp.covariance_ridge_abs = float(self.covariance_ridge_abs)
        cpp.retry_failed_predictions = bool(self.retry_failed_predictions)
        cpp.max_prediction_failures = int(self.max_prediction_failures)
        return cpp


@dataclass
class MCObservableCovariance:
    ids: List[BinnedObservableId] = field(default_factory=list)
    mean: List[float] = field(default_factory=list)
    covariance: Matrix | None = None
    covariance_inv: Matrix | None = None

    @classmethod
    def from_cpp(cls, cpp_obj) -> "MCObservableCovariance":
        return cls(
            ids=[BinnedObservableId.from_cpp(x) for x in cpp_obj.ids],
            mean=[float(x) for x in cpp_obj.mean],
            covariance=Matrix.from_cpp(cpp_obj.covariance),
            covariance_inv=Matrix.from_cpp(cpp_obj.covariance_inv),
        )

    def to_cpp(self):
        cpp = _cpp_stat.MCObservableCovariance()
        cpp.ids = [_cpp_binned_observable_id(x) for x in self.ids]
        cpp.mean = [float(x) for x in self.mean]
        if self.covariance is not None:
            cpp.covariance = self.covariance.to_cpp()
        if self.covariance_inv is not None:
            cpp.covariance_inv = self.covariance_inv.to_cpp()
        return cpp


@dataclass
class MCRealization:
    sampled_obss: ObsSamples = field(default_factory=list)
    sampled_params: NuisanceSamples = field(default_factory=list)

    @classmethod
    def from_cpp(cls, cpp_obj) -> "MCRealization":
        return cls(
            sampled_obss=[_obs_sample_from_cpp(row) for row in cpp_obj.sampled_obss],
            sampled_params=[_nuisance_sample_from_cpp(row) for row in cpp_obj.sampled_params],
        )

    def to_cpp(self):
        cpp = _cpp_stat.MCRealization()
        cpp.sampled_obss = [_obs_sample_to_cpp(row) for row in self.sampled_obss]
        cpp.sampled_params = [_nuisance_sample_to_cpp(row) for row in self.sampled_params]
        return cpp


@dataclass
class MCResult:
    mc_real: MCRealization = field(default_factory=MCRealization)
    summary: List[GaussianSummary] = field(default_factory=list)
    covariance: MCObservableCovariance = field(default_factory=MCObservableCovariance)

    @classmethod
    def from_cpp(cls, cpp_obj) -> "MCResult":
        return cls(
            mc_real=MCRealization.from_cpp(cpp_obj.mc_real),
            summary=[GaussianSummary.from_cpp(gs) for gs in cpp_obj.summary],
            covariance=MCObservableCovariance.from_cpp(cpp_obj.covariance),
        )

    def to_cpp(self):
        cpp = _cpp_stat.MCResult()
        cpp.mc_real = self.mc_real.to_cpp()
        cpp.summary = [gs.to_cpp() for gs in self.summary]
        cpp.covariance = self.covariance.to_cpp()
        return cpp


def covariance_ids_from_first_sample(samples: Sequence[Mapping[BinnedObservableId, float]]) -> List[BinnedObservableId]:
    cpp_samples = [_obs_sample_to_cpp(row) for row in samples]
    return [BinnedObservableId.from_cpp(x) for x in _cpp_stat.covariance_ids_from_first_sample(cpp_samples)]


def covariance_from_obs_samples(
    samples: Sequence[Mapping[BinnedObservableId, float]],
    ids: Sequence[BinnedObservableId],
    ridge_rel: float = 1e-8,
    ridge_abs: float = 1e-12,
) -> MCObservableCovariance:
    cpp_samples = [_obs_sample_to_cpp(row) for row in samples]
    cpp_ids = [_cpp_binned_observable_id(x) for x in ids]
    return MCObservableCovariance.from_cpp(
        _cpp_stat.covariance_from_obs_samples(cpp_samples, cpp_ids, float(ridge_rel), float(ridge_abs))
    )


__all__ = [
    "ObsSample",
    "ObsSamples",
    "NuisanceSample",
    "NuisanceSamples",
    "MCConfig",
    "MCObservableCovariance",
    "MCRealization",
    "MCResult",
    "covariance_ids_from_first_sample",
    "covariance_from_obs_samples",
]
