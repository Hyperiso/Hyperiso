from __future__ import annotations

from typing import List, Optional, Sequence, Union, cast

from pyhyperiso.phyperiso.pyhyperiso import statistic as st
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Statistic.ExperimentObs import ExperimentObs
from pyhyperiso.core.Statistic.MarginalConfig import (
    FlatMarginalConfig,
    GaussianMarginalConfig,
    LikelihoodMarginalConfig,
    MarginalConfig,
    MarginalKind,
    SplitGaussianMarginalConfig,
    _config_from_cpp,
)


def _require(value, typ, name: str):
    if not isinstance(value, typ):
        raise TypeError(f"{name} doit être {typ.__name__}, reçu {type(value)!r}.")
    return value


def _cpp_marginal_kind(kind: MarginalKind):
    return _require(kind, MarginalKind, "MarginalKind").to_cpp()


def _cpp_param_id(pid: ParamId):
    return _require(pid, ParamId, "ParamId").to_cpp()


def _cpp_experiment_obs(obs: ExperimentObs):
    return _require(obs, ExperimentObs, "ExperimentObs").to_cpp()


def _cpp_marginal_config(config: MarginalConfig):
    if not isinstance(config, (FlatMarginalConfig, GaussianMarginalConfig, SplitGaussianMarginalConfig, LikelihoodMarginalConfig)):
        raise TypeError(f"config marginale non supportée : {type(config)!r}.")
    return config.to_cpp()


class MarginalDistribution:
    def __init__(self, cpp_obj) -> None:
        self._cpp_obj = cpp_obj

    @classmethod
    def from_cpp(cls, cpp_obj) -> "MarginalDistribution":
        if isinstance(cpp_obj, st.GaussianMarginal):
            return GaussianMarginalDist(cpp_obj)
        if isinstance(cpp_obj, st.SplitGaussianMarginal):
            return SplitGaussianMarginalDist(cpp_obj)
        if isinstance(cpp_obj, st.FlatMarginal):
            return FlatMarginalDist(cpp_obj)
        if isinstance(cpp_obj, st.LikelihoodMarginal):
            return LikelihoodMarginalDist(cpp_obj)
        return cls(cpp_obj)

    def _to_cpp(self):
        return self._cpp_obj

    def rvs(self, n: int) -> List[float]:
        if int(n) < 0:
            raise ValueError("n doit être >= 0.")
        return [float(v) for v in self._cpp_obj.rvs(int(n))]

    def logpdf(self, x: float) -> float:
        return float(self._cpp_obj.logpdf(float(x)))

    def cdf(self, x: float) -> float:
        return float(self._cpp_obj.cdf(float(x)))

    def ppf(self, p: float) -> float:
        p = float(p)
        if not 0.0 <= p <= 1.0:
            raise ValueError("p doit être dans [0, 1].")
        return float(self._cpp_obj.ppf(p))

    def mean(self) -> float:
        return float(self._cpp_obj.mean())

    def std(self) -> float:
        return float(self._cpp_obj.std())


class GaussianMarginalDist(MarginalDistribution):
    pass


class SplitGaussianMarginalDist(MarginalDistribution):
    pass


class FlatMarginalDist(MarginalDistribution):
    pass


class LikelihoodMarginalDist(MarginalDistribution):
    pass


class DistributionFactoryWrapper:
    @staticmethod
    def create(kind: MarginalKind, config: MarginalConfig, seed: Optional[int] = None) -> MarginalDistribution:
        cpp_dist = st.MarginalFactory.create(_cpp_marginal_kind(kind), _cpp_marginal_config(config), seed)
        return MarginalDistribution.from_cpp(cpp_dist)

    @staticmethod
    def gaussian(mu: float, sigma: float, seed: Optional[int] = None) -> GaussianMarginalDist:
        return cast(GaussianMarginalDist, DistributionFactoryWrapper.create(MarginalKind.GAUSSIAN, GaussianMarginalConfig(mu=mu, sigma=sigma), seed))

    @staticmethod
    def split_gaussian(mu: float, sigma_p: float, sigma_m: float, seed: Optional[int] = None) -> SplitGaussianMarginalDist:
        return cast(
            SplitGaussianMarginalDist,
            DistributionFactoryWrapper.create(
                MarginalKind.HALF_GAUSSIAN,
                SplitGaussianMarginalConfig(mu=mu, sigma_p=sigma_p, sigma_m=sigma_m),
                seed,
            ),
        )

    @staticmethod
    def flat(a: float, b: float, seed: Optional[int] = None) -> FlatMarginalDist:
        return cast(FlatMarginalDist, DistributionFactoryWrapper.create(MarginalKind.FLAT, FlatMarginalConfig(a=a, b=b), seed))

    @staticmethod
    def likelihood(values: Sequence[float], weights: Sequence[float], seed: Optional[int] = None, standardize: bool = False) -> LikelihoodMarginalDist:
        cfg = LikelihoodMarginalConfig(values=values, weights=weights)
        if standardize:
            cpp_obj = st.LikelihoodMarginal(list(map(float, cfg.values)), list(map(float, cfg.weights)), int(0 if seed is None else seed), True)
            return cast(LikelihoodMarginalDist, MarginalDistribution.from_cpp(cpp_obj))
        return cast(LikelihoodMarginalDist, DistributionFactoryWrapper.create(MarginalKind.LIKELIHOOD, cfg, seed))


class MarginalConfigFactoryWrapper:
    def __init__(self) -> None:
        self._cpp_obj = st.MarginalConfigFactory()

    def create_from_param(self, pid: ParamId, kind: MarginalKind) -> MarginalConfig:
        return _config_from_cpp(self._cpp_obj.create(_cpp_param_id(pid), _cpp_marginal_kind(kind)))

    def create_from_observable(self, obs: ExperimentObs, kind: MarginalKind) -> MarginalConfig:
        return _config_from_cpp(self._cpp_obj.create(_cpp_experiment_obs(obs), _cpp_marginal_kind(kind)))


__all__ = [
    "MarginalKind",
    "FlatMarginalConfig",
    "GaussianMarginalConfig",
    "SplitGaussianMarginalConfig",
    "LikelihoodMarginalConfig",
    "MarginalConfig",
    "MarginalDistribution",
    "GaussianMarginalDist",
    "SplitGaussianMarginalDist",
    "FlatMarginalDist",
    "LikelihoodMarginalDist",
    "DistributionFactoryWrapper",
    "MarginalConfigFactoryWrapper",
]


if __name__ == "__main__":
    
    from pyhyperiso.core.Core.HyperisoMaster import HyperisoMaster
    from pyhyperiso.core.Core.HyperisoConfig import HyperisoConfig, ExternalFlag
    from pyhyperiso.core.Common.GeneralEnum import Model, Observables, QCDOrder, ParameterType
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
    
    cfg = GaussianMarginalConfig(mu=1.2, sigma=0.3)
    dist = DistributionFactoryWrapper.create(MarginalKind.GAUSSIAN, cfg, seed=123)

    samples = dist.rvs(5)         # list[float]
    m = dist.mean()               # float
    s = dist.std()                # float

    mcf = MarginalConfigFactoryWrapper()
    pid = ParamId(ParameterType.SM, block = "MASS", code = 1)
    py_cfg = mcf.create_from_param(pid, MarginalKind.FLAT)
    dist2 = DistributionFactoryWrapper.create(MarginalKind.FLAT, py_cfg)