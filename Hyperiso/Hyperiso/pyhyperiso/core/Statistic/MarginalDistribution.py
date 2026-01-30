from pyhyperiso.phyperiso.pyhyperiso import statistic as st
from typing import Any, List, Optional, Sequence, Union, cast
from pyhyperiso.core.Statistic.MarginalConfig import MarginalKind, MarginalConfig, FlatMarginalConfig, GaussianMarginalConfig, LikelihoodMarginalConfig, SplitGaussianMarginalConfig, _config_from_cpp
from pyhyperiso.core.Common.General import ParamId
from pyhyperiso.core.Common.GeneralEnum import ParameterType
class MarginalDistribution:

    def __init__(self, cpp_obj: Any):
        self._cpp_obj = cpp_obj  # st.IMarginalDistribution (shared_ptr)

    @classmethod
    def from_cpp(cls, cpp_obj: Any) -> "MarginalDistribution":
        if isinstance(cpp_obj, st.GaussianMarginal):
            return GaussianMarginalDist(cpp_obj)
        if isinstance(cpp_obj, st.SplitGaussianMarginal):
            return SplitGaussianMarginalDist(cpp_obj)
        if isinstance(cpp_obj, st.FlatMarginal):
            return FlatMarginalDist(cpp_obj)
        if isinstance(cpp_obj, st.LikelihoodMarginal):
            return LikelihoodMarginalDist(cpp_obj)
        return cls(cpp_obj)

    def rvs(self, n: int) -> List[float]:
        if n < 0:
            raise ValueError("n doit être >= 0")
        return list(map(float, self._cpp_obj.rvs(int(n))))

    def logpdf(self, x: float) -> float:
        return float(self._cpp_obj.logpdf(float(x)))

    def cdf(self, x: float) -> float:
        return float(self._cpp_obj.cdf(float(x)))

    def ppf(self, p: float) -> float:
        p = float(p)
        if not (0.0 <= p <= 1.0):
            raise ValueError("p doit être dans [0, 1]")
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
    def create(
        kind: MarginalKind,
        config: MarginalConfig,
        seed: Optional[int] = None,
    ) -> MarginalDistribution:
        cpp_kind = kind.to_cpp()
        cpp_cfg = config.to_cpp()

        cpp_dist = st.DistributionFactory.create(cpp_kind, cpp_cfg, seed)
        return MarginalDistribution.from_cpp(cpp_dist)

    @staticmethod
    def gaussian(mu: float, sigma: float, seed: Optional[int] = None) -> GaussianMarginalDist:
        dist = DistributionFactoryWrapper.create(
            MarginalKind.GAUSSIAN, GaussianMarginalConfig(mu=mu, sigma=sigma), seed=seed
        )
        return cast(GaussianMarginalDist, dist)

    @staticmethod
    def split_gaussian(mu: float, sigma_p: float, sigma_m: float, seed: Optional[int] = None) -> SplitGaussianMarginalDist:
        dist = DistributionFactoryWrapper.create(
            MarginalKind.HALF_GAUSSIAN, SplitGaussianMarginalConfig(mu=mu, sigma_p=sigma_p, sigma_m=sigma_m), seed=seed
        )
        return cast(SplitGaussianMarginalDist, dist)

    @staticmethod
    def flat(a: float, b: float, seed: Optional[int] = None) -> FlatMarginalDist:
        dist = DistributionFactoryWrapper.create(
            MarginalKind.FLAT, FlatMarginalConfig(a=a, b=b), seed=seed
        )
        return cast(FlatMarginalDist, dist)

    @staticmethod
    def likelihood(
        values: Sequence[float],
        weights: Sequence[float],
        seed: Optional[int] = None,
        standardize: bool = False,
    ) -> LikelihoodMarginalDist:
        #bypass factory !! because standardize
        cfg = LikelihoodMarginalConfig(values=values, weights=weights)

        if standardize:
            if seed is None:
                seed = 0
            cpp_obj = st.LikelihoodMarginal(list(cfg.values), list(cfg.weights), int(seed), bool(standardize))
            return cast(LikelihoodMarginalDist, MarginalDistribution.from_cpp(cpp_obj))

        dist = DistributionFactoryWrapper.create(MarginalKind.LIKELIHOOD, cfg, seed=seed)
        return cast(LikelihoodMarginalDist, dist)


class MarginalConfigFactoryWrapper:
    def __init__(self) -> None:
        self._cpp_obj = st.MarginalConfigFactory()

    @staticmethod
    def _maybe_cpp_id(x: Any) -> Any:
        if hasattr(x, "to_cpp") and callable(getattr(x, "to_cpp")):
            return x.to_cpp()
        return x

    def create_from_param(self, pid: Any, kind: MarginalKind) -> MarginalConfig:
        pid_cpp = self._maybe_cpp_id(pid)
        cfg_cpp = self._cpp_obj.create(pid_cpp, kind.to_cpp())
        return _config_from_cpp(cfg_cpp)

    def create_from_observable(self, oid: Any, kind: MarginalKind) -> MarginalConfig:
        oid_cpp = self._maybe_cpp_id(oid)
        cfg_cpp = self._cpp_obj.create(oid_cpp, kind.to_cpp())
        return _config_from_cpp(cfg_cpp)


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
    
    from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster
    from pyhyperiso.core.Core.HyperisoConfig import PyHyperisoConfig, ExternalFlag
    from pyhyperiso.core.Common.GeneralEnum import Model, Observables, QCDOrder
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
    
    cfg = GaussianMarginalConfig(mu=1.2, sigma=0.3)
    dist = DistributionFactoryWrapper.create(MarginalKind.GAUSSIAN, cfg, seed=123)

    samples = dist.rvs(5)         # list[float]
    m = dist.mean()               # float
    s = dist.std()                # float

    mcf = MarginalConfigFactoryWrapper()
    pid = ParamId(ParameterType.SM, block = "MASS", code = 1)
    py_cfg = mcf.create_from_param(pid, MarginalKind.FLAT)
    dist2 = DistributionFactoryWrapper.create(MarginalKind.FLAT, py_cfg)