from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, List, Optional, Sequence, Union, cast


from pyhyperiso.phyperiso.pyhyperiso import statistic as st


class MarginalKind(Enum):
    GAUSSIAN = "GAUSSIAN"
    HALF_GAUSSIAN = "HALF_GAUSSIAN"
    FLAT = "FLAT"
    LIKELIHOOD = "LIKELIHOOD"

    def to_cpp(self) -> Any:
        return getattr(st.MarginalType, self.value)


@dataclass(frozen=True)
class FlatMarginalConfig:
    a: float = 0.0
    b: float = 1.0

    def to_cpp(self) -> Any:
        return st.FlatMarginalCfg(self.a, self.b)

    @classmethod
    def from_cpp(cls, cpp: Any) -> "FlatMarginalConfig":
        return cls(a=float(cpp.a), b=float(cpp.b))


@dataclass(frozen=True)
class GaussianMarginalConfig:
    mu: float = 0.0
    sigma: float = 1.0

    def to_cpp(self) -> Any:
        return st.GaussianMarginalCfg(self.mu, self.sigma)

    @classmethod
    def from_cpp(cls, cpp: Any) -> "GaussianMarginalConfig":
        return cls(mu=float(cpp.mu), sigma=float(cpp.sigma))


@dataclass(frozen=True)
class SplitGaussianMarginalConfig:
    mu: float = 0.0
    sigma_p: float = 1.0
    sigma_m: float = 1.0

    def to_cpp(self) -> Any:
        return st.SplitGaussianMarginalCfg(self.mu, self.sigma_p, self.sigma_m)

    @classmethod
    def from_cpp(cls, cpp: Any) -> "SplitGaussianMarginalConfig":
        return cls(mu=float(cpp.mu), sigma_p=float(cpp.sigma_p), sigma_m=float(cpp.sigma_m))


@dataclass(frozen=True)
class LikelihoodMarginalConfig:
    values: Sequence[float]
    weights: Sequence[float]

    def __post_init__(self) -> None:
        if len(self.values) != len(self.weights):
            raise ValueError("LikelihoodMarginalConfig: values et weights doivent avoir la même taille.")

    def to_cpp(self) -> Any:
        return st.LikelihoodMarginalCfg(list(map(float, self.values)), list(map(float, self.weights)))

    @classmethod
    def from_cpp(cls, cpp: Any) -> "LikelihoodMarginalConfig":
        return cls(values=list(map(float, cpp.values)), weights=list(map(float, cpp.weights)))


MarginalConfig = Union[
    FlatMarginalConfig,
    GaussianMarginalConfig,
    SplitGaussianMarginalConfig,
    LikelihoodMarginalConfig,
]


def _config_from_cpp(cpp_cfg: Any) -> MarginalConfig:
    if isinstance(cpp_cfg, st.FlatMarginalCfg):
        return FlatMarginalConfig.from_cpp(cpp_cfg)
    if isinstance(cpp_cfg, st.GaussianMarginalCfg):
        return GaussianMarginalConfig.from_cpp(cpp_cfg)
    if isinstance(cpp_cfg, st.SplitGaussianMarginalCfg):
        return SplitGaussianMarginalConfig.from_cpp(cpp_cfg)
    if isinstance(cpp_cfg, st.LikelihoodMarginalCfg):
        return LikelihoodMarginalConfig.from_cpp(cpp_cfg)

    raise TypeError(f"Config C++ not known: {type(cpp_cfg)!r}")