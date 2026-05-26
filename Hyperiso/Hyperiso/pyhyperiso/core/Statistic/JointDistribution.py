from __future__ import annotations

from typing import List, Optional, Sequence, Union

from pyhyperiso.phyperiso.pyhyperiso import statistic as st
from pyhyperiso.core.Statistic.Copula import CopulaKind
from pyhyperiso.core.Statistic.CopulaConfig import GaussianCopulaConfigPy as GaussianCopulaConfig, StudentTCopulaConfigPy as StudentTCopulaConfig
from pyhyperiso.core.Statistic.MarginalConfig import (
    FlatMarginalConfig,
    GaussianMarginalConfig,
    LikelihoodMarginalConfig,
    MarginalKind,
    SplitGaussianMarginalConfig,
)


MarginalConfig = Union[FlatMarginalConfig, GaussianMarginalConfig, SplitGaussianMarginalConfig, LikelihoodMarginalConfig]
CopulaConfig = Union[GaussianCopulaConfig, StudentTCopulaConfig]


def _require(value, typ, name: str):
    if not isinstance(value, typ):
        raise TypeError(f"{name} doit être {typ.__name__}, reçu {type(value)!r}.")
    return value


def _cpp_marginal_kind(kind: MarginalKind):
    return _require(kind, MarginalKind, "MarginalKind").to_cpp()


def _cpp_copula_kind(kind: CopulaKind):
    return _require(kind, CopulaKind, "CopulaKind").to_cpp()


def _cpp_marginal_config(cfg: MarginalConfig):
    if not isinstance(cfg, (FlatMarginalConfig, GaussianMarginalConfig, SplitGaussianMarginalConfig, LikelihoodMarginalConfig)):
        raise TypeError(f"config marginale non supportée : {type(cfg)!r}.")
    return cfg.to_cpp()


def _cpp_copula_config(cfg: CopulaConfig):
    if not isinstance(cfg, (GaussianCopulaConfig, StudentTCopulaConfig)):
        raise TypeError(f"config copule non supportée : {type(cfg)!r}.")
    return cfg.to_cpp()


class JointDistribution:
    __slots__ = ("_cpp_obj",)

    def __init__(self, cpp_obj) -> None:
        self._cpp_obj = cpp_obj

    @classmethod
    def from_cpp(cls, cpp_obj) -> "JointDistribution":
        return cls(cpp_obj)

    def _to_cpp(self):
        return self._cpp_obj

    def sample(self, n: Optional[int] = None) -> Union[List[float], List[List[float]]]:
        if n is None:
            return [float(v) for v in self._cpp_obj.sample()]
        if int(n) < 0:
            raise ValueError("n doit être >= 0.")
        return [[float(v) for v in row] for row in self._cpp_obj.sample(int(n))]

    def logpdf(self, x: Sequence[float]) -> float:
        return float(self._cpp_obj.logpdf([float(v) for v in x]))

    def dim(self) -> int:
        return int(self._cpp_obj.dim())

    @property
    def ndim(self) -> int:
        return self.dim()

    def __repr__(self) -> str:
        return f"JointDistribution(dim={self.dim()})"


class JointDistributionFactory:
    @staticmethod
    def create(
        marginal_types: Sequence[MarginalKind],
        marginal_configs: Sequence[MarginalConfig],
        copula_type: CopulaKind,
        copula_config: CopulaConfig,
        *,
        seed: Optional[int] = None,
    ) -> JointDistribution:
        if len(marginal_types) != len(marginal_configs):
            raise ValueError("marginal_types et marginal_configs doivent avoir la même taille.")
        if not marginal_types:
            raise ValueError("Au moins une marginale est requise.")

        cpp = st.JointDistribution.create(
            [_cpp_marginal_kind(t) for t in marginal_types],
            [_cpp_marginal_config(c) for c in marginal_configs],
            _cpp_copula_kind(copula_type),
            _cpp_copula_config(copula_config),
            seed if seed is not None else None,
        )
        return JointDistribution.from_cpp(cpp)

    @staticmethod
    def create_with_seeds(
        marginal_types: Sequence[MarginalKind],
        marginal_configs: Sequence[MarginalConfig],
        marginal_seeds: Sequence[int],
        copula_type: CopulaKind,
        copula_config: CopulaConfig,
        *,
        copula_seed: int,
    ) -> JointDistribution:
        if len(marginal_types) != len(marginal_configs) or len(marginal_types) != len(marginal_seeds):
            raise ValueError("marginal_types, marginal_configs et marginal_seeds doivent avoir la même taille.")
        if not marginal_types:
            raise ValueError("Au moins une marginale est requise.")

        cpp = st.JointDistribution.create_with_seeds(
            [_cpp_marginal_kind(t) for t in marginal_types],
            [_cpp_marginal_config(c) for c in marginal_configs],
            [int(s) for s in marginal_seeds],
            _cpp_copula_kind(copula_type),
            _cpp_copula_config(copula_config),
            int(copula_seed),
        )
        return JointDistribution.from_cpp(cpp)


__all__ = ["JointDistribution", "JointDistributionFactory"]


if __name__ == "__main__":
    from pyhyperiso.core.Statistic.MarginalDistribution import MarginalKind
    from pyhyperiso.core.Statistic.MarginalConfig import GaussianMarginalConfig, FlatMarginalConfig
    from pyhyperiso.core.Statistic.Copula import CopulaKind
    from pyhyperiso.core.Statistic.CopulaConfig import GaussianCopulaConfigPy
    from pyhyperiso.core.Math.real_matrix import Matrix
    m_types = [MarginalKind.GAUSSIAN, MarginalKind.FLAT]
    m_cfgs  = [GaussianMarginalConfig(mu=0.0, sigma=1.0),
            FlatMarginalConfig(a=-2.0, b=2.0)]

    R = Matrix([[1.0, 0.6],[0.6, 1.0]])
    c_cfg = GaussianCopulaConfigPy(R=R)

    jd = JointDistributionFactory.create(m_types, m_cfgs, CopulaKind.GAUSSIAN, c_cfg, seed=123)

    x = jd.sample()
    X = jd.sample(5)
    lp = jd.logpdf(x)
    print(jd, jd.ndim, lp)