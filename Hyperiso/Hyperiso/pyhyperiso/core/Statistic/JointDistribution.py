from __future__ import annotations

from dataclasses import dataclass
from typing import Any, List, Optional, Sequence, Tuple, Union, cast

from pyhyperiso.phyperiso.pyhyperiso import statistic as st


from pyhyperiso.core.Statistic.MarginalDistribution import MarginalKind
from pyhyperiso.core.Statistic.MarginalConfig import (
        FlatMarginalConfig,
        GaussianMarginalConfig,
        SplitGaussianMarginalConfig,
        LikelihoodMarginalConfig,
    )

MarginalConfig = Union[
    FlatMarginalConfig,
    GaussianMarginalConfig,
    SplitGaussianMarginalConfig,
    LikelihoodMarginalConfig,
]


from pyhyperiso.core.Statistic.Copula import CopulaKind  

from pyhyperiso.core.Statistic.CopulaConfig import (
        GaussianCopulaConfigPy as GaussianCopulaConfig,
        StudentTCopulaConfigPy as StudentTCopulaConfig,
    )

CopulaConfig = Union[GaussianCopulaConfig, StudentTCopulaConfig]


def _as_cpp_enum_marginal(x: Any) -> Any:
    """Accept MarginalKind wrapper OR C++ enum st.MarginalType."""
    if hasattr(x, "to_cpp") and callable(getattr(x, "to_cpp")):
        return x.to_cpp()
    return x


def _as_cpp_enum_copula(x: Any) -> Any:
    """Accept CopulaKind wrapper OR C++ enum st.CopulaType."""
    if hasattr(x, "to_cpp") and callable(getattr(x, "to_cpp")):
        return x.to_cpp()
    return x


def _as_cpp_config(x: Any) -> Any:
    """Accept config wrapper with to_cpp() OR already a C++ config object."""
    if hasattr(x, "to_cpp") and callable(getattr(x, "to_cpp")):
        return x.to_cpp()
    return x


class JointDistribution:
    """
    Python wrapper for the C++ JointDistribution.

    - All I/O is Python-native (lists, floats).
    - Underlying C++ object is hidden in _cpp_obj.
    """

    __slots__ = ("_cpp_obj",)

    def __init__(self, cpp_obj: Any):
        self._cpp_obj = cpp_obj

    @classmethod
    def _from_cpp(cls, cpp_obj: Any) -> "JointDistribution":
        inst = cls.__new__(cls)
        inst._cpp_obj = cpp_obj
        return inst

    def sample(self, n: Optional[int] = None) -> Union[List[float], List[List[float]]]:
        """
        If n is None: returns one sample as List[float]
        If n is int: returns n samples as List[List[float]]
        """
        if n is None:
            x = self._cpp_obj.sample()
            return [float(v) for v in x]

        n = int(n)
        if n < 0:
            raise ValueError("n must be >= 0")
        X = self._cpp_obj.sample(n)
        return [[float(v) for v in row] for row in X]

    def logpdf(self, x: Sequence[float]) -> float:
        return float(self._cpp_obj.logpdf([float(v) for v in x]))

    def dim(self) -> int:
        return int(self._cpp_obj.dim())

    @property
    def ndim(self) -> int:
        return self.dim()

    def __repr__(self) -> str:
        return f"JointDistribution(dim={self.dim()})"

    def to_cpp(self) -> Any:
        """Optional: keep private in user-facing docs if you want 'no cpp leak'."""
        return self._cpp_obj



class JointDistributionFactory:
    """
    Pure-Python facade over the bound C++ static constructors:

      - st.JointDistribution.create(...)
      - st.JointDistribution.create_with_seeds(...)

    Users pass wrapper enums/configs, and receive a JointDistribution wrapper.
    """

    @staticmethod
    def create(
        marginal_types: Sequence[Any],
        marginal_configs: Sequence[Any],
        copula_type: Any,
        copula_config: Any,
        *,
        seed: Optional[int] = None,
    ) -> JointDistribution:
        if len(marginal_types) != len(marginal_configs):
            raise ValueError("marginal_types and marginal_configs must have the same length.")
        if len(marginal_types) == 0:
            raise ValueError("At least one marginal is required.")

        m_types_cpp = [_as_cpp_enum_marginal(t) for t in marginal_types]
        m_cfgs_cpp = [_as_cpp_config(c) for c in marginal_configs]
        c_type_cpp = _as_cpp_enum_copula(copula_type)
        c_cfg_cpp = _as_cpp_config(copula_config)

        cpp = st.JointDistribution.create(
            m_types_cpp,
            m_cfgs_cpp,
            c_type_cpp,
            c_cfg_cpp,
            seed if seed is not None else None,
        )
        return JointDistribution._from_cpp(cpp)

    @staticmethod
    def create_with_seeds(
        marginal_types: Sequence[Any],
        marginal_configs: Sequence[Any],
        marginal_seeds: Sequence[int],
        copula_type: Any,
        copula_config: Any,
        *,
        copula_seed: int,
    ) -> JointDistribution:
        if len(marginal_types) != len(marginal_configs) or len(marginal_types) != len(marginal_seeds):
            raise ValueError("marginal_types, marginal_configs, marginal_seeds must have the same length.")
        if len(marginal_types) == 0:
            raise ValueError("At least one marginal is required.")

        m_types_cpp = [_as_cpp_enum_marginal(t) for t in marginal_types]
        m_cfgs_cpp = [_as_cpp_config(c) for c in marginal_configs]
        m_seeds = [int(s) for s in marginal_seeds]
        c_type_cpp = _as_cpp_enum_copula(copula_type)
        c_cfg_cpp = _as_cpp_config(copula_config)

        cpp = st.JointDistribution.create_with_seeds(
            m_types_cpp,
            m_cfgs_cpp,
            m_seeds,
            c_type_cpp,
            c_cfg_cpp,
            int(copula_seed),
        )
        return JointDistribution._from_cpp(cpp)


__all__ = [
    "JointDistribution",
    "JointDistributionFactory",
]

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