"""Joint probability distributions built from marginals and a copula.

The C++ backend represents a joint law as a set of one-dimensional marginal
distributions plus a copula. This Python module exposes wrappers and factory
helpers to build such distributions from Python configuration objects.
"""

from __future__ import annotations

from typing import List, Optional, Sequence, Union

from pyhyperiso.phyperiso.pyhyperiso import statistic as st
from pyhyperiso.core.Statistic.Copula import CopulaKind
from pyhyperiso.core.Statistic.CopulaConfig import (
    GaussianCopulaConfigPy as GaussianCopulaConfig,
    StudentTCopulaConfigPy as StudentTCopulaConfig,
)
from pyhyperiso.core.Statistic.MarginalConfig import (
    FlatMarginalConfig,
    GaussianMarginalConfig,
    LikelihoodMarginalConfig,
    MarginalKind,
    SplitGaussianMarginalConfig,
)


MarginalConfig = Union[
    FlatMarginalConfig,
    GaussianMarginalConfig,
    SplitGaussianMarginalConfig,
    LikelihoodMarginalConfig,
]
CopulaConfig = Union[GaussianCopulaConfig, StudentTCopulaConfig]


def _require(value, typ, name: str):
    """Validate the type of a user-facing argument.

    Args:
        value: Object to validate.
        typ: Expected Python type.
        name: Argument name used in the error message.

    Returns:
        The original value.

    Raises:
        TypeError: If ``value`` is not an instance of ``typ``.
    """
    if not isinstance(value, typ):
        raise TypeError(f"{name} must be {typ.__name__}, received {type(value)!r}.")
    return value


def _cpp_marginal_kind(kind: MarginalKind):
    """Convert a Python marginal kind to the C++ enum value."""
    return _require(kind, MarginalKind, "MarginalKind").to_cpp()


def _cpp_copula_kind(kind: CopulaKind):
    """Convert a Python copula kind to the C++ enum value."""
    return _require(kind, CopulaKind, "CopulaKind").to_cpp()


def _cpp_marginal_config(cfg: MarginalConfig):
    """Convert a Python marginal configuration to its C++ representation.

    Raises:
        TypeError: If the configuration family is not supported.
    """
    if not isinstance(
        cfg,
        (
            FlatMarginalConfig,
            GaussianMarginalConfig,
            SplitGaussianMarginalConfig,
            LikelihoodMarginalConfig,
        ),
    ):
        raise TypeError(f"unsupported marginal configuration: {type(cfg)!r}.")
    return cfg.to_cpp()


def _cpp_copula_config(cfg: CopulaConfig):
    """Convert a Python copula configuration to its C++ representation.

    Raises:
        TypeError: If the configuration family is not supported.
    """
    if not isinstance(cfg, (GaussianCopulaConfig, StudentTCopulaConfig)):
        raise TypeError(f"unsupported copula configuration: {type(cfg)!r}.")
    return cfg.to_cpp()


class JointDistribution:
    """Python wrapper around a C++ joint distribution.

    A joint distribution combines marginal distributions with a dependence
    structure encoded by a copula. It is used by the statistic layer for nuisance
    parameters and experimental observables.

    Args:
        cpp_obj: Bound C++ ``JointDistribution`` instance.
    """

    __slots__ = ("_cpp_obj",)

    def __init__(self, cpp_obj) -> None:
        """Store the bound C++ joint distribution object."""
        self._cpp_obj = cpp_obj

    @classmethod
    def from_cpp(cls, cpp_obj) -> "JointDistribution":
        """Wrap a bound C++ joint distribution.

        Args:
            cpp_obj: Bound C++ ``JointDistribution`` instance.

        Returns:
            A Python wrapper retaining the C++ object.
        """
        return cls(cpp_obj)

    def _to_cpp(self):
        """Return the underlying C++ object for internal binding calls."""
        return self._cpp_obj

    def sample(self, n: Optional[int] = None) -> Union[List[float], List[List[float]]]:
        """Draw one or more samples from the joint distribution.

        Args:
            n: Optional number of samples. When omitted, a single sample vector is
                returned.

        Returns:
            A single sample ``list[float]`` when ``n`` is ``None``; otherwise a
            list of sample vectors with shape ``n x dim``.

        Raises:
            ValueError: If ``n`` is negative.
        """
        if n is None:
            return [float(v) for v in self._cpp_obj.sample()]
        if int(n) < 0:
            raise ValueError("n must be >= 0.")
        return [[float(v) for v in row] for row in self._cpp_obj.sample(int(n))]

    def logpdf(self, x: Sequence[float]) -> float:
        """Evaluate the joint log-density at ``x``.

        Args:
            x: Point in physical variable space. Its length must match the joint
                distribution dimension.

        Returns:
            The scalar log-density ``log f(x)``.
        """
        return float(self._cpp_obj.logpdf([float(v) for v in x]))

    def dim(self) -> int:
        """Return the distribution dimension."""
        return int(self._cpp_obj.dim())

    @property
    def ndim(self) -> int:
        """Alias for :meth:`dim`, following NumPy naming conventions."""
        return self.dim()

    def __repr__(self) -> str:
        """Return a compact representation useful in notebooks and logs."""
        return f"JointDistribution(dim={self.dim()})"


class JointDistributionFactory:
    """Factory helpers for C++ joint distributions.

    The factory receives Python marginal and copula configurations, converts
    them to C++ objects and returns a :class:`JointDistribution` wrapper.
    """

    @staticmethod
    def create(
        marginal_types: Sequence[MarginalKind],
        marginal_configs: Sequence[MarginalConfig],
        copula_type: CopulaKind,
        copula_config: CopulaConfig,
        *,
        seed: Optional[int] = None,
    ) -> JointDistribution:
        """Create a joint distribution with a shared optional seed.

        Args:
            marginal_types: Marginal distribution family for each dimension.
            marginal_configs: Configuration object for each marginal.
            copula_type: Copula family encoding dependence.
            copula_config: Copula configuration, typically a correlation matrix.
            seed: Optional seed passed to the C++ factory.

        Returns:
            A ``JointDistribution`` wrapper.

        Raises:
            ValueError: If marginal type/config lengths differ or no marginal is
                provided.

        Examples:
            >>> jd = JointDistributionFactory.create(
            ...     [MarginalKind.GAUSSIAN, MarginalKind.FLAT],
            ...     [GaussianMarginalConfig(0.0, 1.0), FlatMarginalConfig(-1.0, 1.0)],
            ...     CopulaKind.GAUSSIAN,
            ...     GaussianCopulaConfig(R=[[1.0, 0.2], [0.2, 1.0]]),
            ...     seed=123,
            ... )
            >>> jd.ndim
            2
        """
        if len(marginal_types) != len(marginal_configs):
            raise ValueError("marginal_types et marginal_configs must have the same size.")
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
        """Create a joint distribution with explicit marginal and copula seeds.

        Args:
            marginal_types: Marginal distribution family for each dimension.
            marginal_configs: Configuration object for each marginal.
            marginal_seeds: Seed used for each marginal RNG.
            copula_type: Copula family encoding dependence.
            copula_config: Copula configuration.
            copula_seed: Seed used for the copula RNG.

        Returns:
            A ``JointDistribution`` wrapper.

        Raises:
            ValueError: If input sequence lengths differ or no marginal is
                provided.
        """
        if len(marginal_types) != len(marginal_configs) or len(marginal_types) != len(
            marginal_seeds
        ):
            raise ValueError(
                "marginal_types, marginal_configs et marginal_seeds must have the same size."
            )
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
