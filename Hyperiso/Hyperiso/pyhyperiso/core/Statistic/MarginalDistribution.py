"""Python wrappers and factories for one-dimensional marginal distributions."""

from __future__ import annotations

from typing import List, Optional, Sequence, cast

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
    """Validate a user-facing argument type and return the original value."""
    if not isinstance(value, typ):
        raise TypeError(f"{name} must be {typ.__name__}, received {type(value)!r}.")
    return value


def _cpp_marginal_kind(kind: MarginalKind):
    """Convert a Python marginal kind to the C++ enum value."""
    return _require(kind, MarginalKind, "MarginalKind").to_cpp()


def _cpp_param_id(pid: ParamId):
    """Convert a Python parameter identifier to its C++ representation."""
    return _require(pid, ParamId, "ParamId").to_cpp()


def _cpp_experiment_obs(obs: ExperimentObs):
    """Convert a Python experimental observable wrapper to C++."""
    return _require(obs, ExperimentObs, "ExperimentObs").to_cpp()


def _cpp_marginal_config(config: MarginalConfig):
    """Convert a supported Python marginal config to C++.

    Raises:
        TypeError: If ``config`` is not one of the supported marginal config
            dataclasses.
    """
    if not isinstance(
        config,
        (
            FlatMarginalConfig,
            GaussianMarginalConfig,
            SplitGaussianMarginalConfig,
            LikelihoodMarginalConfig,
        ),
    ):
        raise TypeError(f"unsupported marginal configuration: {type(config)!r}.")
    return config.to_cpp()


class MarginalDistribution:
    """Base wrapper around a bound C++ marginal distribution.

    Args:
        cpp_obj: Bound C++ marginal distribution object.
    """

    def __init__(self, cpp_obj) -> None:
        """Store the bound C++ marginal distribution object."""
        self._cpp_obj = cpp_obj

    @classmethod
    def from_cpp(cls, cpp_obj) -> "MarginalDistribution":
        """Wrap a C++ marginal with the most specific Python subclass.

        Args:
            cpp_obj: Bound C++ marginal distribution object.

        Returns:
            A subclass wrapper when the C++ type is recognized, otherwise a
            generic ``MarginalDistribution``.
        """
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
        """Return the underlying C++ object for internal calls."""
        return self._cpp_obj

    def rvs(self, n: int) -> List[float]:
        """Draw independent random variates from the marginal.

        Args:
            n: Number of variates to draw.

        Returns:
            A list of ``n`` floats.

        Raises:
            ValueError: If ``n`` is negative.
        """
        if int(n) < 0:
            raise ValueError("n must be >= 0.")
        return [float(v) for v in self._cpp_obj.rvs(int(n))]

    def logpdf(self, x: float) -> float:
        """Evaluate the log-density at ``x``."""
        return float(self._cpp_obj.logpdf(float(x)))

    def cdf(self, x: float) -> float:
        """Evaluate the cumulative distribution function at ``x``."""
        return float(self._cpp_obj.cdf(float(x)))

    def ppf(self, p: float) -> float:
        """Evaluate the quantile function.

        Args:
            p: Probability in the closed interval ``[0, 1]``.

        Returns:
            The value ``x`` such that ``cdf(x) ~= p``.

        Raises:
            ValueError: If ``p`` is outside ``[0, 1]``.
        """
        p = float(p)
        if not 0.0 <= p <= 1.0:
            raise ValueError("p must be dans [0, 1].")
        return float(self._cpp_obj.ppf(p))

    def mean(self) -> float:
        """Return the marginal mean."""
        return float(self._cpp_obj.mean())

    def std(self) -> float:
        """Return the marginal standard deviation."""
        return float(self._cpp_obj.std())


class GaussianMarginalDist(MarginalDistribution):
    """Wrapper for a symmetric Gaussian marginal distribution."""

    pass


class SplitGaussianMarginalDist(MarginalDistribution):
    """Wrapper for an asymmetric split-Gaussian marginal distribution."""

    pass


class FlatMarginalDist(MarginalDistribution):
    """Wrapper for a uniform marginal distribution."""

    pass


class LikelihoodMarginalDist(MarginalDistribution):
    """Wrapper for an empirical likelihood marginal distribution."""

    pass


class DistributionFactoryWrapper:
    """Factory helpers for marginal distributions."""

    @staticmethod
    def create(
        kind: MarginalKind, config: MarginalConfig, seed: Optional[int] = None
    ) -> MarginalDistribution:
        """Create a marginal distribution from a kind and configuration.

        Args:
            kind: Marginal family to instantiate.
            config: Configuration object compatible with ``kind``.
            seed: Optional RNG seed.

        Returns:
            A Python wrapper around the bound C++ marginal.
        """
        cpp_dist = st.MarginalFactory.create(
            _cpp_marginal_kind(kind), _cpp_marginal_config(config), seed
        )
        return MarginalDistribution.from_cpp(cpp_dist)

    @staticmethod
    def gaussian(mu: float, sigma: float, seed: Optional[int] = None) -> GaussianMarginalDist:
        """Create a Gaussian marginal.

        Args:
            mu: Central value.
            sigma: Standard deviation.
            seed: Optional RNG seed.

        Returns:
            A ``GaussianMarginalDist`` wrapper.
        """
        return cast(
            GaussianMarginalDist,
            DistributionFactoryWrapper.create(
                MarginalKind.GAUSSIAN, GaussianMarginalConfig(mu=mu, sigma=sigma), seed
            ),
        )

    @staticmethod
    def split_gaussian(
        mu: float, sigma_p: float, sigma_m: float, seed: Optional[int] = None
    ) -> SplitGaussianMarginalDist:
        """Create an asymmetric split-Gaussian marginal.

        Args:
            mu: Central value or mode.
            sigma_p: Right-side standard deviation.
            sigma_m: Left-side standard deviation.
            seed: Optional RNG seed.

        Returns:
            A ``SplitGaussianMarginalDist`` wrapper.
        """
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
        """Create a uniform marginal over ``[a, b]``."""
        return cast(
            FlatMarginalDist,
            DistributionFactoryWrapper.create(
                MarginalKind.FLAT, FlatMarginalConfig(a=a, b=b), seed
            ),
        )

    @staticmethod
    def likelihood(
        values: Sequence[float],
        weights: Sequence[float],
        seed: Optional[int] = None,
        standardize: bool = False,
    ) -> LikelihoodMarginalDist:
        """Create an empirical likelihood marginal.

        Args:
            values: Empirical support values.
            weights: Weights associated with each support value.
            seed: Optional RNG seed.
            standardize: If true, use the backend constructor that standardizes
                the likelihood marginal.

        Returns:
            A ``LikelihoodMarginalDist`` wrapper.

        Raises:
            ValueError: If ``values`` and ``weights`` have different lengths.
        """
        cfg = LikelihoodMarginalConfig(values=values, weights=weights)
        if standardize:
            cpp_obj = st.LikelihoodMarginal(
                list(map(float, cfg.values)),
                list(map(float, cfg.weights)),
                int(0 if seed is None else seed),
                True,
            )
            return cast(LikelihoodMarginalDist, MarginalDistribution.from_cpp(cpp_obj))
        return cast(
            LikelihoodMarginalDist,
            DistributionFactoryWrapper.create(MarginalKind.LIKELIHOOD, cfg, seed),
        )


class MarginalConfigFactoryWrapper:
    """Wrapper around the C++ ``MarginalConfigFactory``.

    The C++ factory can infer a marginal configuration from a model parameter or
    an experimental observable, using the current HyperIso parameter database.
    """

    def __init__(self) -> None:
        """Create the underlying C++ marginal configuration factory."""
        self._cpp_obj = st.MarginalConfigFactory()

    def create_from_param(self, pid: ParamId, kind: MarginalKind) -> MarginalConfig:
        """Infer a marginal config for a parameter.

        Args:
            pid: Parameter identifier.
            kind: Desired marginal family.

        Returns:
            A Python marginal configuration dataclass.
        """
        return _config_from_cpp(self._cpp_obj.create(_cpp_param_id(pid), _cpp_marginal_kind(kind)))

    def create_from_observable(self, obs: ExperimentObs, kind: MarginalKind) -> MarginalConfig:
        """Infer a marginal config for an experimental observable.

        Args:
            obs: Experimental observable key returned by the statistic layer.
            kind: Desired marginal family.

        Returns:
            A Python marginal configuration dataclass.
        """
        return _config_from_cpp(
            self._cpp_obj.create(_cpp_experiment_obs(obs), _cpp_marginal_kind(kind))
        )


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
