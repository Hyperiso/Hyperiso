"""Configuration dataclasses for one-dimensional marginal distributions.

Each dataclass mirrors a C++ marginal configuration struct and exposes a
``to_cpp`` method used by the distribution factories. These configurations are
combined with copulas to build joint distributions in the statistic backend.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, Sequence, Union

from pyhyperiso.phyperiso.pyhyperiso import statistic as st


class MarginalKind(Enum):
    """Supported one-dimensional marginal distribution families.

    Attributes:
        GAUSSIAN: Symmetric Gaussian marginal.
        HALF_GAUSSIAN: Split/asymmetric Gaussian marginal in the backend naming.
        FLAT: Uniform marginal over a finite interval.
        LIKELIHOOD: Empirical marginal represented by values and weights.
    """

    GAUSSIAN = "GAUSSIAN"
    HALF_GAUSSIAN = "HALF_GAUSSIAN"
    FLAT = "FLAT"
    LIKELIHOOD = "LIKELIHOOD"

    def to_cpp(self) -> Any:
        """Convert the Python enum value to the bound C++ ``MarginalType``."""
        return getattr(st.MarginalType, self.value)


@dataclass(frozen=True)
class FlatMarginalConfig:
    """Configuration for a uniform marginal distribution.

    Args:
        a: Lower bound of the support.
        b: Upper bound of the support.

    Examples:
        >>> cfg = FlatMarginalConfig(a=-1.0, b=1.0)
        >>> cfg.a, cfg.b
        (-1.0, 1.0)
    """

    a: float = 0.0
    b: float = 1.0

    def to_cpp(self) -> Any:
        """Build a bound C++ ``FlatMarginalCfg`` instance."""
        return st.FlatMarginalCfg(self.a, self.b)

    @classmethod
    def from_cpp(cls, cpp: Any) -> "FlatMarginalConfig":
        """Create a Python config from a bound C++ flat config."""
        return cls(a=float(cpp.a), b=float(cpp.b))


@dataclass(frozen=True)
class GaussianMarginalConfig:
    """Configuration for a Gaussian marginal distribution.

    Args:
        mu: Central value of the Gaussian.
        sigma: Standard deviation. It should be strictly positive.
    """

    mu: float = 0.0
    sigma: float = 1.0

    def to_cpp(self) -> Any:
        """Build a bound C++ ``GaussianMarginalCfg`` instance."""
        return st.GaussianMarginalCfg(self.mu, self.sigma)

    @classmethod
    def from_cpp(cls, cpp: Any) -> "GaussianMarginalConfig":
        """Create a Python config from a bound C++ Gaussian config."""
        return cls(mu=float(cpp.mu), sigma=float(cpp.sigma))


@dataclass(frozen=True)
class SplitGaussianMarginalConfig:
    """Configuration for an asymmetric split-Gaussian marginal.

    Args:
        mu: Central value or mode of the distribution.
        sigma_p: Right-side standard deviation, used for values above ``mu``.
        sigma_m: Left-side standard deviation, used for values below ``mu``.
    """

    mu: float = 0.0
    sigma_p: float = 1.0
    sigma_m: float = 1.0

    def to_cpp(self) -> Any:
        """Build a bound C++ ``SplitGaussianMarginalCfg`` instance."""
        return st.SplitGaussianMarginalCfg(self.mu, self.sigma_p, self.sigma_m)

    @classmethod
    def from_cpp(cls, cpp: Any) -> "SplitGaussianMarginalConfig":
        """Create a Python config from a bound C++ split-Gaussian config."""
        return cls(mu=float(cpp.mu), sigma_p=float(cpp.sigma_p), sigma_m=float(cpp.sigma_m))


@dataclass(frozen=True)
class LikelihoodMarginalConfig:
    """Configuration for an empirical likelihood marginal.

    Args:
        values: Grid or sample values supporting the empirical likelihood.
        weights: Non-negative weights associated with ``values``. The two
            sequences must have the same length.

    Raises:
        ValueError: If ``values`` and ``weights`` have different lengths.

    Examples:
        >>> cfg = LikelihoodMarginalConfig(values=[0.0, 1.0], weights=[0.25, 0.75])
        >>> len(cfg.values) == len(cfg.weights)
        True
    """

    values: Sequence[float]
    weights: Sequence[float]

    def __post_init__(self) -> None:
        """Validate the empirical likelihood support."""
        if len(self.values) != len(self.weights):
            raise ValueError("LikelihoodMarginalConfig: values et weights must have the same size.")

    def to_cpp(self) -> Any:
        """Build a bound C++ ``LikelihoodMarginalCfg`` instance."""
        return st.LikelihoodMarginalCfg(
            list(map(float, self.values)), list(map(float, self.weights))
        )

    @classmethod
    def from_cpp(cls, cpp: Any) -> "LikelihoodMarginalConfig":
        """Create a Python config from a bound C++ likelihood config."""
        return cls(values=list(map(float, cpp.values)), weights=list(map(float, cpp.weights)))


MarginalConfig = Union[
    FlatMarginalConfig,
    GaussianMarginalConfig,
    SplitGaussianMarginalConfig,
    LikelihoodMarginalConfig,
]


def _config_from_cpp(cpp_cfg: Any) -> MarginalConfig:
    """Dispatch a bound C++ marginal config to the matching Python dataclass.

    Args:
        cpp_cfg: Bound C++ marginal configuration object.

    Returns:
        The corresponding Python configuration dataclass.

    Raises:
        TypeError: If the C++ config type is not recognized.
    """
    if isinstance(cpp_cfg, st.FlatMarginalCfg):
        return FlatMarginalConfig.from_cpp(cpp_cfg)
    if isinstance(cpp_cfg, st.GaussianMarginalCfg):
        return GaussianMarginalConfig.from_cpp(cpp_cfg)
    if isinstance(cpp_cfg, st.SplitGaussianMarginalCfg):
        return SplitGaussianMarginalConfig.from_cpp(cpp_cfg)
    if isinstance(cpp_cfg, st.LikelihoodMarginalCfg):
        return LikelihoodMarginalConfig.from_cpp(cpp_cfg)

    raise TypeError(f"Config C++ not known: {type(cpp_cfg)!r}")
