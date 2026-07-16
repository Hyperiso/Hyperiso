"""Summary statistics for Monte-Carlo observable distributions."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from pyhyperiso.phyperiso.pyhyperiso import statistic as _cpp_stat
from pyhyperiso.core.Common.BinnedObservableId import BinnedObservableId


@dataclass
class GaussianSummary:
    """Gaussian or split-Gaussian approximation of one observable distribution.

    This class mirrors the C++ ``GaussianSummary`` struct returned by the
    Monte-Carlo uncertainty machinery. When the empirical skewness is small, the
    distribution is summarized by ``mu`` and ``sigma``. Otherwise, the backend
    also provides a mode and asymmetric widths ``sigma_p`` and ``sigma_m``.

    Attributes:
        id: Binned observable identifier.
        mu: Population mean of the sampled observable.
        sigma: Unbiased standard deviation used in the symmetric approximation.
        sigma_p: Right-side width for a split-Gaussian approximation.
        sigma_m: Left-side width for a split-Gaussian approximation.
        mode: Estimated population mode.
        skew: Empirical skewness estimator.
        symmetric: Whether the backend classified the sample as symmetric.

    Examples:
        >>> summary = GaussianSummary(id=obs_id, mu=1.0, sigma=0.1, symmetric=True)
        >>> summary.mu, summary.sigma
        (1.0, 0.1)
    """

    id: BinnedObservableId
    mu: float = 0.0
    sigma: float = 0.0
    sigma_p: float = 0.0
    sigma_m: float = 0.0
    mode: float = 0.0
    skew: float = 0.0
    symmetric: bool = False

    @classmethod
    def from_cpp(cls, cpp_obj: Any) -> "GaussianSummary":
        """Create a Python summary from a bound C++ summary.

        Args:
            cpp_obj: Bound C++ ``GaussianSummary`` instance.

        Returns:
            The equivalent Python dataclass.
        """
        return cls(
            id=BinnedObservableId.from_cpp(cpp_obj.id),
            mu=float(cpp_obj.mu),
            sigma=float(cpp_obj.sigma),
            sigma_p=float(cpp_obj.sigma_p),
            sigma_m=float(cpp_obj.sigma_m),
            mode=float(cpp_obj.mode),
            skew=float(cpp_obj.skew),
            symmetric=bool(cpp_obj.symmetric),
        )

    def to_cpp(self):
        """Convert this dataclass to a bound C++ ``GaussianSummary``.

        Returns:
            A newly allocated bound C++ summary object.
        """
        cpp = _cpp_stat.GaussianSummary()
        cpp.id = self.id.to_cpp()
        cpp.mu = float(self.mu)
        cpp.sigma = float(self.sigma)
        cpp.sigma_p = float(self.sigma_p)
        cpp.sigma_m = float(self.sigma_m)
        cpp.mode = float(self.mode)
        cpp.skew = float(self.skew)
        cpp.symmetric = bool(self.symmetric)
        return cpp
