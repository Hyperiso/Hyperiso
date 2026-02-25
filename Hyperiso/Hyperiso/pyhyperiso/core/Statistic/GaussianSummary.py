from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from pyhyperiso.phyperiso.pyhyperiso import statistic as _cpp_stat
from pyhyperiso.core.Common.BinnedObservableId import BinnedObservableId  # adapte l'import


@dataclass
class GaussianSummary:
    """Summary statistics of a (possibly split) Gaussian approximation.

    Python wrapper around the C++ ``GaussianSummary`` struct.

    Attributes:
        id (BinnedObservableId): Binned observable identifier (wrapper).
        mu (float): Population mean.
        sigma (float): Population standard deviation.
        sigma_p (float): Right-side population std (split Gaussian).
        sigma_m (float): Left-side population std (split Gaussian).
        mode (float): Population mode.
        skew (float): Sample skewness.
        symmetric (bool): Whether the distribution is considered symmetric.
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
        """Wrap a bound C++ ``GaussianSummary`` instance."""
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
        """Convert this wrapper into a bound C++ ``GaussianSummary`` instance."""
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