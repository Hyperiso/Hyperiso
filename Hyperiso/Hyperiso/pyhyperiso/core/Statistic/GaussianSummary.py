from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from pyhyperiso.phyperiso.pyhyperiso import statistic as _cpp_stat
from pyhyperiso.core.Common.GeneralEnum import Observables  # adapte si ton wrapper s'appelle autrement


@dataclass
class GaussianSummary:
    """Summary statistics of a (possibly split) Gaussian approximation.

    This is a Python wrapper around the C++ ``GaussianSummary`` struct.

    Attributes:
        id: Observable identifier.
        mu: Population mean.
        sigma: Population standard deviation.
        sigma_p: Right-side population std (for split Gaussian).
        sigma_m: Left-side population std (for split Gaussian).
        mode: Population mode.
        skew: Sample skewness.
        symmetric: Whether the distribution is considered symmetric.
    """

    id: ObservableId
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
            id=cpp_obj.id,
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
        cpp.id = self.id
        cpp.mu = float(self.mu)
        cpp.sigma = float(self.sigma)
        cpp.sigma_p = float(self.sigma_p)
        cpp.sigma_m = float(self.sigma_m)
        cpp.mode = float(self.mode)
        cpp.skew = float(self.skew)
        cpp.symmetric = bool(self.symmetric)
        return cpp