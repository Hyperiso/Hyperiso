from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, List

from pyhyperiso.phyperiso.pyhyperiso import statistic as _cpp_stat
from pyhyperiso.core.Statistic.GaussianSummary import GaussianSummary


ObsSamples = Any
NuisanceSamples = Any


@dataclass
class MCRealization:
    """One Monte Carlo realization of sampled observables and nuisance parameters.

    Python wrapper around the C++ ``MCRealization`` struct.

    Attributes:
        sampled_obss: Samples for observables (Python-native structure).
        sampled_params: Samples for nuisance parameters (Python-native structure).
    """
    sampled_obss: ObsSamples
    sampled_params: NuisanceSamples

    @classmethod
    def from_cpp(cls, cpp_obj: Any) -> "MCRealization":
        """Wrap a bound C++ ``MCRealization`` instance."""
        # pybind11/stl va déjà convertir vector->list si c’est du STL
        return cls(
            sampled_obss=cpp_obj.sampled_obss,
            sampled_params=cpp_obj.sampled_params,
        )

    def to_cpp(self):
        """Convert this wrapper into a bound C++ ``MCRealization`` instance."""
        cpp = _cpp_stat.MCRealization()
        cpp.sampled_obss = self.sampled_obss
        cpp.sampled_params = self.sampled_params
        return cpp


@dataclass
class MCResult:
    """Monte Carlo output: one realization + summaries.

    Python wrapper around the C++ ``MCResult`` struct.

    Attributes:
        mc_real: The underlying Monte Carlo realization.
        summary: List of GaussianSummary objects.
    """
    mc_real: MCRealization
    summary: List[GaussianSummary] = field(default_factory=list)

    @classmethod
    def from_cpp(cls, cpp_obj: Any) -> "MCResult":
        """Wrap a bound C++ ``MCResult`` instance."""
        return cls(
            mc_real=MCRealization.from_cpp(cpp_obj.mc_real),
            summary=[GaussianSummary.from_cpp(gs) for gs in cpp_obj.summary],
        )

    def to_cpp(self):
        """Convert this wrapper into a bound C++ ``MCResult`` instance."""
        cpp = _cpp_stat.MCResult()
        cpp.mc_real = self.mc_real.to_cpp()
        cpp.summary = [gs.to_cpp() for gs in self.summary]
        return cpp