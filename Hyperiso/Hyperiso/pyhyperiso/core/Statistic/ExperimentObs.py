"""Experiment-scoped observable identifiers for statistical selections."""

from __future__ import annotations

from dataclasses import dataclass, field

from pyhyperiso.phyperiso.pyhyperiso import common
from pyhyperiso.core.Common.BinnedObservableId import BinnedObservableId


@dataclass(frozen=True)
class ExperimentObs:
    """Identify one experimental measurement by experiment and observable bin.

    Args:
        experiment: Experiment label as stored in the experimental database.
        obs: Binned observable identifier associated with the measurement.

    Attributes:
        experiment: Experiment label.
        obs: Python :class:`BinnedObservableId` wrapper.

    Examples:
        Build an exact selection entry for a CMS measurement::

            selection = ExperimentObs(
                "CMS",
                BinnedObservableId(Observables.F_L_B0__KSTAR0_MU_MU, (1.1, 2.0)),
            )
    """

    experiment: str
    obs: BinnedObservableId
    _cpp_obj: common.ExperimentObs = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        """Validate fields and construct the bound C++ value object."""
        if not isinstance(self.experiment, str):
            raise TypeError(f"experiment must be str, received {type(self.experiment)!r}.")
        if not self.experiment:
            raise ValueError("experiment must not be empty.")
        if not isinstance(self.obs, BinnedObservableId):
            raise TypeError(
                f"obs must be a Python BinnedObservableId, received {type(self.obs)!r}."
            )

        object.__setattr__(
            self,
            "_cpp_obj",
            common.ExperimentObs(self.experiment, self.obs.to_cpp()),
        )

    @classmethod
    def from_cpp(cls, cpp_obj: common.ExperimentObs) -> "ExperimentObs":
        """Create a Python wrapper from a bound C++ ``ExperimentObs``."""
        inst = cls.__new__(cls)
        object.__setattr__(inst, "experiment", str(cpp_obj.experiment))
        object.__setattr__(inst, "obs", BinnedObservableId.from_cpp(cpp_obj.obs))
        object.__setattr__(inst, "_cpp_obj", cpp_obj)
        return inst

    def to_cpp(self) -> common.ExperimentObs:
        """Return the bound C++ value object."""
        return self._cpp_obj

    def __hash__(self) -> int:
        """Hash both the experiment label and binned observable identifier."""
        return hash((self.experiment, self.obs))

    def __str__(self) -> str:
        """Return a compact human-readable identifier."""
        return f"{self.experiment} :: {self.obs}"

    def __repr__(self) -> str:
        """Return an unambiguous debugging representation."""
        return f"ExperimentObs(experiment={self.experiment!r}, obs={self.obs!r})"


__all__ = ["ExperimentObs"]
