"""Python representation of experimental observable identifiers.

The C++ backend uses ``ExperimentObs`` keys in maps storing experimental central
values and correlations. The current binding exposes these objects when they are
returned by the statistic layer, but does not expose a public pure-Python
constructor for the underlying C++ type. This wrapper therefore keeps the C++
object opaque while exposing the associated binned observable identifier.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from pyhyperiso.core.Common.BinnedObservableId import BinnedObservableId


@dataclass(frozen=True)
class ExperimentObs:
    """Opaque Python wrapper around a C++ ``ExperimentObs`` key.

    Args:
        obs: Binned observable identifier associated with the experimental
            measurement.
        _cpp_obj: Bound C++ object. This is normally set by :meth:`from_cpp` and
            should not be passed by user code.

    Attributes:
        obs: Python ``BinnedObservableId`` exposed for display, hashing and map
            lookups.

    Raises:
        TypeError: If ``obs`` is not a ``BinnedObservableId``.

    Examples:
        ``ExperimentObs`` instances are usually obtained from statistic calls::

            exp_values = stat.get_obs_exp()
            for exp_obs, value in exp_values.items():
                print(exp_obs.obs, value)
    """

    obs: BinnedObservableId
    _cpp_obj: Any = field(default=None, repr=False, compare=False)

    def __post_init__(self) -> None:
        """Validate the public observable identifier field."""
        if not isinstance(self.obs, BinnedObservableId):
            raise TypeError(
                f"obs must be a Python BinnedObservableId, received {type(self.obs)!r}."
            )

    @classmethod
    def from_cpp(cls, cpp_obj) -> "ExperimentObs":
        """Create an ``ExperimentObs`` wrapper from a bound C++ object.

        Args:
            cpp_obj: Bound C++ ``ExperimentObs`` instance.

        Returns:
            A Python wrapper retaining the original C++ object.
        """
        inst = cls.__new__(cls)
        object.__setattr__(inst, "obs", BinnedObservableId.from_cpp(cpp_obj.obs))
        object.__setattr__(inst, "_cpp_obj", cpp_obj)
        return inst

    def to_cpp(self):
        """Return the underlying C++ object.

        Returns:
            The bound C++ ``ExperimentObs`` instance.

        Raises:
            RuntimeError: If the object was built without a C++ counterpart.
        """
        if self._cpp_obj is None:
            raise RuntimeError(
                "ExperimentObs cannot yet be constructed in pure Python: "
                "use an ExperimentObs returned by StatisticInterface.get_obs_exp() "
                "or expose a dedicated C++ constructor in the binding."
            )
        return self._cpp_obj

    def __hash__(self) -> int:
        """Hash the wrapper through its binned observable identifier."""
        return hash(self.obs)

    def __str__(self) -> str:
        """Return a compact string representation of the observable identifier."""
        return str(self.obs)

    def __repr__(self) -> str:
        """Return an unambiguous debugging representation."""
        return f"ExperimentObs(obs={self.obs!r})"


__all__ = ["ExperimentObs"]
