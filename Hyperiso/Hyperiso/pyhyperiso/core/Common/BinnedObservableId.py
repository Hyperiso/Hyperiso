"""Identifiers for observables evaluated over a finite kinematic bin."""

from dataclasses import dataclass, field
from typing import Tuple, Union, Any

from pyhyperiso.phyperiso.pyhyperiso import common
from pyhyperiso.core.Common.SymbolId import ObservableId
from pyhyperiso.core.Common.Mapper import ObservableMapper
from pyhyperiso.core.Common.GeneralEnum import Observables


@dataclass
class BinnedObservableId:
    """Python wrapper around the C++ ``BinnedObservableId`` composite identifier.

    A ``BinnedObservableId`` uniquely identifies a *binned* observable by:
      - an unbinned observable id (``ObservableId``),
      - a bin range ``(low, high)`` stored as a pair of doubles.

    This mirrors the C++ struct::

        struct BinnedObservableId {
            ObservableId s;
            std::pair<double,double> p;
            ...
        };

    Args:
        s (Union[ObservableId, str]): Unbinned observable identifier.
            - If a ``str`` is provided, it is converted to ``ObservableId(str)``.
        p (Tuple[float, float]): Bin range (low, high). Defaults to (0.0, 0.0).

    Attributes:
        s (ObservableId): Unbinned observable identifier (Python wrapper).
        p (Tuple[float, float]): Bin range.
        _cpp_obj (common.BinnedObservableId): Underlying bound C++ object.

    Notes:
        - ``flha()`` returns a bound C++ ``LhaID``.
        - ``from_flha()`` reconstructs a binned id from such an ``LhaID``.
        - Equality/ordering are implemented in C++ (== and <). This wrapper forwards them.
    """

    s: Union[ObservableId, str] = field(default_factory=lambda: ObservableId("NULL"))
    p: Tuple[float, float] = (0.0, 0.0)
    _cpp_obj: common.BinnedObservableId = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """Normalize inputs and construct the underlying C++ object.

        Raises:
            TypeError: If ``p`` is not a pair (low, high).
        """

        if isinstance(self.s, str):
            self.s = ObservableId(self.s)
        elif isinstance(self.s, Observables):
            self.s = ObservableMapper.to_id(self.s)
        elif not isinstance(self.s, ObservableId):
            raise TypeError(f"s must be ObservableId or str, got {type(self.s)}")

        # Normalize p to (float, float)
        if not (isinstance(self.p, (tuple, list, set)) and len(self.p) == 2):
            raise TypeError("p must be a (low, high) pair")
        if isinstance(self.p, set):
            low, high = sorted(float(value) for value in self.p)
        else:
            low = float(self.p[0])
            high = float(self.p[1])

        self.p = (low, high)

        s_cpp = self.s._to_cpp() 

        if low == 0.0 and high == 0.0:
            self._cpp_obj = common.BinnedObservableId(s_cpp)
        else:
            self._cpp_obj = common.BinnedObservableId(s_cpp, (low, high))

    def flha(self):
        """Convert this binned id to its FLHA/LHA representation.

        Returns:
            common.LhaID: Bound C++ ``LhaID`` for the binned observable.

        Raises:
            RuntimeError: If the mapping from ``ObservableId`` to FLHA is unknown.
        """
        return self._cpp_obj.flha()

    @classmethod
    def from_flha(cls, lhaid: Any) -> "BinnedObservableId":
        """Construct a ``BinnedObservableId`` from a binned FLHA/LHA id.

        Args:
            lhaid: A bound C++ ``LhaID`` or a wrapper exposing ``to_cpp()``.

        Returns:
            BinnedObservableId: Wrapped instance reconstructed from the id.

        Raises:
            RuntimeError: If the id cannot be decoded or the mapping is unknown.
        """
        lhaid_cpp = lhaid.to_cpp() if hasattr(lhaid, "to_cpp") else lhaid
        cpp_obj = common.BinnedObservableId.from_flha(lhaid_cpp)
        return cls.from_cpp(cpp_obj)

    @classmethod
    def from_cpp(cls, cpp_obj: common.BinnedObservableId) -> "BinnedObservableId":
        """Wrap an existing bound C++ ``BinnedObservableId`` instance.

        Args:
            cpp_obj (common.BinnedObservableId): C++ object coming from pybind.

        Returns:
            BinnedObservableId: Python wrapper around the provided C++ object.
        """
        inst = cls.__new__(cls)
        inst._cpp_obj = cpp_obj

        inst.s = ObservableId(str(cpp_obj.s))
        inst.p = (float(cpp_obj.p[0]), float(cpp_obj.p[1]))
        return inst

    def to_cpp(self) -> common.BinnedObservableId:
        """Return the underlying bound C++ object."""
        return self._cpp_obj

    def to_dict(self) -> dict:
        """Serialize this identifier into a debug-friendly dict."""
        return {"s": str(self.s), "p": self.p}

    def __repr__(self) -> str:
        """Return a debug representation."""
        return f"BinnedObservableId(s={self.s!s}, p={self.p})"

    def __str__(self) -> str:
        """Return a human-readable representation."""
        return f"{self.s} [{self.p[0]}, {self.p[1]}]"

    def __eq__(self, other: object) -> bool:
        """Value equality using the underlying C++ implementation."""
        if not isinstance(other, BinnedObservableId):
            return NotImplemented
        return self._cpp_obj == other._cpp_obj

    def __lt__(self, other: "BinnedObservableId") -> bool:
        """Strict ordering using the underlying C++ implementation."""
        if not isinstance(other, BinnedObservableId):
            return NotImplemented
        return self._cpp_obj < other._cpp_obj

    def __hash__(self) -> int:
        """Hash using the underlying C++ hash (bound __hash__) if available.

        Falls back to hashing (s, p) in Python if needed.
        """
        try:
            return hash(self._cpp_obj)
        except TypeError:
            return hash((str(self.s), float(self.p[0]), float(self.p[1])))
