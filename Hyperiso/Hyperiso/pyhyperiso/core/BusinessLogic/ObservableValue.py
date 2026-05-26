"""Python representation of computed observable values.

The C++ observable layer returns ``ObservableValue`` objects containing an
observable id, a numerical prediction, and optionally a bin range. This module
keeps the public Python type lightweight and immutable while preserving explicit
conversion helpers to and from the bound C++ object.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

from pyhyperiso.phyperiso.pyhyperiso.observable import ObservableValue as _CppObservableValue
from pyhyperiso.core.Common.GeneralEnum import Observables
from pyhyperiso.core.Common.Mapper import ObservableMapper
from pyhyperiso.core.Common.SymbolId import ObservableId


def _require_observable_id(value: ObservableId, name: str = "observable id") -> ObservableId:
    if not isinstance(value, ObservableId):
        raise TypeError(f"{name} doit être un ObservableId Python, reçu {type(value)!r}.")
    return value


def _cpp_observable_id(value: ObservableId):
    return _require_observable_id(value)._to_cpp()


@dataclass(frozen=True)
class ObservableValue:
    """Value returned by an observable computation.

    Attributes:
        id: Internal observable identifier.
        value: Numerical prediction or experimental value.
        bin: Optional bin range ``(low, high)``. ``None`` means that the
            observable is unbinned or that the C++ value did not carry bin
            information.

    Example:
        >>> from pyhyperiso.core.BusinessLogic.ObservableValue import ObservableValue
        >>> from pyhyperiso.core.Common.GeneralEnum import Observables
        >>> ov = ObservableValue.from_observable(Observables.BR_BS_MUMU, 3.6e-9)
        >>> ov.value
        3.6e-09
    """

    id: ObservableId
    value: float
    bin: Optional[Tuple[float, float]] = None

    def __post_init__(self) -> None:
        """Validate and normalize the immutable dataclass fields.

        Raises:
            TypeError: If ``id`` is not an ``ObservableId`` or if ``bin`` is
                neither ``None`` nor a two-entry tuple.
        """
        _require_observable_id(self.id, "id")
        if self.bin is not None:
            if not (isinstance(self.bin, tuple) and len(self.bin) == 2):
                raise TypeError("bin doit être None ou un tuple (low, high).")
            object.__setattr__(self, "bin", (float(self.bin[0]), float(self.bin[1])))
        object.__setattr__(self, "value", float(self.value))

    @classmethod
    def from_observable(
        cls,
        obs: Observables,
        value: float,
        bin: Optional[Tuple[float, float]] = None,
    ) -> "ObservableValue":
        """Create an observable value from a public observable enum.

        Args:
            obs: Public observable enum value.
            value: Numerical value to store.
            bin: Optional bin range ``(low, high)`` for binned observables.

        Returns:
            A normalized ``ObservableValue`` instance.

        Raises:
            TypeError: If ``obs`` is not an ``Observables`` enum value.
        """
        if not isinstance(obs, Observables):
            raise TypeError(f"obs doit être un Observables Python, reçu {type(obs)!r}.")
        return cls(id=ObservableMapper.to_id(obs), value=float(value), bin=bin)

    @classmethod
    def from_cpp(cls, cpp_obj) -> "ObservableValue":
        """Wrap a bound C++ ``ObservableValue`` object.

        Args:
            cpp_obj: Bound C++ object exposing ``id``, ``value`` and optional
                ``bin`` fields.

        Returns:
            Equivalent Python ``ObservableValue`` instance.
        """
        py_id = ObservableId(str(cpp_obj.id))
        b = cpp_obj.bin
        py_bin = None if b is None else (float(b[0]), float(b[1]))
        return cls(id=py_id, value=float(cpp_obj.value), bin=py_bin)

    def to_cpp(self):
        """Convert this value to the bound C++ representation.

        Returns:
            A C++ ``ObservableValue`` pybind11 object.
        """
        cpp_id = _cpp_observable_id(self.id)
        if self.bin is None:
            return _CppObservableValue(cpp_id, float(self.value))
        return _CppObservableValue(cpp_id, float(self.value), self.bin)


__all__ = [
    "ObservableValue",
]
