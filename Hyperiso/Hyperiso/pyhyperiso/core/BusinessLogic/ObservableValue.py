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
    """Valeur d'observable côté Python.

    L'identifiant public est toujours le wrapper Python ``ObservableId``.
    Les ``Observables`` enum peuvent être convertis explicitement avec
    ``PyObservableValue.from_observable(...)`` ou ``ObservableMapper.to_id(...)``.
    """

    id: ObservableId
    value: float
    bin: Optional[Tuple[float, float]] = None

    def __post_init__(self) -> None:
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
        if not isinstance(obs, Observables):
            raise TypeError(f"obs doit être un Observables Python, reçu {type(obs)!r}.")
        return cls(id=ObservableMapper.to_id(obs), value=float(value), bin=bin)

    @classmethod
    def from_cpp(cls, cpp_obj) -> "ObservableValue":
        # Conversion interne binding -> wrapper Python.
        py_id = ObservableId(str(cpp_obj.id))
        b = cpp_obj.bin
        py_bin = None if b is None else (float(b[0]), float(b[1]))
        return cls(id=py_id, value=float(cpp_obj.value), bin=py_bin)

    def to_cpp(self):
        # Bridge interne wrapper Python -> binding.
        cpp_id = _cpp_observable_id(self.id)
        if self.bin is None:
            return _CppObservableValue(cpp_id, float(self.value))
        return _CppObservableValue(cpp_id, float(self.value), self.bin)


__all__ = [
    "ObservableValue",
]
