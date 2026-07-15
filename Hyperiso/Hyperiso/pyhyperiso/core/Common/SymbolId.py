"""Symbolic identifiers used by dynamic mappers and registries."""

from dataclasses import dataclass
from pyhyperiso.phyperiso.pyhyperiso.common import _CppObservableId
from pyhyperiso.phyperiso.pyhyperiso.common import _CppDecayId
from pyhyperiso.phyperiso.pyhyperiso.common import _CppWGroupId
from pyhyperiso.phyperiso.pyhyperiso.common import _CppWCoefId


@dataclass(frozen=True)
class _SymbolId:
    """Small Python wrapper around Hyperiso's dynamic symbol identifiers.

    Dynamic identifiers are string-backed, strongly typed C++ ids. They allow
    runtime observables, decays, Wilson groups and Wilson coefficients to live
    next to the legacy static enums without forcing a conversion back to an enum.
    """

    _name: str
    _cpp_type = None

    def __str__(self) -> str:
        return self._name

    @property
    def name(self) -> str:
        """Canonical name stored by the C++ id."""
        return self._name

    def _to_cpp(self):
        """Return the bound C++ id object."""
        return self._cpp_type(self._name)

    def to_cpp(self):
        """Public alias for :meth:`_to_cpp`."""
        return self._to_cpp()


class ObservableId(_SymbolId):
    """Dynamic id for builtin or runtime observables."""

    _cpp_type = _CppObservableId


class DecayId(_SymbolId):
    """Dynamic id for builtin or runtime decays."""

    _cpp_type = _CppDecayId


class WGroupId(_SymbolId):
    """Dynamic id for builtin or runtime Wilson groups."""

    _cpp_type = _CppWGroupId


class WCoefId(_SymbolId):
    """Dynamic id for builtin or runtime Wilson coefficients."""

    _cpp_type = _CppWCoefId


def _unwrap_optional(x, what="value"):
    if x is None:
        raise KeyError(f"{what} not found")
    return x
