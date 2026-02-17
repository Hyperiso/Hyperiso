from dataclasses import dataclass
from pyhyperiso.phyperiso.pyhyperiso.common import _CppObservableId
from pyhyperiso.phyperiso.pyhyperiso.common import _CppDecayId
@dataclass(frozen=True)
class _SymbolId:
    _name: str
    _cpp_type = None

    def __str__(self) -> str:
        return self._name

    @property
    def name(self) -> str:
        return self._name

    def _to_cpp(self):
        return self._cpp_type(self._name)
    
class ObservableId(_SymbolId):
    _cpp_type = _CppObservableId

class DecayId(_SymbolId):
    _cpp_type = _CppDecayId
    
def _unwrap_optional(x, what="value"):
    if x is None:
        raise KeyError(f"{what} not found")
    return x