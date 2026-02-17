from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Set, Tuple, Union

from pyhyperiso.phyperiso.pyhyperiso.observable import ObservableValue as _CppObservableValue
from pyhyperiso.phyperiso.pyhyperiso.common import ObservableMapper as _CppObservableMapper
from pyhyperiso.core.Common.GeneralEnum import Observables, Decays, QCDOrder
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Math.scalar import Scalar


ObsLike = Union[Observables, str]

def _key_to_cpp_obs_arg(key: ObsLike):
    """Retourne soit un _CppObservables (pour l'overload Observables),
    soit un _CppObservableId (pour l'overload ObservableId)."""
    print("key : ", key)
    if isinstance(key, Observables):
        return key.value  # -> _CppObservables
    if isinstance(key, str):
        return _CppObservableMapper.id_of(key)  # -> _CppObservableId
    raise TypeError(f"ObservableKey doit être Observables ou str, reçu: {type(key)}")


def _cpp_id_to_python_key(cpp_obs_id) -> ObsLike:
    """Convertit un ObservableId C++ en clé user-friendly:
       - Observables si c'est un builtin ET que le nom matche
       - sinon str canonical."""
    name = str(cpp_obs_id)
    try:
        return Observables[name]
    except KeyError:
        return name
    
    
@dataclass(frozen=True)
class PyObservableValue:
    """Wrapper Python de ObservableValue (C++)."""
    id: Observables
    value: float
    bin: Optional[Tuple[float, float]] = None

    @staticmethod
    def from_cpp(cpp: _CppObservableValue) -> "PyObservableValue":
        # cpp.id est un ObservableId (C++)
        py_id = _cpp_id_to_python_key(cpp.id)

        b = cpp.bin  # None ou (low, high) grâce au binding property
        py_bin = None if b is None else (float(b[0]), float(b[1]))

        return PyObservableValue(
            id=py_id,
            value=float(cpp.value),
            bin=py_bin,
        )

    def to_cpp(self) -> _CppObservableValue:
        """Optionnel: si tu veux reconstruire un ObservableValue C++ (rarement nécessaire)."""
        cpp_id = _key_to_cpp_obs_arg(self.id)
        # Attention: _key_to_cpp_obs_arg renvoie un enum C++ pour Observables,
        # alors que ObservableValue ctor attend ObservableId.
        # Donc pour Observables, on passe par le nom -> id_of.
        if isinstance(self.id, Observables):
            name = self.id.name
            cpp_id = _CppObservableMapper.id_of(name)

        if self.bin is None:
            return _CppObservableValue(cpp_id, float(self.value))
        return _CppObservableValue(cpp_id, float(self.value), (float(self.bin[0]), float(self.bin[1])))
    
if __name__ == "__main__":
    value = PyObservableValue(Observables.R_D, 2, {2,3})
    
    print(value)