from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from pyhyperiso.core.Common.BinnedObservableId import BinnedObservableId


@dataclass(frozen=True)
class ExperimentObs:
    """Wrapper Python opaque de ``ExperimentObs``.

    Le binding actuel retourne des ``ExperimentObs`` dans les maps statistiques,
    mais ne fournit pas de constructeur Python public dédié. On garde donc
    l'objet C++ seulement comme détail interne pour pouvoir le renvoyer au
    binding quand c'est nécessaire, tout en exposant le champ utile ``obs`` sous
    forme de ``BinnedObservableId`` Python.
    """

    obs: BinnedObservableId
    _cpp_obj: Any = field(default=None, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not isinstance(self.obs, BinnedObservableId):
            raise TypeError(f"obs doit être un BinnedObservableId Python, reçu {type(self.obs)!r}.")

    @classmethod
    def from_cpp(cls, cpp_obj) -> "ExperimentObs":
        inst = cls.__new__(cls)
        object.__setattr__(inst, "obs", BinnedObservableId.from_cpp(cpp_obj.obs))
        object.__setattr__(inst, "_cpp_obj", cpp_obj)
        return inst

    def to_cpp(self):
        if self._cpp_obj is None:
            raise RuntimeError(
                "ExperimentObs ne peut pas encore être construit côté Python pur : "
                "utilise un ExperimentObs renvoyé par StatisticInterface.get_obs_exp() "
                "ou expose un constructeur C++ dédié dans le binding."
            )
        return self._cpp_obj

    def __hash__(self) -> int:
        return hash(self.obs)

    def __str__(self) -> str:
        return str(self.obs)

    def __repr__(self) -> str:
        return f"ExperimentObs(obs={self.obs!r})"


__all__ = ["ExperimentObs"]
