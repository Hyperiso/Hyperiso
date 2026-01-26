from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, List, Optional, Sequence, Tuple, Union, cast


from pyhyperiso.phyperiso.pyhyperiso import statistic as st

from pyhyperiso.core.Math.real_matrix import Matrix

NestedList = Sequence[Sequence[float]]
MatrixLike = Union[Matrix, NestedList]


def _to_matrix(obj: MatrixLike) -> Matrix:
    if isinstance(obj, Matrix):
        return obj
    return Matrix(data=obj)


class CopulaKind(Enum):
    GAUSSIAN = "GAUSSIAN"
    STUDENT_T = "STUDENT_T"

    def to_cpp(self) -> Any:
        return getattr(st.CopulaType, self.value)


@dataclass(frozen=True)
class GaussianCopulaConfigPy:
    R: MatrixLike

    def to_cpp(self) -> Any:
        Rm = _to_matrix(self.R)
        return st.GaussianCopulaConfig(Rm.to_cpp())

    @classmethod
    def from_cpp(cls, cpp: Any) -> "GaussianCopulaConfigPy":
        return cls(R=Matrix._from_cpp(cpp.R))


@dataclass(frozen=True)
class StudentTCopulaConfigPy:
    R: MatrixLike
    nu: int = 4

    def __post_init__(self) -> None:
        if int(self.nu) <= 0:
            raise ValueError("StudentTCopulaConfigPy.nu doit être > 0")

    def to_cpp(self) -> Any:
        Rm = _to_matrix(self.R)
        cfg = st.StudentTCopulaConfig()
        cfg.R = Rm.to_cpp()
        cfg.nu = int(self.nu)
        return cfg

    @classmethod
    def from_cpp(cls, cpp: Any) -> "StudentTCopulaConfigPy":
        return cls(R=Matrix._from_cpp(cpp.R), nu=int(cpp.nu))


CopulaConfigPy = Union[GaussianCopulaConfigPy, StudentTCopulaConfigPy]


def _config_from_cpp(cpp_cfg: Any) -> CopulaConfigPy:
    if isinstance(cpp_cfg, st.GaussianCopulaConfig):
        return GaussianCopulaConfigPy.from_cpp(cpp_cfg)
    if isinstance(cpp_cfg, st.StudentTCopulaConfig):
        return StudentTCopulaConfigPy.from_cpp(cpp_cfg)
    raise TypeError(f"Copula config C++ inconnue: {type(cpp_cfg)!r}")