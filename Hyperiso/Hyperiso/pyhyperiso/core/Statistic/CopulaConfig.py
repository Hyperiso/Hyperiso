"""Configuration objects used to build statistical copulas.

This module contains small Python dataclasses mirroring the C++ copula
configuration structs exposed through the ``statistic`` binding. They are meant
for user-facing Python code and convert themselves to the bound C++ objects when
calling factories such as :class:`~pyhyperiso.core.Statistic.Copula.CopulaFactoryWrapper`.

A copula models dependence between marginal distributions in a joint law. In the
C++ backend, the Gaussian and Student-t copulas are parameterized by a
correlation matrix ``R``. The Student-t copula additionally uses a number of
degrees of freedom ``nu``.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, Sequence, Union

from pyhyperiso.phyperiso.pyhyperiso import statistic as st
from pyhyperiso.core.Math.RealMatrix import Matrix

NestedList = Sequence[Sequence[float]]
MatrixLike = Union[Matrix, NestedList]


def _to_matrix(obj: MatrixLike) -> Matrix:
    """Normalize a Python matrix-like object into a ``Matrix`` wrapper.

    Args:
        obj: Either an existing ``Matrix`` instance or a nested sequence of
            numbers with shape ``d x d``.

    Returns:
        A ``Matrix`` instance ready to be converted to C++.

    Raises:
        Exception: Propagates validation errors raised by ``Matrix`` when the
            nested sequence cannot be interpreted as a numeric matrix.
    """
    if isinstance(obj, Matrix):
        return obj
    return Matrix(data=obj)


class CopulaKind(Enum):
    """Available copula families supported by the statistic backend.

    Attributes:
        GAUSSIAN: Gaussian copula based on a multivariate normal latent vector.
        STUDENT_T: Student-t copula, useful when stronger tail dependence is
            desired.
    """

    GAUSSIAN = "GAUSSIAN"
    STUDENT_T = "STUDENT_T"

    def to_cpp(self) -> Any:
        """Convert the Python enum value to the bound C++ ``CopulaType``.

        Returns:
            The corresponding C++ enum value from ``statistic.CopulaType``.
        """
        return getattr(st.CopulaType, self.value)


@dataclass(frozen=True)
class GaussianCopulaConfigPy:
    """Configuration for a Gaussian copula.

    Args:
        R: Correlation matrix of the latent Gaussian vector. The matrix should
            be square, symmetric and have diagonal entries close to one. The C++
            implementation regularizes/projects it before decomposition.

    Examples:
        >>> cfg = GaussianCopulaConfigPy(R=[[1.0, 0.4], [0.4, 1.0]])
        >>> cpp_cfg = cfg.to_cpp()
    """

    R: MatrixLike

    def to_cpp(self) -> Any:
        """Build the bound C++ ``GaussianCopulaConfig`` object.

        Returns:
            A C++ configuration object containing the correlation matrix.
        """
        Rm = _to_matrix(self.R)
        return st.GaussianCopulaConfig(Rm.to_cpp())

    @classmethod
    def from_cpp(cls, cpp: Any) -> "GaussianCopulaConfigPy":
        """Create a Python config from a bound C++ config.

        Args:
            cpp: Bound C++ ``GaussianCopulaConfig`` instance.

        Returns:
            The corresponding Python wrapper.
        """
        return cls(R=Matrix._from_cpp(cpp.R))


@dataclass(frozen=True)
class StudentTCopulaConfigPy:
    """Configuration for a Student-t copula.

    Args:
        R: Correlation matrix of the latent Student-t vector.
        nu: Degrees of freedom. Smaller values increase tail dependence. The C++
            implementation requires a physically meaningful positive value; in
            practice values such as ``4`` or larger are common starting points.

    Raises:
        ValueError: If ``nu`` is not strictly positive.

    Examples:
        >>> cfg = StudentTCopulaConfigPy(R=[[1.0, 0.2], [0.2, 1.0]], nu=5)
        >>> cfg.nu
        5
    """

    R: MatrixLike
    nu: int = 4

    def __post_init__(self) -> None:
        """Validate the number of degrees of freedom."""
        if int(self.nu) <= 0:
            raise ValueError("StudentTCopulaConfigPy.nu must be > 0")

    def to_cpp(self) -> Any:
        """Build the bound C++ ``StudentTCopulaConfig`` object.

        Returns:
            A C++ configuration object containing ``R`` and ``nu``.
        """
        Rm = _to_matrix(self.R)
        cfg = st.StudentTCopulaConfig()
        cfg.R = Rm.to_cpp()
        cfg.nu = int(self.nu)
        return cfg

    @classmethod
    def from_cpp(cls, cpp: Any) -> "StudentTCopulaConfigPy":
        """Create a Python config from a bound C++ Student-t config.

        Args:
            cpp: Bound C++ ``StudentTCopulaConfig`` instance.

        Returns:
            The corresponding Python wrapper.
        """
        return cls(R=Matrix._from_cpp(cpp.R), nu=int(cpp.nu))


CopulaConfigPy = Union[GaussianCopulaConfigPy, StudentTCopulaConfigPy]


def _config_from_cpp(cpp_cfg: Any) -> CopulaConfigPy:
    """Dispatch a bound C++ copula config to the matching Python wrapper.

    Args:
        cpp_cfg: Bound C++ copula configuration object.

    Returns:
        A ``GaussianCopulaConfigPy`` or ``StudentTCopulaConfigPy`` instance.

    Raises:
        TypeError: If the C++ object type is not recognized.
    """
    if isinstance(cpp_cfg, st.GaussianCopulaConfig):
        return GaussianCopulaConfigPy.from_cpp(cpp_cfg)
    if isinstance(cpp_cfg, st.StudentTCopulaConfig):
        return StudentTCopulaConfigPy.from_cpp(cpp_cfg)
    raise TypeError(f"Copula config C++ inconnue: {type(cpp_cfg)!r}")
