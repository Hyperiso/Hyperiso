"""Dense real-matrix utilities backed by the C++ ``RealMatrix`` type.

This module provides a Python-friendly wrapper around the C++ dense matrix
implementation used throughout the statistics and copula stack. Inputs and
outputs are converted to native Python containers, while heavy linear-algebra
operations are delegated to the C++ backend.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple, Union

from pyhyperiso.phyperiso.pyhyperiso import math as ma


Number = Union[int, float]
NestedList = Sequence[Sequence[Number]]


@dataclass(frozen=True)
class SignedLogDet:
    """Signed logarithmic determinant.

    Attributes:
        logdet: Natural logarithm of the determinant absolute value.
        sign: Determinant sign, usually ``-1``, ``0``, or ``1``.
    """

    logdet: float
    sign: int


class Matrix:
    """Python wrapper for the C++ ``RealMatrix`` class.

    The wrapper accepts common Python matrix representations and exposes Python
    arithmetic operators. Matrix-matrix multiplication is also mapped to the
    ``@`` operator.

    Args:
        data: Optional nested sequence representing a dense matrix.
        rows: Number of rows for an empty matrix or flat initialization.
        cols: Number of columns for an empty matrix or flat initialization.
        flat: Optional row-major flat values. Requires ``rows`` and ``cols``.

    Raises:
        ValueError: If constructor arguments are incomplete or inconsistent.

    Example:
        >>> A = Matrix([[1.0, 2.0], [3.0, 4.0]])
        >>> I = eye(2)
        >>> (A @ I).to_list()
        [[1.0, 2.0], [3.0, 4.0]]
    """

    __slots__ = ("_cpp_obj",)

    def __init__(
        self,
        data: Optional[NestedList] = None,
        *,
        rows: Optional[int] = None,
        cols: Optional[int] = None,
        flat: Optional[Sequence[Number]] = None,
    ):
        if data is not None:
            self._cpp_obj = ma.matrix.RealMatrix([[float(x) for x in row] for row in data])
            return

        if flat is not None:
            if rows is None or cols is None:
                raise ValueError("When using flat=..., you must also provide rows=... and cols=...")
            self._cpp_obj = ma.matrix.RealMatrix([float(x) for x in flat], int(rows), int(cols))
            return

        if rows is not None or cols is not None:
            if rows is None or cols is None:
                raise ValueError("Provide both rows and cols.")
            self._cpp_obj = ma.matrix.RealMatrix(int(rows), int(cols))
            return

        self._cpp_obj = ma.matrix.RealMatrix()

    @classmethod
    def _from_cpp(cls, cpp_obj) -> "Matrix":
        """Wrap an existing C++ ``RealMatrix`` object.

        Args:
            cpp_obj: Bound C++ matrix object.

        Returns:
            Matrix: Python wrapper sharing the provided C++ object.
        """
        inst = cls.__new__(cls)
        inst._cpp_obj = cpp_obj
        return inst

    @classmethod
    def from_cpp(cls, cpp_obj) -> "Matrix":
        """Public alias for wrapping a C++ ``RealMatrix`` object."""
        return cls._from_cpp(cpp_obj)

    @property
    def rows(self) -> int:
        """Number of matrix rows."""
        return int(self._cpp_obj.rows)

    @property
    def cols(self) -> int:
        """Number of matrix columns."""
        return int(self._cpp_obj.cols)

    @property
    def shape(self) -> Tuple[int, int]:
        """Matrix shape as ``(rows, cols)``."""
        r, c = self._cpp_obj.shape
        return int(r), int(c)

    def to_list(self) -> List[List[float]]:
        """Return a nested Python list copy of the matrix.

        Returns:
            list[list[float]]: Row-major matrix values.
        """
        return self._cpp_obj.to_list()

    def __repr__(self) -> str:
        """Return a compact representation containing the matrix shape."""
        r, c = self.shape
        return f"Matrix(shape=({r}, {c}))"

    def __getitem__(self, idx: Tuple[int, int]) -> float:
        """Return one matrix entry using ``matrix[i, j]`` indexing."""
        i, j = idx
        return float(self._cpp_obj[(int(i), int(j))])

    def __setitem__(self, idx: Tuple[int, int], value: Number) -> None:
        """Set one matrix entry using ``matrix[i, j] = value`` indexing."""
        i, j = idx
        self._cpp_obj[(int(i), int(j))] = float(value)

    def is_symmetric(self) -> bool:
        """Return whether the C++ backend considers the matrix symmetric."""
        return bool(self._cpp_obj.is_symmetric())

    def T(self) -> "Matrix":
        """Return the transpose of the matrix.

        Returns:
            Matrix: Transposed matrix.
        """
        return Matrix._from_cpp(self._cpp_obj.transpose())

    def inv(self) -> "Matrix":
        """Return the matrix inverse.

        Returns:
            Matrix: Inverse matrix computed by C++.

        Raises:
            RuntimeError: Propagated from the C++ backend if inversion fails.
        """
        return Matrix._from_cpp(self._cpp_obj.inv())

    def slogdet(self) -> SignedLogDet:
        """Compute the signed logarithmic determinant.

        Returns:
            SignedLogDet: Determinant sign and log-absolute determinant.
        """
        s = self._cpp_obj.slogdet()
        return SignedLogDet(logdet=float(s.logdet), sign=int(s.sign))

    def eig(self) -> Tuple["Matrix", "Matrix"]:
        """Compute the eigensystem of the matrix.

        Returns:
            tuple[Matrix, Matrix]: ``(D, P)``, where ``D`` is a diagonal matrix
            of eigenvalues and ``P`` contains eigenvectors as returned by C++.
        """
        es = self._cpp_obj.eig()
        return Matrix._from_cpp(es.D), Matrix._from_cpp(es.P)

    def to_cpp(self):
        """Return the wrapped C++ matrix object."""
        return self._cpp_obj

    @staticmethod
    def _as_cpp(other: Union["Matrix", NestedList]):
        """Convert a Python matrix-like object to a C++ ``RealMatrix``."""
        if isinstance(other, Matrix):
            return other._cpp_obj
        return ma.matrix.RealMatrix([[float(x) for x in row] for row in other])

    def __neg__(self) -> "Matrix":
        """Return ``-self``."""
        return Matrix._from_cpp(-self._cpp_obj)

    def __add__(self, other: Union["Matrix", NestedList]) -> "Matrix":
        """Return matrix addition with another matrix-like object."""
        return Matrix._from_cpp(self._cpp_obj + self._as_cpp(other))

    def __sub__(self, other: Union["Matrix", NestedList]) -> "Matrix":
        """Return matrix subtraction with another matrix-like object."""
        return Matrix._from_cpp(self._cpp_obj - self._as_cpp(other))

    def __mul__(self, other: Union["Matrix", NestedList, Number]) -> "Matrix":
        """Return C++ multiplication by a scalar or matrix-like object."""
        if isinstance(other, (int, float)):
            return Matrix._from_cpp(self._cpp_obj * float(other))
        return Matrix._from_cpp(self._cpp_obj * self._as_cpp(other))

    def __rmul__(self, other: Number) -> "Matrix":
        """Return scalar multiplication from the left."""
        if not isinstance(other, (int, float)):
            return NotImplemented
        return Matrix._from_cpp(float(other) * self._cpp_obj)

    def __truediv__(self, other: Number) -> "Matrix":
        """Return division by a scalar."""
        if not isinstance(other, (int, float)):
            return NotImplemented
        return Matrix._from_cpp(self._cpp_obj / float(other))

    def __matmul__(self, other: Union["Matrix", NestedList]) -> "Matrix":
        """Return matrix multiplication using the Python ``@`` operator."""
        return Matrix._from_cpp(self._cpp_obj * self._as_cpp(other))


def eye(n: int) -> Matrix:
    """Create an identity matrix.

    Args:
        n: Matrix dimension.

    Returns:
        Matrix: ``n x n`` identity matrix.
    """
    return Matrix._from_cpp(ma.matrix.eye(int(n)))


def nearest_psd(R: Union[Matrix, NestedList], thr: float = 1e-12) -> Matrix:
    """Project a matrix to the nearest positive semidefinite matrix.

    This helper is used by copula construction to regularize correlation
    matrices before decomposition.

    Args:
        R: Input matrix or nested list.
        thr: Numerical threshold used by the C++ routine.

    Returns:
        Matrix: Regularized positive-semidefinite matrix.
    """
    cpp_R = (
        R._cpp_obj
        if isinstance(R, Matrix)
        else ma.matrix.RealMatrix([[float(x) for x in row] for row in R])
    )
    return Matrix._from_cpp(ma.matrix.nearest_psd(cpp_R, float(thr)))


def cholesky_L(R: Union[Matrix, NestedList]) -> Matrix:
    """Compute the lower Cholesky factor of a matrix.

    Args:
        R: Input symmetric positive-semidefinite matrix.

    Returns:
        Matrix: Lower triangular Cholesky factor as returned by C++.
    """
    cpp_R = (
        R._cpp_obj
        if isinstance(R, Matrix)
        else ma.matrix.RealMatrix([[float(x) for x in row] for row in R])
    )
    return Matrix._from_cpp(ma.matrix.cholesky_L(cpp_R))


__all__ = ["Matrix", "SignedLogDet", "eye", "nearest_psd", "cholesky_L"]
