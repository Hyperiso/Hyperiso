from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple, Union, Optional, overload

from pyhyperiso.phyperiso.pyhyperiso import math as ma


Number = Union[int, float]
NestedList = Sequence[Sequence[Number]]


@dataclass(frozen=True)
class SignedLogDet:
    """Pure-Python result for slogdet."""
    logdet: float
    sign: int


class Matrix:
    """
    Python wrapper for the C++ RealMatrix class (ma.matrix.RealMatrix).

    - All I/O is Python-native (lists, floats, tuples).
    - The underlying C++ object is kept internal in _cpp_obj.
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
        """
        Construct a Matrix from:
          - data: nested list/sequence (rows x cols), OR
          - rows, cols: empty matrix (initialized as in C++), OR
          - flat + rows + cols: flat row-major data

        Exactly one of (data) or (rows, cols) or (flat, rows, cols) should be used.
        """
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
        """Internal: wrap an existing C++ RealMatrix object."""
        inst = cls.__new__(cls)
        inst._cpp_obj = cpp_obj
        return inst

    @property
    def rows(self) -> int:
        return int(self._cpp_obj.rows)

    @property
    def cols(self) -> int:
        return int(self._cpp_obj.cols)

    @property
    def shape(self) -> Tuple[int, int]:
        r, c = self._cpp_obj.shape  # bound as tuple
        return int(r), int(c)

    def to_list(self) -> List[List[float]]:
        """Return a Python nested list (copy)."""
        return self._cpp_obj.to_list()

    def __repr__(self) -> str:
        r, c = self.shape
        return f"Matrix(shape=({r}, {c}))"

    def __getitem__(self, idx: Tuple[int, int]) -> float:
        i, j = idx
        return float(self._cpp_obj[(int(i), int(j))])

    def __setitem__(self, idx: Tuple[int, int], value: Number) -> None:
        i, j = idx
        self._cpp_obj[(int(i), int(j))] = float(value)

    def is_symmetric(self) -> bool:
        return bool(self._cpp_obj.is_symmetric())

    def T(self) -> "Matrix":
        """Transpose."""
        return Matrix._from_cpp(self._cpp_obj.transpose())

    def inv(self) -> "Matrix":
        return Matrix._from_cpp(self._cpp_obj.inv())

    def slogdet(self) -> SignedLogDet:
        s = self._cpp_obj.slogdet()  # C++ SignedLogDet
        return SignedLogDet(logdet=float(s.logdet), sign=int(s.sign))

    def eig(self) -> Tuple["Matrix", "Matrix"]:
        """
        Returns (D, P) where:
          - D: diagonal matrix of eigenvalues
          - P: eigenvectors matrix
        """
        es = self._cpp_obj.eig()  # C++ EigenSystem with fields D, P
        return Matrix._from_cpp(es.D), Matrix._from_cpp(es.P)

    def to_cpp(self):
        return self._cpp_obj

    @staticmethod
    def _as_cpp(other: Union["Matrix", NestedList]) -> "ma.matrix.RealMatrix":
        if isinstance(other, Matrix):
            return other._cpp_obj
        # allow nested list on the fly
        return ma.matrix.RealMatrix([[float(x) for x in row] for row in other])

    def __neg__(self) -> "Matrix":
        return Matrix._from_cpp(-self._cpp_obj)

    def __add__(self, other: Union["Matrix", NestedList]) -> "Matrix":
        return Matrix._from_cpp(self._cpp_obj + self._as_cpp(other))

    def __sub__(self, other: Union["Matrix", NestedList]) -> "Matrix":
        return Matrix._from_cpp(self._cpp_obj - self._as_cpp(other))

    def __mul__(self, other: Union["Matrix", NestedList, Number]) -> "Matrix":
        if isinstance(other, (int, float)):
            return Matrix._from_cpp(self._cpp_obj * float(other))
        return Matrix._from_cpp(self._cpp_obj * self._as_cpp(other))

    def __rmul__(self, other: Number) -> "Matrix":
        if not isinstance(other, (int, float)):
            return NotImplemented
        return Matrix._from_cpp(float(other) * self._cpp_obj)

    def __truediv__(self, other: Number) -> "Matrix":
        if not isinstance(other, (int, float)):
            return NotImplemented
        return Matrix._from_cpp(self._cpp_obj / float(other))

    def __matmul__(self, other: Union["Matrix", NestedList]) -> "Matrix":
        """Python @ operator mapped to C++ matrix-matrix multiplication."""
        return Matrix._from_cpp(self._cpp_obj * self._as_cpp(other))


def eye(n: int) -> Matrix:
    return Matrix._from_cpp(ma.matrix.eye(int(n)))

def nearest_psd(R: Union[Matrix, NestedList], thr: float = 1e-12) -> Matrix:
    cpp_R = R._cpp_obj if isinstance(R, Matrix) else ma.matrix.RealMatrix([[float(x) for x in row] for row in R])
    return Matrix._from_cpp(ma.matrix.nearest_psd(cpp_R, float(thr)))

def cholesky_L(R: Union[Matrix, NestedList]) -> Matrix:
    cpp_R = R._cpp_obj if isinstance(R, Matrix) else ma.matrix.RealMatrix([[float(x) for x in row] for row in R])
    return Matrix._from_cpp(ma.matrix.cholesky_L(cpp_R))

if __name__ == "__main__":
    
    A = Matrix([[1.0, 2.0],
            [3.0, 4.0]])

    A[0, 1] = 20.0
    print(A.to_list())

    B = eye(2)
    C = A + B
    D = A @ B 
    E = 2.0 * A

    logdet = A.slogdet()
    print(logdet.logdet, logdet.sign)

    Dmat, P = A.eig()
    print(Dmat.to_list())
    print(P.to_list())

    Rpsd = nearest_psd([[1, -2], [-2, 1]])
    L = cholesky_L(Rpsd)