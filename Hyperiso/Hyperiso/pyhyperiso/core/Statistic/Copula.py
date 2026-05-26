"""Python wrappers around C++ copula objects.

The classes in this module provide a small, Pythonic API over the C++ statistical
copula implementations. Copulas operate on the unit hypercube: samples returned
by :meth:`Copula.sample_u` are dependent uniforms, and :meth:`Copula.log_density`
evaluates the copula density contribution ``log c(u)``.
"""

from pyhyperiso.core.Statistic.CopulaConfig import (
    CopulaConfigPy,
    CopulaKind,
    _config_from_cpp,
    GaussianCopulaConfigPy,
    StudentTCopulaConfigPy,
    MatrixLike,
    Matrix,
)
from typing import Any, List, Optional, Sequence, Union, cast
from pyhyperiso.phyperiso.pyhyperiso import statistic as st


class Copula:
    """Base Python wrapper for a bound C++ copula.

    The wrapper intentionally stores the C++ object as an opaque implementation
    detail. Use :class:`CopulaFactoryWrapper` to create instances from Python
    configuration objects.

    Args:
        cpp_obj: Bound C++ copula instance.
    """

    __slots__ = ("_cpp_obj",)

    def __init__(self, cpp_obj: Any):
        """Store the bound C++ copula object."""
        self._cpp_obj = cpp_obj

    @classmethod
    def from_cpp(cls, cpp_obj: Any) -> "Copula":
        """Wrap a bound C++ copula with the most specific Python class.

        Args:
            cpp_obj: Bound C++ copula instance.

        Returns:
            ``GaussianCopula`` for Gaussian C++ objects,
            ``StudentTCopula`` for Student-t C++ objects, otherwise the generic
            ``Copula`` wrapper.
        """
        if isinstance(cpp_obj, st.GaussianCopula):
            return GaussianCopula(cpp_obj)
        if isinstance(cpp_obj, st.StudentTCopula):
            return StudentTCopula(cpp_obj)
        return cls(cpp_obj)

    def sample_u(self, n: Optional[int] = None) -> Union[List[float], List[List[float]]]:
        """Draw one or more dependent uniform samples.

        Args:
            n: Number of samples. When omitted, returns a single vector of
                dependent uniforms. When provided, returns a list of ``n``
                vectors.

        Returns:
            A single ``list[float]`` if ``n`` is ``None``; otherwise a
            ``list[list[float]]`` with shape ``n x dim``.

        Raises:
            ValueError: If ``n`` is negative.

        Examples:
            >>> cop = CopulaFactoryWrapper.gaussian([[1.0, 0.5], [0.5, 1.0]], seed=123)
            >>> u = cop.sample_u()
            >>> len(u)
            2
            >>> U = cop.sample_u(3)
            >>> len(U)
            3
        """
        if n is None:
            u = self._cpp_obj.sample_u()
            return [float(x) for x in u]
        if int(n) < 0:
            raise ValueError("n doit être >= 0")
        U = self._cpp_obj.sample_u(int(n))
        return [[float(x) for x in row] for row in U]

    def log_density(self, u: Sequence[float]) -> float:
        """Evaluate the copula log-density at a point in the unit hypercube.

        Args:
            u: Uniform coordinates. Each component is expected to lie in
                ``[0, 1]`` and the length must match the copula dimension.

        Returns:
            The scalar value ``log c(u)``.
        """
        return float(self._cpp_obj.log_density([float(x) for x in u]))

    def density(self, u: Sequence[float]) -> float:
        """Evaluate the copula density at a point in the unit hypercube.

        Args:
            u: Uniform coordinates passed to :meth:`log_density`.

        Returns:
            The density value ``c(u)``.
        """
        return float(pow(2.718281828459045, self.log_density(u)))


class GaussianCopula(Copula):
    """Wrapper for the Gaussian copula backend.

    The underlying C++ implementation samples a correlated Gaussian latent
    vector and maps it component-wise through the standard normal CDF.
    """

    pass


class StudentTCopula(Copula):
    """Wrapper for the Student-t copula backend.

    The Student-t copula follows the same interface as the Gaussian copula but
    can model stronger tail dependence through its degrees of freedom.
    """

    pass


class CopulaFactoryWrapper:
    """Factory helpers for creating Python copula wrappers.

    This class forwards to the C++ ``CopulaFactory`` and wraps the returned C++
    object in the appropriate Python class.
    """

    @staticmethod
    def create(kind: CopulaKind, config: CopulaConfigPy, seed: Optional[int] = None) -> Copula:
        """Create a copula from an explicit family and configuration.

        Args:
            kind: Copula family to instantiate.
            config: Python configuration object compatible with ``kind``.
            seed: Optional random seed passed to the C++ RNG.

        Returns:
            A Python wrapper around the newly created C++ copula.
        """
        cpp_kind = kind.to_cpp()
        cpp_cfg = config.to_cpp()
        cpp_obj = st.CopulaFactory.create(cpp_kind, cpp_cfg, seed)
        return Copula.from_cpp(cpp_obj)

    @staticmethod
    def gaussian(R: MatrixLike, seed: Optional[int] = None) -> GaussianCopula:
        """Create a Gaussian copula from a correlation matrix.

        Args:
            R: Square correlation matrix.
            seed: Optional random seed.

        Returns:
            A ``GaussianCopula`` instance.

        Examples:
            >>> cop = CopulaFactoryWrapper.gaussian([[1.0, 0.3], [0.3, 1.0]], seed=7)
            >>> isinstance(cop.sample_u(), list)
            True
        """
        cop = CopulaFactoryWrapper.create(CopulaKind.GAUSSIAN, GaussianCopulaConfigPy(R=R), seed=seed)
        return cast(GaussianCopula, cop)

    @staticmethod
    def student_t(R: MatrixLike, nu: int = 4, seed: Optional[int] = None) -> StudentTCopula:
        """Create a Student-t copula from a correlation matrix.

        Args:
            R: Square correlation matrix.
            nu: Degrees of freedom of the latent Student-t vector.
            seed: Optional random seed.

        Returns:
            A ``StudentTCopula`` instance.
        """
        cop = CopulaFactoryWrapper.create(CopulaKind.STUDENT_T, StudentTCopulaConfigPy(R=R, nu=nu), seed=seed)
        return cast(StudentTCopula, cop)


__all__ = [
    "CopulaKind",
    "GaussianCopulaConfigPy",
    "StudentTCopulaConfigPy",
    "CopulaConfigPy",
    "Copula",
    "GaussianCopula",
    "StudentTCopula",
    "CopulaFactoryWrapper",
    "_config_from_cpp",
]
