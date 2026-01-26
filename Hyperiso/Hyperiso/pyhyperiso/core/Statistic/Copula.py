from pyhyperiso.core.Statistic.CopulaConfig import CopulaConfigPy, CopulaKind, _config_from_cpp, GaussianCopulaConfigPy, StudentTCopulaConfigPy, MatrixLike, Matrix
from typing import Any, List, Optional, Sequence, Tuple, Union, cast
from pyhyperiso.phyperiso.pyhyperiso import statistic as st

class Copula:
    __slots__ = ("_cpp_obj",)

    def __init__(self, cpp_obj: Any):
        self._cpp_obj = cpp_obj 

    @classmethod
    def from_cpp(cls, cpp_obj: Any) -> "Copula":
        if isinstance(cpp_obj, st.GaussianCopula):
            return GaussianCopula(cpp_obj)
        if isinstance(cpp_obj, st.StudentTCopula):
            return StudentTCopula(cpp_obj)
        return cls(cpp_obj)

    def sample_u(self, n: Optional[int] = None) -> Union[List[float], List[List[float]]]:
        if n is None:
            u = self._cpp_obj.sample_u()
            return [float(x) for x in u]
        if int(n) < 0:
            raise ValueError("n doit être >= 0")
        U = self._cpp_obj.sample_u(int(n))
        return [[float(x) for x in row] for row in U]

    def log_density(self, u: Sequence[float]) -> float:
        """log c(u)"""
        return float(self._cpp_obj.log_density([float(x) for x in u]))

    def density(self, u: Sequence[float]) -> float:
        return float(pow(2.718281828459045, self.log_density(u)))  # exp without np :p


class GaussianCopula(Copula):
    pass


class StudentTCopula(Copula):
    pass



class CopulaFactoryWrapper:
    @staticmethod
    def create(kind: CopulaKind, config: CopulaConfigPy, seed: Optional[int] = None) -> Copula:
        cpp_kind = kind.to_cpp()
        cpp_cfg = config.to_cpp()
        cpp_obj = st.CopulaFactory.create(cpp_kind, cpp_cfg, seed)
        return Copula.from_cpp(cpp_obj)

    @staticmethod
    def gaussian(R: MatrixLike, seed: Optional[int] = None) -> GaussianCopula:
        cop = CopulaFactoryWrapper.create(CopulaKind.GAUSSIAN, GaussianCopulaConfigPy(R=R), seed=seed)
        return cast(GaussianCopula, cop)

    @staticmethod
    def student_t(R: MatrixLike, nu: int = 4, seed: Optional[int] = None) -> StudentTCopula:
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