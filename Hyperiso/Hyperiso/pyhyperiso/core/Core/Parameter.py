from pyhyperiso.phyperiso.pyhyperiso.core import Parameter as _CppParameter
from pyhyperiso.phyperiso.pyhyperiso.core import ParameterMode as _CppParameterMode

from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Common.GeneralEnum import ParameterType
from pyhyperiso.core.Math.scalar import Scalar, pow_scalar, sqrt, _to_scalar
from pyhyperiso.phyperiso.pyhyperiso.core import DependentParameter as _CppDependentParameter
from typing import Dict, Tuple, Union, List
from enum import Enum

class ParameterMode(Enum):
    FIXED = _CppParameterMode.FIXED
    SHIFTABLE = _CppParameterMode.SHIFTABLE
    
class PyParameter:
    def __init__(
        self,
        pid: ParamId,
        expected: Union[float, complex, Scalar],
        stat: Union[float, complex, Scalar],
        syst: Union[float, complex, Scalar],
        scale: float = None,
        bin: List = None
    ):
        self._cpp_obj = _CppParameter(
            pid._cpp_obj,
            _to_scalar(expected)._cpp_obj,
            _to_scalar(stat)._cpp_obj,
            _to_scalar(syst)._cpp_obj,
        )
        if scale is not None:
            self._cpp_obj.set_scale(scale)
        if bin is not None:
            self._cpp_obj.set_bin(bin)

    @property
    def value(self) -> Scalar:
        return Scalar.from_cpp(self._cpp_obj.get_val())

    @property
    def std(self) -> Tuple[Scalar, Scalar]:
        stat, syst = self._cpp_obj.get_std()
        return Scalar.from_cpp(stat), Scalar.from_cpp(syst)

    @property
    def combined_std(self) -> Scalar:
        return Scalar.from_cpp(self._cpp_obj.get_combined_std())

    @property
    def pid(self) -> ParamId:
        return ParamId.from_cpp(self._cpp_obj.get_id())

    def set_expected(self, val: Union[float, complex, Scalar]):
        self._cpp_obj.set_expected(_to_scalar(val)._cpp_obj)

    def set_std(self, stat: Union[float, complex, Scalar], syst: Union[float, complex, Scalar]):
        self._cpp_obj.set_std(_to_scalar(stat)._cpp_obj, _to_scalar(syst)._cpp_obj)

    def set_shift(self, shift: Union[float, complex, Scalar]):
        self._cpp_obj.set_shift(_to_scalar(shift)._cpp_obj)

    def set_mode(self, mode: ParameterMode):
        self._cpp_obj.set_mode(mode.value)

    def __repr__(self):
        return (
            f"<PyParameter pid={self.pid}, value={self.value}, "
            f"std={self.std}, combined_std={self.combined_std}>"
        )
    
    @staticmethod
    def from_cpp(cpp_obj):
        param = PyParameter(ParamId(), 0, 0, 0)
        param._cpp_obj = cpp_obj
        return param

class PyDependentParameter(PyParameter):
    """Python wrapper around DependentParameter."""

    def __init__(self, cpp_obj: _CppDependentParameter):
        super().__init__(cpp_obj)
        self._cpp_obj = cpp_obj

    def depends_on(self, pid: ParamId) -> bool:
        return self._cpp_obj.dependsOn(pid._cpp_obj)

    def init(self):
        self._cpp_obj.init()

    def update(self):
        self._cpp_obj.update()

    def freeze(self):
        self._cpp_obj.freeze()

    def unfreeze(self):
        self._cpp_obj.unfreeze()
        
        
if __name__ == "__main__":
    print("🧪 Scalar operations")

    a = Scalar(3.0)
    b = Scalar(4.0)
    c = a + b
    print(f"{a} + {b} = {c}")
    print("Real part of c:", c.real())
    print("Imaginary part of c:", c.imag())
    print("As float:", float(c))
    print("As complex:", complex(c))

    d = pow_scalar(a, 2)
    print(f"{a}^2 =", d)

    e = sqrt(b)
    print(f"sqrt({b}) =", e)
    
    pid = ParamId(type=ParameterType.SM, block="MASS", code=25)
    param = PyParameter(pid, 125.0, 0.1, 0.2)

    print("Initial:", param)
    param.set_expected(126.5)
    param.set_std(0.05, 0.15)
    param.set_expected(130)
    print("Updated:", param)
    param.set_mode(ParameterMode.SHIFTABLE)
    param.set_shift(0.1)
    print("Updated by shift:", param)