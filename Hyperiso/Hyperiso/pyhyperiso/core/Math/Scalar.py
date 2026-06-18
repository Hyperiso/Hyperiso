from pyhyperiso.phyperiso.pyhyperiso import math as ma
import math
import cmath
from typing import  Union


class Scalar:
    """Python wrapper for the C++ scalar_t class.

    scalar_t is an enhanced complex number type that supports real projection,
    full arithmetic, and compatibility with Python's numeric types.
    """

    def __init__(self, re: float = 0.0, im: float = 0.0):
        """Initializes a Scalar from real and imaginary parts.

        Args:
            re (float): Real part. Defaults to 0.0.
            im (float): Imaginary part. Defaults to 0.0.
        """
        self._cpp_obj = ma.scalar_t(re, im)

    @classmethod
    def from_complex(cls, z: complex) -> "Scalar":
        """Creates a Scalar from a complex number.

        Args:
            z (complex): A Python complex.

        Returns:
            Scalar: A wrapped scalar_t.
        """
        instance = cls()
        instance._cpp_obj = ma.scalar_t(z)
        return instance

    def real(self) -> float:
        """Returns the real part.

        ``Scalar.from_cpp`` can receive either a bound ``scalar_t`` object or a
        native Python ``complex`` depending on the pybind overload that produced
        the value.  The C++ object exposes ``real()`` as a method, while Python
        complex exposes ``real`` as a float attribute; both forms are supported.

        Returns:
            float: Real component.
        """
        real_attr = getattr(self._cpp_obj, "real")
        return float(real_attr() if callable(real_attr) else real_attr)

    def imag(self) -> float:
        """Returns the imaginary part.

        Supports both bound ``scalar_t.imag()`` and Python ``complex.imag``.

        Returns:
            float: Imaginary component.
        """
        imag_attr = getattr(self._cpp_obj, "imag")
        return float(imag_attr() if callable(imag_attr) else imag_attr)

    def to_double(self) -> float:
        """Casts to float using the real part.

        Bound ``scalar_t`` still uses its native ``to_double()`` implementation.
        Native Python numbers/complex values use ``real()``.

        Returns:
            float: Real part.
        """
        to_double = getattr(self._cpp_obj, "to_double", None)
        if callable(to_double):
            return float(to_double())
        return self.real()

    def __float__(self) -> float:
        """Casts to float, calling to_double.

        Returns:
            float: Real part.
        """
        return self.to_double()

    def __complex__(self) -> complex:
        """Casts to Python complex.

        Returns:
            complex: Native Python complex value.
        """
        return complex(self.real(), self.imag())

    def __add__(self, other: "Scalar") -> "Scalar":
        result = Scalar()
        result._cpp_obj = self._cpp_obj + other._cpp_obj
        return result

    def __sub__(self, other: "Scalar") -> "Scalar":
        result = Scalar()
        result._cpp_obj = self._cpp_obj - other._cpp_obj
        return result

    def __mul__(self, other: "Scalar") -> "Scalar":
        result = Scalar()
        result._cpp_obj = self._cpp_obj * other._cpp_obj
        return result

    def __truediv__(self, other: "Scalar") -> "Scalar":
        result = Scalar()
        result._cpp_obj = self._cpp_obj / other._cpp_obj
        return result

    def __neg__(self) -> "Scalar":
        result = Scalar()
        result._cpp_obj = -self._cpp_obj
        return result

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Scalar):
            return False
        return self.real() == other.real() and self.imag() == other.imag()

    @classmethod
    def from_cpp(cls, cpp_obj) -> "Scalar":
        """Wraps an existing C++ scalar_t object."""
        instance = cls()
        instance._cpp_obj = cpp_obj
        return instance
    
    def __repr__(self) -> str:
        return f"Scalar({self.real()}, {self.imag()})"


def _to_scalar(value: Union[float, complex, Scalar]) -> Scalar:
    """Helper to convert float/complex to Scalar."""
    if isinstance(value, Scalar):
        return value
    elif isinstance(value, (float, int)):
        return Scalar(value)
    elif isinstance(value, complex):
        return Scalar.from_complex(value)
    else:
        raise TypeError(f"Expected float, complex or Scalar, got {type(value)}")
    

# -------- Math function wrappers -------- #

def _wrap_math_func(cpp_func):
    def wrapped(x: Scalar) -> Scalar:
        result = Scalar()
        result._cpp_obj = cpp_func(x._cpp_obj)
        return result
    return wrapped

# Map to C++ functions
sqrt = _wrap_math_func(ma.sqrt)
sin = _wrap_math_func(ma.sin)
cos = _wrap_math_func(ma.cos)
tan = _wrap_math_func(ma.tan)
asin = _wrap_math_func(ma.asin)
acos = _wrap_math_func(ma.acos)
atan = _wrap_math_func(ma.atan)
exp = _wrap_math_func(ma.exp)
log = _wrap_math_func(ma.log)
sinh = _wrap_math_func(ma.sinh)
cosh = _wrap_math_func(ma.cosh)
tanh = _wrap_math_func(ma.tanh)
abs_scalar = _wrap_math_func(ma.abs)
arg = _wrap_math_func(ma.arg)
norm = _wrap_math_func(ma.norm)

def pow_scalar(base: Scalar, exponent: Union[Scalar, float, int]) -> Scalar:
    """Raises a Scalar to a scalar/int/float exponent.

    Args:
        base (Scalar): The base.
        exponent (Union[Scalar, float, int]): The exponent.

    Returns:
        Scalar: Result of exponentiation.
    """
    if isinstance(exponent, Scalar):
        cpp = ma.pow(base._cpp_obj, exponent._cpp_obj)
    elif isinstance(exponent, (float, int)):
        cpp = ma.pow(base._cpp_obj, float(exponent))
    else:
        raise TypeError("Exponent must be Scalar, float, or int")
    result = Scalar()
    result._cpp_obj = cpp
    return result
