"""Python wrappers for C++ parameters."""

from __future__ import annotations

from enum import Enum
from typing import List, Optional, Tuple, Union

from pyhyperiso.phyperiso.pyhyperiso.core import DependentParameter as _CppDependentParameter
from pyhyperiso.phyperiso.pyhyperiso.core import Parameter as _CppParameter
from pyhyperiso.phyperiso.pyhyperiso.core import ParameterMode as _CppParameterMode
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Math.Scalar import Scalar, _to_scalar


class ParameterMode(Enum):
    """Update mode used by shiftable C++ parameters.

    Attributes:
        FIXED: The parameter value is fixed to its expected value.
        SHIFTABLE: The parameter can be shifted through the C++ shift machinery.
    """

    FIXED = _CppParameterMode.FIXED
    SHIFTABLE = _CppParameterMode.SHIFTABLE


class Parameter:
    """Python wrapper around the C++ ``Parameter`` object.

    Parameters store a central value and statistical/systematic uncertainties.
    They may also carry optional metadata such as a scale or bin range. Numeric
    inputs are converted to the package ``Scalar`` type before crossing the
    pybind boundary.

    Args:
        pid: Parameter identifier.
        expected: Central value.
        stat: Statistical uncertainty.
        syst: Systematic uncertainty.
        scale: Optional scale metadata.
        bin: Optional bin metadata forwarded to the C++ object.

    Examples:
        >>> pid = ParamId(type=ParameterType.SM, block="MASS", code=25)
        >>> p = Parameter(pid, expected=125.0, stat=0.1, syst=0.2)
        >>> p.combined_std
        Scalar(...)
    """

    def __init__(
        self,
        pid: ParamId,
        expected: Union[float, complex, Scalar],
        stat: Union[float, complex, Scalar],
        syst: Union[float, complex, Scalar],
        scale: Optional[float] = None,
        bin: Optional[List] = None,
    ) -> None:
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
        """Current parameter value after any active shift."""
        return Scalar.from_cpp(self._cpp_obj.get_val())

    @property
    def std(self) -> Tuple[Scalar, Scalar]:
        """Statistical and systematic uncertainties as ``(stat, syst)``."""
        stat, syst = self._cpp_obj.get_std()
        return Scalar.from_cpp(stat), Scalar.from_cpp(syst)

    @property
    def combined_std(self) -> Scalar:
        """Combined uncertainty according to the C++ parameter rules."""
        return Scalar.from_cpp(self._cpp_obj.get_combined_std())

    @property
    def pid(self) -> ParamId:
        """Identifier of this parameter."""
        return ParamId.from_cpp(self._cpp_obj.get_id())

    def set_expected(self, val: Union[float, complex, Scalar]) -> None:
        """Set the central value.

        Args:
            val: New expected value.
        """
        self._cpp_obj.set_expected(_to_scalar(val)._cpp_obj)

    def set_std(
        self,
        stat: Union[float, complex, Scalar],
        syst: Union[float, complex, Scalar],
    ) -> None:
        """Set statistical and systematic uncertainties.

        Args:
            stat: Statistical uncertainty.
            syst: Systematic uncertainty.
        """
        self._cpp_obj.set_std(_to_scalar(stat)._cpp_obj, _to_scalar(syst)._cpp_obj)

    def set_shift(self, shift: Union[float, complex, Scalar]) -> None:
        """Set the active shift value for a shiftable parameter.

        Args:
            shift: Shift value forwarded to C++.
        """
        self._cpp_obj.set_shift(_to_scalar(shift)._cpp_obj)

    def set_mode(self, mode: ParameterMode) -> None:
        """Set whether the parameter is fixed or shiftable.

        Args:
            mode: New parameter mode.
        """
        self._cpp_obj.set_mode(mode.value)

    def __repr__(self) -> str:
        """Return a readable representation for debugging."""
        return (
            f"<Parameter pid={self.pid}, value={self.value}, "
            f"std={self.std}, combined_std={self.combined_std}>"
        )

    @staticmethod
    def from_cpp(cpp_obj) -> "Parameter":
        """Wrap an existing bound C++ ``Parameter``.

        Args:
            cpp_obj: Bound C++ parameter instance.

        Returns:
            Python wrapper sharing the provided C++ object.
        """
        param = Parameter.__new__(Parameter)
        param._cpp_obj = cpp_obj
        return param


class PyDependentParameter(Parameter):
    """Python wrapper around a C++ ``DependentParameter``.

    Dependent parameters can update their value from other parameters, and can
    be frozen to avoid repeated recomputation in workflows that scan many
    points.
    """

    def __init__(self, cpp_obj: _CppDependentParameter) -> None:
        """Wrap an existing C++ dependent parameter."""
        self._cpp_obj = cpp_obj

    def depends_on(self, pid: ParamId) -> bool:
        """Return whether this parameter depends on ``pid``."""
        return bool(self._cpp_obj.dependsOn(pid._cpp_obj))

    def init(self) -> None:
        """Initialize the underlying dependent parameter."""
        self._cpp_obj.init()

    def update(self) -> None:
        """Recompute the parameter value from its dependencies."""
        self._cpp_obj.update()

    def freeze(self) -> None:
        """Freeze updates for the dependent parameter."""
        self._cpp_obj.freeze()

    def unfreeze(self) -> None:
        """Re-enable updates for the dependent parameter."""
        self._cpp_obj.unfreeze()


__all__ = ["Parameter", "ParameterMode", "PyDependentParameter"]
