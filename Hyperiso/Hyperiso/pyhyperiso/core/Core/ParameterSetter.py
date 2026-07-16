"""Mutable parameter setter for central values and modes."""

from __future__ import annotations

from typing import Union

from pyhyperiso.phyperiso.pyhyperiso.core import ParameterSetter as _CppParameterSetter
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Core.Parameter import ParameterMode
from pyhyperiso.core.Math.Scalar import Scalar, _to_scalar


class ParameterSetter:
    """Python wrapper around the C++ ``ParameterSetter``.

    ``ParameterSetter`` mutates stored parameter values directly. For workflows
    that require recomputing observables after mutation, remember to call the
    relevant reload/update method on the observable or statistic interface.

    Examples:
        >>> setter = ParameterSetter()
        >>> setter.mutate(ParamId(ParameterType.SM, "MASS", 24), 80.379)
    """

    def __init__(self) -> None:
        """Create a setter bound to the current C++ parameter storage."""
        self._cpp_obj = _CppParameterSetter()

    def mutate(self, pid: ParamId, value: Union[float, complex, Scalar]) -> None:
        """Set a parameter value.

        Args:
            pid: Parameter to mutate.
            value: New value. Numeric inputs are converted to ``Scalar``.
        """
        self._cpp_obj.mutate(pid._cpp_obj, _to_scalar(value)._cpp_obj)

    def change_mode(self, pid: ParamId, mode: ParameterMode) -> None:
        """Change a parameter update mode.

        Args:
            pid: Parameter to update.
            mode: New mode, such as ``ParameterMode.FIXED`` or
                ``ParameterMode.SHIFTABLE``.
        """
        self._cpp_obj.change_mode(pid._cpp_obj, mode.value)


__all__ = ["ParameterSetter"]
