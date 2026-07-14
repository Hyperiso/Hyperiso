"""Read-only access to parameters stored in the C++ core.

The filename keeps the historical spelling ``ParamaterProvider.py`` for import
compatibility. The public class is :class:`ParameterProvider`.
"""

from __future__ import annotations

from typing import Optional, Union

from pyhyperiso.phyperiso.pyhyperiso.core import ParameterProvider as _CppParameterProvider
from pyhyperiso.core.Common.GeneralEnum import DataType, ParameterType
from pyhyperiso.core.Common.ParamId import LhaID, ParamId
from pyhyperiso.core.Core.Parameter import Parameter
from pyhyperiso.core.Math.Scalar import Scalar


class ParameterProvider:
    """Read parameter values and metadata from a C++ parameter namespace.

    A provider can be bound to a specific :class:`ParameterType`, or left
    unfiltered to use the C++ provider default. Parameters can be addressed by
    a full ``ParamId`` or by the pair ``(block, code)``.

    Examples:
        >>> provider = ParameterProvider(ParameterType.SM)
        >>> pid = ParamId(type=ParameterType.SM, block="MASS", code=24)
        >>> provider.exists_by_pid(pid)
        True
        >>> provider.get_by_block("MASS", 24)
        80.379
    """

    def __init__(self, param_type: Optional[ParameterType] = None) -> None:
        """Create a provider.

        Args:
            param_type: Optional parameter namespace. When omitted, the C++
                default provider constructor is used.
        """
        self._type = param_type
        self._cpp_obj = (
            _CppParameterProvider(param_type.value)
            if param_type is not None
            else _CppParameterProvider()
        )

    def get_by_pid(self, pid: ParamId, dtype: DataType = DataType.VALUE) -> Scalar:
        """Return a parameter component addressed by ``ParamId``.

        Args:
            pid: Parameter identifier.
            dtype: Component to read, for example central value or uncertainty.

        Returns:
            Requested parameter component as a Scalar (real,complex).
        """
        return Scalar.from_cpp(self._cpp_obj(pid._cpp_obj, dtype.value))

    def get_by_block(
        self,
        block: str,
        code: Union[int, str, list, LhaID],
        dtype: DataType = DataType.VALUE,
    ) -> Scalar:
        """Return a parameter component addressed by block and LHA code.

        Args:
            block: LHA-style block name.
            code: LHA code. Integers, strings, lists and ``LhaID`` are accepted.
            dtype: Component to read.

        Returns:
            Requested parameter component as a Scalar (real,complex).
        """
        if not isinstance(code, LhaID):
            code = LhaID(code)
        return Scalar.from_cpp(self._cpp_obj(block, code._cpp_obj, dtype.value))

    def exists_by_pid(self, pid: ParamId) -> bool:
        """Return whether a parameter exists for a full ``ParamId``."""
        return bool(self._cpp_obj.exists(pid._cpp_obj))

    def exists_by_block(self, block: str, code: Union[int, str, list, LhaID]) -> bool:
        """Return whether a parameter exists for ``(block, code)``."""
        if not isinstance(code, LhaID):
            code = LhaID(code)
        return bool(self._cpp_obj.exists(block, code._cpp_obj))

    def get_parameter(self, pid: ParamId) -> Parameter:
        """Return the full parameter object for ``pid``.

        Args:
            pid: Parameter identifier to retrieve.

        Returns:
            Python ``Parameter`` wrapper around the bound C++ object.
        """
        cpp_param = self._cpp_obj.get_parameter(pid._cpp_obj)
        return Parameter.from_cpp(cpp_param)

    def get_type(self) -> ParameterType:
        """Return the namespace served by this provider."""
        return ParameterType(self._cpp_obj.get_type())

    def __repr__(self) -> str:
        """Return a compact representation including the provider type."""
        return f"<PyParameterProvider type={self.get_type().name}>"


__all__ = ["ParameterProvider"]
