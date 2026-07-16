"""High-level Python wrapper for dependency pruning operations.

This module wraps the bound C++ ``DependencyPruner`` adapter. It provides a
small Python-friendly façade for detaching and reattaching dependent blocks or
individual dependent parameters from their upstream dependencies.

The Hyperiso core must be initialized before using this wrapper, typically via
``HyperisoMaster.init(...)``, because the underlying C++ adapter delegates to
``Parameters::GetInstance(...)``.

Example:
    >>> pruner = DependencyPruner()
    >>> pruner.detach_block(ParameterType.SM, "EW")
    >>> pruner.reattach_block(ParameterType.SM, "EW")
"""

from typing import Any

from pyhyperiso.phyperiso.pyhyperiso.core import DependencyPruner as _CppDependencyPruner
from pyhyperiso.core.Common.GeneralEnum import ParameterType


class DependencyPruner:
    """Detach and reattach dependency links in the Hyperiso parameter graph.

    ``DependencyPruner`` is a Python façade over the C++ adapter with the same
    name. It can operate at two levels:

    - block level, through :meth:`detach_block` and :meth:`reattach_block`,
    - parameter level, through :meth:`detach_parameter` and
      :meth:`reattach_parameter`.

    Args:
        cpp_obj: Optional already-bound C++ ``DependencyPruner`` instance. When
            omitted, a fresh C++ adapter is created.
    """

    def __init__(self, cpp_obj: _CppDependencyPruner | None = None):
        """Create a Python wrapper around the C++ ``DependencyPruner``."""
        self._cpp_obj = cpp_obj if cpp_obj is not None else _CppDependencyPruner()

    @staticmethod
    def _normalize_block_name(block_name: Any) -> str:
        """Convert a Python block-name-like object to the string expected by pybind.

        Args:
            block_name: Block name, usually a ``str``. Objects implementing
                ``__str__`` are also accepted.

        Returns:
            str: Block name passed to the C++ binding, where it is converted to
            ``BlockName``.
        """
        return str(block_name)

    @staticmethod
    def _to_cpp_lha_id(id_: Any) -> Any:
        """Return the bound C++ ``LhaID`` object for a Python LHA id wrapper.

        Args:
            id_: Either a bound C++ ``LhaID`` object or a Python wrapper exposing
                ``to_cpp()``.

        Returns:
            Any: Object forwarded to the pybind layer.
        """
        return id_.to_cpp() if hasattr(id_, "to_cpp") else id_

    def detach_block(self, parameter_type: ParameterType, block_name: Any) -> None:
        """Detach a dependent block from its upstream source blocks.

        After detachment, updates of upstream blocks no longer propagate through
        this dependency link until the block is reattached.

        Args:
            parameter_type: Parameter namespace containing the block, e.g.
                ``ParameterType.SM``.
            block_name: Name of the dependent block to detach.
        """
        self._cpp_obj.detach_block(parameter_type.value, self._normalize_block_name(block_name))

    def reattach_block(self, parameter_type: Any, block_name: Any) -> None:
        """Reattach a previously detached dependent block to its sources.

        Args:
            parameter_type: Parameter namespace containing the block, e.g.
                ``ParameterType.SM``.
            block_name: Name of the dependent block to reattach.
        """
        self._cpp_obj.reattach_block(parameter_type.value, self._normalize_block_name(block_name))

    def detach_parameter(self, parameter_type: ParameterType, block_name: Any, id_: Any) -> None:
        """Detach a dependent parameter from its upstream source parameters.

        Args:
            parameter_type: Parameter namespace containing the parameter.
            block_name: Name of the block containing the parameter.
            id_: LHA identifier of the dependent parameter. This may be a bound
                C++ ``LhaID`` or a Python wrapper exposing ``to_cpp()``.
        """
        self._cpp_obj.detach_parameter(
            parameter_type.value,
            self._normalize_block_name(block_name),
            self._to_cpp_lha_id(id_),
        )

    def reattach_parameter(self, parameter_type: ParameterType, block_name: Any, id_: Any) -> None:
        """Reattach a previously detached dependent parameter to its sources.

        Args:
            parameter_type: Parameter namespace containing the parameter.
            block_name: Name of the block containing the parameter.
            id_: LHA identifier of the dependent parameter. This may be a bound
                C++ ``LhaID`` or a Python wrapper exposing ``to_cpp()``.
        """
        self._cpp_obj.reattach_parameter(
            parameter_type.value,
            self._normalize_block_name(block_name),
            self._to_cpp_lha_id(id_),
        )

    def detach_param(self, parameter_type: ParameterType, block_name: Any, id_: Any) -> None:
        """Alias for :meth:`detach_parameter`.

        Args:
            parameter_type: Parameter namespace containing the parameter.
            block_name: Name of the block containing the parameter.
            id_: LHA identifier of the dependent parameter.
        """
        self.detach_parameter(parameter_type, block_name, id_)

    def reattach_param(self, parameter_type: Any, block_name: Any, id_: Any) -> None:
        """Alias for :meth:`reattach_parameter`.

        Args:
            parameter_type: Parameter namespace containing the parameter.
            block_name: Name of the block containing the parameter.
            id_: LHA identifier of the dependent parameter.
        """
        self.reattach_parameter(parameter_type, block_name, id_)

    @property
    def cpp_obj(self) -> _CppDependencyPruner:
        """Bound C++ ``DependencyPruner`` instance used by this wrapper."""
        return self._cpp_obj


__all__ = ["DependencyPruner"]
