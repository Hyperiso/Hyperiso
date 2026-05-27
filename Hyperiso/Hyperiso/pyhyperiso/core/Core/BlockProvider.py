"""Helpers for inspecting and logging C++ parameter blocks.

The C++ ``BlockProvider`` gives access to HyperISO parameter blocks. It can
check whether a block exists, log blocks through the C++ logging system, and
retrieve block names or block contents as Python objects.
"""

from __future__ import annotations

from typing import Any, TypeAlias

from pyhyperiso.phyperiso.pyhyperiso.core import BlockProvider as _CppBlockProvider
from pyhyperiso.core.Common.GeneralEnum import ParameterType


BlockKey: TypeAlias = Any
BlockContent: TypeAlias = dict[BlockKey, float]


class BlockLogger:
    """Diagnostic and inspection wrapper around the C++ ``BlockProvider``.

    This class exposes block-level utilities for HyperISO parameters. It keeps
    the original C++ logging methods available while also providing Pythonic
    accessors for block names and block contents.

    Examples:
        >>> logger = BlockLogger()
        >>> logger.exists("SMINPUTS", ParameterType.SM)
        True
        >>> logger.get_all_blocks(ParameterType.SM)
        {"SMINPUTS", "MASS", ...}
        >>> logger.get_block(ParameterType.SM, "MASS")
        {25: 125.1, ...}
    """

    def __init__(self) -> None:
        """Create a block logger for the current C++ parameter storage."""
        self._cpp_obj = _CppBlockProvider()

    @staticmethod
    def _to_cpp_param_type(param_type: ParameterType) -> int:
        """Convert a Python ``ParameterType`` enum to its C++ value.

        Args:
            param_type: Python parameter namespace enum.

        Returns:
            Integer value expected by the C++ binding.

        Raises:
            TypeError: If ``param_type`` is not a ``ParameterType`` instance.
        """
        if not isinstance(param_type, ParameterType):
            raise TypeError(
                "param_type must be an instance of ParameterType, "
                f"got {type(param_type).__name__}."
            )

        return param_type.value

    def exists(self, blockname: str, param_type: ParameterType) -> bool:
        """Return whether a block exists in a parameter namespace.

        Args:
            blockname: LHA-style block name, for example ``"MASS"``.
            param_type: Parameter namespace to inspect.

        Returns:
            ``True`` if the block exists, otherwise ``False``.
        """
        return bool(
            self._cpp_obj.exists(
                blockname,
                self._to_cpp_param_type(param_type),
            )
        )

    def log_all_blocks(self, param_type: ParameterType) -> None:
        """Log all blocks from a parameter namespace through the C++ logger.

        Args:
            param_type: Parameter namespace to log.
        """
        self._cpp_obj.log_all_blocks(self._to_cpp_param_type(param_type))

    def log_block(self, param_type: ParameterType, blockname: str) -> None:
        """Log one block through the C++ logger.

        Args:
            param_type: Parameter namespace containing the block.
            blockname: Name of the block to log.
        """
        self._cpp_obj.log_block(
            self._to_cpp_param_type(param_type),
            blockname,
        )

    def get_block(self, param_type: ParameterType, blockname: str) -> BlockContent:
        """Return the content of one block as a Python dictionary.

        Args:
            param_type: Parameter namespace containing the block.
            blockname: Name of the block to retrieve.

        Returns:
            Dictionary mapping LHA identifiers to scalar values.

        Examples:
            >>> logger = BlockLogger()
            >>> logger.get_block(ParameterType.SM, "MASS")
            {25: 125.1, ...}
        """
        return dict(
            self._cpp_obj.get_block(
                self._to_cpp_param_type(param_type),
                blockname,
            )
        )

    def get_all_blocks(self, param_type: ParameterType) -> set[str]:
        """Return all block names available in a parameter namespace.

        Args:
            param_type: Parameter namespace to inspect.

        Returns:
            Set containing all available block names.

        Examples:
            >>> logger = BlockLogger()
            >>> logger.get_all_blocks(ParameterType.SM)
            {"SMINPUTS", "MASS", ...}
        """
        return set(
            self._cpp_obj.get_all_blocks(
                self._to_cpp_param_type(param_type),
            )
        )


__all__ = ["BlockLogger", "BlockContent", "BlockKey"]
        
if __name__ == "__main__":
    from pyhyperiso.core.Core.HyperisoMaster import HyperisoMaster
    from pathlib import Path
    from pyhyperiso.core.Core.HyperisoConfig import HyperisoConfig, ExternalFlag
    from pyhyperiso.core.Common.GeneralEnum import Model
    
    print("🔧 Initializing PyHyperisoMaster with custom PyHyperisoConfig...")

    config = HyperisoConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: True,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            # ExternalFlag.USE_MARTY: False
        },
        model=Model.SM,
        mty_model_name="MSSM_UFO",
        mty_model_path=Path("/my/custom/marty/path")
    )

    print("🔧 PyHyperisoConfig content:")
    print(config)

    hyp = HyperisoMaster()
    lha_file_path = "lha/camilia.flha" 

    print("\n🚀 Calling init with config...")
    hyp.init(lha_file=lha_file_path, config=config)
    
    block_prov = BlockLogger()

    print("block mass : ")
    block_prov.log_block(ParameterType.SM, "MASS")
    print("all blocks in wilson: ")
    block_prov.log_all_blocks(ParameterType.WILSON)
    
    print("does block SMINPUTS exists ? ", block_prov.exists("SMINPUTS", ParameterType.SM))
    
    print("block MASS:")
    print(block_prov.get_block(ParameterType.SM, "MASS"))

    print("all blocks in SM:")
    print(block_prov.get_all_blocks(ParameterType.SM))

    print(
        "does block SMINPUTS exist?",
        block_prov.exists("SMINPUTS", ParameterType.SM),
    )