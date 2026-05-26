"""Logging helpers for C++ parameter blocks.

The C++ ``BlockProvider`` is primarily a diagnostic component: it can check
whether a block exists and print block contents to the C++ logging stream. This
Python wrapper keeps that behavior available to notebooks, scripts and tests.
"""

from __future__ import annotations

from pyhyperiso.phyperiso.pyhyperiso.core import BlockProvider as _CppBlockProvider
from pyhyperiso.core.Common.GeneralEnum import ParameterType


class BlockLogger:
    """Diagnostic wrapper around the C++ ``BlockProvider``.

    The methods mirror C++ logging utilities. They are useful when inspecting
    which parameter blocks were loaded after initializing Hyperiso.

    Examples:
        >>> logger = BlockLogger()
        >>> logger.exists("SMINPUTS", ParameterType.SM)
        True
    """

    def __init__(self) -> None:
        """Create a block logger for the current C++ parameter storage."""
        self._cpp_obj = _CppBlockProvider()

    def exists(self, blockname: str, param_type: ParameterType) -> bool:
        """Return whether a block exists in a parameter namespace.

        Args:
            blockname: LHA-style block name, for example ``"MASS"``.
            param_type: Parameter namespace to inspect.

        Returns:
            ``True`` if the block exists, otherwise ``False``.
        """
        return bool(self._cpp_obj.exists(blockname, param_type.value))

    def log_all_blocks(self, param_type: ParameterType) -> None:
        """Print all blocks from a parameter namespace via the C++ logger.

        Args:
            param_type: Parameter namespace to print.
        """
        return self._cpp_obj.log_all_blocks(param_type.value)

    def log_block(self, param_type: ParameterType, blockname: str) -> None:
        """Print one block via the C++ logger.

        Args:
            param_type: Parameter namespace containing the block.
            blockname: Block name to print.
        """
        return self._cpp_obj.log_block(param_type.value, blockname)


__all__ = ["BlockLogger"]
        
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