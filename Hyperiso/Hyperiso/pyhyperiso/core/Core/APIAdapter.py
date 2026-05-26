"""Python accessors for global Hyperiso API state.

This module exposes a small Pythonic wrapper around the C++ ``APIAdapter``.
The adapter is intended for inspection and diagnostics after the global
Hyperiso runtime has been initialized with :class:`HyperisoMaster`.

Typical uses include checking initialization flags, retrieving configured
paths, listing available parameter blocks and reading block entries from a
specific parameter namespace.
"""

from __future__ import annotations

from enum import Enum
from typing import Set

from pyhyperiso.phyperiso.pyhyperiso.core import APIAdapter as _CppAPIAdapter
from pyhyperiso.phyperiso.pyhyperiso.core import APIPath as _CppAPIPath
from pyhyperiso.core.Common.BlockName import BlockName
from pyhyperiso.core.Common.GeneralEnum import ParameterType
from pyhyperiso.core.Common.LhaID import LhaID
from pyhyperiso.core.Core.HyperisoMaster import ExternalFlag
from pyhyperiso.core.Math.scalar import Scalar


class APIPath(Enum):
    """Known path keys exposed by the C++ API adapter.

    Attributes:
        LHA_PATH: Path to the active LHA input file.
    """

    LHA_PATH = _CppAPIPath.LHA_PATH


class APIAdapter:
    """Inspection wrapper around the C++ ``APIAdapter``.

    The adapter reads process-wide state managed by the C++ core. It should
    normally be used after ``HyperisoMaster.init(...)`` has loaded a model and
    an LHA input file.

    Examples:
        >>> from pyhyperiso.core.Core.HyperisoMaster import HyperisoMaster
        >>> from pyhyperiso.core.Core.HyperisoConfig import HyperisoConfig
        >>> hyp = HyperisoMaster()
        >>> hyp.init("lha/si_input.flha", HyperisoConfig())
        >>> api = APIAdapter()
        >>> api.check_flag(ExternalFlag.IS_LHA_SPECTRUM)
        False
        >>> "MASS" in {str(b) for b in api.get_blocks_list(ParameterType.SM)}
        True
    """

    def __init__(self) -> None:
        """Create an adapter bound to the current C++ runtime state."""
        self._cpp_obj = _CppAPIAdapter()

    def check_flag(self, ext_flag: ExternalFlag) -> bool:
        """Return whether an external initialization flag is enabled.

        Args:
            ext_flag: Python wrapper for a C++ ``ExternalFlag`` value.

        Returns:
            ``True`` if the flag is active in the current Hyperiso session,
            otherwise ``False``.
        """
        return bool(self._cpp_obj.check_flag(ext_flag.value))

    def get_path(self, path_type: APIPath) -> str:
        """Return a configured path from the C++ runtime.

        Args:
            path_type: Path key to query.

        Returns:
            Path string associated with ``path_type``.
        """
        return str(self._cpp_obj.get_path(path_type.value))

    def get_all_blocks(self) -> Set[BlockName]:
        """Return all known block names across parameter namespaces.

        Returns:
            Set of Python ``BlockName`` wrappers.
        """
        raw_blocks = self._cpp_obj.get_all_blocks()
        return {BlockName(x) for x in raw_blocks}

    def get_blocks_list(self, param_type: ParameterType = ParameterType.SM) -> Set[BlockName]:
        """Return the block names available in one parameter namespace.

        Args:
            param_type: Parameter namespace to inspect. For example,
                ``ParameterType.SM`` reads Standard Model input blocks, while
                ``ParameterType.WILSON`` reads Wilson-coefficient blocks.

        Returns:
            Set of available block names.
        """
        raw_blocks = self._cpp_obj.get_blocks_list(param_type.value)
        return {BlockName(x) for x in raw_blocks}

    def get_block_infos(
        self,
        block_name: str,
        param_type: ParameterType = ParameterType.SM,
    ) -> dict[LhaID, Scalar]:
        """Return all entries stored in a parameter block.

        Args:
            block_name: Name of the LHA-style block, such as ``"MASS"`` or
                ``"SMINPUTS"``.
            param_type: Parameter namespace containing the block.

        Returns:
            Mapping from LHA code to scalar value.
        """
        raw_block_infos = self._cpp_obj.get_block_infos(
            BlockName(block_name)._cpp_obj,
            param_type.value,
        )
        return {LhaID(x): Scalar.from_cpp(y) for x, y in raw_block_infos.items()}

    def get_type_of_block(self, block_name: str) -> list[ParameterType]:
        """Return namespaces in which a block name exists.

        Args:
            block_name: Name of the block to locate.

        Returns:
            List of ``ParameterType`` values containing ``block_name``.
        """
        cpp_result = self._cpp_obj.get_type_of_block(BlockName(block_name)._cpp_obj)
        return [ParameterType(item) for item in cpp_result]


__all__ = ["APIAdapter", "APIPath"]
    
if __name__ == "__main__":
    from pyhyperiso.core.Core.HyperisoMaster import HyperisoMaster
    from pathlib import Path
    from pyhyperiso.core.Core.HyperisoConfig import HyperisoConfig, ExternalFlag
    from pyhyperiso.core.Common.GeneralEnum import Model
    
    print("🔧 Initializing PyHyperisoMaster with custom HyperisoConfig...")

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
    
    api_adapter = APIAdapter()

    print("all blocks : ", api_adapter.get_all_blocks())
    
    print("Block mass : ", api_adapter.get_block_infos("MASS"))
    
    print("blocks : ", api_adapter.get_blocks_list(ParameterType.FLAVOR), " are in the flavor part")
    
    print("block SMINPUTS is in : ", api_adapter.get_type_of_block("SMINPUTS"))
    