"""Python wrapper for dependent-block graph introspection.

This module exposes a high-level Python facade over the C++
``DependantBlockInfoProvider``. It lets Python workflows inspect whether a
block is dependent, which blocks it depends on, and which blocks depend on it.

The Hyperiso core must be initialized before calling these methods, typically
through ``HyperisoMaster.init(...)``.
"""

from typing import List, Optional

from pyhyperiso.phyperiso.pyhyperiso.core import (
    DependantBlockInfoProvider as _CppDependantBlockInfoProvider,
)
from pyhyperiso.core.Common.GeneralEnum import ParameterType


class DependantBlockInfoProvider:
    """Inspect block-level dependencies from Python.

    Args:
        cpp_obj: Optional already-bound C++ provider. When omitted, a new C++
            ``DependantBlockInfoProvider`` is created.

    Example:
        >>> provider = DependantBlockInfoProvider()
        >>> provider.is_dependent_block(ParameterType.SM, "EW")
        True
        >>> provider.get_source_blocks(ParameterType.SM, "EW")
        ['QCD', 'SMINPUTS']
    """

    def __init__(self, cpp_obj: Optional[_CppDependantBlockInfoProvider] = None):
        """Create the Python wrapper around the C++ provider."""
        self._cpp_obj = cpp_obj or _CppDependantBlockInfoProvider()

    def is_dependent_block(self, parameter_type: ParameterType, block_name: str) -> bool:
        """Return whether ``block_name`` is a dependent block.

        Args:
            parameter_type: Parameter namespace containing the block.
            block_name: Name or alias of the block to inspect.

        Returns:
            bool: ``True`` when the underlying C++ block is a dependent block.
        """
        return self._cpp_obj.is_dependent_block(parameter_type.value, block_name)

    def get_source_blocks(self, parameter_type: ParameterType, block_name: str) -> List[str]:
        """Return the direct upstream blocks used by ``block_name``.

        This answers: "which blocks does this block depend on directly?"

        Args:
            parameter_type: Parameter namespace containing the block.
            block_name: Name or alias of the block to inspect.

        Returns:
            list[str]: Sorted direct source block names.
        """
        return list(self._cpp_obj.get_source_blocks(parameter_type.value, block_name))

    def get_dependent_blocks(self, parameter_type: ParameterType, block_name: str) -> List[str]:
        """Return the direct blocks depending on ``block_name``.

        This answers: "which blocks directly depend on this block?"

        Args:
            parameter_type: Parameter namespace containing the block.
            block_name: Name or alias of the block to inspect.

        Returns:
            list[str]: Sorted direct dependent/downstream block names.
        """
        return list(self._cpp_obj.get_dependent_blocks(parameter_type.value, block_name))

    def get_all_source_blocks(self, parameter_type: ParameterType, block_name: str) -> List[str]:
        """Return all transitive upstream blocks used by ``block_name``.

        Args:
            parameter_type: Parameter namespace containing the block.
            block_name: Name or alias of the block to inspect.

        Returns:
            list[str]: Sorted deduplicated source block names.
        """
        return list(self._cpp_obj.get_all_source_blocks(parameter_type.value, block_name))

    def get_all_dependent_blocks(self, parameter_type: ParameterType, block_name: str) -> List[str]:
        """Return all transitive blocks depending on ``block_name``.

        Args:
            parameter_type: Parameter namespace containing the block.
            block_name: Name or alias of the block to inspect.

        Returns:
            list[str]: Sorted deduplicated dependent/downstream block names.
        """
        return list(self._cpp_obj.get_all_dependent_blocks(parameter_type.value, block_name))


__all__ = ["DependantBlockInfoProvider"]

if __name__ == "__main__":
    from pathlib import Path
    from pyhyperiso.core.Core.HyperisoMaster import HyperisoMaster, HyperisoConfig, ExternalFlag,Model
    
    print("🔧 Initializing PyHyperisoMaster with custom PyHyperisoConfig...")

    config = HyperisoConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: True,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            ExternalFlag.HYP_AS_SM_MARTY : True
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

    print("✅ Current model:", hyp.model.name)
    print("✅ Flag IS_LHA_SPECTRUM:", hyp.check_flag(ExternalFlag.IS_LHA_SPECTRUM))

    hyp.switch_lha("lha/testinput_thdm.lha",config)
    
    dep = DependantBlockInfoProvider()
    
    print(dep.is_dependent_block(ParameterType.SM, "QCD"))