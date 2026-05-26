"""Python accessors for MARTY-related configuration."""

from __future__ import annotations

from enum import Enum

from pyhyperiso.phyperiso.pyhyperiso.core import MartyAdapter as _CppMartyAdapter
from pyhyperiso.phyperiso.pyhyperiso.core import MartyPath as _CppMartyPath


class MartyPath(Enum):
    """MARTY path keys exposed by the C++ runtime.

    Attributes:
        MODEL_FILE: Path to the generated or configured MARTY model file.
        TEMPLATE_DIR: Path to the template directory.
        PARAM_MAPPING_DIR: Path to the parameter mapping directory.
    """

    MODEL_FILE = _CppMartyPath.MODEL_FILE
    TEMPLATE_DIR = _CppMartyPath.TEMPLATE_DIR
    PARAM_MAPPING_DIR = _CppMartyPath.PARAM_MAPPING_DIR


class MartyAdapter:
    """Inspect MARTY paths and flags used by Hyperiso.

    This adapter is read-only from Python. It is typically useful after
    ``HyperisoMaster.init(...)`` when debugging which MARTY model and support
    files were selected.
    """

    def __init__(self) -> None:
        """Create an adapter bound to the current C++ MARTY configuration."""
        self._cpp_obj = _CppMartyAdapter()

    def get_path(self, path_kind: MartyPath) -> str:
        """Return a MARTY path by kind.

        Args:
            path_kind: Path category to query.

        Returns:
            Path string returned by C++.
        """
        return str(self._cpp_obj.get_path(path_kind.value))

    def check_flag(self, flag) -> bool:
        """Return the value of a MARTY-related C++ flag.

        Args:
            flag: Bound C++ flag accepted by ``MartyAdapter.check_flag``.

        Returns:
            Boolean flag value.
        """
        return bool(self._cpp_obj.check_flag(flag))

    def get_model_name(self) -> str:
        """Return the configured MARTY model name."""
        return str(self._cpp_obj.get_marty_model_name())

    def __repr__(self) -> str:
        """Return a compact representation including the model name."""
        return f"<PyMartyAdapter model='{self.get_model_name()}'>"


__all__ = ["MartyAdapter", "MartyPath"]
    
    
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
    
    marty_adapt = MartyAdapter()
    
    print("model name : ", marty_adapt.get_model_name())
    print("moddel path : ", marty_adapt.get_path(MartyPath.MODEL_FILE))