from pyhyperiso.phyperiso.pyhyperiso.core import MartyAdapter as _CppMartyAdapter
from pyhyperiso.phyperiso.pyhyperiso.core import MartyPath as _CppMartyPath
from pyhyperiso.core.Common.GeneralEnum import Model
from enum import Enum

class MartyPath(Enum):
    MODEL_FILE = _CppMartyPath.MODEL_FILE
    TEMPLATE_DIR = _CppMartyPath.TEMPLATE_DIR
    PARAM_MAPPING_DIR = _CppMartyPath.PARAM_MAPPING_DIR

class PyMartyAdapter:
    """Wrapper for MartyAdapter (path + flags)."""

    def __init__(self):
        self._cpp_obj = _CppMartyAdapter()

    def get_path(self, path_kind : MartyPath):
        return self._cpp_obj.get_path(path_kind.value)

    def check_flag(self, flag):
        return self._cpp_obj.check_flag(flag)

    def get_model_name(self):
        return self._cpp_obj.get_marty_model_name()

    def __repr__(self):
        return f"<PyMartyAdapter model='{self.get_model_name()}'>"
    
    
if __name__ == "__main__":
    from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster
    from pathlib import Path
    from Hyperiso.Hyperiso.pyhyperiso.core.Core.HyperisoConfig import PyHyperisoConfig, ExternalFlag
    print("🔧 Initializing PyHyperisoMaster with custom PyHyperisoConfig...")

    config = PyHyperisoConfig(
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

    hyp = PyHyperisoMaster()
    lha_file_path = "lha/camilia.flha"

    print("\n🚀 Calling init with config...")
    hyp.init(lha_file=lha_file_path, config=config)
    
    marty_adapt = PyMartyAdapter()
    
    print("model name : ", marty_adapt.get_model_name())
    print("moddel path : ", marty_adapt.get_path(MartyPath.MODEL_FILE))