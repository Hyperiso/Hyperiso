

from pyhyperiso.phyperiso.pyhyperiso.core import APIAdapter as _CppAPIAdapter
from pyhyperiso.phyperiso.pyhyperiso.core import APIPath as _CppAPIPath
from pyhyperiso.core.Common.General import PyBlockName, PyLhaID
from pyhyperiso.core.Common.GeneralEnum import Model, ParameterType
from pyhyperiso.core.Common.Configs import PyMassConfig, PyAlphasConfig
from pyhyperiso.core.Core.QCDConstants import PyQCDConstants
from pyhyperiso.core.Core.HyperisoMaster import ExternalFlag
from enum import Enum
from pyhyperiso.core.Math.scalar import Scalar

class APIPath(Enum):
    LHA_PATH = _CppAPIPath.LHA_PATH  
class PyAPIAdapter:
    """Wrapper for C++ PyParameterSetter."""

    def __init__(self):
        self._cpp_obj = _CppAPIAdapter()

    def check_flag(self, ext_flag : ExternalFlag):
        return self._cpp_obj.check_flag(ext_flag.value)

    def get_path(self, path_type : APIPath):
        return self._cpp_obj.get_path(path_type.value)
    
    def get_all_blocks(self):
        raw_blocks = self._cpp_obj.get_all_blocks()
        return {PyBlockName(x) for x in raw_blocks}
    
    def get_blocks_list(self, param_type : ParameterType = ParameterType.SM):
        raw_blocks = self._cpp_obj.get_blocks_list(param_type.value)
        return {PyBlockName(x) for x in raw_blocks}
    
    def get_block_infos(self, block_name : str, param_type : ParameterType = ParameterType.SM):
        raw_block_infos = self._cpp_obj.get_block_infos(PyBlockName(block_name)._cpp_obj, param_type.value)
        return {PyLhaID(x) : Scalar().from_cpp(y) for x,y in raw_block_infos.items()}
    
    def get_type_of_block(self, block_name : str):
        cpp_result = self._cpp_obj.get_type_of_block(PyBlockName(block_name)._cpp_obj)
        return [ParameterType(item) for item in cpp_result]
    
     
if __name__ == "__main__":
    from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster
    from pathlib import Path
    from pyhyperiso.core.Core.Config import PyConfig, ExternalFlag
    from pyhyperiso.core.Core.ParamaterProvider import PyParameterProvider
    print("🔧 Initializing PyHyperisoMaster with custom PyConfig...")

    config = PyConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: True,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            ExternalFlag.USE_MARTY: False
        },
        model=Model.SM,
        mty_model_name="MSSM_UFO",
        mty_model_path=Path("/my/custom/marty/path")
    )

    print("🔧 PyConfig content:")
    print(config)

    hyp = PyHyperisoMaster()
    lha_file_path = "lha/camilia.flha" 

    print("\n🚀 Calling init with config...")
    hyp.init(lha_file=lha_file_path, config=config)
    
    api_adapter = PyAPIAdapter()

    print("all blocks : ", api_adapter.get_all_blocks())
    
    print("Block mass : ", api_adapter.get_block_infos("MASS"))
    
    print("blocks : ", api_adapter.get_blocks_list(ParameterType.FLAVOR), " are in the flavor part")
    
    print("block SMINPUTS is in : ", api_adapter.get_type_of_block("SMINPUTS"))
    