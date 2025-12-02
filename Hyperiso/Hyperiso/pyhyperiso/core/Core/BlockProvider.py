from pyhyperiso.phyperiso.pyhyperiso.core import BlockProvider as _CppBlockProvider
from pyhyperiso.core.Common.General import PyParamId
from pyhyperiso.core.Common.GeneralEnum import Model, ParameterType

class PyBlockLogger:
    """Wrapper for C++ PyParameterSetter."""

    def __init__(self):
        self._cpp_obj = _CppBlockProvider()

    def exists(self, blockname : str, param_type : ParameterType):
        return self._cpp_obj.exists(blockname, param_type.value)

    def log_all_blocks(self, param_type : ParameterType):
        return self._cpp_obj.log_all_blocks(param_type.value)
        
    def log_block(self, param_type : ParameterType, blockname : str):
        return self._cpp_obj.log_block(param_type.value, blockname)
        
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
            # ExternalFlag.USE_MARTY: False
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
    
    block_prov = PyBlockLogger()

    print("block mass : ")
    block_prov.log_block(ParameterType.SM, "MASS")
    print("all blocks in wilson: ")
    block_prov.log_all_blocks(ParameterType.WILSON)
    
    print("does block SMINPUTS exists ? ", block_prov.exists("SMINPUTS", ParameterType.SM))