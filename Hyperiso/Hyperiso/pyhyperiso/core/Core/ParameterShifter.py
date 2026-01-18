from pyhyperiso.phyperiso.pyhyperiso.core import ParameterShifter as _CppParameterShifter
from pyhyperiso.core.Common.General import PyParamId
from pyhyperiso.core.Common.GeneralEnum import Model, ParameterType
from pyhyperiso.core.Math.scalar import Scalar, _to_scalar
from pyhyperiso.core.Core.Parameter import ParameterMode

class PyParameterShifter:
    """Wrapper for C++ ParameterShifter."""

    def __init__(self):
        self._cpp_obj = _CppParameterShifter()

    def mutate(self, pid: PyParamId, value: float):
        self._cpp_obj.mutate(pid._cpp_obj, _to_scalar(value)._cpp_obj)

    def change_mode(self, pid: PyParamId, mode: ParameterMode):
        self._cpp_obj.change_mode(pid._cpp_obj, mode.value)
        
if __name__ == "__main__":
    from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster
    from pathlib import Path
    from Hyperiso.Hyperiso.pyhyperiso.core.Core.HyperisoConfig import PyHyperisoConfig, ExternalFlag
    from pyhyperiso.core.Core.ParamaterProvider import PyParameterProvider
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
    
    param_shifter = PyParameterShifter()
    param_prov = PyParameterProvider(ParameterType.SM)
    
    print("param before mutate : ", param_prov.get_by_block("MASS", 24))
    param_shifter.mutate(PyParamId(ParameterType.SM, "MASS", 24), 82)
    print("param before mutate : ", param_prov.get_by_block("MASS", 24))
    