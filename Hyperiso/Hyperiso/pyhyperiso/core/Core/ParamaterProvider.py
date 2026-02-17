from pyhyperiso.phyperiso.pyhyperiso.core import ParameterProvider as _CppParameterProvider
from pyhyperiso.core.Common.GeneralEnum import ParameterType, DataType, Model
from pyhyperiso.core.Common.ParamId import LhaID, ParamId
from pyhyperiso.core.Core.Parameter import PyParameter
from typing import Dict, Optional, Union



class PyParameterProvider:
    """
    Python wrapper for the C++ ParameterProvider class.
    Provides dual access via ParamId or (block, code).
    """

    def __init__(self, param_type: Optional[ParameterType] = None):
        """Initialize with optional ParameterType filter."""
        self._type = param_type
        self._cpp_obj = (_CppParameterProvider(param_type.value)
                         if param_type is not None else
                         _CppParameterProvider())

    def get_by_pid(self, pid: ParamId, dtype: DataType = DataType.VALUE) -> float:
        """Get parameter value using a ParamId."""
        return self._cpp_obj(pid._cpp_obj, dtype.value)

    def get_by_block(self, block: str, code: Union[int, str, list, LhaID], dtype: DataType = DataType.VALUE) -> float:
        """Get parameter value using block name and code."""
        if not isinstance(code, LhaID):
            code = LhaID(code)
        return self._cpp_obj(block, code._cpp_obj, dtype.value)

    def exists_by_pid(self, pid: ParamId) -> bool:
        """Check if a parameter exists via ParamId."""
        return self._cpp_obj.exists(pid._cpp_obj)

    def exists_by_block(self, block: str, code: Union[int, str, list, LhaID]) -> bool:
        """Check if a parameter exists using block name and code."""
        if not isinstance(code, LhaID):
            code = LhaID(code)
        return self._cpp_obj.exists(block, code._cpp_obj)

    def get_parameter(self, pid: ParamId) -> PyParameter:
        """Calls the `get_parameter` method (different from __call__)."""
        _cpp_param = self._cpp_obj.get_parameter(pid._cpp_obj)
        
        return PyParameter.from_cpp(_cpp_param)

    def get_type(self) -> ParameterType:
        """Returns the ParameterType of this provider."""
        return ParameterType(self._cpp_obj.get_type())

    def __repr__(self):
        return f"<PyParameterProvider type={self.get_type().name}>"
    
    
if __name__ == "__main__" :
    
    from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster
    from pathlib import Path
    from pyhyperiso.core.Core.HyperisoConfig import PyHyperisoConfig, ExternalFlag
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
    
    provider = PyParameterProvider(ParameterType.SM)
    pid = ParamId(type=ParameterType.SM, block="MASS", code=24)

    print("🔍 By ParamId")
    print("exists:", provider.exists_by_pid(pid))
    print("value:", provider.get_by_pid(pid))
    print("stored (get_parameter):", provider.get_parameter(pid))

    print("\n🔍 By (block, code)")
    print("exists:", provider.exists_by_block("MASS", [24]))
    print("value:", provider.get_by_block("MASS", "24"))