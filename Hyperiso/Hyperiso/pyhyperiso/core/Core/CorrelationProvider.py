from pyhyperiso.phyperiso.pyhyperiso.core import CorrelationProvider as _CppCorrelationProvider
from pyhyperiso.phyperiso.pyhyperiso.core import CorrelationType as _CppCorrelationType
from pyhyperiso.core.Common.General import PyParamId
from pyhyperiso.core.Common.GeneralEnum import Model, ParameterType, Observables
from enum import Enum

class CorrelationType(Enum):
    STAT = _CppCorrelationType.STAT
    SYST = _CppCorrelationType.SYST
    COMBINED = _CppCorrelationType.COMBINED
    
class PyCorrelationProvider:
    """Wrapper for C++ PyCorrelationProvider."""

    def __init__(self):
        self._cpp_obj = _CppCorrelationProvider()

    def correlation_from_paramid(self, pid_1: PyParamId, pid_2 : PyParamId, corr_type : CorrelationType):
        return self._cpp_obj.correlation_from_paramid(pid_1._cpp_obj, pid_2._cpp_obj, corr_type.value)

    def correlation_from_observable(self, obs_1 : Observables, obs_2 : Observables, corr_type : CorrelationType):
        return self._cpp_obj.correlation_from_observable(obs_1.value, obs_2.value, corr_type.value)
        
if __name__ == "__main__":
    from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster
    from pathlib import Path
    from pyhyperiso.core.Core.HyperisoConfig import PyHyperisoConfig, ExternalFlag
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
    
    corr_provider = PyCorrelationProvider()
    
    print("correlation between two obs (same obs)", corr_provider.correlation_from_observable(Observables.BR_B_XS_GAMMA, Observables.BR_B_XS_GAMMA, CorrelationType.COMBINED))