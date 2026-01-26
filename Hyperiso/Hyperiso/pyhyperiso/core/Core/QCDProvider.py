from pyhyperiso.phyperiso.pyhyperiso.core import QCDProvider as _CppQCDProvider
from pyhyperiso.core.Common.General import PyParamId
from pyhyperiso.core.Common.GeneralEnum import Model, MassType
from pyhyperiso.core.Common.Configs import PyMassConfig, PyAlphasConfig
from pyhyperiso.core.Core.QCDConstants import PyQCDConstants
    
class PyQCDProvider:
    """Wrapper for C++ PyParameterSetter."""

    def __init__(self):
        self._cpp_obj = _CppQCDProvider()

    def get_alphas(self, alpha_config : PyAlphasConfig):
        return self._cpp_obj.compute_alphas(alpha_config.to_cpp())

    def get_qcd_masses(self, mass_config : PyMassConfig):
        return self._cpp_obj.compute_mass(mass_config.to_cpp())
    
    def get_qcd_constants(self):
        cpp_constants = self._cpp_obj.get_constants()
        return PyQCDConstants(cpp_constants)
        
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
    
    qcd_params = PyQCDProvider()

    print("alphas_s at 42 GeV : ", qcd_params.get_alphas(PyAlphasConfig(42)))
    
    print("mass at 5 GeV : ", qcd_params.get_qcd_masses(PyMassConfig(5, pdg_id=6)))
    
    consts = qcd_params.get_qcd_constants()

    print(consts.Nc)      # 3
    print(consts.C_F)     # 4/3
    print(consts.beta[0]) # par ex., [31/3, 134/3, ...]