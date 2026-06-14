from pathlib import Path

from pyhyperiso.Common import (
    Model
)
from pyhyperiso.Core import HyperisoConfig, HyperisoMaster, ExternalFlag

if __name__ == "__main__":

    #Configuration for HyperisoMaster, with global flags and model's informations.
    config = HyperisoConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: False, #If True, Hyperiso will consider that all parameters are inside the lha, and doesn't need a calculator. Useful only with THDM or SUSY models.
            ExternalFlag.HAS_WILSON_INPUT: False, #If True, Hyperiso will search for wilson coefficients inside the flha. If only BSM are provided, SM contributions will be calculated within Hyperiso.
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False, #If True, Hyperiso will search for theoretical values for observables inside the flha.
            ExternalFlag.HYP_AS_SM_MARTY: True, #If True, Hyperiso will calculate wilson coefficients up to NNLO for the SM and let MARTY calculate only the BSM contributions.
        },
        model=Model.SM, #Model (SM, THDM, SUSY, MARTY), if MARTY, the name of the model and the path to the marty model file is necessary.
        mty_model_name="marty_model_name", #Only use in Model.MARTY mode. Name of the model (class name in MARTY).
        mty_model_path=Path("/my/custom/marty/path") #Only use in Model.MARTY mode. Path to the model file (class name in MARTY).
    )

    #Creation of the HyperisoMaster, base of all other interfaces and singleton, need to exist for all others interfaces.
    hyp = HyperisoMaster()
    
    #Path to the lha file for inputs. If relative, start in the Assets/ directory.
    lha_file_path = "lha/si_input.flha" 

    #Initialization of Hyperiso, with configuration and path to a lha for inputs.
    hyp.init(lha_file=lha_file_path, config=config)

    #Some flag checks.
    print("Current model:", hyp.model.name)
    print("Flag IS_LHA_SPECTRUM:", hyp.check_flag(ExternalFlag.IS_LHA_SPECTRUM))

    #If you want to switch LHA you can use this API, this will reload all parameters from the LHA.
    #Be careful if you have already done calculation, wilson coefficients will be removed from the database and need to be recalculated.
    hyp.switch_lha("lha/testinput_thdm.lha",config)
    