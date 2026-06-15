from pathlib import Path
import os

from pyhyperiso.Common import ParameterType, ParamId, DataType

from pyhyperiso.Core import HyperisoConfig, HyperisoMaster, APIPath, ParameterProvider

if __name__ == "__main__":

    #All informations for the creation of the HyperisoMaster are provided in the base_core_example.py example.
    config = HyperisoConfig()
    hyp = HyperisoMaster()
    parameters_values = os.path.join(__file__, "..", "inputs", "sm_param_test.yaml")
    
    #This API allows to set special paths for Hyperiso's inputs. One can change parameters through this API, Observables, Observables correlations, etc.
    #We recommend not to change DEFAULT values, even if this API allows it.
    #User inputs take the form of .yml/.yaml files, where you can specify the value and uncertainties of your parameters.
    hyp.pre_init_set_paths({
        APIPath.USER_SM_PARAMS: parameters_values
    })
    
    #One can also add specials blocks through a pre_init. This is useful for BSM models with new groups not in the flha convention.
    #If its Matrix block, with multiple ID, use the options available within the API.
    hyp.pre_init_add_block("SPECIAL_BLOCk")
    
    #If you already have MARTY installed on your machine, please specify the path through the following API.
    # hyp.pre_init_set_marty_path("Path/To/Marty")
    
    lha_file_path = "lha/si_input.flha" 
    hyp.init(lha_file=lha_file_path, config=config)
    
    print(ParameterProvider().get_by_pid(ParamId(ParameterType.SM, "MASS", 18)))
    
    
