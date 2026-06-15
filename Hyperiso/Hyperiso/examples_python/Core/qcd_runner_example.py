from pathlib import Path
import os

from pyhyperiso.Common import ParameterType, ParamId, DataType

from pyhyperiso.Core import HyperisoConfig, HyperisoMaster, QCDProvider, AlphasConfig, MassConfig

if __name__ == "__main__":

    #All informations for the creation of the HyperisoMaster are provided in the base_core_example.py example.
    config = HyperisoConfig()
    hyp = HyperisoMaster()
    parameters_values = os.path.join(__file__, "..", "inputs", "sm_param_test.yaml")
    
    lha_file_path = "lha/si_input.flha" 
    hyp.init(lha_file=lha_file_path, config=config)
    
    #This API allows to run alpha_s and masses.
    qcd_runner = QCDProvider()
    
    #Define the configuration, with the scale at which you want to get alpha_s, and the mass scheme.
    ac = AlphasConfig(81)
    
    #The get_alphas API return alpha_s at 81 GeV.
    print(qcd_runner.get_alphas(ac))

    #Configuration to run masses, please specify the pdg_id of the mass to run. Here 6 is the top quark (pdg id)
    mc = MassConfig(81, pdg_id = 6)
    
    #The get_alphas API return the top mass at 81 GeV.
    print(qcd_runner.get_qcd_masses(mc))