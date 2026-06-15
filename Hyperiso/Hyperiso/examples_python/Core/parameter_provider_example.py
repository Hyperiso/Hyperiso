from pyhyperiso.Common import ParameterType, ParamId

from pyhyperiso.Core import HyperisoConfig, HyperisoMaster, ParameterProvider

if __name__ == "__main__":

    #All informations for the creation of the HyperisoMaster are provided in the base_core_example.py example.
    config = HyperisoConfig()
    hyp = HyperisoMaster()
    lha_file_path = "lha/si_input.flha" 
    hyp.init(lha_file=lha_file_path, config=config)

    #This create a parameterprovider for SM parameters. This allows to get parameters without specify each time the parameterType
    sm_param_provider = ParameterProvider(ParameterType.SM)
    
    #If you want parameters from differents sections, use a no-typed ParameterProvider:
    no_pt_param_provider = ParameterProvider()
    
    #Using the sm_param_provider, one can test if a parameter exist using the following API:
    print(sm_param_provider.exists_by_block("SMINPUTS", 2))
    
    #Using the no_pt_param_provider, one can test if a parameter exist using the following API:
    print(no_pt_param_provider.exists_by_pid(ParamId(ParameterType.SM, "SMINPUTS", 2)))