from pyhyperiso.Common import ParamId, ParameterType

from pyhyperiso.Core import HyperisoConfig, HyperisoMaster, ParameterSetter, ParameterProvider

if __name__ == "__main__":

    #All informations for the creation of the HyperisoMaster are provided in the base_core_example.py example.
    config = HyperisoConfig()
    hyp = HyperisoMaster()
    lha_file_path = "lha/si_input.flha" 
    hyp.init(lha_file=lha_file_path, config=config)

    #All informations for the usage of the ParameterProvider are provided in the parameter_provider_example.py example.
    param_provider = ParameterProvider()
    
    #Creation of the ParameterSetter proxy, which is meant for changing value of parameters within Hyperiso.
    #Careful with this API, some parameters are DependentParameters, or inside DependentBlocks, which cannot be modified
    #If you want to modify these kind of parameters, please see the dependency_pruner_example.py example
    param_setter = ParameterSetter()
    
    print("Value before the changes of the Electroweak scale within Hyperiso :", param_provider.get_by_pid(ParamId(ParameterType.WILSON, "EW_SCALE", 1)))
    
    #The following API can change the value of a Parameter. One can use a float (this will change the real part of the parameter) or a Scalar.
    param_setter.mutate(ParamId(ParameterType.WILSON, "EW_SCALE", 1), 160)
    
    print("Value after the changes of the Electroweak scale within Hyperiso :", param_provider.get_by_pid(ParamId(ParameterType.WILSON, "EW_SCALE", 1)))