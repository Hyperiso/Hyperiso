from pyhyperiso.Common import ParameterType, ParamId, DataType

from pyhyperiso.Core import HyperisoConfig, HyperisoMaster, DependencyPruner, DependantBlockInfoProvider, ParameterProvider, ParameterSetter

if __name__ == "__main__":

    #All informations for the creation of the HyperisoMaster are provided in the base_core_example.py example.
    config = HyperisoConfig()
    hyp = HyperisoMaster()
    lha_file_path = "lha/si_input.flha" 
    hyp.init(lha_file=lha_file_path, config=config)
    
    #The DependencyPruner Proxy allows to unlink dependentblocks or parameters from their dependencies.
    #Useful when wanting to scan over a dependentParameter (or parameter inside a DependentBlock).
    dep_pruner = DependencyPruner()
    
    #All informations on the ParameterProvider proxy are given in the parameter_provider_example.py example
    pp = ParameterProvider()
    
    #All informations on the ParameterSetter proxy are given in the parameter_setter_example.py example
    ps = ParameterSetter()
    
    #The DependantBlockInfoProvider class allows to get information on depentblocks (get the source blocks, get the depend blocks of a block, know if a block is a dependentblock).
    print("Is VCKM a dependentblock : ", DependantBlockInfoProvider().is_dependent_block(ParameterType.SM, "VCKM"))
    
    #get_source_blocks will return all blocks used in the calculation of the parameter.
    print("Sources of VCKM : ", DependantBlockInfoProvider().get_source_blocks(ParameterType.SM, "VCKM"))
    
    print("Current VCKM[1,1] value: ", pp.get_by_pid(ParamId(ParameterType.SM, "VCKM", [1,1])))
    
    
    
    dep_pruner.detach_block(ParameterType.SM, "VCKM")
    
    #VCKM is still a dependentblock (can be reatach later).
    print("Is VCKM a dependentblock : ", DependantBlockInfoProvider().is_dependent_block(ParameterType.SM, "VCKM"))
    
    #But it has no source for the moment.
    print("Sources of VCKM after detach : ", DependantBlockInfoProvider().get_source_blocks(ParameterType.SM, "VCKM"))
    ps.mutate(ParamId(ParameterType.SM, "VCKM", [1,1]), 1)
    
    print("Changed VCKM[1,1] value: ", pp.get_by_pid(ParamId(ParameterType.SM, "VCKM", [1,1])))
    
    #This allows to reattach a block to its sources, and will get back its original value.
    dep_pruner.reattach_block(ParameterType.SM, "VCKM")
    
    print("Changed VCKM[1,1] value after reatach: ", pp.get_by_pid(ParamId(ParameterType.SM, "VCKM", [1,1])))