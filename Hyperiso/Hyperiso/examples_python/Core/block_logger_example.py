from pyhyperiso.Common import ParameterType

from pyhyperiso.Core import HyperisoConfig, HyperisoMaster
from pyhyperiso.Core import BlockLogger

if __name__ == "__main__":

    #All informations for the creation of the HyperisoMaster are provided in the base_core_example.py example.
    config = HyperisoConfig()
    hyp = HyperisoMaster()
    lha_file_path = "lha/si_input.flha" 
    hyp.init(lha_file=lha_file_path, config=config)

    #This API allows to search within Hyperiso to get informations on the differents parameters.
    #Can be very helpful when facing a strange result and need to know which parameters (SM inputs, form factors, etc.) are used inside the pipeline.
    block_logger = BlockLogger()
    
    #The parameters are divided in multiple sections : 
    #   -ParameterType.SM : Standard Model inputs, within the lha, such as GAUGE, MASS, SMINPUTS, etc.
    #   -ParameterType.BSM : All Model inputs which are not SM. Custom block go in this category.
    #   -ParameterType.DECAY : Decays's inputs, such as form factors, etc.
    #   -ParameterType.FLAVOR : Flavor's inputs, such as FCONST, FLIFE, etc.
    #   -ParameterType.WILSON : Scales inputs, Wilson Coefficients (after calculating them inside the pipeline)
    #   -ParameterType.OBSERVABLE : Experimental Observables Values.
    
    print("\nSMINPUTS : \n")
    #For example, one can get the content of the block using:
    print(block_logger.get_block(ParameterType.SM, "SMINPUTS"))
    
    print("\nWILSON BLOCKS : \n")
    #One can also get all the content within a section:
    print(block_logger.get_all_blocks(ParameterType.WILSON))