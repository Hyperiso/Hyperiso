from Hyperiso.Phyperiso import MemoryManager, Parameters, ParameterType
from Hyperiso.API.utils.infos import db_path
import os

class MemoryManagerCache:
    def __init__(self, lha_input, model, use_marty = False, is_spectrum = False, has_wilson = False, has_obs = False):
        self.mm = MemoryManager()
        self.mm.init(os.path.join(db_path, lha_input), model, use_marty, is_spectrum, has_wilson, has_obs)
        # self.lha = lha_input
        # self.model = model
        # self.use_marty = use_marty
        # self.is_spectrum = is_spectrum
        # self.has_wilson = has_wilson
        # self.has_obs = has_obs
        self.infos = {"lha" : os.path.join(db_path, lha_input), "model" : model, "use_marty" : use_marty, "is_spectrum" : is_spectrum, "has_wilson" : has_wilson,
                      "has_obs" : has_obs}
        self.model_available = {"SM"}
    def switch_info(self, inputs):
        new_infos = {}
        for elem in ["lha", "model", "use_marty", "is_spectrum", "has_wilson", "has_obs"]:
            if inputs[elem] is not None:
                new_infos[elem] = inputs[elem]
            else:
                new_infos[elem] = self.infos[elem]

        self.mm.switch_lha(os.path.join(db_path,new_infos["lha"]), new_infos["model"], new_infos["use_marty"], new_infos["is_spectrum"],
                            new_infos["has_wilson"], new_infos["has_obs"])
        self.infos = new_infos.copy()
        print("I was here")

    def get_lha(self):
        with open(self.infos["lha"], 'r') as f:
            data = f.read()
        return data
    
    def get_blocks_list(self, param_type = ParameterType.SM):
        return self.mm.get_blocks_list(param_type)
    
    def get_block_infos(self, block, param_type = ParameterType.SM):
        return self.mm.get_block_infos(block, param_type)
    
    def get_parameters_types(self):
        return self.mm.get_parameters_types()
    
    def get_type_of_block(self, block : str) ->list:
        return self.mm.get_type_of_block(block)

class ParametersCache:
    def __init__(self, param_type):
        self.param = Parameters(param_type)

    def exists(self, block : str, code : int):
        print(block, code)
        return Parameters().exists(block, code)
    
    def __call__(self, block : str, code : int):
        return self.param(block, code)
    
    def switch_param(self, param_type):
        self.param = Parameters(param_type)
    
    def alpha_s(self, q):
        return Parameters().alpha_s(q)
    
    def running_mass(self, quark_mass, q_init, q_end, option_massb = "running", option_masst = "pole"):
        return Parameters().running_mass(quark_mass, q_init,q_end, option_massb, option_masst)
    
    def set_block_value(self, block, pdgcode, value):
        self.param.set_block_value(block, pdgcode, value)

    def get_qcd_mass(self, masstype):
        return Parameters().get_qcd_mass(masstype)


    





    