from Python.Phyperiso import MemoryManager, Parameters

class MemoryManagerCache:
    def __init__(self, lha_input, model, use_marty = False, is_spectrum = False, has_wilson = False, has_obs = False):
        self.mm = MemoryManager()
        self.mm.init(lha_input, model, use_marty, is_spectrum, has_wilson, has_obs)
        self.lha = lha_input

    def switch_lha(self, lha_input, model = False, use_marty = False, is_spectrum = False, has_wilsons = False, has_obs = False):
        self.lha = lha_input
        self.mm.switch_lha(lha_input, model, use_marty, is_spectrum, has_wilsons, has_obs)

    def get_lha(self):
        with open(self.lha, 'r') as f:
            data = f.read()
        return data
    
    def get_blocks_list(self):
        return self.mm.get_blocks_list()
    
    def get_block_infos(self, block):
        return self.mm.get_block_infos(block)
    

class ParametersCache:
    def __init__(self, param_type):
        self.param = Parameters(param_type)

    def exists(self, block : str, code : int):
        return self.param.exists(block, code)
    
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
    





    