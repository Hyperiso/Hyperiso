from pyhyperiso.phyperiso.pyhyperiso import common
from pyhyperiso.core.Common.GeneralEnum import ParameterType

class ParamId:
    def __init__(self, type : ParameterType, block : str, code : int):
        self.Paramid = common.ParamId(type.value, block, code)
        # self.Paramid.type = type
        # self.Paramid.block = block
        # self.Paramid.type = type