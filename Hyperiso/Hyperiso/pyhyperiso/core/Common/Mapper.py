from pyhyperiso.phyperiso.pyhyperiso.common import GroupMapper as _CppGroupMapper
from pyhyperiso.core.Common.GeneralEnum import WGroup
class GroupMapper:
    def __init__(self):
        pass
    
    def str(self, group_id : WGroup):
        return _CppGroupMapper.str(group_id.value)
    
    def id_of(self, group_id: WGroup):
        return _CppGroupMapper.id_of(self.str(group_id))
if __name__ == "__main__":
    gm = GroupMapper()
    
    print(gm.str(WGroup.B))
    
    print(gm.id_of(WGroup.B))
    print(type(gm.id_of(WGroup.B)))