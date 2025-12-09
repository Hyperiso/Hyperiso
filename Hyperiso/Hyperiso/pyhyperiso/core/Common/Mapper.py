from pyhyperiso.phyperiso.pyhyperiso.common import GroupMapper as _CppGroupMapper
from pyhyperiso.phyperiso.pyhyperiso.common import ObservableMapper as _CppObservableMapper
from pyhyperiso.core.Common.GeneralEnum import WGroup, Observables

class GroupMapper:
    def __init__(self):
        pass
    
    def str(self, group_id : WGroup):
        return _CppGroupMapper.str(group_id.value)
    
    def id_of(self, group_id: WGroup):
        return _CppGroupMapper.id_of(self.str(group_id))
    
    
class ObservableMapper:
    def __init__(self):
        pass
    
    def str(self, obs_id : Observables):
        return _CppObservableMapper.str(obs_id.value)
    
    def id_of(self, obs_id: Observables):
        return _CppObservableMapper.id_of(self.str(obs_id))
    
    
if __name__ == "__main__":
    gm = GroupMapper()
    
    print(gm.str(WGroup.B))
    
    print(gm.id_of(WGroup.B))
    print(type(gm.id_of(WGroup.B)))
    
    om = ObservableMapper()
    
    print(om.str(Observables.A_FB_B__D_TAU_NU))
    
    print(om.id_of(Observables.A_FB_B__D_TAU_NU))
    print(type(om.id_of(Observables.A_FB_B__D_TAU_NU)))