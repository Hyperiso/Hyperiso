from pyhyperiso.phyperiso.pyhyperiso.common import OrderMapper as _CppOrderMapper
from pyhyperiso.phyperiso.pyhyperiso.common import ParameterTypeMapper as _CppParameterTypeMapper
from pyhyperiso.phyperiso.pyhyperiso.common import ModelMapper as _CppModelMapper
from pyhyperiso.phyperiso.pyhyperiso.common import WilsonBasisMapper as _CppWilsonBasisMapper
from pyhyperiso.phyperiso.pyhyperiso.common import ContributionTypeMapper as _CppContributionTypeMapper
from pyhyperiso.phyperiso.pyhyperiso.common import MassTypeMapper as _CppMassTypeMapper
from pyhyperiso.phyperiso.pyhyperiso.common import ScaleTypeMapper as _CppScaleTypeMapper
from pyhyperiso.phyperiso.pyhyperiso.common import GroupMapper as _CppGroupMapper
from pyhyperiso.phyperiso.pyhyperiso.common import WCoefMapper as _CppWCoefMapper
from pyhyperiso.phyperiso.pyhyperiso.common import ObservableMapper as _CppObservableMapper
from pyhyperiso.phyperiso.pyhyperiso.common import DecayMapper as _CppDecayMapper
from pyhyperiso.core.Common.GeneralEnum import WGroup, Observables, QCDOrder, ParameterType, Model, WilsonBasis, ContributionType, MassType, ScaleType, Decays
from pyhyperiso.core.Common.LhaID import LhaID
from pyhyperiso.core.Common.SymbolId import ObservableId, _unwrap_optional
from typing import AnyStr

class OrderMapper:
    def __init__(self):
        pass
    
    def str(self, obs_id : QCDOrder):
        return _CppOrderMapper.str(obs_id.value)
    
    def id_of(self, obs_id: QCDOrder):
        return _CppOrderMapper.id_of(self.str(obs_id))
    
    def get_str(self):
        return _CppOrderMapper.get_str()
    
    def get_str_all(self):
        return _CppOrderMapper.get_str_all()
    
    def get_enum(self):
        return [x for x in _CppOrderMapper.get_enum()]
    
class ParameterTypeMapper:
    def __init__(self):
        pass
    
    def str(self, obs_id : ParameterType):
        return _CppParameterTypeMapper.str(obs_id.value)
    
    def id_of(self, obs_id: ParameterType):
        return _CppParameterTypeMapper.id_of(self.str(obs_id))
    
    def get_str(self):
        return _CppParameterTypeMapper.get_str()
    
    def get_str_all(self):
        return _CppParameterTypeMapper.get_str_all()
    
    def get_enum(self):
        return [x for x in _CppParameterTypeMapper.get_enum()]
    
class ModelMapper:
    def __init__(self):
        pass
    
    def str(self, obs_id : Model):
        return _CppModelMapper.str(obs_id.value)
    
    def id_of(self, obs_id: Model):
        return _CppModelMapper.id_of(self.str(obs_id))
    
    def get_str(self):
        return _CppModelMapper.get_str()
    
    def get_str_all(self):
        return _CppModelMapper.get_str_all()
    
    def get_enum(self):
        return [x for x in _CppModelMapper.get_enum()]
    
class WilsonBasisMapper:
    def __init__(self):
        pass
    
    def str(self, obs_id : WilsonBasis):
        return _CppWilsonBasisMapper.str(obs_id.value)
    
    def id_of(self, obs_id: WilsonBasis):
        return _CppWilsonBasisMapper.id_of(self.str(obs_id))
    
    def get_str(self):
        return _CppWilsonBasisMapper.get_str()
    
    def get_str_all(self):
        return _CppWilsonBasisMapper.get_str_all()
    
    def get_enum(self):
        return [x for x in _CppWilsonBasisMapper.get_enum()]
    
class ContributionTypeMapper:
    def __init__(self):
        pass
    
    def str(self, obs_id : ContributionType):
        return _CppContributionTypeMapper.str(obs_id.value)
    
    def id_of(self, obs_id: ContributionType):
        return _CppContributionTypeMapper.id_of(self.str(obs_id))
    
    def get_str(self):
        return _CppContributionTypeMapper.get_str()
    
    def get_str_all(self):
        return _CppContributionTypeMapper.get_str_all()
    
    def get_enum(self):
        return [x for x in _CppContributionTypeMapper.get_enum()]
    
class MassTypeMapper:
    def __init__(self):
        pass
    
    def str(self, obs_id : MassType):
        return _CppMassTypeMapper.str(obs_id.value)
    
    def id_of(self, obs_id: MassType):
        return _CppMassTypeMapper.id_of(self.str(obs_id))
    
    def get_str(self):
        return _CppMassTypeMapper.get_str()
    
    def get_str_all(self):
        return _CppMassTypeMapper.get_str_all()
    
    def get_enum(self):
        return [x for x in _CppMassTypeMapper.get_enum()]
    
class ScaleTypeMapper:
    def __init__(self):
        pass
    
    def str(self, obs_id : ScaleType):
        return _CppScaleTypeMapper.str(obs_id.value)
    
    def id_of(self, obs_id: ScaleType):
        return _CppScaleTypeMapper.id_of(self.str(obs_id))
    
    def get_str(self):
        return _CppScaleTypeMapper.get_str()
    
    def get_str_all(self):
        return _CppScaleTypeMapper.get_str_all()
    
    def get_enum(self):
        return [x for x in _CppScaleTypeMapper.get_enum()]
    

    
    
class GroupMapper:
    def __init__(self):
        pass
    
    def str(self, group_id : WGroup):
        return _CppGroupMapper.str(group_id.value)
    
    def id_of(self, group_id: WGroup):
        return _CppGroupMapper.id_of(self.str(group_id))
    
    def get_str(self):
        return _CppGroupMapper.get_str()
    
    def get_str_all(self):
        return _CppGroupMapper.get_str_all()
    
    def get_enum(self):
        return [x for x in _CppGroupMapper.get_enum()]
    
    
class WCoefMapper:
    def __init__(self):
        pass
    
    def str(self, obs_id : Observables):
        return _CppWCoefMapper.str(obs_id.value)
    
    def id_of(self, obs_id: Observables):
        return _CppWCoefMapper.id_of(self.str(obs_id))
    
    def get_str(self):
        return _CppWCoefMapper.get_str()
    
    def get_str_all(self):
        return _CppWCoefMapper.get_str_all()
    
    def get_enum(self):
        return [x for x in _CppWCoefMapper.get_enum()]
    
class ObservableMapper:
    def __init__(self):
        pass
    
    @staticmethod
    def str(obs_id : Observables):
        return _CppObservableMapper.str(obs_id.value)
    
    @staticmethod
    def id_of(name: AnyStr) -> ObservableId:
        # return _CppObservableMapper.id_of(self.str(obs_id))
        cpp_id = _CppObservableMapper.id_of(name)  # retourne _CppObservableId
        return ObservableId(str(cpp_id))           # wrap python pur
    
    @staticmethod
    def to_id(obs : Observables) -> ObservableId:
        cpp_id = _CppObservableMapper.to_id(obs.value)
        return ObservableId(str(cpp_id))
    
    @staticmethod
    def canonical(obs_id: ObservableId) -> str:
        return _CppObservableMapper.canonical(obs_id._to_cpp())
    
    @staticmethod
    def get_str():
        return _CppObservableMapper.get_str()
    
    @staticmethod
    def get_str_all():
        return _CppObservableMapper.get_str_all()
    
    @staticmethod
    def get_enum():
        # si _CppObservableMapper.get_enum() renvoie une liste d'enums C++ (ou d'int)
        return [Observables(x) for x in _CppObservableMapper.get_enum()]
    
    @staticmethod
    def from_flha(flha_id: LhaID) -> ObservableId:
        # binding C++: from_flha(LhaID) -> optional<ObservableId> -> None ou _CppObservableId
        cpp_opt = _CppObservableMapper.from_flha(flha_id._cpp_obj)
        cpp_id = _unwrap_optional(cpp_opt, "ObservableId from FLHA")
        return ObservableId(str(cpp_id))
    
    @staticmethod
    def flha(obs) -> LhaID:
        """
        obs peut être:
          - Observables  -> appelle flha(Observables)
          - ObservableId -> appelle flha(ObservableId) (celle qui throw si pas d'ext)
        """
        if isinstance(obs, Observables):
            cpp_lha = _CppObservableMapper.flha(obs.value)  # returns cpp LhaID
            return LhaID(cpp_lha)                  # à adapter à ton wrapper LhaID
        elif isinstance(obs, ObservableId):
            cpp_lha = _CppObservableMapper.flha(obs._to_cpp())
            return LhaID(cpp_lha)
        else:
            raise TypeError("flha() expects Observables or ObservableId")
    
class DecayMapper:
    def __init__(self):
        pass
    
    def str(self, obs_id : Decays):
        return _CppDecayMapper.str(obs_id.value)
    
    def id_of(self, obs_id: Decays):
        return _CppDecayMapper.id_of(self.str(obs_id))
    
    def get_str(self):
        return _CppDecayMapper.get_str()
    
    def get_str_all(self):
        return _CppDecayMapper.get_str_all()
    
    def get_enum(self):
        return [x for x in _CppDecayMapper.get_enum()]
    
    def get_observables(self, decay : Decays):
        cpp_list = _CppDecayMapper.get_observables(decay.value)
        return [Observables(o) for o in cpp_list]
    
    def get_decay(self, obs : Observables):
        cpp_decay = _CppDecayMapper.get_decay(obs.value) 
        return Decays(cpp_decay)
    
    
    

if __name__ == "__main__":
    
    ordermapper = OrderMapper()
    
    print(ordermapper.get_str())
    print(ordermapper.get_str_all())
    print(ordermapper.get_enum())
    
    gm = GroupMapper()
    
    print(gm.str(WGroup.B))
    
    print(gm.id_of(WGroup.B))
    print(type(gm.id_of(WGroup.B)))
    
    om = ObservableMapper()
    
    print(om.str(Observables.A_FB_B__D_TAU_NU))
    
    print(om.flha(Observables.BR_BU_TAU_NU))
    print(type(om.flha(Observables.BR_BU_TAU_NU)))
    
    decaymapper = DecayMapper()
    

    print(decaymapper.get_observables(Decays.B__D_l_nu))
    print(decaymapper.get_decay(Observables.BR_D__MU_NU))
    
    wcoefmapper = WCoefMapper()
    
