from enum import Enum
from pyhyperiso.phyperiso.pyhyperiso.common import Model as _CppModel
from pyhyperiso.phyperiso.pyhyperiso.common import ParameterType as _CppParameterType
from pyhyperiso.phyperiso.pyhyperiso.common import Observables as _CppObservables
from pyhyperiso.phyperiso.pyhyperiso.common import QCDOrder as _CppQCDOrder
from pyhyperiso.phyperiso.pyhyperiso.common import WGroup as _CppWGroup
from pyhyperiso.phyperiso.pyhyperiso.common import BWilsonBasis as _CppBWilsonBasis
from pyhyperiso.phyperiso.pyhyperiso.common import WCoef as _CppWCoef
from pyhyperiso.phyperiso.pyhyperiso.common import Decays as _CppDecays
from pyhyperiso.phyperiso.pyhyperiso.common import MassType as _CppMassType
from pyhyperiso.phyperiso.pyhyperiso.common import ScaleType as _CppScaleType
from pyhyperiso.phyperiso.pyhyperiso.common import ContributionType  as _CppContributionType 

class Model(Enum):
    SM = _CppModel.SM
    SUSY = _CppModel.SUSY
    THDM = _CppModel.THDM
    CUSTOM = _CppModel.CUSTOM


class ParameterType(Enum):
    SM = _CppParameterType.SM
    BSM = _CppParameterType.BSM
    FLAVOR = _CppParameterType.FLAVOR
    WILSON = _CppParameterType.WILSON
    DECAY = _CppParameterType.DECAY
    PASSTHROUGH = _CppParameterType.PASSTHROUGH
    OBSERVABLE = _CppParameterType.OBSERVABLE

class QCDOrder(Enum):
    NONE = _CppQCDOrder.NONE
    LO = _CppQCDOrder.LO
    NLO = _CppQCDOrder.NLO
    NNLO = _CppQCDOrder.NNLO


class WCoeff(Enum):
    C1 = _CppWCoef.C1
    C2 = _CppWCoef.C2
    C3 = _CppWCoef.C3
    C4 = _CppWCoef.C4
    C5 = _CppWCoef.C5
    C6 = _CppWCoef.C6
    C7 = _CppWCoef.C7
    C8 = _CppWCoef.C8
    C9 = _CppWCoef.C9
    C10 = _CppWCoef.C10
    CQ1 = _CppWCoef.CQ1
    CQ2 = _CppWCoef.CQ2
    CP1 = _CppWCoef.CP1
    CP2 = _CppWCoef.CP2
    CP3 = _CppWCoef.CP3
    CP4 = _CppWCoef.CP4
    CP5 = _CppWCoef.CP5
    CP6 = _CppWCoef.CP6
    CP7 = _CppWCoef.CP7
    CP8 = _CppWCoef.CP8
    CP9 = _CppWCoef.CP9
    CP10 = _CppWCoef.CP10
    CPQ1 = _CppWCoef.CPQ1
    CPQ2 = _CppWCoef.CPQ2
    CBlnu_A = _CppWCoef.CBlnu_A
    CBlnu_P = _CppWCoef.CBlnu_P
    C_V1 = _CppWCoef.C_V1
    C_V2 = _CppWCoef.C_V2
    C_S1 = _CppWCoef.C_S1
    C_S2 = _CppWCoef.C_S2
    C_T = _CppWCoef.C_T


class WGroup(Enum):
    B = _CppWGroup.B
    BPrime = _CppWGroup.BPrime
    BScalar = _CppWGroup.BScalar
    Blnu = _CppWGroup.Blnu
    BCLNU = _CppWGroup.BCLNU
    
class BWilsonBasis(Enum):
    STANDARD = _CppBWilsonBasis.STANDARD
    TRADITIONAL = _CppBWilsonBasis.TRADITIONAL

class Observables(Enum):
    BR_BS_MUMU = _CppObservables.BR_BS_MUMU
    BR_BS_MUMU_UNTAG = _CppObservables.BR_BS_MUMU_UNTAG
    BR_BD_MUMU = _CppObservables.BR_BD_MUMU
    R_TAU_NU = _CppObservables.R_TAU_NU
    BR_BU_TAU_NU = _CppObservables.BR_BU_TAU_NU
    ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA = _CppObservables.ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA
    BR_B_XS_GAMMA = _CppObservables.BR_B_XS_GAMMA
    BR_B__D_TAU_NU = _CppObservables.BR_B__D_TAU_NU
    A_FB_B__D_TAU_NU = _CppObservables.A_FB_B__D_TAU_NU
    P_TAU_B__D_TAU_NU = _CppObservables.P_TAU_B__D_TAU_NU
    R_D = _CppObservables.R_D
    BR_B__DSTAR_TAU_NU = _CppObservables.BR_B__DSTAR_TAU_NU
    A_FB_B__DSTAR_TAU_NU = _CppObservables.A_FB_B__DSTAR_TAU_NU
    P_TAU_B__DSTAR_TAU_NU = _CppObservables.P_TAU_B__DSTAR_TAU_NU
    P_D_B__DSTAR_TAU_NU = _CppObservables.P_D_B__DSTAR_TAU_NU
    R_DSTAR = _CppObservables.R_DSTAR
    
class Decays(Enum):
    B__D_l_nu = _CppDecays.B__D_l_nu
    B__Dstar_l_nu = _CppDecays.B__Dstar_l_nu
    B__Kstar = _CppDecays.B__Kstar
    B__l_l = _CppDecays.B__l_l
    B__l_nu = _CppDecays.B__l_nu
    B__Xs = _CppDecays.B__Xs
    
class MassType(Enum):
    POLE = _CppMassType.POLE
    MSBAR = _CppMassType.MSBAR
    
class ScaleType(Enum):
    MATCHING = _CppScaleType.MATCHING
    HADRONIC = _CppScaleType.HADRONIC
    
class ContributionType(Enum):
    SM = _CppContributionType.SM
    BSM = _CppContributionType.BSM
    TOTAL = _CppContributionType.TOTAL
    