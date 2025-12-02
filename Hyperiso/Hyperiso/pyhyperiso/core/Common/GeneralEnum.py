from enum import Enum
from pyhyperiso.phyperiso.pyhyperiso.common import Model as _CppModel
from pyhyperiso.phyperiso.pyhyperiso.common import ParameterType as _CppParameterType
from pyhyperiso.phyperiso.pyhyperiso.common import Observables as _CppObservables
from pyhyperiso.phyperiso.pyhyperiso.common import QCDOrder as _CppQCDOrder
from pyhyperiso.phyperiso.pyhyperiso.common import WGroup as _CppWGroup
from pyhyperiso.phyperiso.pyhyperiso.common import WilsonBasis as _CppWilsonBasis
from pyhyperiso.phyperiso.pyhyperiso.common import WCoef as _CppWCoef
from pyhyperiso.phyperiso.pyhyperiso.common import Decays as _CppDecays
from pyhyperiso.phyperiso.pyhyperiso.common import MassType as _CppMassType
from pyhyperiso.phyperiso.pyhyperiso.common import ScaleType as _CppScaleType
from pyhyperiso.phyperiso.pyhyperiso.common import ContributionType  as _CppContributionType 
from pyhyperiso.phyperiso.pyhyperiso.common import DataType  as _CppDataType
from pyhyperiso.phyperiso.pyhyperiso.common import UncertaintyType as _CppUncertaintyType

class Model(Enum):
    SM = _CppModel.SM
    SUSY = _CppModel.SUSY
    THDM = _CppModel.THDM
    MARTY = _CppModel.MARTY


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
    C_V1_bc = _CppWCoef.C_V1_bc
    C_V2_bc = _CppWCoef.C_V2_bc
    C_S1_bc = _CppWCoef.C_S1_bc
    C_S2_bc = _CppWCoef.C_S2_bc
    C_T_bc = _CppWCoef.C_T_bc
    C_V1_bu = _CppWCoef.C_V1_bu
    C_V2_bu = _CppWCoef.C_V2_bu
    C_S1_bu = _CppWCoef.C_S1_bu
    C_S2_bu = _CppWCoef.C_S2_bu
    C_T_bu = _CppWCoef.C_T_bu
    C_V1_cs = _CppWCoef.C_V1_cs
    C_V2_cs = _CppWCoef.C_V2_cs
    C_S1_cs = _CppWCoef.C_S1_cs
    C_S2_cs = _CppWCoef.C_S2_cs
    C_T_cs = _CppWCoef.C_T_cs
    C_V1_cd = _CppWCoef.C_V1_cd
    C_V2_cd = _CppWCoef.C_V2_cd
    C_S1_cd = _CppWCoef.C_S1_cd
    C_S2_cd = _CppWCoef.C_S2_cd
    C_T_cd = _CppWCoef.C_T_cd
    C_V1_su = _CppWCoef.C_V1_su
    C_V2_su = _CppWCoef.C_V2_su
    C_S1_su = _CppWCoef.C_S1_su
    C_S2_su = _CppWCoef.C_S2_su
    C_T_su = _CppWCoef.C_T_su
    C_V1_du = _CppWCoef.C_V1_du
    C_V2_du = _CppWCoef.C_V2_du
    C_S1_du = _CppWCoef.C_S1_du
    C_S2_du = _CppWCoef.C_S2_du
    C_T_du = _CppWCoef.C_T_du
    C_BD_1 = _CppWCoef.C_BD_1
    CT_BD_1 = _CppWCoef.CT_BD_1
    C_BD_2 = _CppWCoef.C_BD_2
    CT_BD_2 = _CppWCoef.CT_BD_2
    C_BD_3 = _CppWCoef.C_BD_3
    CT_BD_3 = _CppWCoef.CT_BD_3
    C_BD_4 = _CppWCoef.C_BD_4
    C_BD_5 = _CppWCoef.C_BD_5
    
    C_BS_1 = _CppWCoef.C_BS_1
    CT_BS_1 = _CppWCoef.CT_BS_1
    C_BS_2 = _CppWCoef.C_BS_2
    CT_BS_2 = _CppWCoef.CT_BS_2
    C_BS_3 = _CppWCoef.C_BS_3
    CT_BS_3 = _CppWCoef.CT_BS_3
    C_BS_4 = _CppWCoef.C_BS_4
    C_BS_5 = _CppWCoef.C_BS_5
    
    C_SD_1 = _CppWCoef.C_SD_1
    CT_SD_1 = _CppWCoef.CT_SD_1
    C_SD_2 = _CppWCoef.C_SD_2
    CT_SD_2 = _CppWCoef.CT_SD_2
    C_SD_3 = _CppWCoef.C_SD_3
    CT_SD_3 = _CppWCoef.CT_SD_3
    C_SD_4 = _CppWCoef.C_SD_4
    C_SD_5 = _CppWCoef.C_SD_5
    
    C_CU_1 = _CppWCoef.C_CU_1
    CT_CU_1 = _CppWCoef.CT_CU_1
    C_CU_2 = _CppWCoef.C_CU_2
    CT_CU_2 = _CppWCoef.CT_CU_2
    C_CU_3 = _CppWCoef.C_CU_3
    CT_CU_3 = _CppWCoef.CT_CU_3
    C_CU_4 = _CppWCoef.C_CU_4
    C_CU_5 = _CppWCoef.C_CU_5
    
    CK9 = _CppWCoef.CK9
    CPK9 = _CppWCoef.CPK9
    CK10 = _CppWCoef.CK10
    CPK10 = _CppWCoef.CPK10
    CKQ1 = _CppWCoef.CKQ1
    CKQ2 = _CppWCoef.CKQ2
    CPKQ1 = _CppWCoef.CPKQ1
    CPKQ2 = _CppWCoef.CPKQ2
    CK_L = _CppWCoef.CK_L


class WGroup(Enum):
    B = _CppWGroup.B
    BPrime = _CppWGroup.BPrime
    BScalar = _CppWGroup.BScalar
    CC_bc = _CppWGroup.CC_bc
    
class WilsonBasis(Enum):
    STANDARD = _CppWilsonBasis.STANDARD
    TRADITIONAL = _CppWilsonBasis.TRADITIONAL

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
    
class DataType(Enum):
    VALUE = _CppDataType.VALUE
    STD_STAT = _CppDataType.STD_STAT
    STD_SYST = _CppDataType.STD_SYST
    STD_COMBINED = _CppDataType.STD_COMBINED
    
class UncertaintyType(Enum):
    STAT = _CppUncertaintyType.STAT
    SYST = _CppUncertaintyType.SYST
    COMBINED = _CppUncertaintyType.COMBINED