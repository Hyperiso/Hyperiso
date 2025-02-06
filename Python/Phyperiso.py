from enum import Enum
from phyperiso.pyhyperiso import core, wilson, observable
from phyperiso.pyhyperiso.core import Model as _CppModel
from phyperiso.pyhyperiso.core import ParameterType as _CppParameterType
from phyperiso.pyhyperiso.core import Observables as _CppObservables
from phyperiso.pyhyperiso.core import QCDOrder as _CppQCDOrder
from phyperiso.pyhyperiso.core import WGroup as _CppWGroup
from phyperiso.pyhyperiso.core import BWilsonBasis as _CppBWilsonBasis
from phyperiso.pyhyperiso.core import WCoef as _CppWCoef
from phyperiso.pyhyperiso.core import ModelMapper as _CppModelMapper
from phyperiso.pyhyperiso.core import ParameterTypeMapper as _CppParameterTypeMapper
from phyperiso.pyhyperiso.core import WCoefMapper as _CppWCoefMapper
from phyperiso.pyhyperiso.core import GroupMapper as _CppGroupMapper
from phyperiso.pyhyperiso.core import OrderMapper as _CppOrderMapper
from phyperiso.pyhyperiso.core import ObservableMapper as _CppObservableMapper


class Model(Enum):
    SM = _CppModel.SM
    SUSY = _CppModel.SUSY
    THDM = _CppModel.THDM
    CUSTOM = _CppModel.CUSTOM


class ParameterType(Enum):
    SM = _CppParameterType.SM
    SUSY = _CppParameterType.SUSY
    THDM = _CppParameterType.THDM
    CUSTOM = _CppParameterType.CUSTOM
    FLAVOR = _CppParameterType.FLAVOR
    WILSON = _CppParameterType.WILSON
    FF = _CppParameterType.FF

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



class WGroup(Enum):
    B = _CppWGroup.B
    BPrime = _CppWGroup.BPrime
    BScalar = _CppWGroup.BScalar
    
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
    XI__D_L_NU = _CppObservables.XI__D_L_NU

class ParamId:
    def __init__(self, type : ParameterType, block : str, code : int):
        self.Paramid = core.ParamId(type.value, block, code)
        # self.Paramid.type = type
        # self.Paramid.block = block
        # self.Paramid.type = type
 

class ModelMapper:
    @staticmethod
    def get_model_str_list():
        return _CppModelMapper.get_str()
    
    @staticmethod
    def get_model_enum_list():
        return _CppModelMapper.get_enum()
    
    @staticmethod
    def to_str(model : Model):
        return _CppModelMapper.str(model.value)
    
    @staticmethod
    def enum_elt(model : str):
        return _CppModelMapper.enum_elt(model)
    
class ParameterTypeMapper:
    @staticmethod
    def get_model_str_list():
        return _CppParameterTypeMapper.get_str()
    
    @staticmethod
    def get_model_enum_list():
        return _CppParameterTypeMapper.get_enum()
    
    @staticmethod
    def to_str(model : ParameterType):
        return _CppParameterTypeMapper.str(model.value)
    
    @staticmethod
    def enum_elt(model : str):
        return _CppParameterTypeMapper.enum_elt(model)
    
class WCoefMapper:
    @staticmethod
    def get_model_str_list():
        return _CppWCoefMapper.get_str()
    
    @staticmethod
    def get_model_enum_list():
        return _CppWCoefMapper.get_enum()
    
    @staticmethod
    def to_str(model : WCoeff):
        return _CppWCoefMapper.str(model.value)
    
    @staticmethod
    def enum_elt(model : str):
        return _CppWCoefMapper.enum_elt(model)
    
    @staticmethod
    def get_group(group : WGroup):
        return _CppWCoefMapper.get_group(group.value)
    
class GroupMapper:
    @staticmethod
    def get_model_str_list():
        return _CppGroupMapper.get_str()
    
    @staticmethod
    def get_model_enum_list():
        return _CppGroupMapper.get_enum()
    
    @staticmethod
    def to_str(model : WGroup):
        return _CppGroupMapper.str(model.value)
    
    @staticmethod
    def enum_elt(model : str):
        return _CppGroupMapper.enum_elt(model)
    
class OrderMapper:
    @staticmethod
    def get_model_str_list():
        return _CppOrderMapper.get_str()
    
    @staticmethod
    def get_model_enum_list():
        return _CppOrderMapper.get_enum()
    
    @staticmethod
    def to_str(model : QCDOrder):
        return _CppOrderMapper.str(model.value)
    
    @staticmethod
    def enum_elt(model : str):
        return _CppOrderMapper.enum_elt(model)
    
class ObservableMapper:
    @staticmethod
    def get_model_str_list():
        return _CppObservableMapper.get_str()
    
    @staticmethod
    def get_model_enum_list():
        return _CppObservableMapper.get_enum()
    
    @staticmethod
    def to_str(model : Observables):
        return _CppObservableMapper.str(model.value)
    
    @staticmethod
    def enum_elt(model : str):
        return _CppObservableMapper.enum_elt(model)
    
class MemoryManager:
    """Interface for managing memory and caching."""

    def __init__(self):
        self._manager = core.MemoryManager.get_instance()

    def init(self, lha_file: str, model: Model = Model.SM, use_marty : bool = False, is_spectrum: bool = False, has_wilsons: bool = False, has_obs: bool = False):
        """
        Initialize the memory manager.

        :param lha_file: Path to the LHA file.
        :param model: The model to use (default is Model.SM).
        :param is_spectrum: Whether the user want to use Marty for the wilson calculation.
        :param is_spectrum: Whether the data includes a spectrum.
        :param has_wilsons: Whether the data includes Wilson coefficients.
        :param has_obs: Whether the data includes observables.
        """
        self._manager.init(lha_file, model.value, use_marty, is_spectrum, has_wilsons, has_obs)

    def get_input_lha_path(self) -> str:
        """Get the path to the input LHA file."""
        return self._manager.get_input_lha_path()

    def get_data(self):
        """Get the underlying LhaReader instance (DO NOT USE RIGHT NOW)."""
        return self._manager.get_data()

    def switch_lha(self, lhaFile : str, model : Model, use_marty : bool = False, is_spectrum : bool = False, has_wilsons : bool = False, has_obs : bool = False):
        self._manager.switch_lha(lhaFile, model.value, use_marty, is_spectrum, has_wilsons, has_obs)

    def switch_model(self, model : Model, use_marty : bool = False):
        self._manager.switch_model(model.value, use_marty)

    def get_blocks_list(self, paramtype : ParameterType = ParameterType.SM):
        return self._manager.get_blocks_list(paramtype.value)
    
    def get_block_infos(self, block :str, paramtype : ParameterType = ParameterType.SM):
        return self._manager.get_block_infos(block, paramtype.value)
    
class Parameters:
    """Interface for handling parameters."""

    def __init__(self, param_type: ParameterType = ParameterType.SM):
        """
        Initialize a Parameters object for a specific parameter type.

        :param param_type: The type of parameters (default is ParameterType.SM).
        """
        self._parameters = core.Parameters.get_instance(param_type.value)

    def alpha_s(self, q: float) -> float:
        """
        Compute the strong coupling constant at a given scale.

        :param q: The scale in GeV.
        :return: The value of alpha_s at the given scale.
        """
        return self._parameters.alpha_s(q)

    def running_mass(self, quark_mass: float, q_init: float, q_end: float,
                     option_massb: str = "running", option_masst: str = "pole") -> float:
        """
        Compute the running mass of a quark.

        :param quark_mass: Initial quark mass.
        :param q_init: Initial scale.
        :param q_end: Final scale.
        :param option_massb: Option for bottom mass (default "running").
        :param option_masst: Option for top mass (default "pole").
        :return: The running mass at the final scale.
        """
        return self._parameters.running_mass(quark_mass, q_init, q_end, option_massb, option_masst)

    def set_block_value(self, block: str, code: int, value: float, force: bool = True):
        """
        Set a parameter value in a specific block.

        :param block: The block name.
        :param code: The parameter code.
        :param value: The parameter value.
        :param force: Whether to force the update.
        """
        self._parameters.set_block_value(block, code, value, force)

    def get_qcd_mass(self, masstype: str) -> float:
        """
        Retrieve a QCD mass of a specific type.

        :param masstype: The type of mass to retrieve (e.g., "mt_mt").
        :return: The QCD mass of the specified type.
        """
        return self._parameters.get_qcd_masse(masstype)

    def exists(self, block: str, code: int) -> bool:
        """
        Check if a parameter exists in a block.

        :param block: The block name.
        :param code: The parameter code.
        :return: True if the parameter exists, otherwise False.
        """
        return self._parameters.exists(block, code)

    def shift_parameter(self, param_id : int, shift_value: float) -> None:
        """
        Shift a parameter value by a specified amount.

        :param param_id: The parameter ID.
        :param shift_value: The value by which to shift the parameter.
        """
        self._parameters.shift_parameter(param_id, shift_value)

    def __call__(self, block: str, code: int) -> float:
        """
        Retrieve the value of a parameter in a block using the callable interface.

        :param block: The block name.
        :param code: The parameter code.
        :return: The parameter value.
        """
        return self._parameters(block, code)

class QCDHelper:
    @staticmethod
    def mass_b_1S():
        return core.QCDHelper.mass_b_1S()
    

class BCoefficientGroup:
    def __init__(self):
        self.group = wilson.coefficient_groups.BCoefficientGroup()

class BPrimeCoefficientGroup:
    def __init__(self):
        self.group = wilson.coefficient_groups.BPrimeCoefficientGroup()

class BScalarCoefficientGroup:
    def __init__(self):
        self.group = wilson.coefficient_groups.BScalarCoefficientGroup()

class BCoefficientSUSY:
    def __init__(self):
        self.group = wilson.coefficient_groups.BCoefficientGroup_susy()

class BCoefficientTHDM:
    def __init__(self):
        self.group = wilson.coefficient_groups.BCoefficientGroup_THDM()

class BPrimeCoefficientSUSY:
    def __init__(self):
        self.group = wilson.coefficient_groups.BPrimeCoefficientGroup_susy()

class BPrimeCoefficientTHDM:
    def __init__(self):
        self.group = wilson.coefficient_groups.BPrimeCoefficientGroup_THDM()

class BScalarCoefficientSUSY:
    def __init__(self):
        self.group = wilson.coefficient_groups.BScalarCoefficientGroup_susy()

class BScalarCoefficientTHDM:
    def __init__(self):
        self.group = wilson.coefficient_groups.BScalarCoefficientGroup_THDM()


class WilsonManager:
    """Interface for managing Wilson coefficients."""

    def __init__(self, model_name: str = "SM"):
        self._manager = wilson.coefficient_manager.CoefficientManager.get_instance()

    def initialize(self, lha_file: str, model: Model = Model.SM, use_marty : bool = False, is_spectrum : bool = False, has_wilsons : bool = False, has_obs : bool = False):
        """
        Initialize the Wilson coefficient manager.

        :param lha_file: Path to the LHA file.
        :param model: The model to use.
        :use marty: Whether the user want to use Marty or not.
        :is spectrum: Is the lhaFile containing spectrum or not.
        :has wilson: Does the lhaFile already contains wilson coefficients.
        :has obs: Does the lhaFile already contains observables.
        """
        self._manager.initialize(lha_file, model.value, use_marty, is_spectrum, has_wilsons, has_obs)

    def register_coefficient_group(self, groupName : str, coeff_group) -> None:
        self._manager.register_coefficient_group(groupName, coeff_group.group)

    def get_state(self) -> str:
        """Get the current state of the Wilson coefficient manager."""
        return self._manager.get_state("BPrimeCoefficientGroup")

    def set_q_match(self, groupName : str, value: float) -> None:
        """Set the matching scale."""
        self._manager.set_q_match(groupName, value)

    def set_params(self, group : str, block : str, pdgcode : int, value : float) -> None:
        """
            Set a params (SM only) 
            :block: lhablock for the element the user want to change
            :pdgcode: pdgcode of the element the user want to change
            :value: value of the element after the change
        """
        self._manager.set_params(group, block, pdgcode, value)

    def get_params(self, block : str, pdgcode : int):
        """
            Get a param (SM only) 
            :block: lhablock for the element the user want to get
            :pdgcode: pdgcode of the element the user want to get
        """
        return self._manager.get_params(block, pdgcode)
    
    def set_group_scale(self, groupName : str, Q : float) -> None:
        """
            Set The running scale for a wilson coefficients group.

            :groupName: the name of the group.
            :Q: The running energy.
        """
        self._manager.set_group_scale(groupName, Q)
    
    def set_matching_coefficient(self, groupName : str, QCDorder : str):
        self._manager.set_matching_coefficient(groupName, QCDorder)

    def set_run_coefficient(self, groupName :str, QCDorder : str):
        self._manager.set_run_coefficient(groupName, QCDorder)

    def get_matching_coefficient(self, groupName : str, coeffName : str, QCDorder : str):
        return self._manager.get_matching_coefficient(groupName, coeffName, QCDorder, False)
    
    def get_run_coefficient(self, groupName : str, coeffName : str, QCDorder : str):
        return self._manager.get_run_coefficient(groupName, coeffName, QCDorder, False)
    
    def get_coefficient_group(self, groupname : str):
        return self._manager.get_coefficient_group(groupname)
    
class WilsonInterface:
    def __init__(self):
        self._wilsoninterface = wilson.wilson_interface.WilsonInterface()

    def set_q_match(self):
        self._wilsoninterface.set_q_match()

    def set_matching_coefficient(self):
        self._wilsoninterface.set_matching_coefficient()

    def set_group_scale(self):
        self._wilsoninterface.set_group_scale()
    
    def set_run_coefficient(self):
        self._wilsoninterface.set_run_coefficient()

    def get_matching_coefficient(self):
        return self._wilsoninterface.get_matching_coefficient()
    
    def get_full_matching_coefficient(self, group : WGroup, coeff : WCoeff, order : QCDOrder):
        return self._wilsoninterface.get_full_matching_coefficient(group.value, coeff.value, order.value, False)
    
    def get_full_run_coefficient(self, group : WGroup, coeff : WCoeff, order : QCDOrder):
        return self._wilsoninterface.get_full_run_coefficient(group.value, coeff.value, order.value, False)
    
    def build(self, group : list, Q_match : float, Q : float, order : QCDOrder):
        self._wilsoninterface.build(group, Q_match, Q, order.value)



class ObservableInterface:
    """Interface for computing observables."""

    def __init__(self):
        self._interface = observable.ObservableInterface()

    def add_observable(self, obs : Observables, qcd_order : QCDOrder) -> None:
        self._interface.add_observable(obs.value, qcd_order.value, False)
    
    def add_observables(self, obss : list) -> None:
        self._interface.add_observables(obss, False)

    def add_observable_parameter(self, obs : Observables, pid) -> None:
        self._interface.add_observable_parameter(obs.value, pid)

    def add_observable_parameters(self, obs : Observables, pids : list) -> None:
        self._interface.add_observable_parameters(obs.value, pids)

    def compute_observable(self, observable: Observables) -> float:
        """
        Compute the value of a given observable.

        :param observable: The observable to compute.
        :return: The computed value.
        """
        return self._interface.compute_observable(observable.value)

    def compute_all_observables(self) -> dict:
        return self._interface.compute_all_observables()
    
    def compute_uncertainty(self, obs : Observables) -> float:
        return self._interface.compute_uncertainty(obs.value)
    
    def compute_leading_uncertainties(self,obs : Observables, n : int) -> dict:
        return self._interface.compute_leading_uncertainties(obs.value, n)
    
    def compute_all_uncertainties(self):
        return self._interface.compute_all_uncertainties()
    
    def compute_chi2(self):
        return self._interface.compute_chi2()
    
    def set_param(self, block: str, code: int, value: float, param_type: ParameterType):
        """
        Set a parameter value.

        :param block: The block name.
        :param code: The parameter code.
        :param value: The parameter value.
        :param param_type: The parameter type.
        """
        self._interface.set_param(block, code, value, param_type.value)

    def get_param(self, block : str, code : int) -> float:
        return self._interface.get_param(block, code)
    
    def get_current_obss(self) -> list:
        return self._interface.get_current_obss()