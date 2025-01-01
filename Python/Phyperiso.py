from enum import Enum
from phyperiso.pyhyperiso import core, wilson, observable
from phyperiso.pyhyperiso.core import Model as _CppModel
from phyperiso.pyhyperiso.core import ParameterType as _CppParameterType
from phyperiso.pyhyperiso.observable import Observables as _CppObservables
class Model(Enum):
    SM = _CppModel.SM
    CUSTOM = _CppModel.CUSTOM


class ParameterType(Enum):
    SM = _CppParameterType.SM
    FLAVOR = _CppParameterType.FLAVOR


class Observables(Enum):
    BR_BS_MUMU = _CppObservables.BR_BS_MUMU
    BR_BS_MUMU_UNTAG = _CppObservables.BR_BS_MUMU_UNTAG
    BR_BD_MUMU = _CppObservables.BR_BD_MUMU
    BR_BU_TAUNU = _CppObservables.BR_BU_TAUNU
    ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA = _CppObservables.ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA


class MemoryManager:
    """Interface for managing memory and caching."""

    def __init__(self):
        self._manager = core.MemoryManager.get_instance()

    def init(self, lha_file: str, model: Model = Model.SM, is_spectrum: bool = False, has_wilsons: bool = False, has_obs: bool = False):
        """
        Initialize the memory manager.

        :param lha_file: Path to the LHA file.
        :param model: The model to use (default is Model.SM).
        :param is_spectrum: Whether the data includes a spectrum.
        :param has_wilsons: Whether the data includes Wilson coefficients.
        :param has_obs: Whether the data includes observables.
        """
        self._manager.init(lha_file, model.value, is_spectrum, has_wilsons, has_obs)

    def get_input_lha_path(self) -> str:
        """Get the path to the input LHA file."""
        return self._manager.get_input_lha_path()

    def get_data(self):
        """Get the underlying LhaReader instance."""
        return self._manager.get_data()


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

    def set_block_value(self, block: str, code: int, value: float, force: bool = False):
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

    def shift_parameter(self, param_id, shift_value: float):
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

class WilsonManager:
    """Interface for managing Wilson coefficients."""

    def __init__(self, model_name: str = "SM"):
        self._manager = wilson.coefficient_manager.CoefficientManager.get_instance(model_name)

    def initialize(self, lha_file: str, model: Model = Model.SM):
        """
        Initialize the Wilson coefficient manager.

        :param lha_file: Path to the LHA file.
        :param model: The model to use.
        """
        self._manager.initialize(lha_file, model.value)

    def get_state(self) -> str:
        """Get the current state of the Wilson coefficient manager."""
        return self._manager.get_state()

    def set_q_match(self, value: float):
        """Set the matching scale."""
        self._manager.set_q_match(value)


class ObservableInterface:
    """Interface for computing observables."""

    def __init__(self):
        self._interface = observable.ObservableInterface()

    def compute_observable(self, observable: Observables) -> float:
        """
        Compute the value of a given observable.

        :param observable: The observable to compute.
        :return: The computed value.
        """
        return self._interface.compute_observable(observable.value)

    def set_param(self, block: str, code: int, value: float, param_type: ParameterType):
        """
        Set a parameter value.

        :param block: The block name.
        :param code: The parameter code.
        :param value: The parameter value.
        :param param_type: The parameter type.
        """
        self._interface.set_param(block, code, value, param_type.value)